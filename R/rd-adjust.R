# R/rd-adjust.R
#
# Flexible covariate adjustment for regression discontinuity designs with
# cross-fitted machine learning (Noack, Olma, and Rothe 2024).

#' Flexible covariate adjustment for regression discontinuity designs
#'
#' Subtracts a cross-fitted prediction of the outcome from the covariates,
#' `eta_hat(Z)`, before the local polynomial step. Because the RD estimand
#' is the jump at the cutoff of the conditional mean given `x`, and
#' `E[eta(Z) | x]` is continuous at the cutoff for any function `eta` of
#' the covariates alone, the jump of `Y - eta(Z)` equals the jump of `Y`.
#' The adjustment therefore needs no rate condition on the learner; it only
#' removes outcome variation explained by the covariates and shrinks the
#' variance (Noack, Olma, and Rothe 2024). With a linear learner it is the
#' local linear adjustment of Calonico et al. (2019) fitted globally; with
#' `regr.cv_glmnet` it is the local post-lasso of Kreiss and Rothe (2023).
#'
#' The prediction is cross-fitted, so no observation's outcome enters the
#' model that adjusts it. Two variants: `pooled = TRUE` fits one model of
#' `Y` on `Z` on both sides of the cutoff; `pooled = FALSE` fits one model
#' per side and averages the two predictions, which is the authors'
#' recommended version when the covariate effects differ across sides
#' (both are functions of `Z` only). In fuzzy designs the treatment is
#' adjusted the same way.
#'
#' When `rdrobust` is installed and `fit = TRUE`, the function also runs
#' `rdrobust::rdrobust()` on the raw and the adjusted outcome and reports
#' both, so the precision gain is visible in one call.
#'
#' @param data A data frame.
#' @param y,x Outcome and running variable column names.
#' @param cutoff Cutoff value.
#' @param covariates Character vector of covariate column names.
#' @param d Optional treatment column name (fuzzy design).
#' @param learner `mlr3` regression learner (default linear regression;
#'   `lrn("regr.ranger")`, `lrn("regr.cv_glmnet")`, ... for flexible
#'   adjustment).
#' @param h_train Optional half-width of the window used to train the
#'   learner (`NULL` uses all rows). Predictions are produced for all rows.
#' @param weighted Train with kernel weights `K((x - c) / h_train)`
#'   (requires `h_train`).
#' @param kernel Kernel for the training weights.
#' @param pooled One model on both sides (`TRUE`) or the average of two
#'   side-specific models (`FALSE`).
#' @param folds,seed Cross-fitting folds and seed.
#' @param fit Run `rdrobust` on the raw and adjusted outcomes.
#' @param ... Arguments passed to `rdrobust::rdrobust()` (for instance
#'   `p`, `kernel`, `bwselect`, `cluster`, `vce`, `h`).
#'
#' @return A list of class `cm_rd_adjust` with `data` (the input rows plus
#'   `y_adj` and, for fuzzy designs, `d_adj`), `eta` (the predictions),
#'   `fits` (raw and adjusted `rdrobust` objects when computed),
#'   `comparison` (tidy rows of both fits), `r2` (share of outcome variance
#'   explained by the adjustment on the estimation window), the learner
#'   label and settings.
#' @references
#' Noack, C., Olma, T., and Rothe, C. (2024). Flexible covariate adjustments
#' in regression discontinuity designs. arXiv:2107.07942.
#'
#' Calonico, S., Cattaneo, M. D., Farrell, M. H., and Titiunik, R. (2019).
#' Regression discontinuity designs using covariates. *Review of Economics
#' and Statistics*, 101(3), 442-451.
#'
#' Kreiss, A. and Rothe, C. (2023). Inference in regression discontinuity
#' designs with high-dimensional covariates. *The Econometrics Journal*,
#' 26(2), 105-123.
#' @examples
#' dat <- sim_rd(2000, "covariates", seed = 1)
#' adj <- rd_adjust(dat, "y", "x", covariates = paste0("z", 1:4), seed = 1)
#' adj
#' @export
rd_adjust <- function(data, y, x, cutoff = 0, covariates, d = NULL, learner = NULL,
                      h_train = NULL, weighted = FALSE, kernel = "triangular", pooled = TRUE,
                      folds = 5L, seed = NULL, fit = TRUE, ...) {
  .cm_require_mlr3()
  if (is.null(learner)) learner <- .cm_default_learner("regr")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, covariates = covariates)
  work <- pr$data
  n <- nrow(work)
  work$.cm_y <- pr$y
  work$.cm_x <- pr$x
  if (!is.null(d)) work$.cm_d <- pr$d
  train_rows <- rep(TRUE, n)
  wts <- NULL
  if (!is.null(h_train)) {
    train_rows <- abs(pr$x - cutoff) <= h_train
    if (weighted) wts <- .cm_rd_kernel((pr$x - cutoff) / h_train, kernel)
  }
  fold_id <- .cm_with_seed(seed, .cm_draw_folds(n, .cm_validate_folds(folds, n)))
  right <- pr$x >= cutoff

  predict_target <- function(target) {
    if (pooled) {
      .cm_crossfit_weighted(work, target, covariates, learner, fold_id,
                            subset = train_rows, weights = wts, task_hint = "rd_adjust")
    } else {
      p_r <- .cm_crossfit_weighted(work, target, covariates, learner, fold_id,
                                   subset = train_rows & right, weights = wts, task_hint = "rd_adjust_right")
      p_l <- .cm_crossfit_weighted(work, target, covariates, learner, fold_id,
                                   subset = train_rows & !right, weights = wts, task_hint = "rd_adjust_left")
      0.5 * (p_r + p_l)
    }
  }
  eta_y <- predict_target(".cm_y")
  out <- pr$data
  out$y_adj <- pr$y - eta_y
  eta <- data.frame(eta_y = eta_y)
  if (!is.null(d)) {
    eta_d <- predict_target(".cm_d")
    out$d_adj <- pr$d - eta_d
    eta$eta_d <- eta_d
  }
  r2 <- {
    v <- stats::var(pr$y[train_rows])
    if (is.finite(v) && v > 0) 1 - stats::var(out$y_adj[train_rows]) / v else NA_real_
  }
  fits <- NULL
  comparison <- NULL
  if (fit && requireNamespace("rdrobust", quietly = TRUE)) {
    dots <- list(...)
    raw <- do.call(.cm_rd_fit, c(list(y = pr$y, x = pr$x, cutoff = cutoff, d = pr$d), dots))
    adj <- do.call(.cm_rd_fit, c(list(y = out$y_adj, x = pr$x, cutoff = cutoff,
                                      d = if (is.null(d)) NULL else out$d_adj), dots))
    fits <- list(raw = raw, adjusted = adj)
    comparison <- rbind(cbind(outcome = "raw", .cm_rd_tidy_rdrobust(raw)),
                        cbind(outcome = "adjusted", .cm_rd_tidy_rdrobust(adj)))
  }
  structure(list(
    data = out, eta = eta, fits = fits, comparison = comparison, r2 = r2,
    learner = .cm_learner_label(learner), pooled = pooled, h_train = h_train, weighted = weighted,
    covariates = covariates, y = y, x = x, d = d, cutoff = cutoff, fold_id = fold_id,
    n = n, n_dropped = pr$n_dropped, call = match.call()
  ), class = "cm_rd_adjust")
}

#' @export
print.cm_rd_adjust <- function(x, ...) {
  cat("Flexible covariate adjustment for RD (", x$learner, if (x$pooled) ", pooled" else ", side-averaged", ")\n", sep = "")
  cat("  n = ", x$n, ", covariates: ", paste(x$covariates, collapse = ", "), "\n", sep = "")
  cat("  share of outcome variance removed by the adjustment: ", format(round(x$r2, 3)), "\n", sep = "")
  if (!is.null(x$comparison)) {
    cat("  rdrobust, robust rows:\n")
    print(x$comparison[x$comparison$method == "robust", c("outcome", "estimate", "std.error", "conf.low", "conf.high", "h_left", "n_left", "n_right")],
          digits = 4, row.names = FALSE)
  }
  invisible(x)
}
