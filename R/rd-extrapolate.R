# R/rd-extrapolate.R
#
# Treatment effects away from the cutoff under conditional independence
# given covariates (Angrist and Rokkanen 2015).

#' Extrapolate regression discontinuity effects away from the cutoff
#'
#' The RD estimand is the effect at the cutoff. Angrist and Rokkanen (2015)
#' extrapolate under a conditional independence assumption (CIA): within a
#' window around the cutoff, potential outcomes are mean independent of the
#' running variable given covariates `Z`. Two consequences follow.
#'
#' * **A test.** On each side of the cutoff, the outcome should not depend
#'   on the running variable once `Z` is controlled for; the function
#'   reports the Wald test of the running-variable terms in a regression of
#'   `y` on `x`, `x^2`, and `Z` per side, with HC1 standard errors.
#' * **Extrapolation.** The untreated mean of units above the cutoff is
#'   `E[mu_0(Z) | above]` with `mu_0` fitted below, so the average effect on
#'   the treated side is `mean(y - mu_0_hat(Z))` over units above the
#'   cutoff, and symmetrically for the average effect on the untreated
#'   side. `method = "aipw"` estimates the same objects doubly robustly
#'   with [est_aipw()] using the side indicator as the treatment; `method =
#'   "regression"` uses outcome regressions only. A curve of effects along
#'   the running variable is obtained by a kernel-weighted local linear
#'   regression of the residual `y - mu_other(Z)` on `x` at a grid of
#'   points.
#'
#' @param data A data frame.
#' @param y,x Outcome and running variable column names.
#' @param cutoff Cutoff value.
#' @param covariates Character vector of covariate names (`Z`).
#' @param window Half-width of the window around the cutoff; `NULL` uses the
#'   MSE-optimal bandwidth of `rdrobust::rdbwselect()` when installed and
#'   the interquartile range of `x - cutoff` otherwise.
#' @param method `"regression"` (default) or `"aipw"`.
#' @param learner Optional `mlr3` regression learner for the outcome
#'   regressions (default linear regression).
#' @param learner_p Optional `mlr3` classification learner for the side
#'   propensity (`method = "aipw"`).
#' @param target Grid of running-variable values for the effect curve;
#'   `NULL` uses ten points per side inside the window.
#' @param curve_bandwidth Bandwidth of the local linear curve (default half
#'   the window).
#' @param folds,seed Cross-fitting settings for `method = "aipw"`.
#' @param conf_level Confidence level.
#'
#' @return A list of class `cm_rd_extrapolate` with `test` (per side: Wald
#'   statistic, df, p-value), `effects` (rows `above` and `below`: estimate,
#'   std.error, conf.low, conf.high, n), `curve` (x, estimate, std.error,
#'   conf.low, conf.high, side), `window`, `method`.
#' @references
#' Angrist, J. D. and Rokkanen, M. (2015). Wanna get away? Regression
#' discontinuity estimation of exam school effects away from the cutoff.
#' *Journal of the American Statistical Association*, 110(512), 1331-1344.
#' @examples
#' dat <- sim_rd(3000, "cia", seed = 1)
#' ex <- rd_extrapolate(dat, "y", "x", covariates = "z", window = 0.5)
#' ex
#' @export
rd_extrapolate <- function(data, y, x, cutoff = 0, covariates, window = NULL,
                           method = c("regression", "aipw"), learner = NULL, learner_p = NULL,
                           target = NULL, curve_bandwidth = NULL, folds = 5L, seed = NULL, conf_level = 0.95) {
  method <- match.arg(method)
  pr <- .cm_rd_prepare(data, y, x, cutoff, covariates = covariates)
  if (is.null(window)) {
    window <- if (requireNamespace("rdrobust", quietly = TRUE)) .cm_rd_bandwidth(pr$y, pr$x, cutoff)$h else
      stats::IQR(pr$x - cutoff) / 2
  }
  keep <- abs(pr$x - cutoff) <= window
  w <- pr$data[keep, , drop = FALSE]
  w$.cm_y <- pr$y[keep]
  w$.cm_xc <- pr$x[keep] - cutoff
  w$.cm_above <- as.integer(w$.cm_xc >= 0)
  n_w <- nrow(w)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # 1. CIA test per side ------------------------------------------------------
  test <- do.call(rbind, lapply(c("below", "above"), function(s) {
    ws <- w[if (s == "above") w$.cm_above == 1L else w$.cm_above == 0L, , drop = FALSE]
    ws$.cm_xc2 <- ws$.cm_xc^2
    fml <- stats::reformulate(c(".cm_xc", ".cm_xc2", covariates), ".cm_y")
    m <- estimatr::lm_robust(fml, data = ws, se_type = "HC1")
    b <- stats::coef(m)[c(".cm_xc", ".cm_xc2")]
    V <- stats::vcov(m)[c(".cm_xc", ".cm_xc2"), c(".cm_xc", ".cm_xc2")]
    stat <- tryCatch(as.numeric(t(b) %*% solve(V, b)), error = function(e) NA_real_)
    data.frame(side = s, n = nrow(ws), coef_x = unname(b[1]), se_x = sqrt(V[1, 1]),
               statistic = stat, df = 2L, p.value = stats::pchisq(stat, 2, lower.tail = FALSE))
  }))

  # 2. Outcome regressions by side and extrapolated effects ------------------
  fit_side <- function(rows, newrows) {
    if (is.null(learner)) {
      m <- stats::lm(stats::reformulate(covariates, ".cm_y"), data = w[rows, , drop = FALSE])
      list(pred = as.numeric(stats::predict(m, newdata = w[newrows, , drop = FALSE])),
           pred_all = as.numeric(stats::predict(m, newdata = w)),
           vcov = stats::vcov(m), X_new = stats::model.matrix(stats::reformulate(covariates), w[newrows, , drop = FALSE]))
    } else {
      .cm_require_mlr3()
      f <- .cm_fit_mlr3(learner, w[rows, , drop = FALSE], ".cm_y", covariates, task_hint = "rd_extrapolate")
      list(pred = .cm_predict_fit(f, w[newrows, , drop = FALSE]), pred_all = .cm_predict_fit(f, w), vcov = NULL, X_new = NULL)
    }
  }
  above <- w$.cm_above == 1L
  effects <- NULL
  if (method == "regression") {
    m0 <- fit_side(!above, above)      # mu_0 fitted below, predicted above
    m1 <- fit_side(above, !above)      # mu_1 fitted above, predicted below
    eff_above <- w$.cm_y[above] - m0$pred
    eff_below <- m1$pred - w$.cm_y[!above]
    se_of <- function(res, m) {
      v <- stats::var(res) / length(res)
      if (!is.null(m$vcov)) { xb <- colMeans(m$X_new); v <- v + as.numeric(t(xb) %*% m$vcov %*% xb) }
      sqrt(v)
    }
    effects <- data.frame(
      side = c("above", "below"), estimate = c(mean(eff_above), mean(eff_below)),
      std.error = c(se_of(eff_above, m0), se_of(eff_below, m1)), n = c(sum(above), sum(!above)),
      stringsAsFactors = FALSE)
    resid_above <- eff_above; resid_below <- eff_below
  } else {
    w$.cm_d_side <- above
    a1 <- est_aipw(w, ".cm_y", ".cm_above", covariates, estimand = "ATT", learner_p = learner_p,
                   learner_mu0 = learner, learner_mu1 = learner, folds = folds, seed = seed)
    w$.cm_below <- 1L - w$.cm_above
    a0 <- est_aipw(w, ".cm_y", ".cm_below", covariates, estimand = "ATT", learner_p = learner_p,
                   learner_mu0 = learner, learner_mu1 = learner, folds = folds, seed = seed)
    effects <- data.frame(side = c("above", "below"), estimate = c(a1$estimate, -a0$estimate),
                          std.error = c(a1$std.error, a0$std.error), n = c(sum(above), sum(!above)),
                          stringsAsFactors = FALSE)
    m0 <- fit_side(!above, above); m1 <- fit_side(above, !above)
    resid_above <- w$.cm_y[above] - m0$pred; resid_below <- m1$pred - w$.cm_y[!above]
  }
  effects$conf.low <- effects$estimate - crit * effects$std.error
  effects$conf.high <- effects$estimate + crit * effects$std.error

  # 3. Effect curve along the running variable --------------------------------
  if (is.null(curve_bandwidth)) curve_bandwidth <- window / 2
  if (is.null(target)) target <- c(seq(cutoff - window, cutoff, length.out = 11L)[-11L],
                                   seq(cutoff, cutoff + window, length.out = 11L))
  curve <- do.call(rbind, lapply(target, function(x0) {
    s <- if (x0 >= cutoff) "above" else "below"
    xs <- w$.cm_xc[if (s == "above") above else !above] + cutoff
    r <- if (s == "above") resid_above else resid_below
    kw <- .cm_rd_kernel((xs - x0) / curve_bandwidth, "triangular")
    if (sum(kw > 0) < 5L) return(NULL)
    m <- estimatr::lm_robust(r ~ I(xs - x0), weights = kw, se_type = "HC1")
    data.frame(x = x0, side = s, estimate = unname(stats::coef(m)[1]), std.error = unname(m$std.error[1]))
  }))
  if (!is.null(curve)) {
    curve$conf.low <- curve$estimate - crit * curve$std.error
    curve$conf.high <- curve$estimate + crit * curve$std.error
  }
  structure(list(test = test, effects = effects, curve = curve, window = window, method = method,
                 cutoff = cutoff, covariates = covariates, n_window = n_w, conf_level = conf_level,
                 call = match.call()), class = "cm_rd_extrapolate")
}

#' @export
print.cm_rd_extrapolate <- function(x, ...) {
  cat("RD extrapolation under conditional independence (", x$method, "), window = ",
      format(round(x$window, 4)), ", n = ", x$n_window, "\n", sep = "")
  cat("\nCIA test (running-variable terms given covariates, per side):\n")
  print(x$test, digits = 4, row.names = FALSE)
  cat("\nAverage effects away from the cutoff:\n")
  print(x$effects, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Plot extrapolated RD effects along the running variable
#'
#' @param x A `cm_rd_extrapolate` object.
#' @param truth Optional true effect to draw as a dashed line.
#' @return A ggplot object.
#' @export
plot_rd_extrapolate <- function(x, truth = NULL) {
  cv <- x$curve
  g <- ggplot2::ggplot(cv, ggplot2::aes(x = .data$x, y = .data$estimate, colour = .data$side)) +
    ggplot2::geom_vline(xintercept = x$cutoff, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$side), alpha = 0.2, colour = NA) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_colour_manual(values = c(below = "#0072B2", above = "#D55E00"), guide = "none") +
    ggplot2::scale_fill_manual(values = c(below = "#0072B2", above = "#D55E00"), guide = "none") +
    ggplot2::labs(x = "Running variable", y = "Extrapolated effect with 95% CI") +
    ggplot2::theme_minimal(base_size = 11)
  if (!is.null(truth)) g <- g + ggplot2::geom_hline(yintercept = truth, linetype = "dotted")
  g
}
