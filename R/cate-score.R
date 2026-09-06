# R/cate-score.R
#
# Out-of-sample scoring of CATE models with the doubly robust loss, and
# stacking of several models into one (best single model, convex weights,
# Q-aggregation, or unconstrained least squares).

#' Compare CATE models out of sample with the doubly robust loss
#'
#' For each candidate CATE model `tau_m`, computes the doubly robust loss
#' `L(tau_m) = mean((Y(eta) - tau_m(X))^2)` on the rows of `scores`, which
#' should be a held-out scoring sample. Differences of losses between two
#' models estimate differences of their mean squared errors against the true
#' CATE, are Neyman orthogonal to the nuisances, and are asymptotically normal
#' (Chernozhukov et al. 2026, Theorem 15.2.1), so each model is reported with
#' its loss difference from the baseline and a confidence interval. The
#' normalized score `1 - L(tau_m) / L(constant)` measures the improvement over
#' the constant-effect model.
#'
#' @param scores A `cm_scores` object on the scoring sample.
#' @param ... Named candidate models: `cm_cate` objects (predicted on
#'   `scores$data`), numeric vectors of predictions, or column names of
#'   `scores$data`. Avoid names that partially match `scores` (such as `s`
#'   or `sc`); R would assign them to that argument.
#' @param baseline `"constant"` (default), which uses the constant model
#'   `ate`, or the name of one of the candidates.
#' @param ate Constant effect used by the constant model; by default the mean
#'   of the score on the scoring sample. Pass the training-sample ATE when
#'   the scoring sample is small.
#' @param conf_level Confidence level.
#'
#' @return A list of class `cm_cate_score` with `table` (model, loss,
#'   diff (baseline loss minus model loss; positive means better than the
#'   baseline), std.error, conf.low, conf.high, score_norm), the
#'   `predictions` matrix, the `baseline`, and `ate`.
#' @examples
#' dat <- sim_hte(1500, dgp = "smooth", seed = 1)
#' x <- paste0("x", 1:5)
#' train <- dat[1:1000, ]; test <- dat[1001:1500, ]
#' m_dr <- cate_learner(train, "y", "d", x, x_het = c("x1", "x2"), method = "dr", seed = 1)
#' m_t <- cate_learner(train, "y", "d", x, x_het = c("x1", "x2"), method = "t", seed = 1)
#' sc_test <- dr_scores(test, "y", "d", x, seed = 2)
#' cate_score(sc_test, dr = m_dr, t = m_t)
#' @seealso [cate_ensemble()], [cate_validate()]
#' @export
cate_score <- function(scores, ..., baseline = "constant", ate = NULL, conf_level = 0.95) {
  .cm_check_scores(scores)
  models <- list(...)
  if (length(models) == 0L) stop("Supply at least one candidate model.", call. = FALSE)
  if (is.null(names(models)) || any(names(models) == "")) {
    names(models) <- ifelse(is.null(names(models)) | names(models) == "",
                            paste0("model", seq_along(models)), names(models))
  }
  y <- scores$score
  n <- scores$n
  if (is.null(ate)) ate <- mean(y)
  P <- vapply(names(models), function(nm) .cm_cate_predictions(models[[nm]], scores, nm), numeric(n))
  P <- matrix(P, nrow = n, dimnames = list(NULL, names(models)))
  loss_const <- mean((y - ate)^2)
  losses <- colMeans((y - P)^2)
  ref <- if (identical(baseline, "constant")) (y - ate)^2 else {
    if (!baseline %in% names(models)) stop("`baseline` must be \"constant\" or a candidate name.", call. = FALSE)
    (y - P[, baseline])^2
  }
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  rows <- lapply(names(models), function(nm) {
    diff_i <- ref - (y - P[, nm])^2
    dlt <- mean(diff_i)
    se <- stats::sd(diff_i) / sqrt(n)
    data.frame(model = nm, loss = losses[[nm]], diff = dlt, std.error = se,
               conf.low = dlt - crit * se, conf.high = dlt + crit * se,
               score_norm = 1 - losses[[nm]] / loss_const, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  structure(list(table = tab, predictions = P, baseline = baseline, ate = ate,
                 loss_constant = loss_const, n = n, conf_level = conf_level),
            class = "cm_cate_score")
}

#' @export
print.cm_cate_score <- function(x, ...) {
  cat("Doubly robust loss of CATE models (n = ", x$n, ")\n", sep = "")
  cat("  constant model loss = ", format(round(x$loss_constant, 4)),
      "; diff = baseline loss - model loss (positive is better)\n", sep = "")
  print(x$table, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Stack CATE models on a scoring sample
#'
#' Combines several fitted CATE models into one by minimizing the doubly
#' robust loss on the scoring sample. `"best"` picks the single model with the
#' smallest loss; `"convex"` finds weights on the simplex; `"q"`
#' (Q-aggregation, Lecué and Rigollet 2014; Lan and Syrgkanis 2023) adds the
#' penalty `sum_m w_m L(tau_m)`, which yields the optimal `log(M) / n` rate;
#' `"ols"` is unconstrained least squares. Base predictions are centered
#' and the intercept is fixed at `ate` (Chernozhukov et al. 2026,
#' Remark 15.2.3).
#'
#' @inheritParams cate_score
#' @param ... Named `cm_cate` objects (they must be able to predict on new
#'   data).
#' @param method `"q"`, `"convex"`, `"best"`, or `"ols"`.
#' @param ate Intercept of the stacked model (default the mean score on the
#'   scoring sample; the training-sample ATE is preferable when the scoring
#'   sample is small).
#' @param max_iter,tol Frank-Wolfe settings for the simplex problems.
#'
#' @return A `cm_cate` object with `method = "ensemble"`, whose `ensemble`
#'   element records the rule, weights, base models, centering means, and
#'   losses. [predict()] combines the base models on new data.
#' @examples
#' dat <- sim_hte(1500, dgp = "smooth", seed = 1)
#' x <- paste0("x", 1:5)
#' train <- dat[1:1000, ]; test <- dat[1001:1500, ]
#' m_dr <- cate_learner(train, "y", "d", x, x_het = c("x1", "x2"), method = "dr", seed = 1)
#' m_r <- cate_learner(train, "y", "d", x, x_het = c("x1", "x2"), method = "r", seed = 1)
#' sc_test <- dr_scores(test, "y", "d", x, seed = 2)
#' ens <- cate_ensemble(sc_test, dr = m_dr, r = m_r, method = "q")
#' ens
#' @references
#' Lecué, G. and Rigollet, P. (2014). Optimal learning with Q-aggregation.
#' *Annals of Statistics*, 42(1), 211-224.
#'
#' Lan, H. and Syrgkanis, V. (2023). Causal Q-aggregation for CATE model
#' selection. arXiv:2310.16945.
#' @export
cate_ensemble <- function(scores, ..., method = c("q", "convex", "best", "ols"),
                          ate = NULL, max_iter = 5000L, tol = 1e-10) {
  method <- match.arg(method)
  .cm_check_scores(scores)
  models <- list(...)
  if (length(models) == 0L) stop("Supply at least one `cm_cate` model.", call. = FALSE)
  if (!all(vapply(models, inherits, logical(1), "cm_cate"))) {
    stop("All base models must be `cm_cate` objects so the ensemble can predict on new data.", call. = FALSE)
  }
  if (is.null(names(models)) || any(names(models) == "")) names(models) <- paste0("model", seq_along(models))
  y <- scores$score
  n <- scores$n
  if (is.null(ate)) ate <- mean(y)
  P <- vapply(names(models), function(nm) .cm_cate_predictions(models[[nm]], scores, nm), numeric(n))
  P <- matrix(P, nrow = n, dimnames = list(NULL, names(models)))
  means <- colMeans(P)
  Pc <- sweep(P, 2, means)
  target <- y - ate
  losses <- colMeans((target - Pc)^2)
  M <- ncol(P)
  intercept <- ate
  if (method == "best") {
    w <- as.numeric(seq_len(M) == which.min(losses))
  } else if (method == "convex") {
    w <- .cm_simplex_qp(Pc, target, max_iter = max_iter, tol = tol)$weights
  } else if (method == "q") {
    w <- .cm_simplex_qp(Pc, target, lin = losses, max_iter = max_iter, tol = tol)$weights
  } else {
    fit <- stats::lm.fit(cbind(1, Pc), target)
    cf <- stats::coef(fit)
    cf[is.na(cf)] <- 0
    intercept <- ate + cf[1]
    w <- as.numeric(cf[-1])
  }
  names(w) <- names(models)
  tau_hat <- as.numeric(intercept + Pc %*% w)
  structure(list(
    tau_hat = tau_hat, label = y, weights = NULL, final = NULL,
    x_het = unique(unlist(lapply(models, `[[`, "x_het"))), method = "ensemble",
    scores = scores,
    ensemble = list(rule = method, weights = w, models = models, means = means,
                    intercept = intercept, losses = losses,
                    loss = mean((target - Pc %*% w)^2)),
    learners = list(final = paste0("stack:", method)), n = n, call = match.call()
  ), class = "cm_cate")
}

.cm_predict_ensemble <- function(object, newdata) {
  e <- object$ensemble
  P <- vapply(names(e$models), function(nm) as.numeric(stats::predict(e$models[[nm]], newdata)),
              numeric(nrow(as.data.frame(newdata))))
  P <- matrix(P, ncol = length(e$models))
  as.numeric(e$intercept + sweep(P, 2, e$means) %*% e$weights)
}
