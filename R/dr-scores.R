# R/dr-scores.R
#
# Cross-fitted doubly robust pseudo-outcomes: the one object that every
# heterogeneous-effect and policy-learning function in the package consumes.

#' Cross-fitted doubly robust pseudo-outcomes
#'
#' `dr_scores()` builds, for every observation, the doubly robust (AIPW)
#' signal
#' \deqn{Y_i(\hat\eta) = \hat\mu_1(Z_i) - \hat\mu_0(Z_i) + H_i(\hat p)\,(Y_i - \hat\mu_{D_i}(Z_i)),
#' \qquad H_i(p) = \frac{D_i}{p(Z_i)} - \frac{1 - D_i}{1 - p(Z_i)},}
#' whose conditional mean given any function of the covariates is the
#' conditional average treatment effect, `E[Y(eta0) | X] = tau(X)`, under
#' conditional exogeneity. Its sample mean is the AIPW estimate of the ATE.
#' The same vector is the regression label of the DR-learner
#' ([cate_learner()]), the outcome of best linear predictors
#' ([cate_blp()], [cate_gate()]), the loss of model selection
#' ([cate_score()], [cate_ensemble()]), the input of every validation
#' statistic ([cate_validate()]), and the reward of policy evaluation and
#' learning ([policy_value()], [policy_learn()]).
#'
#' Nuisances are cross-fitted with `mlr3` learners (as in [est_dml()]) or
#' supplied as out-of-fold predictions. The object also stores the residuals
#' `Y - l(Z)` and `D - p(Z)`, with `l(Z) = p mu1 + (1 - p) mu0`, for the
#' R-learner.
#'
#' @param data A data frame.
#' @param y,d Outcome and binary treatment column names.
#' @param x Character vector of covariate names (the controls `Z`). Needed
#'   when the nuisances are estimated internally, and used as the default
#'   heterogeneity variables downstream.
#' @param p_hat,mu0_hat,mu1_hat Optional out-of-fold predictions of
#'   `P(D = 1 | Z)`, `E[Y | D = 0, Z]`, `E[Y | D = 1, Z]`: numeric vectors of
#'   length `nrow(data)` or column names. When all three are supplied no
#'   learner is fitted.
#' @param learner_p `mlr3` classification learner for the propensity score
#'   (default logistic regression).
#' @param learner_mu `mlr3` regression learner for the two outcome regressions
#'   (default linear regression); `learner_mu0`, `learner_mu1` override it
#'   arm by arm.
#' @param learner_mu0,learner_mu1 Optional arm-specific learners.
#' @param type `"dr"` (default) for the doubly robust signal, `"ipw"` for
#'   `Y H(p)`, `"reg"` for `mu1 - mu0`. The alternatives exist to show why
#'   the doubly robust signal is preferred.
#' @param folds Number of cross-fitting folds.
#' @param fold_id Optional fold identifier (vector or column name) recording
#'   the folds of supplied nuisances or imposing a partition.
#' @param seed Optional seed for the fold partition.
#' @param p_clip Propensity clipping bounds.
#' @param outcome_type `"auto"`, `"continuous"`, or `"binary"`; matters only
#'   when a classification learner predicts the outcome.
#' @param na_action `"fail"` or `"omit"` for rows with missing values.
#'
#' @return A list of class `cm_scores`: `score` (the pseudo-outcome),
#'   `nuisance` (data frame with `p`, `mu0`, `mu1`, `l`), `residuals`
#'   (`y_tilde`, `d_tilde`), `y`, `d`, `x` (covariate names), `data` (the
#'   retained rows), `fold_id`, `ate` (estimate and standard error),
#'   `diagnostics`, and learner labels. [tidy()] returns the ATE row.
#'
#' @references
#' Chernozhukov, V., Hansen, C., Kallus, N., Spindler, M., and Syrgkanis, V.
#' (2026). *Applied Causal Inference Powered by ML and AI*, chapters 14-15.
#'
#' Kennedy, E. H. (2023). Towards optimal doubly robust estimation of
#' heterogeneous causal effects. *Electronic Journal of Statistics*, 17(2).
#'
#' @examples
#' dat <- sim_hte(n = 600, dgp = "smooth", seed = 1)
#' sc <- dr_scores(dat, y = "y", d = "d", x = paste0("x", 1:5), seed = 1)
#' sc
#' # the mean of the score is the AIPW estimate of the ATE
#' tidy(sc)
#' @seealso [cate_learner()], [cate_blp()], [cate_validate()], [policy_learn()]
#' @export
dr_scores <- function(data, y, d, x = NULL,
                      p_hat = NULL, mu0_hat = NULL, mu1_hat = NULL,
                      learner_p = NULL, learner_mu = NULL,
                      learner_mu0 = NULL, learner_mu1 = NULL,
                      type = c("dr", "ipw", "reg"),
                      folds = 5L, fold_id = NULL, seed = NULL,
                      p_clip = c(0.01, 0.99),
                      outcome_type = c("auto", "continuous", "binary"),
                      na_action = c("fail", "omit")) {
  type <- match.arg(type)
  outcome_type <- match.arg(outcome_type)
  na_action <- match.arg(na_action)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data)
  if (!.cm_is_string(y) || !.cm_is_string(d)) stop("`y` and `d` must be column names.", call. = FALSE)
  .cm_check_column(y, data)
  .cm_check_column(d, data)
  if (!is.null(x)) {
    if (!is.character(x)) stop("`x` must be a character vector of column names.", call. = FALSE)
    for (v in x) .cm_check_column(v, data)
  }
  n_all <- nrow(data)
  p_in <- .cm_get_optional_numeric(p_hat, data, n_all, "p_hat")
  mu0_in <- .cm_get_optional_numeric(mu0_hat, data, n_all, "mu0_hat")
  mu1_in <- .cm_get_optional_numeric(mu1_hat, data, n_all, "mu1_hat")
  fold_in <- .cm_get_optional_vector(fold_id, data, n_all, "fold_id")

  keep <- stats::complete.cases(data[, c(y, d, x), drop = FALSE])
  for (v in list(p_in, mu0_in, mu1_in)) if (!is.null(v)) keep <- keep & is.finite(v)
  if (!all(keep)) {
    if (na_action == "fail") stop("Missing values in the required columns; use na_action = \"omit\".", call. = FALSE)
    data <- data[keep, , drop = FALSE]
    p_in <- p_in[keep]; mu0_in <- mu0_in[keep]; mu1_in <- mu1_in[keep]
    fold_in <- fold_in[keep]
  }
  n <- nrow(data)
  yv <- as.numeric(data[[y]])
  dv <- .cm_as_binary(data[[d]], d)
  if (sum(dv == 1L) < 2L || sum(dv == 0L) < 2L) stop("Both treatment arms need at least two observations.", call. = FALSE)

  need_p <- is.null(p_in) && type != "reg"
  need_mu0 <- is.null(mu0_in) && type != "ipw"
  need_mu1 <- is.null(mu1_in) && type != "ipw"
  need_any <- need_p || need_mu0 || need_mu1
  if (need_any && is.null(x)) stop("`x` is required when nuisances are estimated internally.", call. = FALSE)

  if (!is.null(fold_in)) {
    fold_id <- as.integer(factor(fold_in))
  } else if (need_any) {
    fold_id <- .cm_make_folds(dv, folds, seed)
  } else {
    fold_id <- rep(1L, n)
  }

  if (outcome_type == "auto") {
    outcome_type <- if (all(yv %in% c(0, 1))) "binary" else "continuous"
  }
  work <- data
  work$.cm_y <- yv
  work$.cm_d <- dv
  labels <- list(p = "supplied", mu0 = "supplied", mu1 = "supplied")
  if (need_any) {
    .cm_require_mlr3()
    if (is.null(learner_mu0)) learner_mu0 <- learner_mu
    if (is.null(learner_mu1)) learner_mu1 <- learner_mu
    fit <- .cm_learn_nuisances(work, x, fold_id, cross_fit = TRUE,
                               need_p = need_p, need_mu0 = need_mu0, need_mu1 = need_mu1,
                               learner_p = learner_p, learner_mu0 = learner_mu0,
                               learner_mu1 = learner_mu1, outcome_type = outcome_type)
    if (need_p) { p_in <- fit$p_hat; labels$p <- .cm_learner_label(learner_p %||% .cm_default_learner("classif")) }
    if (need_mu0) { mu0_in <- fit$mu0_hat; labels$mu0 <- .cm_learner_label(learner_mu0 %||% .cm_default_learner("regr")) }
    if (need_mu1) { mu1_in <- fit$mu1_hat; labels$mu1 <- .cm_learner_label(learner_mu1 %||% .cm_default_learner("regr")) }
  }
  if (type == "reg") p_in <- rep(mean(dv), n)
  if (type == "ipw") { mu0_in <- rep(0, n); mu1_in <- rep(0, n) }
  p_raw <- p_in
  p <- .cm_clip(p_in, p_clip[1], p_clip[2])
  mu0 <- as.numeric(mu0_in)
  mu1 <- as.numeric(mu1_in)
  H <- dv / p - (1 - dv) / (1 - p)
  mu_d <- ifelse(dv == 1L, mu1, mu0)
  score <- switch(type,
    dr = (mu1 - mu0) + H * (yv - mu_d),
    ipw = H * yv,
    reg = mu1 - mu0
  )
  l <- p * mu1 + (1 - p) * mu0
  ate <- list(estimate = mean(score), std.error = stats::sd(score) / sqrt(n))

  diagnostics <- list(
    propensity = .cm_summary(p_raw),
    clipped_share = mean(p_raw < p_clip[1] | p_raw > p_clip[2]),
    common_support = if (type != "reg") .cm_common_support(p_raw, dv) else NULL,
    nuisance = rbind(
      p = if (type != "reg") .cm_fit_quality(dv, p_raw) else c(rmse = NA, r2 = NA),
      mu0 = if (type != "ipw") .cm_fit_quality(yv[dv == 0L], mu0[dv == 0L]) else c(rmse = NA, r2 = NA),
      mu1 = if (type != "ipw") .cm_fit_quality(yv[dv == 1L], mu1[dv == 1L]) else c(rmse = NA, r2 = NA)
    ),
    score = .cm_summary(score)
  )

  structure(list(
    score = as.numeric(score),
    nuisance = data.frame(p = p, mu0 = mu0, mu1 = mu1, l = l),
    residuals = data.frame(y_tilde = yv - l, d_tilde = dv - p),
    y = yv, d = dv, x = x, y_name = y, d_name = d,
    data = data, fold_id = fold_id, n = n, type = type,
    ate = ate, learners = labels, outcome_type = outcome_type,
    p_clip = p_clip, diagnostics = diagnostics, call = match.call()
  ), class = "cm_scores")
}

#' @export
print.cm_scores <- function(x, ...) {
  cat("Doubly robust pseudo-outcomes (", x$type, ")\n", sep = "")
  cat("  n = ", x$n, ", treated share = ", round(mean(x$d), 3),
      ", folds = ", length(unique(x$fold_id)), "\n", sep = "")
  cat("  nuisances: p = ", x$learners$p, ", mu0 = ", x$learners$mu0,
      ", mu1 = ", x$learners$mu1, "\n", sep = "")
  cat("  ATE (mean of the score) = ", format(round(x$ate$estimate, 4)),
      " (SE ", format(round(x$ate$std.error, 4)), ")\n", sep = "")
  cat("  score quantiles (1%, 50%, 99%): ",
      paste(format(round(x$diagnostics$score[c("q01", "median", "q99")], 3)), collapse = ", "), "\n", sep = "")
  invisible(x)
}
