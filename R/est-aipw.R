#' Augmented inverse propensity weighting for binary-treatment ATE and ATT
#'
#' `est_aipw()` estimates the average treatment effect (ATE) or the average
#' treatment effect on the treated (ATT) of a binary treatment under selection
#' on observables. The estimator combines outcome regression predictions,
#' propensity-score predictions, and inverse-propensity residual corrections.
#' It can use user-supplied nuisance predictions or fit nuisance functions with
#' `mlr3` learners, including cross-fitting.
#'
#' This function deliberately owns only the causal-estimation layer. Machine
#' learning is used only to produce nuisance predictions. To use Python or any
#' external prediction engine, generate out-of-fold predictions externally and
#' pass them through `p_hat`, `mu0_hat`, and `mu1_hat`.
#'
#' `est_aipw()` is the interactive-regression-model special case of
#' [est_dml()]: `est_dml(model = "irm")` with the same nuisances and folds
#' returns the same estimate and standard error.
#'
#' @param data A data frame or `data.table`.
#' @param y Character scalar. Outcome column name.
#' @param d Character scalar. Binary treatment column name. Values must be 0/1,
#'   logical, or character/factor values coercible to 0/1.
#' @param x Character vector of pre-treatment covariate column names used when
#'   nuisance predictions are estimated internally. May be `NULL` when all of
#'   the required nuisance predictions are supplied.
#' @param estimand Character scalar: `"ATE"` (default) or `"ATT"`. The ATT
#'   score uses only the propensity score and the untreated outcome regression,
#'   so `mu1_hat` is neither required nor estimated for it. Overlap-weighted
#'   targets are a planned extension.
#' @param p_hat Optional propensity-score predictions. Either a numeric vector
#'   of length `nrow(data)` or a character scalar naming a column in `data`.
#' @param mu0_hat Optional predictions for `E[Y | D = 0, X]`. Either a numeric
#'   vector of length `nrow(data)` or a character scalar naming a column in
#'   `data`.
#' @param mu1_hat Optional predictions for `E[Y | D = 1, X]`. Either a numeric
#'   vector of length `nrow(data)` or a character scalar naming a column in
#'   `data`. Ignored when `estimand = "ATT"`.
#' @param learner_p Optional `mlr3` learner for the propensity score. If omitted
#'   and `p_hat` is not supplied, `mlr3::lrn("classif.log_reg")` is used when
#'   `mlr3` and `mlr3learners` are installed.
#' @param learner_mu0 Optional `mlr3` learner for the untreated outcome
#'   regression. If omitted and `mu0_hat` is not supplied,
#'   `mlr3::lrn("regr.lm")` is used when `mlr3` and `mlr3learners` are installed.
#' @param learner_mu1 Optional `mlr3` learner for the treated outcome
#'   regression. If omitted and `mu1_hat` is not supplied,
#'   `mlr3::lrn("regr.lm")` is used when `mlr3` and `mlr3learners` are installed.
#' @param folds Number of folds for cross-fitting when internal learners are
#'   used. Ignored when `cross_fit = FALSE`.
#' @param fold_id Optional fold identifier. Either a vector of length
#'   `nrow(data)` or a character scalar naming a column in `data`.
#' @param cross_fit Logical. If `TRUE`, internal nuisance learners are trained on
#'   all folds except the held-out fold and predict only held-out observations.
#' @param p_clip Numeric length-two vector giving lower and upper clipping bounds
#'   for propensity scores. Use `NULL` for no clipping, though this is usually
#'   not recommended.
#' @param trim Optional numeric length-two vector. If supplied, observations with
#'   raw propensity scores outside `[trim[1], trim[2]]` are dropped before the
#'   score is computed. Trimming changes the target population.
#' @param outcome_type One of `"auto"`, `"continuous"`, or `"binary"`. The
#'   value matters only when a classification learner is used for outcome
#'   nuisance estimation.
#' @param seed Optional integer seed for fold creation.
#' @param conf_level Confidence level for the Wald interval.
#' @param na_action One of `"fail"` or `"omit"`. If `"omit"`, rows with missing
#'   required variables or supplied nuisance predictions are dropped.
#'
#' @return A list with class `cm_aipw` containing `estimate`, `std.error`,
#'   confidence limits, nuisance predictions, the score, residual weights,
#'   fold ids, diagnostics, and the matched call. [tidy()] and [glance()]
#'   methods are available, so the object works with `modelsummary`.
#'
#' @details
#' For observations `i = 1, ..., n`, the ATE score is
#'
#' `mu1_hat(X_i) - mu0_hat(X_i) + D_i * (Y_i - mu1_hat(X_i)) / p_hat(X_i) -
#' (1 - D_i) * (Y_i - mu0_hat(X_i)) / (1 - p_hat(X_i))`,
#'
#' and the ATT score is
#'
#' `(D_i - (1 - D_i) * p_hat(X_i) / (1 - p_hat(X_i))) * (Y_i - mu0_hat(X_i)) /
#' mean(D)`.
#'
#' In both cases the estimate is the sample mean of the score. The reported
#' standard error is the influence-function standard error with an `n - 1`
#' denominator; for the ATE it equals `sd(score) / sqrt(n)`. Double robustness
#' is about nuisance-model misspecification: consistency can survive if either
#' the propensity score or the outcome regressions are correctly specified. It
#' is not robustness to unobserved confounding, bad adjustment sets, or lack of
#' overlap.
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.5)
#' p <- plogis(-0.2 + 0.6 * x1 - 0.4 * x2)
#' d <- rbinom(n, 1, p)
#' mu0 <- 1 + x1 + x2
#' tau <- 2
#' mu1 <- mu0 + tau
#' y <- mu0 + tau * d + rnorm(n)
#' dat <- data.frame(y = y, d = d, x1 = x1, x2 = x2,
#'                   p = p, mu0 = mu0, mu1 = mu1)
#' est_aipw(dat, y = "y", d = "d", p_hat = "p",
#'          mu0_hat = "mu0", mu1_hat = "mu1")
#'
#' # The ATT needs only the propensity score and the untreated regression.
#' est_aipw(dat, y = "y", d = "d", estimand = "ATT",
#'          p_hat = "p", mu0_hat = "mu0")
#'
#' # Internal nuisance estimation with cross-fitted mlr3 learners.
#' if (requireNamespace("mlr3", quietly = TRUE) &&
#'     requireNamespace("mlr3learners", quietly = TRUE)) {
#'   est_aipw(
#'     dat,
#'     y = "y",
#'     d = "d",
#'     x = c("x1", "x2"),
#'     learner_p = mlr3::lrn("classif.log_reg", predict_type = "prob"),
#'     learner_mu0 = mlr3::lrn("regr.lm"),
#'     learner_mu1 = mlr3::lrn("regr.lm"),
#'     folds = 5,
#'     seed = 1
#'   )
#' }
#' @seealso [est_dml()] for the partially linear model and the general
#'   double machine learning interface.
#' @export
est_aipw <- function(data,
                     y,
                     d,
                     x = NULL,
                     estimand = "ATE",
                     p_hat = NULL,
                     mu0_hat = NULL,
                     mu1_hat = NULL,
                     learner_p = NULL,
                     learner_mu0 = NULL,
                     learner_mu1 = NULL,
                     folds = 5L,
                     fold_id = NULL,
                     cross_fit = TRUE,
                     p_clip = c(0.01, 0.99),
                     trim = NULL,
                     outcome_type = c("auto", "continuous", "binary"),
                     seed = NULL,
                     conf_level = 0.95,
                     na_action = c("fail", "omit")) {
  call <- match.call()
  outcome_type <- match.arg(outcome_type)
  na_action <- match.arg(na_action)

  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("`data` must be a data frame or data.table.", call. = FALSE)
  }
  if (!.cm_is_string(y) || !.cm_is_string(d)) {
    stop("`y` and `d` must be character scalars naming columns in `data`.", call. = FALSE)
  }
  if (!.cm_is_string(estimand)) {
    stop("`estimand` must be \"ATE\" or \"ATT\".", call. = FALSE)
  }
  estimand <- toupper(estimand)
  if (!estimand %in% c("ATE", "ATT")) {
    stop("`estimand` must be \"ATE\" or \"ATT\". Overlap-weighted targets are a planned extension.", call. = FALSE)
  }
  if (!is.logical(cross_fit) || length(cross_fit) != 1L || is.na(cross_fit)) {
    stop("`cross_fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }

  dt0 <- as.data.frame(data)
  n0 <- nrow(dt0)
  if (n0 < 2L) {
    stop("`data` must contain at least two rows.", call. = FALSE)
  }
  .cm_check_column(y, dt0)
  .cm_check_column(d, dt0)

  if (is.null(x)) {
    x <- character(0L)
  }
  if (!is.character(x)) {
    stop("`x` must be a character vector of covariate column names.", call. = FALSE)
  }
  if (length(x) > 0L) {
    if (anyDuplicated(x)) {
      stop("`x` contains duplicated column names.", call. = FALSE)
    }
    if (any(x %in% c(y, d))) {
      stop("`x` must not include the outcome or treatment column.", call. = FALSE)
    }
    reserved <- c(".cm_y", ".cm_d", ".cm_row_id")
    if (any(x %in% reserved)) {
      stop("`x` contains names reserved for internal use.", call. = FALSE)
    }
    invisible(lapply(x, .cm_check_column, data = dt0))
  }

  y_vec <- dt0[[y]]
  if (is.logical(y_vec)) {
    y_vec <- as.integer(y_vec)
  }
  if (!is.numeric(y_vec)) {
    stop("`y` must be numeric or logical.", call. = FALSE)
  }
  d_vec <- .cm_as_binary(dt0[[d]], "d")

  p_supplied <- .cm_get_optional_numeric(p_hat, dt0, n0, "p_hat")
  mu0_supplied <- .cm_get_optional_numeric(mu0_hat, dt0, n0, "mu0_hat")
  mu1_supplied <- if (estimand == "ATE") .cm_get_optional_numeric(mu1_hat, dt0, n0, "mu1_hat") else NULL
  fid_supplied <- .cm_get_optional_vector(fold_id, dt0, n0, "fold_id")

  need_p <- is.null(p_supplied)
  need_mu0 <- is.null(mu0_supplied)
  need_mu1 <- estimand == "ATE" && is.null(mu1_supplied)
  need_learners <- need_p || need_mu0 || need_mu1

  if (need_learners && length(x) == 0L) {
    stop("`x` must be supplied when any nuisance prediction is estimated internally.", call. = FALSE)
  }

  keep <- is.finite(y_vec) & !is.na(d_vec)
  if (length(x) > 0L) {
    for (xj in x) {
      keep <- keep & !is.na(dt0[[xj]])
    }
  }
  if (!is.null(p_supplied)) keep <- keep & is.finite(p_supplied)
  if (!is.null(mu0_supplied)) keep <- keep & is.finite(mu0_supplied)
  if (!is.null(mu1_supplied)) keep <- keep & is.finite(mu1_supplied)
  if (!is.null(fid_supplied)) keep <- keep & !is.na(fid_supplied)

  n_missing <- sum(!keep)
  if (n_missing > 0L && identical(na_action, "fail")) {
    stop("Missing or non-finite values found in required variables or supplied nuisance predictions. Use na_action = 'omit' to drop them.", call. = FALSE)
  }
  if (n_missing > 0L && identical(na_action, "omit")) {
    warning(n_missing, " row(s) omitted because of missing or non-finite required values.", call. = FALSE)
  }

  dt <- dt0[keep, , drop = FALSE]
  y_vec <- as.numeric(y_vec[keep])
  d_vec <- d_vec[keep]
  if (!is.null(p_supplied)) p_supplied <- p_supplied[keep]
  if (!is.null(mu0_supplied)) mu0_supplied <- mu0_supplied[keep]
  if (!is.null(mu1_supplied)) mu1_supplied <- mu1_supplied[keep]
  if (!is.null(fid_supplied)) fid_supplied <- fid_supplied[keep]

  if (outcome_type == "auto") {
    y_unique <- unique(y_vec[is.finite(y_vec)])
    outcome_type <- if (length(y_unique) <= 2L && all(y_unique %in% c(0, 1))) "binary" else "continuous"
  }

  n_pre_trim <- length(y_vec)
  n_treated_pre <- sum(d_vec == 1L)
  n_control_pre <- sum(d_vec == 0L)
  if (n_treated_pre == 0L || n_control_pre == 0L) {
    stop("Both treated and control observations are required.", call. = FALSE)
  }

  work <- data.frame(
    .cm_y = y_vec,
    .cm_d = d_vec,
    .cm_row_id = seq_along(y_vec),
    check.names = FALSE
  )
  for (xj in x) {
    work[[xj]] <- dt[[xj]]
  }

  if (!is.null(fid_supplied)) {
    fold_vec <- as.integer(as.factor(fid_supplied))
    if (need_learners && cross_fit && length(unique(fold_vec)) < 2L) {
      stop("`fold_id` must contain at least two folds when cross_fit = TRUE and learners are used.", call. = FALSE)
    }
  } else if (need_learners && cross_fit) {
    fold_vec <- .cm_make_folds(d_vec, folds, seed)
  } else {
    fold_vec <- rep.int(1L, length(y_vec))
  }

  fold_summary <- data.frame(
    fold = integer(0L),
    train_n = integer(0L),
    test_n = integer(0L),
    training_excludes_test = logical(0L)
  )

  nuisance_source <- list(
    p_hat = if (need_p) NA_character_ else "supplied",
    mu0_hat = if (need_mu0) NA_character_ else "supplied",
    mu1_hat = if (estimand == "ATT") "not needed for ATT" else if (need_mu1) NA_character_ else "supplied"
  )

  if (need_learners) {
    nuis <- .cm_learn_nuisances(
      work = work,
      features = x,
      fold_id = fold_vec,
      cross_fit = cross_fit,
      need_p = need_p,
      need_mu0 = need_mu0,
      need_mu1 = need_mu1,
      learner_p = learner_p,
      learner_mu0 = learner_mu0,
      learner_mu1 = learner_mu1,
      outcome_type = outcome_type
    )
    fold_summary <- nuis$fold_summary
    if (need_p) {
      p_supplied <- nuis$p_hat
      nuisance_source$p_hat <- nuis$source$p_hat
    }
    if (need_mu0) {
      mu0_supplied <- nuis$mu0_hat
      nuisance_source$mu0_hat <- nuis$source$mu0_hat
    }
    if (need_mu1) {
      mu1_supplied <- nuis$mu1_hat
      nuisance_source$mu1_hat <- nuis$source$mu1_hat
    }
  }

  p_raw <- as.numeric(p_supplied)
  mu0 <- as.numeric(mu0_supplied)
  mu1 <- if (estimand == "ATE") as.numeric(mu1_supplied) else NULL

  .cm_check_finite(p_raw, "p_hat")
  .cm_check_finite(mu0, "mu0_hat")
  if (estimand == "ATE") .cm_check_finite(mu1, "mu1_hat")
  if (any(p_raw < 0 | p_raw > 1)) {
    stop("`p_hat` must be between 0 and 1 before clipping.", call. = FALSE)
  }

  if (!is.null(trim)) {
    trim <- .cm_check_bounds(trim, "trim", strict = FALSE)
    trim_keep <- p_raw >= trim[1L] & p_raw <= trim[2L]
  } else {
    trim_keep <- rep.int(TRUE, length(p_raw))
  }

  n_trimmed <- sum(!trim_keep)
  if (n_trimmed > 0L) {
    y_vec <- y_vec[trim_keep]
    d_vec <- d_vec[trim_keep]
    p_raw <- p_raw[trim_keep]
    mu0 <- mu0[trim_keep]
    if (!is.null(mu1)) mu1 <- mu1[trim_keep]
    fold_vec <- fold_vec[trim_keep]
  }
  if (length(y_vec) < 2L || sum(d_vec == 1L) == 0L || sum(d_vec == 0L) == 0L) {
    stop("Trimming removed too many observations; both treatment groups must remain.", call. = FALSE)
  }

  if (is.null(p_clip)) {
    p_bounds <- c(0, 1)
  } else {
    p_bounds <- .cm_check_bounds(p_clip, "p_clip", strict = TRUE)
  }
  p <- pmin(pmax(p_raw, p_bounds[1L]), p_bounds[2L])
  n_clip_low <- sum(p_raw < p_bounds[1L])
  n_clip_high <- sum(p_raw > p_bounds[2L])
  if (any(p <= 0 | p >= 1)) {
    stop("Propensity scores must be strictly between 0 and 1 after clipping. Use nonzero clipping bounds.", call. = FALSE)
  }

  psi <- .cm_score_irm(y_vec, d_vec, p, mu0, mu1, estimand)
  solution <- .cm_solve_linear_score(psi$a, psi$b, fold_vec, solve = "pooled")
  score <- psi$b
  estimate <- solution$estimate
  std_error <- solution$std.error
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  conf_low <- estimate - z * std_error
  conf_high <- estimate + z * std_error

  if (estimand == "ATE") {
    w_treated <- d_vec / p
    w_control <- (1 - d_vec) / (1 - p)
  } else {
    w_treated <- d_vec
    w_control <- (1 - d_vec) * p / (1 - p)
  }
  weights <- data.table::data.table(
    treated = w_treated,
    control = w_control,
    ipw = w_treated + w_control
  )

  diagnostics <- list(
    call = list(cross_fit = cross_fit, folds = length(unique(fold_vec))),
    sample = list(
      n_before_missing = n0,
      n_after_missing = n_pre_trim,
      n_analysis = length(y_vec),
      n_treated = sum(d_vec == 1L),
      n_control = sum(d_vec == 0L),
      omitted_missing = n_missing
    ),
    propensity = list(
      raw_summary = .cm_summary(p_raw),
      clipped_summary = .cm_summary(p),
      p_clip = p_bounds,
      n_clipped_low = n_clip_low,
      n_clipped_high = n_clip_high,
      common_support = .cm_common_support(p, d_vec)
    ),
    trimming = list(
      trim = trim,
      n_before_trim = n_pre_trim,
      n_trimmed = n_trimmed,
      n_after_trim = length(y_vec)
    ),
    weights = list(
      ipw_summary = .cm_summary(weights$ipw),
      treated_residual_weight_summary = .cm_summary(w_treated[d_vec == 1L]),
      control_residual_weight_summary = .cm_summary(w_control[d_vec == 0L]),
      ess_treated = .cm_ess(w_treated[d_vec == 1L]),
      ess_control = .cm_ess(w_control[d_vec == 0L])
    ),
    nuisance = list(
      source = nuisance_source,
      prediction_mode = if (need_learners && cross_fit) "out_of_fold" else if (need_learners) "full_sample" else "supplied",
      outcome_type = outcome_type,
      mu0_summary = .cm_summary(mu0),
      mu1_summary = if (!is.null(mu1)) .cm_summary(mu1) else NULL,
      fold_summary = fold_summary
    ),
    fold_estimates = solution$fold_estimates
  )

  out <- list(
    estimate = unname(estimate),
    std.error = unname(std_error),
    conf.low = unname(conf_low),
    conf.high = unname(conf_high),
    conf.level = conf_level,
    estimand = estimand,
    outcome = y,
    treatment = d,
    n = length(y_vec),
    n_treated = sum(d_vec == 1L),
    n_control = sum(d_vec == 0L),
    p_hat = p,
    mu0_hat = mu0,
    mu1_hat = mu1,
    score = score,
    weights = weights,
    fold_id = fold_vec,
    diagnostics = diagnostics,
    call = call
  )
  class(out) <- "cm_aipw"
  out
}

#' Print method for `cm_aipw` objects
#'
#' @param x A `cm_aipw` object returned by [est_aipw()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.cm_aipw <- function(x, ...) {
  cat("AIPW estimate (", x$estimand, ")\n", sep = "")
  cat("  Estimate:   ", formatC(x$estimate, digits = 4, format = "f"), "\n", sep = "")
  cat("  Std. Error: ", formatC(x$std.error, digits = 4, format = "f"), "\n", sep = "")
  ci_label <- if (!is.null(x$conf.level)) round(100 * x$conf.level, 1) else 95
  cat("  ", ci_label, "% CI:     [",
      formatC(x$conf.low, digits = 4, format = "f"), ", ",
      formatC(x$conf.high, digits = 4, format = "f"), "]\n", sep = "")
  cat("  N:          ", x$n, " (treated: ", x$n_treated,
      ", control: ", x$n_control, ")\n", sep = "")
  invisible(x)
}

# Cross-fitted (or full-sample) nuisance predictions for est_aipw().
.cm_learn_nuisances <- function(work,
                                features,
                                fold_id,
                                cross_fit,
                                need_p,
                                need_mu0,
                                need_mu1,
                                learner_p,
                                learner_mu0,
                                learner_mu1,
                                outcome_type) {
  .cm_require_mlr3()
  if (need_p && is.null(learner_p)) learner_p <- .cm_default_learner("classif")
  if (need_mu0 && is.null(learner_mu0)) learner_mu0 <- .cm_default_learner("regr")
  if (need_mu1 && is.null(learner_mu1)) learner_mu1 <- .cm_default_learner("regr")
  positive_y <- if (outcome_type == "binary") "1" else NULL
  mode <- if (cross_fit) "mlr3_out_of_fold" else "mlr3_full_sample"

  p_hat <- NULL
  mu0_hat <- NULL
  mu1_hat <- NULL
  fold_summary <- NULL

  if (need_p) {
    fit <- .cm_crossfit_predict(
      work, ".cm_d", features, learner_p, fold_id, cross_fit,
      positive = "1", task_hint = "propensity",
      what = "the propensity score"
    )
    p_hat <- fit$pred
    fold_summary <- fit$fold_summary
  }
  if (need_mu0) {
    fit <- .cm_crossfit_predict(
      work, ".cm_y", features, learner_mu0, fold_id, cross_fit,
      subset = work$.cm_d == 0L, positive = positive_y, task_hint = "outcome",
      what = "mu0 estimation (no controls in a training fold)"
    )
    mu0_hat <- fit$pred
    if (is.null(fold_summary)) fold_summary <- fit$fold_summary
  }
  if (need_mu1) {
    fit <- .cm_crossfit_predict(
      work, ".cm_y", features, learner_mu1, fold_id, cross_fit,
      subset = work$.cm_d == 1L, positive = positive_y, task_hint = "outcome",
      what = "mu1 estimation (no treated observations in a training fold)"
    )
    mu1_hat <- fit$pred
    if (is.null(fold_summary)) fold_summary <- fit$fold_summary
  }

  list(
    p_hat = p_hat,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    source = list(
      p_hat = if (need_p) mode else "supplied",
      mu0_hat = if (need_mu0) mode else "supplied",
      mu1_hat = if (need_mu1) mode else "supplied"
    ),
    fold_summary = fold_summary
  )
}
