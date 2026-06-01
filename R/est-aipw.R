#' Augmented inverse propensity weighting for binary-treatment ATEs
#'
#' `est_aipw()` estimates the average treatment effect (ATE) for a binary
#' treatment under selection on observables. The estimator combines outcome
#' regression predictions, propensity-score predictions, and inverse-propensity
#' residual corrections. It can use user-supplied nuisance predictions or fit
#' nuisance functions with `mlr3` learners, including cross-fitting.
#'
#' This function deliberately owns only the causal-estimation layer. Machine
#' learning is used only to produce nuisance predictions. To use Python or any
#' external prediction engine, generate out-of-fold predictions externally and
#' pass them through `p_hat`, `mu0_hat`, and `mu1_hat`.
#'
#' @param data A data frame or `data.table`.
#' @param y Character scalar. Outcome column name.
#' @param d Character scalar. Binary treatment column name. Values must be 0/1,
#'   logical, or character/factor values coercible to 0/1.
#' @param x Character vector of pre-treatment covariate column names used when
#'   nuisance predictions are estimated internally. May be `NULL` when all of
#'   `p_hat`, `mu0_hat`, and `mu1_hat` are supplied.
#' @param estimand Character scalar. Currently only `"ATE"` is implemented.
#'   ATT and overlap-weighted targets are planned extensions.
#' @param p_hat Optional propensity-score predictions. Either a numeric vector
#'   of length `nrow(data)` or a character scalar naming a column in `data`.
#' @param mu0_hat Optional predictions for `E[Y | D = 0, X]`. Either a numeric
#'   vector of length `nrow(data)` or a character scalar naming a column in
#'   `data`.
#' @param mu1_hat Optional predictions for `E[Y | D = 1, X]`. Either a numeric
#'   vector of length `nrow(data)` or a character scalar naming a column in
#'   `data`.
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
#'   AIPW score is computed. Trimming changes the target population.
#' @param outcome_type One of `"auto"`, `"continuous"`, or `"binary"`. The
#'   value matters only when a classification learner is used for outcome
#'   nuisance estimation.
#' @param seed Optional integer seed for fold creation.
#' @param conf_level Confidence level for the Wald interval.
#' @param na_action One of `"fail"` or `"omit"`. If `"omit"`, rows with missing
#'   required variables or supplied nuisance predictions are dropped.
#'
#' @return A list with class `cm_aipw` containing `estimate`, `std.error`,
#'   confidence limits, nuisance predictions, AIPW scores, residual weights,
#'   fold ids, diagnostics, and the matched call.
#'
#' @details
#' For observations `i = 1, ..., n`, the AIPW score is
#'
#' `mu1_hat(X_i) - mu0_hat(X_i) + D_i * (Y_i - mu1_hat(X_i)) / p_hat(X_i) -
#' (1 - D_i) * (Y_i - mu0_hat(X_i)) / (1 - p_hat(X_i))`.
#'
#' The estimate is the sample mean of this score. The reported standard error is
#' the influence-function standard error, `sd(score) / sqrt(n)`. Double
#' robustness is about nuisance-model misspecification: consistency can survive
#' if either the propensity score or the outcome regressions are correctly
#' specified. It is not robustness to unobserved confounding, bad adjustment
#' sets, or lack of overlap.
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
#' \dontrun{
#' library(mlr3)
#' library(mlr3learners)
#' est_aipw(
#'   dat,
#'   y = "y",
#'   d = "d",
#'   x = c("x1", "x2"),
#'   learner_p = lrn("classif.log_reg", predict_type = "prob"),
#'   learner_mu0 = lrn("regr.lm"),
#'   learner_mu1 = lrn("regr.lm"),
#'   folds = 5,
#'   seed = 1
#' )
#' }
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
  estimand <- toupper(estimand)
  if (!identical(estimand, "ATE")) {
    stop("Only estimand = 'ATE' is implemented. ATT and overlap-weighted targets are planned extensions.", call. = FALSE)
  }
  if (!is.logical(cross_fit) || length(cross_fit) != 1L || is.na(cross_fit)) {
    stop("`cross_fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }

  dt0 <- data.table::copy(data.table::as.data.table(data))
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
  mu1_supplied <- .cm_get_optional_numeric(mu1_hat, dt0, n0, "mu1_hat")
  fid_supplied <- .cm_get_optional_vector(fold_id, dt0, n0, "fold_id")

  need_p <- is.null(p_supplied)
  need_mu0 <- is.null(mu0_supplied)
  need_mu1 <- is.null(mu1_supplied)
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

  dt <- dt0[keep]
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

  work <- data.table::data.table(.cm_y = y_vec, .cm_d = d_vec, .cm_row_id = seq_along(y_vec))
  if (length(x) > 0L) {
    for (xj in x) {
      work[, (xj) := dt[[xj]]]
    }
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
    mu1_hat = if (need_mu1) NA_character_ else "supplied"
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
  mu1 <- as.numeric(mu1_supplied)

  .cm_check_finite(p_raw, "p_hat")
  .cm_check_finite(mu0, "mu0_hat")
  .cm_check_finite(mu1, "mu1_hat")
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
    mu1 <- mu1[trim_keep]
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

  score <- (mu1 - mu0) +
    d_vec * (y_vec - mu1) / p -
    (1 - d_vec) * (y_vec - mu0) / (1 - p)
  estimate <- mean(score)
  std_error <- stats::sd(score) / sqrt(length(score))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  conf_low <- estimate - z * std_error
  conf_high <- estimate + z * std_error

  w_treated <- d_vec / p
  w_control <- (1 - d_vec) / (1 - p)
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
      mu1_summary = .cm_summary(mu1),
      fold_summary = fold_summary
    )
  )

  out <- list(
    estimate = unname(estimate),
    std.error = unname(std_error),
    conf.low = unname(conf_low),
    conf.high = unname(conf_high),
    conf.level = conf_level,
    estimand = estimand,
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

.cm_is_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.cm_check_column <- function(column, data) {
  if (!column %in% names(data)) {
    stop("Column `", column, "` not found in `data`.", call. = FALSE)
  }
  invisible(TRUE)
}

.cm_as_binary <- function(x, nm) {
  if (is.logical(x)) {
    out <- as.integer(x)
  } else if (is.factor(x)) {
    out <- suppressWarnings(as.integer(as.character(x)))
  } else if (is.character(x)) {
    out <- suppressWarnings(as.integer(x))
  } else if (is.numeric(x) || is.integer(x)) {
    if (any(!is.na(x) & !(x %in% c(0, 1)))) {
      stop("`", nm, "` must be binary with values 0 and 1.", call. = FALSE)
    }
    out <- as.integer(x)
  } else {
    stop("`", nm, "` must be binary with values 0 and 1.", call. = FALSE)
  }
  if (any(is.na(out) & !is.na(x)) || !all(out[!is.na(out)] %in% c(0L, 1L))) {
    stop("`", nm, "` must be binary with values 0 and 1.", call. = FALSE)
  }
  out
}

.cm_get_optional_vector <- function(arg, data, n, nm) {
  if (is.null(arg)) return(NULL)
  if (.cm_is_string(arg)) {
    .cm_check_column(arg, data)
    return(data[[arg]])
  }
  if (length(arg) != n) {
    stop("`", nm, "` must be length nrow(data) or a column name.", call. = FALSE)
  }
  arg
}

.cm_get_optional_numeric <- function(arg, data, n, nm) {
  vec <- .cm_get_optional_vector(arg, data, n, nm)
  if (is.null(vec)) return(NULL)
  if (!is.numeric(vec) && !is.integer(vec)) {
    stop("`", nm, "` must be numeric or a numeric column name.", call. = FALSE)
  }
  as.numeric(vec)
}

.cm_check_finite <- function(x, nm) {
  if (length(x) == 0L || any(!is.finite(x))) {
    stop("`", nm, "` contains missing or non-finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

.cm_check_bounds <- function(x, nm, strict) {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
    stop("`", nm, "` must be a numeric vector of length two.", call. = FALSE)
  }
  if (x[1L] >= x[2L] || x[1L] < 0 || x[2L] > 1) {
    stop("`", nm, "` must satisfy 0 <= lower < upper <= 1.", call. = FALSE)
  }
  if (strict && (x[1L] <= 0 || x[2L] >= 1)) {
    stop("`", nm, "` must use bounds strictly inside (0, 1).", call. = FALSE)
  }
  x
}

.cm_make_folds <- function(d, folds, seed) {
  if (!is.numeric(folds) || length(folds) != 1L || !is.finite(folds)) {
    stop("`folds` must be a positive integer.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (folds < 2L) {
    stop("`folds` must be at least 2 when cross_fit = TRUE.", call. = FALSE)
  }
  if (folds > length(d)) {
    stop("`folds` cannot exceed the number of complete observations.", call. = FALSE)
  }
  if (sum(d == 1L) < 2L || sum(d == 0L) < 2L) {
    stop("At least two treated and two control observations are required for cross-fitting.", call. = FALSE)
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  out <- integer(length(d))
  for (val in c(0L, 1L)) {
    ids <- which(d == val)
    out[ids] <- sample(rep(seq_len(folds), length.out = length(ids)))
  }
  out
}

.cm_default_learner <- function(which) {
  if (!requireNamespace("mlr3", quietly = TRUE) ||
      !requireNamespace("mlr3learners", quietly = TRUE)) {
    stop("Internal nuisance estimation requires packages `mlr3` and `mlr3learners`, or supply p_hat, mu0_hat, and mu1_hat directly.", call. = FALSE)
  }
  if (identical(which, "p")) {
    mlr3::lrn("classif.log_reg", predict_type = "prob")
  } else {
    mlr3::lrn("regr.lm")
  }
}

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
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    stop("Internal nuisance estimation requires package `mlr3`.", call. = FALSE)
  }
  if (need_p && is.null(learner_p)) learner_p <- .cm_default_learner("p")
  if (need_mu0 && is.null(learner_mu0)) learner_mu0 <- .cm_default_learner("mu")
  if (need_mu1 && is.null(learner_mu1)) learner_mu1 <- .cm_default_learner("mu")

  n <- nrow(work)
  p_hat <- rep(NA_real_, n)
  mu0_hat <- rep(NA_real_, n)
  mu1_hat <- rep(NA_real_, n)
  rows <- list()

  if (cross_fit) {
    folds <- sort(unique(fold_id))
    for (k in folds) {
      test_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k)
      train <- work[train_idx]
      test <- work[test_idx]
      if (need_p) {
        .cm_check_two_groups(train$.cm_d, "propensity learner training fold")
        p_hat[test_idx] <- .cm_predict_mlr3(
          learner = learner_p,
          train = train,
          test = test,
          target = ".cm_d",
          features = features,
          positive = "1",
          task_hint = "propensity"
        )
      }
      if (need_mu0) {
        train0 <- train[.cm_d == 0L]
        if (nrow(train0) == 0L) stop("A training fold has no controls for mu0 estimation.", call. = FALSE)
        mu0_hat[test_idx] <- .cm_predict_mlr3(
          learner = learner_mu0,
          train = train0,
          test = test,
          target = ".cm_y",
          features = features,
          positive = if (outcome_type == "binary") "1" else NULL,
          task_hint = "outcome"
        )
      }
      if (need_mu1) {
        train1 <- train[.cm_d == 1L]
        if (nrow(train1) == 0L) stop("A training fold has no treated observations for mu1 estimation.", call. = FALSE)
        mu1_hat[test_idx] <- .cm_predict_mlr3(
          learner = learner_mu1,
          train = train1,
          test = test,
          target = ".cm_y",
          features = features,
          positive = if (outcome_type == "binary") "1" else NULL,
          task_hint = "outcome"
        )
      }
      rows[[length(rows) + 1L]] <- data.frame(
        fold = as.integer(k),
        train_n = length(train_idx),
        test_n = length(test_idx),
        training_excludes_test = TRUE
      )
    }
  } else {
    train <- work
    test <- work
    if (need_p) {
      .cm_check_two_groups(train$.cm_d, "propensity learner training sample")
      p_hat <- .cm_predict_mlr3(learner_p, train, test, ".cm_d", features, "1", "propensity")
    }
    if (need_mu0) {
      train0 <- train[.cm_d == 0L]
      mu0_hat <- .cm_predict_mlr3(learner_mu0, train0, test, ".cm_y", features,
                                  if (outcome_type == "binary") "1" else NULL, "outcome")
    }
    if (need_mu1) {
      train1 <- train[.cm_d == 1L]
      mu1_hat <- .cm_predict_mlr3(learner_mu1, train1, test, ".cm_y", features,
                                  if (outcome_type == "binary") "1" else NULL, "outcome")
    }
    rows[[1L]] <- data.frame(
      fold = 1L,
      train_n = n,
      test_n = n,
      training_excludes_test = FALSE
    )
  }

  list(
    p_hat = p_hat,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    source = list(
      p_hat = if (need_p) if (cross_fit) "mlr3_out_of_fold" else "mlr3_full_sample" else "supplied",
      mu0_hat = if (need_mu0) if (cross_fit) "mlr3_out_of_fold" else "mlr3_full_sample" else "supplied",
      mu1_hat = if (need_mu1) if (cross_fit) "mlr3_out_of_fold" else "mlr3_full_sample" else "supplied"
    ),
    fold_summary = do.call(rbind, rows)
  )
}

.cm_check_two_groups <- function(d, where) {
  if (sum(d == 1L) == 0L || sum(d == 0L) == 0L) {
    stop("The ", where, " must contain both treatment groups.", call. = FALSE)
  }
  invisible(TRUE)
}

.cm_predict_mlr3 <- function(learner,
                             train,
                             test,
                             target,
                             features,
                             positive = NULL,
                             task_hint = c("propensity", "outcome")) {
  task_hint <- match.arg(task_hint)
  if (is.null(learner) || is.null(learner$task_type)) {
    stop("Learners must be valid `mlr3` learner objects.", call. = FALSE)
  }
  learner_i <- learner$clone(deep = TRUE)
  cols <- c(target, features)
  train_df <- as.data.frame(train[, cols, with = FALSE])
  test_df <- as.data.frame(test[, features, with = FALSE])

  if (identical(learner_i$task_type, "classif")) {
    if (is.null(positive)) {
      stop("Classification outcome learners require binary `outcome_type`.", call. = FALSE)
    }
    train_df[[target]] <- factor(as.character(train_df[[target]]), levels = c("0", "1"))
    if (length(unique(stats::na.omit(train_df[[target]]))) < 2L) {
      stop("Classification learner training data must contain both outcome classes.", call. = FALSE)
    }
    if (!("prob" %in% learner_i$predict_types)) {
      stop("Classification learners must support predict_type = 'prob'.", call. = FALSE)
    }
    learner_i$predict_type <- "prob"
    task <- mlr3::TaskClassif$new(
      id = paste0("cm_", task_hint),
      backend = train_df,
      target = target,
      positive = positive
    )
    learner_i$train(task)
    pred <- learner_i$predict_newdata(test_df)
    probs <- pred$prob
    if (is.null(probs) || !(positive %in% colnames(probs))) {
      stop("Could not extract positive-class probabilities from learner predictions.", call. = FALSE)
    }
    return(as.numeric(probs[, positive]))
  }

  if (identical(learner_i$task_type, "regr")) {
    train_df[[target]] <- as.numeric(train_df[[target]])
    task <- mlr3::TaskRegr$new(
      id = paste0("cm_", task_hint),
      backend = train_df,
      target = target
    )
    learner_i$train(task)
    pred <- learner_i$predict_newdata(test_df)
    return(as.numeric(pred$response))
  }

  stop("Learner task_type must be 'classif' or 'regr'.", call. = FALSE)
}

.cm_summary <- function(x) {
  x <- as.numeric(x)
  c(
    min = min(x),
    q01 = stats::quantile(x, 0.01, names = FALSE, type = 7),
    q05 = stats::quantile(x, 0.05, names = FALSE, type = 7),
    median = stats::median(x),
    mean = mean(x),
    q95 = stats::quantile(x, 0.95, names = FALSE, type = 7),
    q99 = stats::quantile(x, 0.99, names = FALSE, type = 7),
    max = max(x)
  )
}

.cm_ess <- function(w) {
  w <- as.numeric(w)
  if (length(w) == 0L || sum(w^2) == 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

.cm_common_support <- function(p, d) {
  p1 <- p[d == 1L]
  p0 <- p[d == 0L]
  low <- max(min(p1), min(p0))
  high <- min(max(p1), max(p0))
  outside <- p < low | p > high
  list(
    low = low,
    high = high,
    share_outside = mean(outside),
    treated_share_outside = mean(outside[d == 1L]),
    control_share_outside = mean(outside[d == 0L])
  )
}
