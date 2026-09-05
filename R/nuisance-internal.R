# R/nuisance-internal.R
#
# Internals shared by est_aipw() and est_dml(): argument checks, random folds,
# mlr3 prediction, cross-fitting, Neyman-orthogonal scores, and summaries.
# Nothing in this file is exported.

# Argument checks -------------------------------------------------------------

.cm_is_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.cm_check_column <- function(column, data) {
  if (!column %in% names(data)) {
    stop("Column `", column, "` not found in `data`.", call. = FALSE)
  }
  invisible(TRUE)
}

.cm_check_count <- function(x, nm, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min || x != round(x)) {
    stop("`", nm, "` must be a single integer of at least ", min, ".", call. = FALSE)
  }
  as.integer(x)
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

# Random folds ----------------------------------------------------------------

# Evaluate `code` with the RNG seeded by `seed`, then restore the caller's
# random stream. With `seed = NULL` the code simply consumes the current stream.
.cm_with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
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
  code
}

.cm_validate_folds <- function(folds, n, strata = NULL) {
  if (!is.numeric(folds) || length(folds) != 1L || !is.finite(folds)) {
    stop("`folds` must be a positive integer.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (folds < 2L) {
    stop("`folds` must be at least 2 when cross_fit = TRUE.", call. = FALSE)
  }
  if (folds > n) {
    stop("`folds` cannot exceed the number of complete observations.", call. = FALSE)
  }
  if (!is.null(strata) && any(table(strata) < 2L)) {
    stop("At least two treated and two control observations are required for cross-fitting.", call. = FALSE)
  }
  folds
}

# One random partition into `folds` groups, balanced within each stratum.
.cm_draw_folds <- function(n, folds, strata = NULL) {
  out <- integer(n)
  if (is.null(strata)) {
    out[] <- sample(rep(seq_len(folds), length.out = n))
  } else {
    for (val in sort(unique(strata))) {
      ids <- which(strata == val)
      out[ids] <- sample(rep(seq_len(folds), length.out = length(ids)))
    }
  }
  out
}

# Folds stratified by a binary treatment (used by est_aipw()).
.cm_make_folds <- function(d, folds, seed) {
  folds <- .cm_validate_folds(folds, length(d), strata = d)
  .cm_with_seed(seed, .cm_draw_folds(length(d), folds, strata = d))
}

# `n_rep` independent partitions drawn from one seeded stream, so repetition r
# is reproducible given `seed`.
.cm_make_fold_sets <- function(n, folds, n_rep, seed, strata = NULL) {
  folds <- .cm_validate_folds(folds, n, strata)
  .cm_with_seed(seed, lapply(seq_len(n_rep), function(r) .cm_draw_folds(n, folds, strata)))
}

# mlr3 learners ---------------------------------------------------------------

.cm_require_mlr3 <- function() {
  if (!requireNamespace("mlr3", quietly = TRUE) ||
      !requireNamespace("mlr3learners", quietly = TRUE)) {
    stop("Internal nuisance estimation requires packages `mlr3` and `mlr3learners`, or supply the nuisance predictions directly.", call. = FALSE)
  }
  invisible(TRUE)
}

.cm_default_learner <- function(type = c("regr", "classif")) {
  type <- match.arg(type)
  .cm_require_mlr3()
  if (type == "classif") {
    mlr3::lrn("classif.log_reg", predict_type = "prob")
  } else {
    mlr3::lrn("regr.lm")
  }
}

.cm_learner_label <- function(learner) {
  if (is.null(learner)) return(NA_character_)
  if (!is.null(learner$id)) return(as.character(learner$id))
  class(learner)[1L]
}

# Train `learner` on `train`, predict on `test`. Classification learners must
# return positive-class probabilities; regression learners return responses.
.cm_predict_mlr3 <- function(learner,
                             train,
                             test,
                             target,
                             features,
                             positive = NULL,
                             task_hint = "nuisance") {
  if (is.null(learner) || is.null(learner$task_type)) {
    stop("Learners must be valid `mlr3` learner objects.", call. = FALSE)
  }
  learner_i <- learner$clone(deep = TRUE)
  cols <- c(target, features)
  train_df <- train[, cols, drop = FALSE]
  test_df <- test[, features, drop = FALSE]

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

# Out-of-fold predictions of `target` from `features`. Training rows can be
# restricted with `subset` (e.g. controls only); predictions are always made
# for every row of the held-out fold, or for every row when `cross_fit = FALSE`.
.cm_crossfit_predict <- function(work,
                                 target,
                                 features,
                                 learner,
                                 fold_id,
                                 cross_fit,
                                 subset = NULL,
                                 positive = NULL,
                                 task_hint = "nuisance",
                                 what = task_hint) {
  n <- nrow(work)
  if (is.null(subset)) subset <- rep(TRUE, n)
  pred <- rep(NA_real_, n)
  rows <- list()

  if (cross_fit) {
    for (k in sort(unique(fold_id))) {
      test_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k & subset)
      if (length(train_idx) == 0L) {
        stop("A training fold has no observations for ", what, ".", call. = FALSE)
      }
      pred[test_idx] <- .cm_predict_mlr3(
        learner = learner,
        train = work[train_idx, , drop = FALSE],
        test = work[test_idx, , drop = FALSE],
        target = target,
        features = features,
        positive = positive,
        task_hint = task_hint
      )
      rows[[length(rows) + 1L]] <- data.frame(
        fold = as.integer(k),
        train_n = length(train_idx),
        test_n = length(test_idx),
        training_excludes_test = TRUE
      )
    }
  } else {
    train_idx <- which(subset)
    if (length(train_idx) == 0L) {
      stop("There are no observations to train ", what, ".", call. = FALSE)
    }
    pred <- .cm_predict_mlr3(
      learner = learner,
      train = work[train_idx, , drop = FALSE],
      test = work,
      target = target,
      features = features,
      positive = positive,
      task_hint = task_hint
    )
    rows[[1L]] <- data.frame(
      fold = 1L,
      train_n = length(train_idx),
      test_n = n,
      training_excludes_test = FALSE
    )
  }

  list(pred = pred, fold_summary = do.call(rbind, rows))
}

# Neyman-orthogonal scores ----------------------------------------------------
#
# Every score is linear in the target parameter: psi = psi_b - psi_a * theta,
# so theta_hat = sum(psi_b) / sum(psi_a) and the influence function is
# psi / E[psi_a].

# Partially linear model, partialling-out score (Robinson 1988):
# psi = (Y - l(X) - theta (D - m(X))) (D - m(X)).
.cm_score_plr <- function(y_tilde, d_tilde) {
  list(a = d_tilde^2, b = y_tilde * d_tilde)
}

# Interactive regression model. ATE: the AIPW score. ATT: the score of
# Chernozhukov et al. (2018), which needs only mu0 and the propensity score.
.cm_score_irm <- function(y, d, p, mu0, mu1, estimand = c("ATE", "ATT")) {
  estimand <- match.arg(estimand)
  if (estimand == "ATE") {
    b <- (mu1 - mu0) + d * (y - mu1) / p - (1 - d) * (y - mu0) / (1 - p)
    return(list(a = rep(1, length(y)), b = b))
  }
  p_bar <- mean(d)
  b <- (d - (1 - d) * p / (1 - p)) * (y - mu0) / p_bar
  list(a = d / p_bar, b = b)
}

# Solve the linear score for theta. "pooled" stacks the scores across folds and
# solves once (DML2 in Chernozhukov et al. 2018); "fold_average" solves fold by
# fold and averages (DML1). Both use the influence-function variance evaluated
# at the final estimate, with an n - 1 denominator (HC1-type).
.cm_solve_linear_score <- function(psi_a, psi_b, fold_id, solve = c("pooled", "fold_average")) {
  solve <- match.arg(solve)
  n <- length(psi_a)
  jac <- mean(psi_a)
  if (!is.finite(jac) || jac == 0) {
    stop("The score Jacobian is zero: the target parameter is not identified (no residual variation in the treatment).", call. = FALSE)
  }
  fold_levels <- sort(unique(fold_id))
  fold_est <- vapply(fold_levels, function(k) {
    i <- fold_id == k
    sum(psi_b[i]) / sum(psi_a[i])
  }, numeric(1))
  fold_n <- vapply(fold_levels, function(k) sum(fold_id == k), integer(1))

  estimate <- if (solve == "fold_average") mean(fold_est) else sum(psi_b) / sum(psi_a)
  psi_hat <- psi_b - psi_a * estimate
  variance <- sum(psi_hat^2) / (n - 1) / jac^2

  list(
    estimate = unname(estimate),
    std.error = sqrt(variance / n),
    fold_estimates = data.frame(fold = as.integer(fold_levels), estimate = unname(fold_est), n = fold_n),
    psi = psi_hat,
    jacobian = unname(jac)
  )
}

# Median aggregation over repeated cross-fitting splits
# (Chernozhukov et al. 2018, Section 3.4).
.cm_aggregate_reps <- function(estimates, std_errors) {
  est <- stats::median(estimates)
  se <- sqrt(stats::median(std_errors^2 + (estimates - est)^2))
  list(estimate = est, std.error = se)
}

# Summaries -------------------------------------------------------------------

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

# Root mean squared error and an R^2-style fit measure for a nuisance
# prediction of `target`.
.cm_fit_quality <- function(target, prediction) {
  resid <- target - prediction
  mse <- mean(resid^2)
  v <- stats::var(target)
  c(rmse = sqrt(mse), r2 = if (is.finite(v) && v > 0) 1 - mse / v else NA_real_)
}
