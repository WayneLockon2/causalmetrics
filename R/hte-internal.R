# R/hte-internal.R
#
# Internals shared by the heterogeneous-treatment-effect and policy-learning
# functions: fitted mlr3 models that can predict on new data (with observation
# weights), cross-fitted counterfactual predictions for the S-learner, design
# matrices for best linear predictors, simultaneous bands from influence
# functions, and a small quadratic program on the simplex for stacking.
# Nothing in this file is exported.

utils::globalVariables(c(
  "estimate", "conf.low", "conf.high", "band.low", "band.high", "group",
  "gate", "mean_tau", "q", "curve", "value", "x_var", "selected", "y0_hat",
  "tau_hat", "label", "xend", "yend", "type"
))

#' @importFrom stats predict
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

# Fitted mlr3 models ----------------------------------------------------------

# Train `learner` on `train` (a data.frame holding `target`, `features`, and
# optionally weights) and return a small object that predicts on new data.
# Classification learners return positive-class probabilities.
.cm_fit_mlr3 <- function(learner, train, target, features, weights = NULL,
                         positive = NULL, task_hint = "final") {
  if (is.null(learner) || is.null(learner$task_type)) {
    stop("Learners must be valid `mlr3` learner objects.", call. = FALSE)
  }
  learner_i <- learner$clone(deep = TRUE)
  df <- train[, c(target, features), drop = FALSE]
  if (!is.null(weights)) {
    if (!("weights" %in% learner_i$properties)) {
      stop("Learner `", .cm_learner_label(learner),
           "` does not support observation weights, which this stage requires.",
           call. = FALSE)
    }
    df$.cm_w <- as.numeric(weights)
  }
  if (identical(learner_i$task_type, "classif")) {
    positive <- positive %||% "1"
    df[[target]] <- factor(as.character(df[[target]]), levels = c("0", "1"))
    if (length(unique(stats::na.omit(df[[target]]))) < 2L) {
      stop("Classification training data must contain both classes.", call. = FALSE)
    }
    learner_i$predict_type <- "prob"
    task <- mlr3::TaskClassif$new(id = paste0("cm_", task_hint), backend = df,
                                  target = target, positive = positive)
  } else if (identical(learner_i$task_type, "regr")) {
    df[[target]] <- as.numeric(df[[target]])
    task <- mlr3::TaskRegr$new(id = paste0("cm_", task_hint), backend = df, target = target)
  } else {
    stop("Learner task_type must be 'classif' or 'regr'.", call. = FALSE)
  }
  if (!is.null(weights)) task$set_col_roles(".cm_w", roles = "weights_learner")
  learner_i$train(task)
  structure(list(learner = learner_i, features = features,
                 type = learner_i$task_type, positive = positive,
                 label = .cm_learner_label(learner)),
            class = "cm_mlr3_fit")
}

.cm_predict_fit <- function(fit, newdata) {
  miss <- setdiff(fit$features, names(newdata))
  if (length(miss)) {
    stop("`newdata` lacks columns: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  df <- as.data.frame(newdata)[, fit$features, drop = FALSE]
  pred <- fit$learner$predict_newdata(df)
  if (fit$type == "classif") {
    as.numeric(pred$prob[, fit$positive])
  } else {
    as.numeric(pred$response)
  }
}

# Cross-fitted counterfactual predictions g(1, Z) and g(0, Z) from a single
# model of Y on (D, Z), for the S-learner. `d_col` is set to 1 and 0 on the
# held-out fold before predicting.
.cm_crossfit_counterfactual <- function(work, target, features, d_col, learner,
                                        fold_id, interactions = FALSE, positive = NULL) {
  n <- nrow(work)
  g1 <- g0 <- rep(NA_real_, n)
  z <- setdiff(features, d_col)
  build <- function(df, d_value = NULL) {
    if (!is.null(d_value)) df[[d_col]] <- d_value
    if (interactions) {
      for (v in z) df[[paste0(".cm_dx_", v)]] <- df[[d_col]] * df[[v]]
    }
    df
  }
  feats <- c(features, if (interactions) paste0(".cm_dx_", z))
  for (k in sort(unique(fold_id))) {
    te <- which(fold_id == k)
    tr <- which(fold_id != k)
    fit <- .cm_fit_mlr3(learner, build(work[tr, , drop = FALSE]), target, feats,
                        positive = positive, task_hint = "s_learner")
    g1[te] <- .cm_predict_fit(fit, build(work[te, , drop = FALSE], 1))
    g0[te] <- .cm_predict_fit(fit, build(work[te, , drop = FALSE], 0))
  }
  list(g1 = g1, g0 = g0)
}

# Cross-fitted regression of `target` on `features` restricted to `subset`
# rows for training, with optional weights, predicting every held-out row.
.cm_crossfit_weighted <- function(work, target, features, learner, fold_id,
                                  subset = NULL, weights = NULL, positive = NULL,
                                  task_hint = "stage2") {
  n <- nrow(work)
  if (is.null(subset)) subset <- rep(TRUE, n)
  pred <- rep(NA_real_, n)
  for (k in sort(unique(fold_id))) {
    te <- which(fold_id == k)
    tr <- which(fold_id != k & subset)
    if (length(tr) < 2L) stop("A training fold has too few observations for ", task_hint, ".", call. = FALSE)
    fit <- .cm_fit_mlr3(learner, work[tr, , drop = FALSE], target, features,
                        weights = if (is.null(weights)) NULL else weights[tr],
                        positive = positive, task_hint = task_hint)
    pred[te] <- .cm_predict_fit(fit, work[te, , drop = FALSE])
  }
  pred
}

# Scores objects --------------------------------------------------------------

.cm_check_scores <- function(scores) {
  if (!inherits(scores, "cm_scores")) {
    hint <- if (inherits(scores, "cm_cate")) {
      " A `cm_cate` model arrived here instead: R matches argument names partially, so a model passed as `s = ` or `sc = ` is taken as `scores`. Use longer model names such as `s_learner`."
    } else ""
    stop("`scores` must be a `cm_scores` object from dr_scores().", hint, call. = FALSE)
  }
  invisible(TRUE)
}

# Predictions of a CATE model (or a numeric vector) on the rows of `scores`.
.cm_cate_predictions <- function(model, scores, nm = "model") {
  n <- length(scores$score)
  if (inherits(model, "cm_cate")) {
    return(as.numeric(stats::predict(model, scores$data)))
  }
  if (is.numeric(model)) {
    if (length(model) != n) stop("`", nm, "` must have length ", n, ".", call. = FALSE)
    return(as.numeric(model))
  }
  if (.cm_is_string(model)) {
    .cm_check_column(model, scores$data)
    return(as.numeric(scores$data[[model]]))
  }
  if (is.function(model)) return(as.numeric(model(scores$data)))
  stop("`", nm, "` must be a `cm_cate` object, a numeric vector, a column name, or a function of the data.", call. = FALSE)
}

# Best linear predictors ------------------------------------------------------

# OLS of `y` on the design `X` with HC1 sandwich, influence functions, and
# an optional simultaneous critical value from the multiplier bootstrap.
.cm_ols_if <- function(X, y, conf_level = 0.95, uniform = TRUE, n_boot = 999L, seed = NULL) {
  n <- nrow(X)
  d <- ncol(X)
  Q <- crossprod(X) / n
  if (.cm_rcond(Q) < 1e-12) {
    stop("The design matrix is rank deficient; drop collinear terms.", call. = FALSE)
  }
  Qinv <- solve(Q)
  beta <- as.numeric(Qinv %*% crossprod(X, y) / n)
  e <- as.numeric(y - X %*% beta)
  inff <- (X * e) %*% Qinv
  vcov <- crossprod(inff) / n^2 * n / max(n - d, 1)
  se <- sqrt(diag(vcov))
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  crit_unif <- NA_real_
  if (uniform && d > 1L) {
    bt <- .cm_multiplier_bootstrap(inff, n, n_boot = n_boot, conf_level = conf_level, seed = seed)
    crit_unif <- bt$crit_val
  } else if (uniform) {
    crit_unif <- crit
  }
  list(beta = beta, se = se, vcov = vcov, inffunc = inff, Qinv = Qinv,
       residuals = e, crit = crit, crit_unif = crit_unif, n = n)
}

.cm_rcond <- function(M) {
  tryCatch(rcond(M), error = function(e) 0)
}

# Simultaneous critical value from an n x G influence-function matrix.
# Two-sided: quantile of max_g |Z_g / sigma_g|. One-sided: quantile of
# max_g Z_g / sigma_g. `sigma` are the analytic standard deviations of the
# columns (sqrt(mean(IF^2))).
.cm_sup_crit <- function(inffunc, n, n_boot = 999L, conf_level = 0.95,
                         one_sided = FALSE, seed = NULL) {
  bt <- .cm_multiplier_bootstrap(inffunc, n, n_boot = n_boot, conf_level = conf_level, seed = seed)
  sigma <- sqrt(colMeans(as.matrix(inffunc)^2))
  ok <- is.finite(sigma) & sigma > 0
  if (!any(ok)) return(list(crit = NA_real_, boot = bt$boot))
  z <- sweep(bt$boot[, ok, drop = FALSE], 2, sigma[ok], "/")
  stat <- if (one_sided) apply(z, 1, max) else apply(abs(z), 1, max)
  list(crit = as.numeric(stats::quantile(stat, conf_level, type = 1, na.rm = TRUE)),
       boot = bt$boot)
}

# Simplex quadratic program ---------------------------------------------------

# Minimize ||A w - b||^2 / n + lin' w over the probability simplex by
# Frank-Wolfe with exact line search. Used by the convex and Q-aggregation
# stacking rules; the problem is convex, so the iterates converge to a global
# minimizer.
.cm_simplex_qp <- function(A, b, lin = NULL, max_iter = 5000L, tol = 1e-10, w0 = NULL) {
  n <- nrow(A)
  M <- ncol(A)
  if (is.null(lin)) lin <- rep(0, M)
  w <- if (is.null(w0)) rep(1 / M, M) else w0
  obj <- function(w) sum((A %*% w - b)^2) / n + sum(lin * w)
  f <- obj(w)
  for (it in seq_len(max_iter)) {
    r <- as.numeric(A %*% w - b)
    grad <- 2 * as.numeric(crossprod(A, r)) / n + lin
    i <- which.min(grad)
    dw <- -w
    dw[i] <- dw[i] + 1
    gap <- -sum(grad * dw)
    if (gap <= tol) break
    Ad <- as.numeric(A %*% dw)
    denom <- 2 * sum(Ad^2) / n
    step <- if (denom > 0) min(1, max(0, gap / denom)) else 1
    w <- w + step * dw
    f_new <- obj(w)
    if (f - f_new < tol * max(1, abs(f))) {
      f <- f_new
      break
    }
    f <- f_new
  }
  w[w < 1e-12] <- 0
  w <- w / sum(w)
  list(weights = w, value = obj(w), iterations = it)
}

# Utilities -------------------------------------------------------------------

.cm_quantile_threshold <- function(tau, q) {
  # (1 - q) quantile of the CATE predictions, type 1 so that the threshold is
  # an attained value and tie handling is well defined.
  as.numeric(stats::quantile(tau, probs = 1 - q, type = 1, names = FALSE))
}

.cm_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
