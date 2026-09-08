# R/mediation-internal.R
#
# Internals shared by the mediation functions: robust and clustered covariance
# matrices for lm/glm fits, a generic cluster bootstrap, coefficient draws for
# quasi-Bayesian inference, and summaries. Nothing here is exported.

utils::globalVariables(c("effect", "contribution", "covariate", "rho", "nie", ".cm_w"))

# HC1 (or cluster-robust) covariance of an lm or glm fit.
.cm_vcov_robust <- function(fit, cluster = NULL) {
  X <- stats::model.matrix(fit)
  if (inherits(fit, "glm")) {
    pw <- if (is.null(fit$prior.weights)) rep(1, nrow(X)) else fit$prior.weights
    bread <- solve(crossprod(X * sqrt(fit$weights)))
    score <- X * (pw * (fit$y - fit$fitted.values))
  } else {
    pw <- if (is.null(fit$weights)) rep(1, nrow(X)) else fit$weights
    bread <- solve(crossprod(X * sqrt(pw)))
    score <- X * (pw * stats::residuals(fit))
  }
  n <- nrow(X)
  k <- ncol(X)
  if (is.null(cluster)) {
    meat <- crossprod(score)
    adj <- n / (n - k)
  } else {
    cl <- as.integer(factor(cluster))
    G <- max(cl)
    meat <- crossprod(rowsum(score, cl))
    adj <- G / (G - 1) * (n - 1) / (n - k)
  }
  V <- adj * bread %*% meat %*% bread
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# Draw coefficient vectors from N(coef, vcov); returns an n_sim x k matrix.
.cm_draw_coefs <- function(coef, vcov, n_sim) {
  k <- length(coef)
  Z <- matrix(stats::rnorm(n_sim * k), n_sim, k)
  R <- tryCatch(chol(vcov), error = function(e) NULL)
  draws <- if (!is.null(R)) {
    Z %*% R
  } else {
    ev <- eigen((vcov + t(vcov)) / 2, symmetric = TRUE)
    Z %*% t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), k))
  }
  colnames(draws) <- names(coef)
  sweep(draws, 2, coef, "+")
}

# Generic cluster bootstrap: `stat(data)` returns a named numeric vector.
.cm_boot <- function(data, stat, n_boot, cluster = NULL, seed = NULL, conf_level = 0.95) {
  est <- stat(data)
  n <- nrow(data)
  if (n_boot < 1) {
    return(data.frame(term = names(est), estimate = as.numeric(est), std.error = NA_real_,
                      conf.low = NA_real_, conf.high = NA_real_, stringsAsFactors = FALSE, row.names = NULL))
  }
  # A resample in which a regression cannot be fitted (a covariate constant
  # in the drawn clusters, say) is dropped from the replications.
  safe_stat <- function(dd) tryCatch(stat(dd), error = function(e) rep(NA_real_, length(est)))
  draws <- .cm_with_seed(seed, {
    if (is.null(cluster)) {
      replicate(n_boot, safe_stat(data[sample.int(n, n, replace = TRUE), , drop = FALSE]))
    } else {
      cl <- factor(cluster)
      ids <- split(seq_len(n), cl)
      G <- length(ids)
      replicate(n_boot, {
        pick <- sample.int(G, G, replace = TRUE)
        safe_stat(data[unlist(ids[pick], use.names = FALSE), , drop = FALSE])
      })
    }
  })
  draws <- matrix(draws, nrow = length(est), dimnames = list(names(est), NULL))
  n_failed <- sum(colSums(is.na(draws)) == nrow(draws))
  if (n_failed > 0) warning(n_failed, " of ", n_boot, " bootstrap replications failed and were dropped.", call. = FALSE)
  alpha <- (1 - conf_level) / 2
  data.frame(term = names(est), estimate = as.numeric(est),
             std.error = apply(draws, 1, stats::sd, na.rm = TRUE),
             conf.low = apply(draws, 1, stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE),
             conf.high = apply(draws, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE),
             stringsAsFactors = FALSE, row.names = NULL)
}

# Summary table from a matrix of simulated effects (n_sim x K).
# Per-observation influence contributions of the coefficients of an lm/glm
# fit (bread %*% score), so that several regressions can be stacked and
# their joint covariance computed as in Stata's `suest`. `rows` maps the
# fit's observations to rows of the full sample (for models fitted on a
# subset); other rows contribute zero.
.cm_infl <- function(fit, n_full, rows = seq_len(n_full)) {
  X <- stats::model.matrix(fit)
  if (inherits(fit, "glm")) {
    pw <- if (is.null(fit$prior.weights)) rep(1, nrow(X)) else fit$prior.weights
    bread <- solve(crossprod(X * sqrt(fit$weights)))
    score <- X * (pw * (fit$y - fit$fitted.values))
  } else {
    pw <- if (is.null(fit$weights)) rep(1, nrow(X)) else fit$weights
    bread <- solve(crossprod(X * sqrt(pw)))
    score <- X * (pw * stats::residuals(fit))
  }
  out <- matrix(0, n_full, ncol(X), dimnames = list(NULL, colnames(X)))
  out[rows, ] <- score %*% bread
  out
}

# Joint robust (or cluster-robust) covariance of stacked influence matrices.
.cm_stacked_vcov <- function(infl, cluster = NULL) {
  n <- nrow(infl); k <- ncol(infl)
  if (is.null(cluster)) {
    V <- crossprod(infl) * n / (n - 1)
  } else {
    cl <- as.integer(factor(cluster))
    G <- max(cl)
    V <- crossprod(rowsum(infl, cl)) * G / (G - 1)
  }
  V
}

.cm_sim_table <- function(est, draws, conf_level = 0.95) {
  alpha <- (1 - conf_level) / 2
  data.frame(term = names(est), estimate = as.numeric(est),
             std.error = apply(draws, 2, stats::sd, na.rm = TRUE),
             conf.low = apply(draws, 2, stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE),
             conf.high = apply(draws, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE),
             stringsAsFactors = FALSE, row.names = NULL)
}

.cm_check_binary_d <- function(data, d) {
  .cm_check_column(d, data)
  .cm_as_binary(data[[d]], d)
}

.cm_med_formula <- function(lhs, rhs) {
  stats::as.formula(paste(lhs, "~", if (length(rhs)) paste(rhs, collapse = " + ") else "1"))
}
