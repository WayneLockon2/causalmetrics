# R/iv-internal.R
#
# Internals for the instrumental-variables tools: two-stage least squares on
# arbitrary matrices with robust or clustered sandwich variances,
# residualization on controls, and the Anderson-Rubin statistic. Nothing in
# this file is exported.

utils::globalVariables(c("theta", "statistic", "u", "mte", "weight", "alpha", "beta_k",
                         "r2_d", "r2_y", "bias", "term", "component"))

# Residualize the columns of `M` on `X` (with intercept) by weighted least
# squares. `X = NULL` demeans.
.cm_residualize <- function(M, X = NULL, weights = NULL) {
  M <- as.matrix(M)
  n <- nrow(M)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  Xd <- if (is.null(X)) matrix(1, n, 1) else cbind(1, as.matrix(X))
  fit <- stats::lm.wfit(Xd, M, w)
  R <- as.matrix(fit$residuals)
  colnames(R) <- colnames(M)
  R
}

# Two-stage least squares of y on d (one endogenous column) plus exogenous
# columns X, with instruments Z. Robust (HC1) or cluster-robust variance.
# Everything is done on residualized variables when X is given, which is
# numerically identical to including X in both stages.
.cm_tsls <- function(y, d, Z, X = NULL, weights = NULL, cluster = NULL) {
  n <- length(y)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  yt <- .cm_residualize(y, X, w)[, 1]
  dt <- .cm_residualize(d, X, w)[, 1]
  Zt <- .cm_residualize(Z, X, w)
  # first stage
  pi_hat <- stats::lm.wfit(Zt, dt, w)$coefficients
  pi_hat[is.na(pi_hat)] <- 0
  d_hat <- as.numeric(Zt %*% pi_hat)
  jac <- sum(w * d_hat * dt) / n
  if (!is.finite(jac) || abs(jac) < 1e-12) stop("The instrument has no explanatory power for the treatment.", call. = FALSE)
  theta <- sum(w * d_hat * yt) / sum(w * d_hat * dt)
  e <- yt - theta * dt
  psi <- w * d_hat * e
  k_adj <- ncol(Zt) + if (is.null(X)) 1L else ncol(as.matrix(X)) + 1L
  if (is.null(cluster)) {
    v <- sum(psi^2) / (n - k_adj) / jac^2 / n
    n_cl <- n
  } else {
    cs <- rowsum(psi, as.integer(factor(cluster)))
    n_cl <- nrow(cs)
    v <- sum(cs^2) * n_cl / (n_cl - 1) / jac^2 / n^2
  }
  list(estimate = unname(theta), std.error = sqrt(v), first_stage = pi_hat, d_hat = d_hat,
       residuals = list(y = yt, d = dt, z = Zt, e = e), n = n, n_clusters = n_cl,
       jacobian = jac)
}

# First-stage strength on residualized variables: homoskedastic F, robust
# (HC1 or clustered) Wald F, and the effective F of Montiel Olea and
# Pflueger (2013) for one endogenous regressor.
.cm_first_stage_strength <- function(dt, Zt, weights = NULL, cluster = NULL) {
  n <- length(dt)
  k <- ncol(Zt)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  fit <- stats::lm.wfit(Zt, dt, w)
  pi_hat <- fit$coefficients
  pi_hat[is.na(pi_hat)] <- 0
  u <- fit$residuals
  ZtZ <- crossprod(Zt * sqrt(w))
  ZtZ_inv <- tryCatch(solve(ZtZ), error = function(e) MASS_ginv(ZtZ))
  # homoskedastic F
  sigma2 <- sum(w * u^2) / (n - k)
  F_hom <- as.numeric(t(pi_hat) %*% ZtZ %*% pi_hat) / (k * sigma2)
  # robust variance of pi_hat
  scores <- Zt * (w * u)
  if (is.null(cluster)) {
    meat <- crossprod(scores) * n / (n - k)
  } else {
    cs <- rowsum(scores, as.integer(factor(cluster)))
    G <- nrow(cs)
    meat <- crossprod(cs) * G / (G - 1)
  }
  V <- ZtZ_inv %*% meat %*% ZtZ_inv
  F_rob <- as.numeric(t(pi_hat) %*% solve(V, pi_hat)) / k
  Q <- ZtZ / n
  F_eff <- as.numeric(t(pi_hat) %*% Q %*% pi_hat) / sum(diag(V %*% Q))
  list(F = F_hom, F_robust = F_rob, F_effective = F_eff, coefficients = pi_hat,
       vcov = V, t = if (k == 1L) unname(pi_hat[1] / sqrt(V[1, 1])) else NA_real_, k = k)
}

MASS_ginv <- function(M) {
  s <- svd(M)
  tol <- max(dim(M)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# Anderson-Rubin / C(alpha) statistic on a grid, from moment contributions
# m_i(theta) = (y_i - theta d_i) z_i (z_i a k-vector), with robust or
# clustered variance. Returns the statistic at each grid value and, for k = 1,
# the analytic confidence set.
.cm_ar_statistic <- function(yt, dt, Zt, theta_grid, cluster = NULL, weights = NULL) {
  n <- length(yt)
  k <- ncol(Zt)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  my <- Zt * (w * yt)
  md <- Zt * (w * dt)
  if (!is.null(cluster)) {
    cl <- as.integer(factor(cluster))
    my_c <- rowsum(my, cl); md_c <- rowsum(md, cl)
    G <- nrow(my_c)
  } else {
    my_c <- my; md_c <- md; G <- n
  }
  My <- colSums(my) / n; Md <- colSums(md) / n
  # centered cluster sums scaled so that Omega = (1/n) sum_c (sum_i m_i)(sum_i m_i)' approximately
  Cy <- sweep(my_c, 2, colSums(my_c) / G); Cd <- sweep(md_c, 2, colSums(md_c) / G)
  Syy <- crossprod(Cy) / n * G / (G - 1); Sdd <- crossprod(Cd) / n * G / (G - 1)
  Syd <- crossprod(Cy, Cd) / n * G / (G - 1)
  stat <- vapply(theta_grid, function(th) {
    M <- My - th * Md
    Om <- Syy - th * (Syd + t(Syd)) + th^2 * Sdd
    as.numeric(n * t(M) %*% solve(Om, M))
  }, numeric(1))
  list(statistic = stat, k = k, My = My, Md = Md, Syy = Syy, Sdd = Sdd, Syd = Syd, n = n)
}

# Analytic Anderson-Rubin set for one instrument: solve the quadratic
# n (My - th Md)^2 <= c (Syy - 2 th Syd + th^2 Sdd).
.cm_ar_set_analytic <- function(ar, crit) {
  n <- ar$n
  My <- as.numeric(ar$My); Md <- as.numeric(ar$Md)
  Syy <- as.numeric(ar$Syy); Sdd <- as.numeric(ar$Sdd); Syd <- as.numeric(ar$Syd)
  A <- n * Md^2 - crit * Sdd
  B <- -2 * (n * My * Md - crit * Syd)
  C <- n * My^2 - crit * Syy
  disc <- B^2 - 4 * A * C
  if (abs(A) < 1e-14) {
    if (abs(B) < 1e-14) return(if (C <= 0) list(type = "whole line", intervals = data.frame(lower = -Inf, upper = Inf)) else list(type = "empty", intervals = data.frame(lower = numeric(0), upper = numeric(0))))
    root <- -C / B
    iv <- if (B > 0) data.frame(lower = -Inf, upper = root) else data.frame(lower = root, upper = Inf)
    return(list(type = "unbounded", intervals = iv))
  }
  if (disc < 0) {
    if (A < 0) return(list(type = "whole line", intervals = data.frame(lower = -Inf, upper = Inf)))
    return(list(type = "empty", intervals = data.frame(lower = numeric(0), upper = numeric(0))))
  }
  r <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
  if (A > 0) {
    list(type = "bounded", intervals = data.frame(lower = r[1], upper = r[2]))
  } else {
    list(type = "unbounded", intervals = data.frame(lower = c(-Inf, r[2]), upper = c(r[1], Inf)))
  }
}

# Intervals of a grid where `keep` is TRUE.
.cm_grid_intervals <- function(grid, keep) {
  if (!any(keep)) return(data.frame(lower = numeric(0), upper = numeric(0)))
  r <- rle(keep)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  ok <- r$values
  data.frame(lower = grid[starts[ok]], upper = grid[ends[ok]])
}
