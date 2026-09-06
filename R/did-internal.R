# R/did-internal.R
#
# Internals shared by att_gt(), aggregate_att(), did_imputation(), att_dose()
# and the other difference-in-differences tools: 2x2 cell estimators that
# reproduce Sant'Anna and Zhao (2020) / DRDID exactly, learner-based cells with
# the plain orthogonal score, the multiplier bootstrap on influence functions,
# and small design helpers. Nothing in this file is exported.

utils::globalVariables(c(
  ".N", ".cl", ".cm_id", ".d", ".data", ".dose", ".e", ".g", ".id", ".lead", ".t", ".tau",
  ".treated", ".w", ".y", ".y0", ":=", "N", "cl", "dtil", "eps", "g", "weight", "ybar", ".g0"
))

# Design matrices ------------------------------------------------------------

# Build the covariate design for a cell. `x` is NULL (intercept only), a
# character vector of column names, or a one-sided formula. Always includes an
# intercept, as DRDID does.
.cm_did_design <- function(data, x) {
  if (is.null(x)) {
    return(matrix(1, nrow = nrow(data), ncol = 1L, dimnames = list(NULL, "(Intercept)")))
  }
  if (inherits(x, "formula")) {
    fml <- x
  } else if (is.character(x)) {
    for (v in x) .cm_check_column(v, data)
    fml <- stats::reformulate(x)
  } else {
    stop("`x` must be NULL, a character vector of column names, or a one-sided formula.", call. = FALSE)
  }
  mm <- stats::model.matrix(fml, data = as.data.frame(data), na.action = stats::na.pass)
  if (!"(Intercept)" %in% colnames(mm)) mm <- cbind("(Intercept)" = 1, mm)
  storage.mode(mm) <- "double"
  mm
}

.cm_did_x_vars <- function(x) {
  if (is.null(x)) return(character(0))
  if (inherits(x, "formula")) return(all.vars(x))
  x
}

# Weighted logit and weighted least squares ----------------------------------

.cm_did_logit <- function(X, d, w) {
  fit <- suppressWarnings(stats::glm.fit(
    X, d, weights = w, family = stats::binomial(),
    control = stats::glm.control(epsilon = 1e-10, maxit = 100L)
  ))
  if (!isTRUE(fit$converged)) {
    warning("Propensity score estimation did not converge in a 2x2 cell.", call. = FALSE)
  }
  if (anyNA(fit$coefficients)) {
    stop("Propensity score coefficients contain NA: covariates are collinear or lack variation.", call. = FALSE)
  }
  as.numeric(fit$fitted.values)
}

.cm_did_wols <- function(X, y, w) {
  fit <- stats::lm.wfit(X, y, w)
  if (anyNA(fit$coefficients)) {
    stop("Outcome regression coefficients contain NA: covariates are collinear or lack variation.", call. = FALSE)
  }
  as.numeric(fit$coefficients)
}

.cm_did_check_rcond <- function(M, what) {
  if (!is.finite(rcond(M)) || rcond(M) < .Machine$double.eps) {
    stop("The ", what, " design matrix is singular in a 2x2 cell; remove collinear covariates.", call. = FALSE)
  }
  invisible(TRUE)
}

# Panel 2x2 cell (Sant'Anna and Zhao 2020; DRDID::drdid_panel, reg_did_panel,
# std_ipw_did_panel) ---------------------------------------------------------
#
# y1, y0: outcome in the current and base period; d: 1 for the treated group,
# 0 for the comparison group; X: design with intercept; w: sampling weights.
# Returns the ATT, the influence function on the cell sample, and the fitted
# propensity score. The influence functions include the first-step
# corrections for the logit and OLS coefficients, so they coincide with DRDID.
.cm_did_panel_cell <- function(y1, y0, d, X, w, method = c("dr", "reg", "ipw"), p_trim = 0.995) {
  method <- match.arg(method)
  n <- length(d)
  dy <- as.numeric(y1 - y0)
  w <- w / mean(w)
  ps <- NULL

  if (method %in% c("dr", "ipw")) {
    ps <- pmin(.cm_did_logit(X, d, w), 1 - 1e-6)
    trim <- rep(1, n)
    trim[d == 0] <- as.numeric(ps[d == 0] < p_trim)
    W_ps <- ps * (1 - ps) * w
    XtWX <- crossprod(X, W_ps * X)
    .cm_did_check_rcond(XtWX, "propensity score")
    H_ps <- chol2inv(chol(XtWX)) * n
    alr_ps <- (w * (d - ps) * X) %*% H_ps
  }
  if (method %in% c("dr", "reg")) {
    ctrl <- d == 0
    coef_or <- .cm_did_wols(X[ctrl, , drop = FALSE], dy[ctrl], w[ctrl])
    m_hat <- as.numeric(X %*% coef_or)
    w_ols <- w * (1 - d)
    XpX <- crossprod(w_ols * X, X) / n
    .cm_did_check_rcond(XpX, "outcome regression")
    alr_or <- (w_ols * (dy - m_hat) * X) %*% solve(XpX)
  }

  if (method == "dr") {
    w_t <- trim * w * d
    w_c <- trim * w * ps * (1 - d) / (1 - ps)
    mw_t <- mean(w_t)
    mw_c <- mean(w_c)
    s_t <- w_t * (dy - m_hat)
    s_c <- w_c * (dy - m_hat)
    eta_t <- mean(s_t) / mw_t
    eta_c <- mean(s_c) / mw_c
    att <- eta_t - eta_c
    inf_t1 <- s_t - w_t * eta_t
    inf_t2 <- alr_or %*% (colSums(w_t * X) / n)
    inf_c1 <- s_c - w_c * eta_c
    inf_c2 <- alr_ps %*% (colSums(w_c * (dy - m_hat - eta_c) * X) / n)
    inf_c3 <- alr_or %*% (colSums(w_c * X) / n)
    inf <- (inf_t1 - inf_t2) / mw_t - (inf_c1 + inf_c2 - inf_c3) / mw_c
  } else if (method == "reg") {
    w_t <- w * d
    mw_t <- mean(w_t)
    eta_t <- mean(w_t * dy) / mw_t
    eta_c <- mean(w_t * m_hat) / mw_t
    att <- eta_t - eta_c
    inf_t <- (w_t * dy - w_t * eta_t) / mw_t
    inf_c1 <- w_t * m_hat - w_t * eta_c
    inf_c2 <- alr_or %*% (colSums(w_t * X) / n)
    inf <- inf_t - (inf_c1 + inf_c2) / mw_t
  } else {
    w_t <- trim * w * d
    w_c <- trim * w * ps * (1 - d) / (1 - ps)
    mw_t <- mean(w_t)
    mw_c <- mean(w_c)
    eta_t <- mean(w_t * dy) / mw_t
    eta_c <- mean(w_c * dy) / mw_c
    att <- eta_t - eta_c
    inf_t <- (w_t * dy - w_t * eta_t) / mw_t
    inf_c1 <- w_c * dy - w_c * eta_c
    inf_c2 <- alr_ps %*% (colSums(w_c * (dy - eta_c) * X) / n)
    inf <- inf_t - (inf_c1 + inf_c2) / mw_c
  }
  list(att = unname(att), inffunc = as.numeric(inf), ps = ps)
}

# Repeated cross-section 2x2 cell (DRDID::drdid_rc, reg_did_rc,
# std_ipw_did_rc). `post` marks rows in the current period; `d` marks the
# treated group. The DR version is the locally efficient one of Sant'Anna and
# Zhao (2020), which also models the treated outcomes.
.cm_did_rc_cell <- function(y, post, d, X, w, method = c("dr", "reg", "ipw"), p_trim = 0.995) {
  method <- match.arg(method)
  n <- length(d)
  y <- as.numeric(y)
  w <- w / mean(w)
  ps <- NULL

  if (method %in% c("dr", "ipw")) {
    ps <- pmin(.cm_did_logit(X, d, w), 1 - 1e-6)
    trim <- rep(1, n)
    trim[d == 0] <- as.numeric(ps[d == 0] < p_trim)
    W_ps <- ps * (1 - ps) * w
    XtWX <- crossprod(X, W_ps * X)
    .cm_did_check_rcond(XtWX, "propensity score")
    H_ps <- chol2inv(chol(XtWX)) * n
    alr_ps <- (w * (d - ps) * X) %*% H_ps
  }
  fit_block <- function(mask, target, what) {
    coef <- .cm_did_wols(X[mask, , drop = FALSE], target[mask], w[mask])
    pred <- as.numeric(X %*% coef)
    w_ols <- w * mask
    XpX <- crossprod(w_ols * X, X) / n
    .cm_did_check_rcond(XpX, what)
    list(pred = pred, alr = (w_ols * (target - pred) * X) %*% solve(XpX))
  }
  if (method %in% c("dr", "reg")) {
    c_pre <- fit_block(d == 0 & post == 0, y, "pre-period control regression")
    c_post <- fit_block(d == 0 & post == 1, y, "post-period control regression")
    m_c <- post * c_post$pred + (1 - post) * c_pre$pred
  }

  if (method == "dr") {
    t_pre <- fit_block(d == 1 & post == 0, y, "pre-period treated regression")
    t_post <- fit_block(d == 1 & post == 1, y, "post-period treated regression")
    w_tpre <- trim * w * d * (1 - post)
    w_tpost <- trim * w * d * post
    w_cpre <- trim * w * ps * (1 - d) * (1 - post) / (1 - ps)
    w_cpost <- trim * w * ps * (1 - d) * post / (1 - ps)
    w_d <- trim * w * d
    w_dt1 <- trim * w * d * post
    w_dt0 <- trim * w * d * (1 - post)
    mw <- function(v) mean(v)
    eta_tpre <- w_tpre * (y - m_c) / mw(w_tpre)
    eta_tpost <- w_tpost * (y - m_c) / mw(w_tpost)
    eta_cpre <- w_cpre * (y - m_c) / mw(w_cpre)
    eta_cpost <- w_cpost * (y - m_c) / mw(w_cpost)
    eta_d_post <- w_d * (t_post$pred - c_post$pred) / mw(w_d)
    eta_dt1_post <- w_dt1 * (t_post$pred - c_post$pred) / mw(w_dt1)
    eta_d_pre <- w_d * (t_pre$pred - c_pre$pred) / mw(w_d)
    eta_dt0_pre <- w_dt0 * (t_pre$pred - c_pre$pred) / mw(w_dt0)
    a_tpre <- mean(eta_tpre); a_tpost <- mean(eta_tpost)
    a_cpre <- mean(eta_cpre); a_cpost <- mean(eta_cpost)
    a_d_post <- mean(eta_d_post); a_dt1_post <- mean(eta_dt1_post)
    a_d_pre <- mean(eta_d_pre); a_dt0_pre <- mean(eta_dt0_pre)
    att <- (a_tpost - a_tpre) - (a_cpost - a_cpre) + (a_d_post - a_dt1_post) - (a_d_pre - a_dt0_pre)

    inf_tpre <- eta_tpre - w_tpre * a_tpre / mw(w_tpre)
    inf_tpost <- eta_tpost - w_tpost * a_tpost / mw(w_tpost)
    M1_post <- -colSums(w_tpost * X) / n / mw(w_tpost)
    M1_pre <- -colSums(w_tpre * X) / n / mw(w_tpre)
    inf_t_or <- c_post$alr %*% M1_post + c_pre$alr %*% M1_pre
    inf_cpre <- eta_cpre - w_cpre * a_cpre / mw(w_cpre)
    inf_cpost <- eta_cpost - w_cpost * a_cpost / mw(w_cpost)
    M2_pre <- colSums(w_cpre * (y - m_c - a_cpre) * X) / n / mw(w_cpre)
    M2_post <- colSums(w_cpost * (y - m_c - a_cpost) * X) / n / mw(w_cpost)
    inf_c_ps <- alr_ps %*% (M2_post - M2_pre)
    M3_post <- -colSums(w_cpost * X) / n / mw(w_cpost)
    M3_pre <- -colSums(w_cpre * X) / n / mw(w_cpre)
    inf_c_or <- c_post$alr %*% M3_post + c_pre$alr %*% M3_pre
    inf_eff <- (eta_d_post - w_d * a_d_post / mw(w_d)) - (eta_dt1_post - w_dt1 * a_dt1_post / mw(w_dt1)) -
      ((eta_d_pre - w_d * a_d_pre / mw(w_d)) - (eta_dt0_pre - w_dt0 * a_dt0_pre / mw(w_dt0)))
    mom_post <- colSums((w_d / mw(w_d) - w_dt1 / mw(w_dt1)) * X) / n
    mom_pre <- colSums((w_d / mw(w_d) - w_dt0 / mw(w_dt0)) * X) / n
    inf_or <- (t_post$alr - c_post$alr) %*% mom_post - (t_pre$alr - c_pre$alr) %*% mom_pre
    inf_t <- inf_tpost - inf_tpre + inf_t_or
    inf_c <- inf_cpost - inf_cpre + inf_c_ps + inf_c_or
    inf <- (inf_t - inf_c) + inf_eff + inf_or
  } else if (method == "reg") {
    w_tpre <- w * d * (1 - post)
    w_tpost <- w * d * post
    w_c <- w * d
    eta_tpre <- mean(w_tpre * y) / mean(w_tpre)
    eta_tpost <- mean(w_tpost * y) / mean(w_tpost)
    eta_c <- mean(w_c * (c_post$pred - c_pre$pred)) / mean(w_c)
    att <- (eta_tpost - eta_tpre) - eta_c
    inf_t <- (w_tpost * y - w_tpost * eta_tpost) / mean(w_tpost) -
      (w_tpre * y - w_tpre * eta_tpre) / mean(w_tpre)
    inf_c1 <- w_c * (c_post$pred - c_pre$pred) - w_c * eta_c
    M1 <- colSums(w_c * X) / n
    inf_c <- (inf_c1 + c_post$alr %*% M1 - c_pre$alr %*% M1) / mean(w_c)
    inf <- inf_t - inf_c
  } else {
    w_tpre <- trim * w * d * (1 - post)
    w_tpost <- trim * w * d * post
    w_cpre <- trim * w * ps * (1 - d) * (1 - post) / (1 - ps)
    w_cpost <- trim * w * ps * (1 - d) * post / (1 - ps)
    eta_tpre <- w_tpre * y / mean(w_tpre)
    eta_tpost <- w_tpost * y / mean(w_tpost)
    eta_cpre <- w_cpre * y / mean(w_cpre)
    eta_cpost <- w_cpost * y / mean(w_cpost)
    a_tpre <- mean(eta_tpre); a_tpost <- mean(eta_tpost)
    a_cpre <- mean(eta_cpre); a_cpost <- mean(eta_cpost)
    att <- (a_tpost - a_tpre) - (a_cpost - a_cpre)
    inf_t <- (eta_tpost - w_tpost * a_tpost / mean(w_tpost)) - (eta_tpre - w_tpre * a_tpre / mean(w_tpre))
    inf_c <- (eta_cpost - w_cpost * a_cpost / mean(w_cpost)) - (eta_cpre - w_cpre * a_cpre / mean(w_cpre))
    M2_pre <- colSums(w_cpre * (y - a_cpre) * X) / n / mean(w_cpre)
    M2_post <- colSums(w_cpost * (y - a_cpost) * X) / n / mean(w_cpost)
    inf_c <- inf_c + alr_ps %*% (M2_post - M2_pre)
    inf <- inf_t - inf_c
  }
  list(att = unname(att), inffunc = as.numeric(inf), ps = ps)
}

# Panel 2x2 cell with cross-fitted machine-learning nuisances. The score is
# the Hajek-normalized doubly robust score with the outcome regression fitted
# on the comparison group only; the influence function is the plain
# orthogonal score (no first-step correction).
.cm_did_panel_cell_learner <- function(dy, d, work, features, learner_p, learner_or,
                                       folds, seed, p_trim = 0.995, p_clip = 0.995) {
  n <- length(d)
  work <- as.data.frame(work)
  work$.cm_dy <- dy
  work$.cm_d <- d
  fold_id <- .cm_make_folds(d, folds, seed)
  if (is.null(learner_p)) learner_p <- .cm_default_learner("classif")
  if (is.null(learner_or)) learner_or <- .cm_default_learner("regr")
  ps <- .cm_crossfit_predict(work, ".cm_d", features, learner_p, fold_id, TRUE,
                             positive = "1", task_hint = "propensity", what = "the propensity score")$pred
  ps <- pmin(pmax(ps, 1e-6), p_clip)
  m_hat <- .cm_crossfit_predict(work, ".cm_dy", features, learner_or, fold_id, TRUE,
                                subset = d == 0, task_hint = "outcome", what = "the outcome regression")$pred
  trim <- rep(1, n)
  trim[d == 0] <- as.numeric(ps[d == 0] < p_trim)
  w_t <- trim * d
  w_c <- trim * ps * (1 - d) / (1 - ps)
  mw_t <- mean(w_t)
  mw_c <- mean(w_c)
  s_t <- w_t * (dy - m_hat)
  s_c <- w_c * (dy - m_hat)
  eta_t <- mean(s_t) / mw_t
  eta_c <- mean(s_c) / mw_c
  att <- eta_t - eta_c
  inf <- (s_t - w_t * eta_t) / mw_t - (s_c - w_c * eta_c) / mw_c
  list(att = unname(att), inffunc = as.numeric(inf), ps = ps, m_hat = m_hat, fold_id = fold_id)
}

# Influence-function inference ------------------------------------------------

# Collapse an n x K influence-function matrix (dense or sparse) to cluster
# sums. Returns the matrix and the number of clusters.
.cm_if_cluster_sums <- function(inffunc, cluster = NULL) {
  if (is.null(cluster)) {
    return(list(mat = inffunc, n_clusters = nrow(inffunc)))
  }
  cl <- as.integer(factor(cluster))
  if (inherits(inffunc, "sparseMatrix")) {
    agg <- Matrix::sparseMatrix(i = cl, j = seq_along(cl), x = 1, dims = c(max(cl), length(cl)))
    mat <- agg %*% inffunc
  } else {
    mat <- rowsum(as.matrix(inffunc), cl, reorder = TRUE)
  }
  list(mat = mat, n_clusters = max(cl))
}

# Analytic standard errors: sqrt(sum_c (sum_{i in c} IF_i)^2) / n.
.cm_if_analytic_se <- function(inffunc, n, cluster = NULL) {
  cs <- .cm_if_cluster_sums(inffunc, cluster)
  mat <- cs$mat
  ss <- if (inherits(mat, "sparseMatrix")) Matrix::colSums(mat^2) else colSums(as.matrix(mat)^2)
  sqrt(as.numeric(ss)) / n
}

# Multiplier bootstrap of Callaway and Sant'Anna (2021) on an n x K
# influence-function matrix: draws multiplier weights at the cluster level,
# returns bootstrap standard errors (interquartile-range based) and the
# sup-t critical value for a simultaneous band. Draws are processed in chunks
# so a sparse IF with a million rows never becomes dense.
.cm_multiplier_bootstrap <- function(inffunc, n, n_boot = 999L, cluster = NULL,
                                     boot_weights = c("mammen", "rademacher"),
                                     conf_level = 0.95, seed = NULL, chunk_cells = 2e7) {
  boot_weights <- match.arg(boot_weights)
  cs <- .cm_if_cluster_sums(inffunc, cluster)
  mat <- cs$mat
  n_c <- cs$n_clusters
  K <- ncol(mat)
  chunk <- max(1L, min(n_boot, floor(chunk_cells / max(n_c, 1))))
  draw <- function(m) {
    if (boot_weights == "mammen") {
      k1 <- 0.5 * (1 - sqrt(5))
      k2 <- 0.5 * (1 + sqrt(5))
      pk <- 0.5 * (1 + sqrt(5)) / sqrt(5)
      u <- stats::runif(n_c * m)
      matrix(ifelse(u < pk, k1, k2), n_c, m)
    } else {
      matrix(sample(c(-1, 1), n_c * m, replace = TRUE), n_c, m)
    }
  }
  bres <- .cm_with_seed(seed, {
    out <- matrix(NA_real_, n_boot, K)
    done <- 0L
    while (done < n_boot) {
      m <- min(chunk, n_boot - done)
      U <- draw(m)
      prod <- if (inherits(mat, "sparseMatrix")) Matrix::crossprod(mat, U) else crossprod(as.matrix(mat), U)
      out[(done + 1L):(done + m), ] <- t(as.matrix(prod)) / n_c
      done <- done + m
    }
    out
  })
  bres <- sqrt(n_c) * bres
  ok <- !is.na(colSums(bres)) & colSums(bres^2) > sqrt(.Machine$double.eps) * 10
  iqr_norm <- stats::qnorm(0.75) - stats::qnorm(0.25)
  b_sigma <- apply(bres, 2, function(b) {
    b <- sort.int(b)
    nb <- length(b)
    (b[ceiling(0.75 * nb)] - b[ceiling(0.25 * nb)]) / iqr_norm
  })
  se <- rep(NA_real_, K)
  se[ok] <- b_sigma[ok] * sqrt(n_c) / n
  crit_val <- NA_real_
  if (any(ok)) {
    scaled <- abs(sweep(bres[, ok, drop = FALSE], 2, b_sigma[ok], "/"))
    bT <- apply(scaled, 1, max)
    bT <- bT[is.finite(bT)]
    crit_val <- as.numeric(stats::quantile(bT, conf_level, type = 1, na.rm = TRUE))
  }
  list(se = se, crit_val = crit_val, boot = bres, n_clusters = n_c)
}

# Wald test that a set of estimates is jointly zero, using the influence
# functions for the covariance (clustered when `cluster` is given).
.cm_if_wald <- function(estimates, inffunc, n, cluster = NULL) {
  keep <- which(is.finite(estimates))
  if (length(keep) == 0L) return(list(statistic = NA_real_, df = 0L, p.value = NA_real_))
  cs <- .cm_if_cluster_sums(inffunc[, keep, drop = FALSE], cluster)
  mat <- as.matrix(cs$mat)
  V <- crossprod(mat) / n^2
  est <- estimates[keep]
  stat <- tryCatch(as.numeric(t(est) %*% solve(V, est)), error = function(e) NA_real_)
  list(statistic = stat, df = length(keep),
       p.value = if (is.finite(stat)) stats::pchisq(stat, df = length(keep), lower.tail = FALSE) else NA_real_)
}
