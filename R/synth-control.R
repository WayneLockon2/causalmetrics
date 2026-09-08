# R/synth-control.R
#
# Synthetic control from a long panel: predictors (outcome lags and covariate
# averages), V weights, simplex weights with optional intercept (demeaned SC),
# ridge regularization, ridge augmentation, and a per-unit mode for staggered
# adoption. Inference and robustness tools live in R/synth-inference.R.

#' Synthetic control
#'
#' `synth_control()` builds the synthetic control of Abadie and Gardeazabal
#' (2003) and Abadie, Diamond, and Hainmueller (2010) from a long panel: a
#' weighted average of untreated donor units, with weights on the simplex,
#' that reproduces the treated unit's pre-treatment predictors. The
#' predictors are outcome lags and pre-treatment averages of covariates,
#' each weighted by an importance weight `V`. Options give the demeaned
#' synthetic control of Ferman and Pinto (2021) (`demean = TRUE`), relaxed
#' constraints (Doudchenko and Imbens 2016), a ridge penalty on the weights,
#' the ridge-augmented estimator of Ben-Michael, Feller, and Rothstein (2021)
#' (`augment = "ridge"`), and one synthetic control per treated unit for
#' staggered adoption (`treated_units = "separate"`).
#'
#' @param data A balanced long panel.
#' @param id,time,y Column names of the unit identifier, the period, and the
#'   outcome.
#' @param d Column name of the treatment indicator (1 for treated units in
#'   periods at or after adoption, 0 otherwise).
#' @param x Character vector of covariate columns used as predictors; each
#'   is averaged over `pre_window`.
#' @param lags `"all"` (every pre-treatment period of the outcome is a
#'   predictor; the Ferman-Pinto and `synthdid` choice) or a vector of
#'   periods whose outcomes are predictors (Andersson 2019 uses 1989, 1980,
#'   and 1970).
#' @param pre_window Periods over which covariates are averaged. Default:
#'   all pre-treatment periods.
#' @param v Predictor importance weights. `"equal"` (default), `"mspe"`
#'   (nested optimization: outer search over `V` minimizing the pre-treatment
#'   mean squared prediction error of the outcome, inner simplex problem;
#'   the default of the `Synth` package), `"regression"` (regression-based
#'   weights, `Synth`'s starting value), or a numeric vector with one entry
#'   per predictor. Predictors are standardized by their standard deviation
#'   across units before weighting.
#' @param demean If `TRUE`, outcome lags enter as deviations from each
#'   unit's pre-treatment mean and the counterfactual adds back an
#'   intercept: the demeaned synthetic control of Ferman and Pinto (2021),
#'   equivalent to the intercept of Doudchenko and Imbens (2016).
#' @param constraints `"simplex"` (weights nonnegative and summing to one),
#'   `"nonnegative"`, or `"none"` (unconstrained least squares).
#' @param zeta Ridge penalty on the weights, on the scale of
#'   [sdid_weights()]: the objective adds `zeta^2 * K * sum(w^2)` with `K`
#'   predictors.
#' @param augment `"none"` or `"ridge"`: add the bias correction
#'   `m(z_0) - sum_j w_j m(z_j)` from a ridge regression of each period's
#'   outcome on the predictors, fitted on the donors.
#' @param lambda Ridge penalty of the augmentation model; chosen by
#'   leave-one-donor-out cross-validation on post-treatment outcomes when
#'   `NULL`.
#' @param treated_units `"average"` fits one synthetic control to the average
#'   of the treated units (they must adopt in the same period);
#'   `"separate"` fits one synthetic control per treated unit, using as
#'   donors the units untreated throughout that unit's window, and averages
#'   the gaps in event time. Staggered adoption is allowed in the second
#'   mode.
#' @param horizon In the separate mode, the number of post-treatment periods
#'   kept for each unit (default: all available).
#' @param weights An earlier `cm_synth` object whose donor weights are
#'   applied to the outcome `y` without re-estimation (for example, the
#'   synthetic GDP check of Andersson 2019).
#' @param v_control List with `maxit` (outer iterations, default 500) and
#'   `starts` (number of random starting values in addition to the equal
#'   and regression-based starts, default 0) for `v = "mspe"`.
#' @param standardize Divide each predictor by its standard deviation across
#'   units before weighting, as the `Synth` package does. Default `TRUE`
#'   when covariates are present and `FALSE` for outcome lags alone, so that
#'   the plain synthetic control is the least-squares problem on raw
#'   outcomes.
#'
#' @return An object of class `cm_synth` with `estimate` (average
#'   post-treatment gap), `effects` (a data frame with the treated and
#'   synthetic paths and the gap in every period), `weights` (donor
#'   weights), `v`, `balance` (treated, synthetic, and donor-mean value of
#'   each predictor, with a range check), `intercept`, `augmentation`
#'   (the per-period correction when `augment = "ridge"`), `lambda`,
#'   `pre_rmspe`, `post_rmspe`, `concentration` (`l2`, the squared norm of
#'   the weights of Ferman and Pinto, and `effective_donors = 1 / l2`), the
#'   block dimensions, the block matrices, and the settings. In the separate
#'   mode the object holds `units` (one `cm_synth` per treated unit),
#'   `unit_table`, `event_gap`, and a jackknife `std.error` over treated
#'   units.
#'
#' @details With `lags = "all"`, `demean = FALSE`, and `zeta = 0` the weights
#'   equal those of `sdid_weights(estimator = "sc")`; with `demean = TRUE`
#'   the average post-treatment gap equals the `"difp"` estimate of
#'   [sdid_weights()]. Weights are solved with `quadprog` when installed and
#'   with the Frank-Wolfe routine of [sdid_weights()] otherwise.
#'
#' @references Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic
#'   control methods for comparative case studies. *JASA*, 105(490), 493-505.
#'   Ferman, B. and Pinto, C. (2021). Synthetic controls with imperfect
#'   pretreatment fit. *Quantitative Economics*, 12(4), 1197-1221.
#'   Ben-Michael, E., Feller, A., and Rothstein, J. (2021). The augmented
#'   synthetic control method. *JASA*, 116(536), 1789-1803.
#'   Doudchenko, N. and Imbens, G. W. (2016). Balancing, regression,
#'   difference-in-differences and synthetic control methods: a synthesis.
#'   NBER Working Paper 22791.
#' @examples
#' dat <- sim_synth_panel(n_donors = 15, t_pre = 15, t_post = 5, effect = 2, seed = 1)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
#' fit
#' fit_dm <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
#' c(sc = fit$estimate, demeaned = fit_dm$estimate)
#' @export
synth_control <- function(data, id, time, y, d, x = NULL, lags = "all", pre_window = NULL,
                          v = c("equal", "mspe", "regression"), demean = FALSE,
                          constraints = c("simplex", "nonnegative", "none"), zeta = 0,
                          augment = c("none", "ridge"), lambda = NULL,
                          treated_units = c("average", "separate"), horizon = NULL,
                          weights = NULL, v_control = list(), standardize = NULL) {
  if (!is.numeric(v)) v <- match.arg(v)
  constraints <- match.arg(constraints)
  augment <- match.arg(augment)
  treated_units <- match.arg(treated_units)
  for (nm in c(id, time, y, d, x)) .cm_check_column(nm, data)
  dt <- data.table::as.data.table(data)[, c(id, time, y, d, x), with = FALSE]
  data.table::setnames(dt, c(id, time, y, d), c(".id", ".t", ".y", ".d"))
  dt[, .d := .cm_as_binary(.d, "d")]
  if (!is.null(weights) && !inherits(weights, "cm_synth")) stop("`weights` must be a `cm_synth` object.", call. = FALSE)
  if (is.null(standardize)) standardize <- length(x) > 0L
  settings <- list(x = x, lags = lags, pre_window = pre_window, v = v, demean = demean, standardize = standardize, id_name = id,
                   constraints = constraints, zeta = zeta, augment = augment, lambda = lambda,
                   treated_units = treated_units, horizon = horizon, y_name = y,
                   v_control = utils::modifyList(list(maxit = 500L, starts = 0L), v_control))
  design <- .cm_synth_design(dt)
  columns <- c(id = id, time = time, y = y, d = d)
  if (treated_units == "average") {
    treated <- design$units[!is.na(design$adoption)]
    g <- unique(design$adoption[!is.na(design$adoption)])
    if (length(g) > 1L) stop("Treated units adopt in different periods; use `treated_units = \"separate\"`.", call. = FALSE)
    donors <- design$units[is.na(design$adoption)]
    if (length(donors) < 2L) stop("At least two donor units are required.", call. = FALSE)
    blk <- .cm_synth_block(dt, x, treated, donors, design$times, g)
    fit <- .cm_synth_fit(blk, settings, weights = weights)
    fit$mode <- "average"
    fit$columns <- columns
    fit$call <- match.call()
    class(fit) <- "cm_synth"
    return(fit)
  }
  # separate mode: one synthetic control per treated unit
  treated <- design$units[!is.na(design$adoption)]
  if (!is.null(weights)) stop("`weights` can only be reused in the average mode.", call. = FALSE)
  units <- lapply(treated, function(u) {
    g <- design$adoption[match(u, design$units)]
    pos_g <- match(g, design$times)
    t_end <- if (is.null(horizon)) max(design$times) else design$times[min(length(design$times), pos_g + horizon - 1L)]
    window <- design$times[design$times <= t_end]
    donors <- design$units[design$units != u & (is.na(design$adoption) | design$adoption > t_end)]
    if (length(donors) < 2L) stop("Fewer than two donors are untreated over the window of unit ", u, ".", call. = FALSE)
    blk <- .cm_synth_block(dt, x, u, donors, window, g)
    f <- .cm_synth_fit(blk, settings)
    f$mode <- "average"
    f$columns <- columns
    f$adoption <- g
    f$effects$event_time <- match(f$effects[[1L]], design$times) - pos_g
    class(f) <- "cm_synth"
    f
  })
  names(units) <- as.character(treated)
  unit_table <- data.frame(
    id = treated,
    adoption = design$adoption[match(treated, design$units)],
    n_donors = vapply(units, function(f) f$N0, numeric(1)),
    pre_rmspe = vapply(units, function(f) f$pre_rmspe, numeric(1)),
    estimate = vapply(units, function(f) f$estimate, numeric(1)),
    stringsAsFactors = FALSE
  )
  names(unit_table)[1L] <- id
  ev <- do.call(rbind, lapply(units, function(f) f$effects[, c("event_time", "gap")]))
  event_gap <- stats::aggregate(gap ~ event_time, data = ev, FUN = mean)
  event_gap$n_units <- stats::aggregate(gap ~ event_time, data = ev, FUN = length)$gap
  ests <- unit_table$estimate
  n <- length(ests)
  se <- if (n >= 2L) {
    loo <- vapply(seq_len(n), function(i) mean(ests[-i]), numeric(1))
    sqrt((n - 1) / n * sum((loo - mean(loo))^2))
  } else NA_real_
  out <- list(
    estimate = mean(ests), std.error = se, units = units, unit_table = unit_table,
    event_gap = event_gap, N1 = n, mode = "separate", settings = settings, columns = columns,
    call = match.call()
  )
  class(out) <- "cm_synth"
  out
}

# Balanced-panel check and adoption periods (NA for never treated).
.cm_synth_design <- function(dt) {
  tlist <- sort(unique(dt$.t))
  units <- unique(dt$.id)
  cnt <- dt[, .N, by = .id]
  if (any(cnt$N != length(tlist)) || anyDuplicated(dt[, list(.id, .t)])) stop("A balanced panel is required.", call. = FALSE)
  first <- dt[.d == 1L, list(g = min(.t)), by = .id]
  adoption <- first$g[match(units, first$.id)]
  if (all(is.na(adoption))) stop("No treated unit found.", call. = FALSE)
  if (!any(is.na(adoption))) stop("No untreated unit found.", call. = FALSE)
  g_all <- adoption[match(dt$.id, units)]
  expected <- as.integer(!is.na(g_all) & dt$.t >= g_all)
  if (any(dt$.d != expected)) stop("`d` must equal 1 exactly for treated units in periods at or after adoption.", call. = FALSE)
  if (any(vapply(adoption[!is.na(adoption)], function(g) sum(tlist < g) < 2L, logical(1)))) stop("At least two pre-treatment periods are required.", call. = FALSE)
  list(units = units, times = tlist, adoption = adoption)
}

# Block matrices for one set of treated units and donors over a time window.
.cm_synth_block <- function(dt, x, treated, donors, times, t_adopt) {
  ord <- c(donors, treated)
  sub <- dt[.id %in% ord & .t %in% times]
  idx <- cbind(match(sub$.id, ord), match(sub$.t, times))
  Y <- matrix(NA_real_, length(ord), length(times))
  Y[idx] <- sub$.y
  if (anyNA(Y)) stop("The outcome has missing values inside the block.", call. = FALSE)
  X <- lapply(x, function(v) { M <- matrix(NA_real_, length(ord), length(times)); M[idx] <- sub[[v]]; M })
  names(X) <- x
  T0 <- sum(times < t_adopt)
  list(Y = Y, X = X, N0 = length(donors), N1 = length(treated), T0 = T0, T1 = length(times) - T0,
       unit_names = ord, time_values = times, treated_names = treated)
}

# Predictor matrix (rows: covariates then outcome lags; columns: units).
.cm_synth_predictors <- function(blk, s) {
  pre_t <- blk$time_values[seq_len(blk$T0)]
  if (identical(s$lags, "all")) lag_t <- pre_t else {
    lag_t <- s$lags
    if (!all(lag_t %in% pre_t)) stop("`lags` must be pre-treatment periods.", call. = FALSE)
  }
  win <- if (is.null(s$pre_window)) pre_t else s$pre_window
  if (!all(win %in% pre_t)) stop("`pre_window` must lie in the pre-treatment periods.", call. = FALSE)
  ybar <- rowMeans(blk$Y[, seq_len(blk$T0), drop = FALSE])
  Zlag <- t(blk$Y[, match(lag_t, blk$time_values), drop = FALSE])
  if (s$demean) Zlag <- sweep(Zlag, 2, ybar)
  rownames(Zlag) <- paste0(s$y_name, "(", lag_t, ")")
  Zcov <- NULL
  if (length(blk$X)) {
    Zcov <- do.call(rbind, lapply(blk$X, function(M) rowMeans(M[, match(win, blk$time_values), drop = FALSE], na.rm = TRUE)))
    rownames(Zcov) <- names(blk$X)
  }
  Z <- rbind(Zcov, Zlag)
  colnames(Z) <- NULL
  list(Z = Z, ybar = ybar, type = c(rep("covariate", length(blk$X)), rep("lag", length(lag_t))), lag_times = lag_t)
}

# Simplex / nonnegative / unconstrained weighted least squares with a ridge term.
.cm_synth_qp <- function(Z0s, z0s, V, constraints, zeta) {
  K <- nrow(Z0s); N0 <- ncol(Z0s)
  A <- Z0s * sqrt(V)
  b <- as.numeric(z0s * sqrt(V))
  ridge <- zeta^2 * K
  D <- crossprod(A) + diag(ridge, N0)
  dvec <- as.numeric(crossprod(A, b))
  if (constraints == "none") {
    jit <- 1e-10 * max(mean(diag(D)), 1e-12)
    return(as.numeric(solve(D + diag(jit, N0), dvec)))
  }
  if (requireNamespace("quadprog", quietly = TRUE)) {
    if (constraints == "simplex") { Amat <- cbind(rep(1, N0), diag(N0)); bvec <- c(1, rep(0, N0)); meq <- 1L }
    else { Amat <- diag(N0); bvec <- rep(0, N0); meq <- 0L }
    jit <- 1e-8 * max(mean(diag(D)), 1e-12)
    res <- NULL
    for (k in 1:6) {
      res <- tryCatch(quadprog::solve.QP(D + diag(jit, N0), dvec, Amat, bvec, meq = meq), error = function(e) NULL)
      if (!is.null(res)) break
      jit <- jit * 100
    }
    if (!is.null(res)) {
      w <- res$solution
      w[w < 0] <- 0
      if (constraints == "simplex") w <- w / sum(w)
      return(w)
    }
  }
  if (constraints == "nonnegative") stop("`constraints = \"nonnegative\"` needs the quadprog package.", call. = FALSE)
  .cm_simplex_fw(cbind(A, b), zeta = zeta, intercept = FALSE, min_decrease = 1e-8, max_iter = 20000L)$weights
}

# Moore-Penrose inverse (small matrices).
.cm_pinv <- function(M, tol = 1e-10) {
  s <- svd(M)
  keep <- s$d > tol * max(s$d)
  if (!any(keep)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

# Regression-based V (Synth's starting value): squared coefficients of the
# pre-treatment outcomes regressed on the standardized predictors.
.cm_synth_v_regression <- function(Zs, Ypre) {
  X <- cbind(1, t(Zs))
  Beta <- .cm_pinv(X) %*% t(Ypre)
  v <- rowSums(Beta[-1L, , drop = FALSE]^2)
  if (sum(v) <= 0) return(rep(1 / nrow(Zs), nrow(Zs)))
  v / sum(v)
}

# Outer search over V minimizing the pre-treatment MSPE of the outcome.
.cm_synth_v_mspe <- function(Z0s, z0s, V_start, mspe_fun, constraints, zeta, control) {
  K <- nrow(Z0s)
  if (K == 1L) return(1)
  obj <- function(theta) {
    V <- exp(theta - max(theta)); V <- V / sum(V)
    mspe_fun(.cm_synth_qp(Z0s, z0s, V, constraints, zeta))
  }
  starts <- list(rep(0, K), log(pmax(V_start, 1e-8)))
  if (control$starts > 0L) starts <- c(starts, lapply(seq_len(control$starts), function(i) stats::rnorm(K)))
  best <- NULL
  for (th in starts) {
    o1 <- stats::optim(th, obj, method = "Nelder-Mead", control = list(maxit = control$maxit))
    o2 <- tryCatch(stats::optim(o1$par, obj, method = "BFGS", control = list(maxit = 100L)), error = function(e) o1)
    o <- if (o2$value <= o1$value) o2 else o1
    if (is.null(best) || o$value < best$value) best <- o
  }
  V <- exp(best$par - max(best$par))
  V / sum(V)
}

# Ridge augmentation (Ben-Michael, Feller, Rothstein 2021) on standardized predictors.
.cm_synth_ridge <- function(Z0s, z0s, Y0, w, lambda, post_idx) {
  F <- t(Z0s)
  N0 <- nrow(F); K <- ncol(F)
  fbar <- colMeans(F)
  Fc <- sweep(F, 2, fbar)
  ybar <- colMeans(Y0)
  Yc <- sweep(Y0, 2, ybar)
  G <- crossprod(Fc)
  if (is.null(lambda)) {
    scale <- max(mean(diag(G)), 1e-12)
    grid <- scale * 10^seq(-3, 3, length.out = 25)
    cv_idx <- if (length(post_idx)) post_idx else seq_len(ncol(Y0))
    cv <- vapply(grid, function(l) {
      H <- Fc %*% solve(G + diag(l, K), t(Fc))
      R <- Yc[, cv_idx, drop = FALSE] - H %*% Yc[, cv_idx, drop = FALSE]
      mean((R / (1 - diag(H)))^2)
    }, numeric(1))
    lambda <- grid[which.min(cv)]
  }
  Beta <- solve(G + diag(lambda, K), crossprod(Fc, Yc))
  m0 <- ybar + as.numeric((z0s - fbar) %*% Beta)
  m_sc <- sum(w) * ybar + as.numeric(crossprod(w, Fc) %*% Beta)
  list(correction = m0 - m_sc, lambda = lambda)
}

# Fit on a block with settings s; `weights` reuses donor weights and V.
.cm_synth_fit <- function(blk, s, weights = NULL) {
  N0 <- blk$N0; N1 <- blk$N1; T0 <- blk$T0; Tn <- ncol(blk$Y)
  post_idx <- if (T0 < Tn) (T0 + 1L):Tn else integer(0)
  P <- .cm_synth_predictors(blk, s)
  Z <- P$Z
  K <- nrow(Z)
  donors_idx <- seq_len(N0); treated_idx <- N0 + seq_len(N1)
  z0 <- rowMeans(Z[, treated_idx, drop = FALSE])
  Z0 <- Z[, donors_idx, drop = FALSE]
  sdv <- if (isTRUE(s$standardize)) apply(cbind(Z0, z0), 1, stats::sd) else rep(1, K)
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  Z0s <- Z0 / sdv; z0s <- z0 / sdv
  Y0 <- blk$Y[donors_idx, , drop = FALSE]
  y_tr <- colMeans(blk$Y[treated_idx, , drop = FALSE])
  ybar0 <- mean(P$ybar[treated_idx]); ybar_d <- P$ybar[donors_idx]
  pre_idx <- seq_len(T0)
  # pre-treatment outcome paths used by the MSPE objective (demeaned if demean)
  y_tr_fit <- if (s$demean) y_tr[pre_idx] - ybar0 else y_tr[pre_idx]
  Y0_fit <- if (s$demean) sweep(Y0[, pre_idx, drop = FALSE], 1, ybar_d) else Y0[, pre_idx, drop = FALSE]
  mspe_fun <- function(w) mean((y_tr_fit - as.numeric(crossprod(Y0_fit, w)))^2)
  if (!is.null(weights)) {
    if (!identical(as.character(weights$weights[[1L]]), as.character(blk$unit_names[donors_idx]))) stop("The donor units of `weights` differ from those in `data`.", call. = FALSE)
    w <- weights$weights$weight
    V <- weights$v$v
    if (length(V) != K) V <- rep(1 / K, K)
  } else {
    V <- if (is.numeric(s$v)) {
      if (length(s$v) != K) stop("`v` must have one entry per predictor (", K, ").", call. = FALSE)
      s$v / sum(s$v)
    } else if (s$v == "equal") rep(1 / K, K)
    else {
      Ypre_all <- blk$Y[, pre_idx, drop = FALSE]
      if (s$demean) Ypre_all <- sweep(Ypre_all, 1, P$ybar)
      V_reg <- .cm_synth_v_regression(cbind(Z0s, z0s), t(Ypre_all))
      if (s$v == "regression") V_reg else .cm_synth_v_mspe(Z0s, z0s, V_reg, mspe_fun, s$constraints, s$zeta, s$v_control)
    }
    w <- .cm_synth_qp(Z0s, z0s, V, s$constraints, s$zeta)
  }
  intercept <- if (s$demean) ybar0 - sum(w * ybar_d) else 0
  y_syn <- as.numeric(crossprod(Y0, w)) + intercept
  aug <- NULL; lambda <- NULL
  if (s$augment == "ridge" && is.null(weights)) {
    Y0_aug <- if (s$demean) sweep(Y0, 1, ybar_d) else Y0
    r <- .cm_synth_ridge(Z0s, z0s, Y0_aug, w, s$lambda, post_idx)
    y_syn <- y_syn + r$correction
    aug <- data.frame(time = blk$time_values, correction = r$correction)
    lambda <- r$lambda
  }
  gap <- y_tr - y_syn
  effects <- data.frame(time = blk$time_values, treated = y_tr, synthetic = y_syn, gap = gap,
                        post = seq_len(Tn) > T0)
  # balance on raw levels: lag rows add back the demeaning intercept
  raw_lag <- t(blk$Y[, match(P$lag_times, blk$time_values), drop = FALSE])
  Zraw <- if (length(blk$X)) rbind(Z[P$type == "covariate", , drop = FALSE], raw_lag) else raw_lag
  z0_raw <- rowMeans(Zraw[, treated_idx, drop = FALSE])
  Z0_raw <- Zraw[, donors_idx, drop = FALSE]
  syn_raw <- as.numeric(Z0_raw %*% w) + ifelse(P$type == "lag", intercept, 0)
  balance <- data.frame(predictor = rownames(Z), type = P$type, treated = z0_raw, synthetic = syn_raw,
                        donor_mean = rowMeans(Z0_raw), v = V,
                        in_range = z0_raw >= apply(Z0_raw, 1, min) & z0_raw <= apply(Z0_raw, 1, max),
                        stringsAsFactors = FALSE)
  rownames(balance) <- NULL
  l2 <- sum(w^2)
  out <- list(
    estimate = if (length(post_idx)) mean(gap[post_idx]) else NA_real_,
    effects = effects,
    weights = data.frame(id = blk$unit_names[donors_idx], weight = w, stringsAsFactors = FALSE),
    v = data.frame(predictor = rownames(Z), v = V, stringsAsFactors = FALSE),
    balance = balance, intercept = intercept, augmentation = aug, lambda = lambda,
    pre_rmspe = sqrt(mean(gap[pre_idx]^2)),
    post_rmspe = if (length(post_idx)) sqrt(mean(gap[post_idx]^2)) else NA_real_,
    concentration = c(l2 = l2, effective_donors = 1 / l2),
    N0 = N0, N1 = N1, T0 = T0, T1 = Tn - T0,
    unit_names = blk$unit_names, time_values = blk$time_values, treated_names = blk$treated_names,
    block = blk, Z = Z, settings = s
  )
  names(out$weights)[1L] <- if (!is.null(s$id_name)) s$id_name else "id"
  names(out$effects)[1L] <- "time"
  out
}

#' @export
print.cm_synth <- function(x, digits = 4, ...) {
  s <- x$settings
  fm <- function(v, d = digits) trimws(formatC(v, digits = d, format = "g"))
  label <- paste0(if (s$demean) "Demeaned synthetic control" else "Synthetic control",
                  if (s$augment == "ridge") ", ridge-augmented" else "")
  if (x$mode == "separate") {
    cat(label, " fitted separately to ", x$N1, " treated unit(s)\n", sep = "")
    cat("  Average post-treatment effect: ", fm(x$estimate),
        " (jackknife SE over units ", fm(x$std.error, 3), ")\n", sep = "")
    cat("  Per-unit estimates in `unit_table`, event-time gaps in `event_gap`, fits in `units`\n")
    return(invisible(x))
  }
  cat(label, ": average post-treatment effect ", fm(x$estimate), "\n", sep = "")
  cat("  ", x$N1, " treated and ", x$N0, " donor units; ", x$T0, " pre- and ", x$T1, " post-treatment periods; ",
      nrow(x$v), " predictors (V: ", if (is.numeric(s$v)) "user" else s$v, ")\n", sep = "")
  cat("  Pre-treatment RMSPE ", fm(x$pre_rmspe, 3),
      "; weight concentration ||w||^2 = ", fm(x$concentration[["l2"]], 3),
      " (", fm(x$concentration[["effective_donors"]], 3), " effective donors)\n", sep = "")
  top <- x$weights[order(-x$weights$weight), ]
  top <- top[top$weight > 1e-3, , drop = FALSE]
  cat("  Weights above 0.001: ", paste0(top[[1L]], " (", formatC(top$weight, digits = 3, format = "f"), ")", collapse = ", "), "\n", sep = "")
  if (!all(x$balance$in_range)) cat("  Note: ", sum(!x$balance$in_range), " predictor(s) of the treated unit lie outside the donors' range\n", sep = "")
  invisible(x)
}

#' Plot a synthetic control fit or its inference objects
#'
#' @param x A `cm_synth`, `cm_synth_placebo`, `cm_synth_loo`, or
#'   `cm_synth_spec` object.
#' @param type For `cm_synth`: `"path"` (treated and synthetic series),
#'   `"gap"`, `"demeaned_path"` (both series after subtracting the
#'   cross-sectional mean of all units in each period, the fit picture of
#'   Ferman and Pinto 2021), `"weights"`, `"balance"` (standardized
#'   predictor values of the treated unit, the synthetic unit, and the donor
#'   mean), or `"event"` (separate mode: average gap by event time). For
#'   `cm_synth_placebo`: `"placebos"` (gap paths) or `"ratio"` (post/pre
#'   MSPE ratios). For `cm_synth_loo`: `"paths"`. For `cm_synth_spec`:
#'   `"paths"` (demeaned synthetic control against DiD).
#' @param mspe_limit For `"placebos"`, drop placebos whose pre-treatment MSPE
#'   exceeds this multiple of the treated unit's (default: the limit stored
#'   in the object).
#' @return A `ggplot` object.
#' @export
plot_synth <- function(x, type = NULL, mspe_limit = NULL) {
  if (inherits(x, "cm_synth_placebo")) return(.cm_plot_synth_placebo(x, if (is.null(type)) "placebos" else type, mspe_limit))
  if (inherits(x, "cm_synth_loo")) return(.cm_plot_synth_loo(x))
  if (inherits(x, "cm_synth_spec")) return(.cm_plot_synth_spec(x))
  if (!inherits(x, "cm_synth")) stop("`x` must come from synth_control() or its inference tools.", call. = FALSE)
  if (x$mode == "separate") {
    if (is.null(type)) type <- "event"
    if (type != "event") stop("In the separate mode use `type = \"event\"` or plot the fits in `x$units`.", call. = FALSE)
    d <- x$event_gap
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$event_time, y = .data$gap)) +
             ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
             ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted") +
             ggplot2::geom_line() + ggplot2::geom_point(ggplot2::aes(size = .data$n_units)) +
             ggplot2::labs(x = "Periods since adoption", y = "Average gap", size = "Units") +
             ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom"))
  }
  type <- match.arg(if (is.null(type)) "path" else type, c("path", "gap", "demeaned_path", "weights", "balance"))
  cut <- x$time_values[x$T0] + 0.5 * (if (x$T1 > 0) x$time_values[x$T0 + 1L] - x$time_values[x$T0] else 1)
  base <- ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  if (type == "path" || type == "demeaned_path") {
    e <- x$effects
    tr <- e$treated; syn <- e$synthetic
    if (type == "demeaned_path") {
      m <- colMeans(x$block$Y)
      tr <- tr - m; syn <- syn - m
    }
    d <- data.frame(time = rep(e$time, 2), y = c(tr, syn),
                    series = rep(c("Treated", "Synthetic"), each = nrow(e)))
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$y, colour = .data$series, linetype = .data$series)) +
             ggplot2::geom_line() + ggplot2::geom_vline(xintercept = cut, linetype = "dotted") +
             ggplot2::scale_colour_manual(values = c(Treated = "black", Synthetic = "grey45")) +
             ggplot2::scale_linetype_manual(values = c(Treated = "solid", Synthetic = "longdash")) +
             ggplot2::labs(x = NULL, colour = NULL, linetype = NULL,
                           y = if (type == "path") x$settings$y_name else paste(x$settings$y_name, "minus period mean of all units")) + base)
  }
  if (type == "gap") {
    return(ggplot2::ggplot(x$effects, ggplot2::aes(x = .data$time, y = .data$gap)) +
             ggplot2::geom_hline(yintercept = 0, linetype = "dotted") + ggplot2::geom_vline(xintercept = cut, linetype = "dotted") +
             ggplot2::geom_line() + ggplot2::labs(x = NULL, y = "Gap: treated minus synthetic") + base)
  }
  if (type == "weights") {
    d <- x$weights; names(d) <- c("id", "weight")
    d <- d[d$weight > 1e-6, ]
    d$id <- factor(d$id, levels = d$id[order(d$weight)])
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$id, y = .data$weight)) + ggplot2::geom_col(fill = "grey40") +
             ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "Donor weight") + base)
  }
  b <- x$balance
  sdv <- apply(cbind(x$Z), 1, stats::sd); sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  d <- data.frame(predictor = rep(b$predictor, 3),
                  value = c(b$treated, b$synthetic, b$donor_mean) / rep(sdv, 3),
                  series = rep(c("Treated", "Synthetic", "Donor mean"), each = nrow(b)))
  d$predictor <- factor(d$predictor, levels = rev(b$predictor))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$predictor, shape = .data$series, colour = .data$series)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(Treated = "black", Synthetic = "grey45", `Donor mean` = "#B22222")) +
    ggplot2::labs(x = "Predictor value (in donor standard deviations)", y = NULL, shape = NULL, colour = NULL) + base
}
