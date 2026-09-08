# R/rd-internal.R
#
# Internals for the regression discontinuity helpers: kernel weights, a small
# local polynomial fit with influence functions (used for joint balance tests
# and the extrapolation curve), argument preparation, and thin calls into
# rdrobust when it is installed. Nothing in this file is exported.

utils::globalVariables(c(
  "bin_x", "bin_y", "side", "conf.low", "conf.high", "estimate", "std.error",
  "radius", "h", "cutoff", "covariate", "term", "method", "n_eff", "x_grid", "fit"
))

.cm_rd_kernel <- function(u, kernel = c("triangular", "uniform", "epanechnikov")) {
  kernel <- match.arg(kernel)
  a <- abs(u)
  switch(kernel,
    triangular = pmax(0, 1 - a),
    uniform = 0.5 * as.numeric(a <= 1),
    epanechnikov = 0.75 * pmax(0, 1 - u^2)
  )
}

# Complete rows and numeric vectors for an RD call.
.cm_rd_prepare <- function(data, y, x, cutoff, d = NULL, covariates = NULL, cluster = NULL,
                           weights = NULL, cols_extra = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data)
  for (v in c(y, x, d, covariates, cluster, weights, cols_extra)) .cm_check_column(v, data)
  if (!is.numeric(cutoff) || length(cutoff) != 1L) stop("`cutoff` must be one number.", call. = FALSE)
  keep <- stats::complete.cases(data[, c(y, x, d, covariates, cluster, weights, cols_extra), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  list(
    data = data,
    y = if (is.null(y)) NULL else as.numeric(data[[y]]),
    x = as.numeric(data[[x]]),
    d = if (is.null(d)) NULL else as.numeric(data[[d]]),
    covs = if (is.null(covariates)) NULL else as.matrix(data[, covariates, drop = FALSE]),
    cluster = if (is.null(cluster)) NULL else data[[cluster]],
    weights = if (is.null(weights)) NULL else as.numeric(data[[weights]]),
    n_dropped = sum(!keep)
  )
}

# Bandwidths from rdrobust::rdbwselect (h and b, left and right).
.cm_rd_bandwidth <- function(y, x, cutoff, d = NULL, covs = NULL, cluster = NULL, weights = NULL,
                             p = 1, kernel = "triangular", bwselect = "mserd", deriv = 0) {
  .cm_check_package("rdrobust")
  bw <- suppressWarnings(rdrobust::rdbwselect(
    y = y, x = x, c = cutoff, fuzzy = d, covs = covs, cluster = cluster, weights = weights,
    p = p, deriv = deriv, kernel = kernel, bwselect = bwselect
  ))
  b <- bw$bws
  list(h_left = b[1, 1], h_right = b[1, 2], b_left = b[1, 3], b_right = b[1, 4],
       h = b[1, 1], b = b[1, 3], bwselect = bwselect)
}

# One rdrobust fit with the package's argument names.
.cm_rd_fit <- function(y, x, cutoff, d = NULL, covs = NULL, cluster = NULL, weights = NULL,
                       h = NULL, b = NULL, p = 1, q = NULL, kernel = "triangular", vce = "nn",
                       deriv = 0, bwselect = "mserd", conf_level = 0.95, subset = NULL, masspoints = "adjust") {
  .cm_check_package("rdrobust")
  args <- list(y = y, x = x, c = cutoff, fuzzy = d, covs = covs, cluster = cluster, weights = weights,
               p = p, deriv = deriv, kernel = kernel, vce = vce, bwselect = bwselect,
               level = 100 * conf_level, masspoints = masspoints, subset = subset)
  if (!is.null(h)) args$h <- h
  if (!is.null(b)) args$b <- b
  if (!is.null(q)) args$q <- q
  args <- args[!vapply(args, is.null, logical(1))]
  suppressWarnings(do.call(rdrobust::rdrobust, args))
}

# Tidy rows of an rdrobust object.
.cm_rd_tidy_rdrobust <- function(fit, term = "RD effect") {
  est <- as.numeric(fit$coef)
  se <- as.numeric(fit$se)
  ci <- fit$ci
  data.frame(
    term = term,
    method = c("conventional", "bias_corrected", "robust"),
    estimate = est, std.error = se, statistic = as.numeric(fit$z), p.value = as.numeric(fit$pv),
    conf.low = ci[, 1], conf.high = ci[, 2],
    h_left = fit$bws[1, 1], h_right = fit$bws[1, 2], b_left = fit$bws[2, 1], b_right = fit$bws[2, 2],
    n_left = fit$N_h[1], n_right = fit$N_h[2],
    stringsAsFactors = FALSE
  )
}

# Local polynomial jump estimates for one or several outcomes at a common
# bandwidth, with influence functions. Y is an n x K matrix. Returns the
# jumps (right minus left intercepts) and an n x K influence-function
# matrix, so that joint Wald tests can be formed with .cm_if_wald().
.cm_rd_local_jump <- function(Y, x, cutoff, h, kernel = "triangular", p = 1, cluster = NULL) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  K <- ncol(Y)
  u <- (x - cutoff) / h
  w <- .cm_rd_kernel(u, kernel)
  right <- x >= cutoff
  jumps <- numeric(K)
  inff <- matrix(0, n, K)
  levels <- matrix(NA_real_, 2, K, dimnames = list(c("left", "right"), colnames(Y)))
  design <- function(idx) {
    xc <- x[idx] - cutoff
    X <- cbind(1, stats::poly(xc, degree = p, raw = TRUE))
    X
  }
  for (side in c("left", "right")) {
    idx <- which(if (side == "right") right & w > 0 else !right & w > 0)
    if (length(idx) <= p + 1L) stop("Too few observations on the ", side, " side within the bandwidth.", call. = FALSE)
    X <- design(idx)
    W <- w[idx]
    XtWX <- crossprod(X * W, X)
    if (.cm_rcond(XtWX) < 1e-12) stop("Local design matrix is singular on the ", side, " side.", call. = FALSE)
    Q <- solve(XtWX)
    beta <- Q %*% crossprod(X * W, Y[idx, , drop = FALSE])
    resid <- Y[idx, , drop = FALSE] - X %*% beta
    lev <- as.numeric(beta[1, ])
    # influence function of the intercept: n * e1' Q x_i w_i e_i
    e1Q <- Q[1, ]
    a <- as.numeric(X %*% e1Q) * W
    contrib <- a * resid * n
    sgn <- if (side == "right") 1 else -1
    inff[idx, ] <- inff[idx, ] + sgn * contrib
    jumps <- jumps + sgn * lev
    levels[side, ] <- lev
  }
  se <- .cm_if_analytic_se(inff, n, cluster)
  list(jumps = jumps, std.error = se, inffunc = inff, levels = levels,
       n_left = sum(!right & w > 0), n_right = sum(right & w > 0))
}

# Bin means on each side of the cutoff (evenly spaced or quantile spaced).
.cm_rd_bins <- function(y, x, cutoff, n_left, n_right, bins = c("qs", "es"), weights = NULL) {
  bins <- match.arg(bins)
  one_side <- function(idx, J, side) {
    xs <- x[idx]; ys <- y[idx]
    ws <- if (is.null(weights)) rep(1, length(idx)) else weights[idx]
    if (length(xs) == 0L || J < 1L) return(NULL)
    edges <- if (bins == "es") seq(min(xs), max(xs), length.out = J + 1L) else
      unique(stats::quantile(xs, probs = seq(0, 1, length.out = J + 1L), type = 7, names = FALSE))
    if (length(edges) < 2L) edges <- range(xs)
    g <- cut(xs, breaks = edges, include.lowest = TRUE, labels = FALSE)
    out <- lapply(sort(unique(g)), function(k) {
      i <- g == k
      m <- sum(ws[i] * ys[i]) / sum(ws[i])
      se <- if (sum(i) > 1L) stats::sd(ys[i]) / sqrt(sum(i)) else NA_real_
      data.frame(side = side, bin = k, bin_x = sum(ws[i] * xs[i]) / sum(ws[i]), bin_y = m,
                 std.error = se, n = sum(i), x_min = edges[k], x_max = edges[k + 1L])
    })
    do.call(rbind, out)
  }
  rbind(one_side(which(x < cutoff), n_left, "left"), one_side(which(x >= cutoff), n_right, "right"))
}
