#' Synthetic difference-in-differences weights
#'
#' `sdid_weights()` computes the unit and time weights of synthetic
#' difference-in-differences (Arkhangelsky et al. 2021) for a block design in
#' which a set of units adopts treatment in the same period. It returns the
#' weights, the point estimate, and a cell-weight column so that the
#' difference-in-differences step can be run as a weighted two-way fixed
#' effects regression with `fixest::feols(y ~ d | id + time, weights = w)`.
#' The same function produces the weights of the synthetic control estimator
#' (`"sc"`), the plain difference-in-differences (`"did"`), and the
#' difference-in-differences with unit weights only (`"difp"`), so the four
#' estimators can be compared on one footing.
#'
#' @inheritParams att_gt
#' @param d Column name of the treatment indicator (0/1; 1 for treated units
#'   in post-treatment periods).
#' @param estimator `"sdid"` (default), `"sc"`, `"did"`, or `"difp"`.
#' @param zeta_omega,zeta_lambda Regularization for the unit and time
#'   weights. Defaults follow Arkhangelsky et al. (2021):
#'   `zeta_omega = (N1 * T1)^(1/4) * sigma` and `zeta_lambda = 1e-6 * sigma`
#'   for `"sdid"`, where `sigma` is the standard deviation of first
#'   differences of untreated outcomes before treatment; `"sc"` and `"difp"`
#'   use `1e-6 * sigma` for the unit weights.
#' @param sparsify Zero out unit and time weights below a quarter of the
#'   largest weight and re-optimize, as `synthdid` does.
#' @param max_iter,min_decrease Frank-Wolfe iteration limit and stopping rule.
#'
#' @return An object of class `cm_sdid` with `estimate`, `omega` (unit
#'   weights of the untreated units), `lambda` (time weights of the
#'   pre-treatment periods), `cell_weights` (a data frame with `id`, `time`,
#'   and `weight` for the weighted regression), `sigma`, `zeta`, the
#'   dimensions `N0`, `N1`, `T0`, `T1`, the outcome matrix, and the call.
#'
#' @details Unit weights solve
#'   `min ||omega_0 + Y_pre' omega - Ybar_treated_pre||^2 + zeta^2 T0 ||omega||^2`
#'   over the simplex, and time weights solve the analogous problem that
#'   matches pre-treatment to post-treatment outcomes of the untreated units.
#'   The estimate is the weighted double difference, which equals the
#'   coefficient on `d` in the two-way fixed effects regression weighted by
#'   `omega_i * lambda_t` (treated cells weight one).
#'
#' @references Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W.,
#'   and Wager, S. (2021). Synthetic difference-in-differences. *American
#'   Economic Review*, 111(12), 4088-4118.
#' @examples
#' dat <- sim_did_panel(n_units = 40, n_periods = 12, groups = 9, never_share = 0.8, seed = 1)
#' w <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "treated")
#' w
#' fit <- fixest::feols(y ~ treated | id + time, data = merge(dat, w$cell_weights), weights = ~weight)
#' coef(fit)
#' @export
sdid_weights <- function(data, id, time, y, d, estimator = c("sdid", "sc", "did", "difp"),
                         zeta_omega = NULL, zeta_lambda = NULL, sparsify = TRUE,
                         max_iter = 10000L, min_decrease = NULL) {
  estimator <- match.arg(estimator)
  for (v in c(id, time, y, d)) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, c(id, time, y, d), with = FALSE]
  data.table::setnames(dt, c(".id", ".t", ".y", ".d"))
  dt[, .d := .cm_as_binary(.d, "d")]
  setup <- .cm_sdid_setup(dt)
  fit <- .cm_sdid_fit(setup$Y, setup$N0, setup$T0, estimator, zeta_omega, zeta_lambda, sparsify, max_iter, min_decrease)
  N <- nrow(setup$Y)
  Tn <- ncol(setup$Y)
  unit_w <- c(fit$omega, rep(1 / setup$N1, setup$N1))
  time_w <- c(fit$lambda, rep(1 / setup$T1, setup$T1))
  cw <- expand.grid(.row = seq_len(N), .col = seq_len(Tn))
  cw$weight <- unit_w[cw$.row] * time_w[cw$.col]
  cell_weights <- data.frame(setup$unit_names[cw$.row], setup$time_values[cw$.col], weight = cw$weight, stringsAsFactors = FALSE)
  names(cell_weights)[1:2] <- c(id, time)
  out <- list(
    estimate = fit$estimate, estimator = estimator,
    omega = data.frame(id = setup$unit_names[seq_len(setup$N0)], weight = fit$omega, stringsAsFactors = FALSE),
    lambda = data.frame(time = setup$time_values[seq_len(setup$T0)], weight = fit$lambda),
    cell_weights = cell_weights,
    sigma = fit$sigma, zeta = c(omega = fit$zeta_omega, lambda = fit$zeta_lambda),
    N0 = setup$N0, N1 = setup$N1, T0 = setup$T0, T1 = setup$T1,
    Y = setup$Y, unit_names = setup$unit_names, time_values = setup$time_values,
    settings = list(sparsify = sparsify, max_iter = max_iter, min_decrease = min_decrease,
                    zeta_omega = zeta_omega, zeta_lambda = zeta_lambda),
    columns = c(id = id, time = time, y = y, d = d),
    call = match.call()
  )
  names(out$omega)[1L] <- id
  names(out$lambda)[1L] <- time
  class(out) <- "cm_sdid"
  out
}

# Arrange a long panel into the block matrix Y (controls first, then treated;
# pre-periods first, then post) and check the block design.
.cm_sdid_setup <- function(dt) {
  tlist <- sort(unique(dt$.t))
  units <- unique(dt$.id)
  cnt <- dt[, .N, by = .id]
  if (any(cnt$N != length(tlist)) || anyDuplicated(dt[, list(.id, .t)])) stop("A balanced panel is required.", call. = FALSE)
  first <- dt[.d == 1L, list(g = min(.t)), by = .id]
  ever <- units %in% first$.id
  if (!any(ever) || all(ever)) stop("Both treated and untreated units are required.", call. = FALSE)
  if (length(unique(first$g)) > 1L) stop("All treated units must adopt in the same period (block design). Apply the estimator cohort by cohort otherwise.", call. = FALSE)
  g <- first$g[1L]
  chk <- dt[, all(.d == as.integer(.t >= g & .id %in% first$.id)), by = .id]
  if (!all(chk$V1)) stop("`d` must equal 1 exactly for treated units in periods at or after adoption.", call. = FALSE)
  T0 <- sum(tlist < g)
  if (T0 < 2L) stop("At least two pre-treatment periods are required.", call. = FALSE)
  ord_units <- c(units[!ever], units[ever])
  Y <- matrix(NA_real_, length(units), length(tlist))
  Y[cbind(match(dt$.id, ord_units), match(dt$.t, tlist))] <- dt$.y
  list(Y = Y, N0 = sum(!ever), N1 = sum(ever), T0 = T0, T1 = length(tlist) - T0,
       unit_names = ord_units, time_values = tlist)
}

# Frank-Wolfe step and solver for the regularized simplex regression
# (following synthdid).
.cm_fw_step <- function(A, x, b, eta) {
  Ax <- A %*% x
  half_grad <- as.numeric(t(Ax - b) %*% A) + eta * x
  i <- which.min(half_grad)
  dx <- -x
  dx[i] <- 1 - x[i]
  if (all(dx == 0)) return(x)
  d_err <- A[, i] - Ax
  step <- -sum(half_grad * dx) / (sum(d_err^2) + eta * sum(dx^2))
  x + min(1, max(0, step)) * dx
}

.cm_simplex_fw <- function(Y, zeta, intercept = TRUE, x0 = NULL, min_decrease = 1e-3, max_iter = 1000L) {
  T0 <- ncol(Y) - 1L
  N0 <- nrow(Y)
  x <- if (is.null(x0)) rep(1 / T0, T0) else x0
  if (intercept) Y <- sweep(Y, 2, colMeans(Y))
  A <- Y[, seq_len(T0), drop = FALSE]
  b <- Y[, T0 + 1L]
  eta <- N0 * zeta^2
  vals <- numeric(max_iter)
  t <- 0L
  while (t < max_iter && (t < 2L || vals[t - 1L] - vals[t] > min_decrease^2)) {
    t <- t + 1L
    x <- .cm_fw_step(A, x, b, eta)
    err <- A %*% x - b
    vals[t] <- zeta^2 * sum(x^2) + sum(err^2) / N0
  }
  list(weights = as.numeric(x), iterations = t, value = vals[t])
}

.cm_sdid_sparsify <- function(v) {
  v[v <= max(v) / 4] <- 0
  v / sum(v)
}

.cm_sdid_fit <- function(Y, N0, T0, estimator, zeta_omega = NULL, zeta_lambda = NULL,
                         sparsify = TRUE, max_iter = 10000L, min_decrease = NULL) {
  N1 <- nrow(Y) - N0
  T1 <- ncol(Y) - T0
  sigma <- stats::sd(apply(Y[seq_len(N0), seq_len(T0), drop = FALSE], 1, diff))
  if (is.null(min_decrease)) min_decrease <- 1e-5 * sigma
  # collapsed form: control rows, treated mean row; pre columns, post mean column
  Yc <- rbind(cbind(Y[seq_len(N0), seq_len(T0), drop = FALSE], rowMeans(Y[seq_len(N0), T0 + seq_len(T1), drop = FALSE])),
              c(colMeans(Y[N0 + seq_len(N1), seq_len(T0), drop = FALSE]), mean(Y[N0 + seq_len(N1), T0 + seq_len(T1)])))
  if (estimator == "did") {
    omega <- rep(1 / N0, N0)
    lambda <- rep(1 / T0, T0)
    zo <- NA_real_; zl <- NA_real_
  } else {
    zo <- if (!is.null(zeta_omega)) zeta_omega else if (estimator == "sdid") (N1 * T1)^(1 / 4) * sigma else 1e-6 * sigma
    zl <- if (!is.null(zeta_lambda)) zeta_lambda else 1e-6 * sigma
    intercept_omega <- estimator != "sc"
    solve_w <- function(M, zeta, intercept) {
      if (zeta == 0 && !sparsify) {
        return(.cm_simplex_fw(M, zeta, intercept, max_iter = max_iter, min_decrease = min_decrease)$weights)
      }
      w <- .cm_simplex_fw(M, zeta, intercept, min_decrease = min_decrease, max_iter = if (sparsify) 100L else max_iter)$weights
      if (sparsify) {
        w <- .cm_simplex_fw(M, zeta, intercept, x0 = .cm_sdid_sparsify(w), min_decrease = min_decrease, max_iter = max_iter)$weights
      }
      w
    }
    omega <- solve_w(t(Yc[, seq_len(T0), drop = FALSE]), zo, intercept_omega)
    if (estimator %in% c("sc", "difp")) {
      lambda <- if (estimator == "sc") rep(0, T0) else rep(1 / T0, T0)
      zl <- NA_real_
    } else {
      lambda <- solve_w(Yc[seq_len(N0), , drop = FALSE], zl, TRUE)
    }
  }
  unit_contrast <- c(-omega, rep(1 / N1, N1))
  time_contrast <- c(-lambda, rep(1 / T1, T1))
  estimate <- as.numeric(t(unit_contrast) %*% Y %*% time_contrast)
  list(estimate = estimate, omega = omega, lambda = lambda, sigma = sigma, zeta_omega = zo, zeta_lambda = zl)
}

#' @export
print.cm_sdid <- function(x, digits = 4, ...) {
  label <- switch(x$estimator, sdid = "Synthetic difference-in-differences", sc = "Synthetic control",
                  did = "Difference-in-differences", difp = "Difference-in-differences with unit weights (DIFP)")
  cat(label, ": estimate ", formatC(x$estimate, digits = digits, format = "g"), "\n", sep = "")
  cat("  ", x$N1, " treated and ", x$N0, " untreated units; ", x$T0, " pre- and ", x$T1, " post-treatment periods\n", sep = "")
  if (x$estimator != "did") {
    cat("  Regularization zeta: ", formatC(x$zeta[["omega"]], digits = 3, format = "g"),
        " (sigma = ", formatC(x$sigma, digits = 3, format = "g"), "); nonzero unit weights: ",
        sum(x$omega$weight > 1e-8), "; nonzero time weights: ", sum(x$lambda$weight > 1e-8), "\n", sep = "")
  }
  cat("  Use `cell_weights` as regression weights: feols(y ~ d | id + time, weights = ~weight)\n")
  invisible(x)
}

#' Standard errors for synthetic difference-in-differences
#'
#' Recomputes the weights and the estimate under placebo, bootstrap, or
#' jackknife resampling of units (Arkhangelsky et al. 2021, Section 5) and
#' returns the standard error of the estimate in `x`.
#'
#' @param x An object from [sdid_weights()].
#' @param method `"placebo"` (default; treats random sets of untreated units
#'   as if treated, valid with few treated units), `"bootstrap"` (resample
#'   units with replacement), or `"jackknife"` (leave one unit out; needs at
#'   least two treated units).
#' @param n_reps Number of replications for placebo and bootstrap.
#' @param seed Optional seed.
#' @return A list with `std.error`, `estimate`, `conf.low`, `conf.high`,
#'   `replicates`, and `method`.
#' @export
sdid_se <- function(x, method = c("placebo", "bootstrap", "jackknife"), n_reps = 200L, seed = NULL) {
  method <- match.arg(method)
  if (!inherits(x, "cm_sdid")) stop("`x` must come from sdid_weights().", call. = FALSE)
  Y <- x$Y; N0 <- x$N0; N1 <- x$N1; T0 <- x$T0
  s <- x$settings
  refit <- function(Ym, n0) {
    .cm_sdid_fit(Ym, n0, T0, x$estimator, s$zeta_omega, s$zeta_lambda, s$sparsify, s$max_iter, s$min_decrease)$estimate
  }
  reps <- .cm_with_seed(seed, {
    if (method == "placebo") {
      if (N0 <= N1) stop("Placebo standard errors need more untreated than treated units.", call. = FALSE)
      vapply(seq_len(n_reps), function(b) {
        pseudo <- sample(N0, N1)
        Ym <- rbind(Y[setdiff(seq_len(N0), pseudo), , drop = FALSE], Y[pseudo, , drop = FALSE])
        refit(Ym, N0 - N1)
      }, numeric(1))
    } else if (method == "bootstrap") {
      vapply(seq_len(n_reps), function(b) {
        repeat {
          idx <- sample(nrow(Y), replace = TRUE)
          n0 <- sum(idx <= N0)
          if (n0 >= 2L && n0 < nrow(Y)) break
        }
        Ym <- Y[c(idx[idx <= N0], idx[idx > N0]), , drop = FALSE]
        refit(Ym, n0)
      }, numeric(1))
    } else {
      if (N1 < 2L) stop("The jackknife needs at least two treated units.", call. = FALSE)
      vapply(seq_len(nrow(Y)), function(i) refit(Y[-i, , drop = FALSE], N0 - as.integer(i <= N0)), numeric(1))
    }
  })
  se <- if (method == "jackknife") {
    n <- length(reps)
    sqrt((n - 1) / n * sum((reps - mean(reps))^2))
  } else {
    sqrt((length(reps) - 1) / length(reps)) * stats::sd(reps)
  }
  z <- stats::qnorm(0.975)
  list(estimate = x$estimate, std.error = se, conf.low = x$estimate - z * se, conf.high = x$estimate + z * se,
       replicates = reps, method = method, n_reps = length(reps))
}

#' Plot synthetic difference-in-differences trajectories or weights
#'
#' @param x An object from [sdid_weights()].
#' @param type `"trends"` (treated average against the weighted untreated
#'   average, with the time weights shown as bars) or `"weights"` (unit
#'   weights).
#' @return A `ggplot` object.
#' @export
plot_sdid <- function(x, type = c("trends", "weights")) {
  type <- match.arg(type)
  if (type == "weights") {
    d <- x$omega
    names(d) <- c("id", "weight")
    d$id <- factor(d$id, levels = d$id[order(d$weight)])
    return(ggplot2::ggplot(d[d$weight > 0, ], ggplot2::aes(x = .data$id, y = .data$weight)) +
             ggplot2::geom_col(fill = "grey40") + ggplot2::coord_flip() +
             ggplot2::labs(x = NULL, y = "Unit weight") + ggplot2::theme_minimal(base_size = 11))
  }
  Y <- x$Y; N0 <- x$N0
  treated <- colMeans(Y[(N0 + 1):nrow(Y), , drop = FALSE])
  control <- as.numeric(t(x$omega$weight) %*% Y[seq_len(N0), , drop = FALSE])
  if (x$estimator == "sc") {
    control <- control
  }
  d <- data.frame(time = rep(x$time_values, 2), y = c(treated, control),
                  series = rep(c("Treated", "Weighted untreated"), each = length(x$time_values)))
  lam <- data.frame(time = x$lambda[[1L]], weight = x$lambda$weight)
  rng <- range(d$y)
  lam$height <- rng[1] + lam$weight / max(c(lam$weight, 1e-12)) * diff(rng) * 0.25
  ggplot2::ggplot() +
    ggplot2::geom_col(data = lam, ggplot2::aes(x = .data$time, y = .data$height), fill = "grey80", width = 0.8, alpha = 0.7) +
    ggplot2::geom_line(data = d, ggplot2::aes(x = .data$time, y = .data$y, colour = .data$series)) +
    ggplot2::geom_vline(xintercept = x$time_values[x$T0] + 0.5, linetype = "dotted") +
    ggplot2::labs(x = NULL, y = x$columns[["y"]], colour = NULL,
                  subtitle = paste0(toupper(x$estimator), " estimate ", formatC(x$estimate, digits = 3, format = "g"),
                                    "; bars show time weights")) +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}
