# R/iv-weak.R
#
# Weak-instrument tools that two-stage least squares packages do not give in
# one place: first-stage strength including the effective F, the
# Anderson-Rubin confidence set by test inversion, and the plausibly-exogenous
# union of intervals.

#' First-stage strength for one endogenous regressor
#'
#' Reports the homoskedastic first-stage F, the heteroskedasticity- or
#' cluster-robust Wald F, and the effective F of Montiel Olea and Pflueger
#' (2013), after partialling the exogenous controls out of the treatment and
#' the instruments. With one instrument the three coincide up to the variance
#' estimator; with several the effective F is the statistic to compare with
#' the Stock-Yogo style thresholds under heteroskedasticity.
#'
#' @param data A data frame.
#' @param d Endogenous regressor column name.
#' @param z Character vector of instrument column names.
#' @param x Optional character vector of exogenous control column names.
#' @param cluster Optional cluster column name.
#' @param weights Optional weights column name.
#'
#' @return A list of class `cm_first_stage` with `F`, `F_robust`,
#'   `F_effective`, the first-stage `coefficients` and robust `vcov`, the
#'   `t` statistic (one instrument), and the numbers of instruments,
#'   observations, and clusters.
#' @references
#' Montiel Olea, J. L. and Pflueger, C. (2013). A robust test for weak
#' instruments. *Journal of Business & Economic Statistics*, 31(3), 358-369.
#' @examples
#' dat <- sim_iv(1000, dgp = "weak", seed = 1)
#' iv_first_stage(dat, d = "d", z = "z", x = paste0("x", 1:3))
#' @export
iv_first_stage <- function(data, d, z, x = NULL, cluster = NULL, weights = NULL) {
  data <- as.data.frame(data)
  for (v in c(d, z, x, cluster, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(d, z, x, cluster, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  w <- if (is.null(weights)) NULL else data[[weights]]
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  dt <- .cm_residualize(as.numeric(data[[d]]), X, w)[, 1]
  Zt <- .cm_residualize(as.matrix(data[, z, drop = FALSE]), X, w)
  fs <- .cm_first_stage_strength(dt, Zt, w, cl)
  structure(c(fs, list(n = nrow(data), n_clusters = if (is.null(cl)) nrow(data) else length(unique(cl)),
                       instruments = z)), class = "cm_first_stage")
}

#' @export
print.cm_first_stage <- function(x, ...) {
  cat("First stage: ", x$k, " instrument(s), n = ", x$n, if (x$n_clusters < x$n) paste0(", clusters = ", x$n_clusters), "\n", sep = "")
  cat("  F (homoskedastic) = ", format(round(x$F, 2)), "; F (robust) = ", format(round(x$F_robust, 2)),
      "; effective F = ", format(round(x$F_effective, 2)), "\n", sep = "")
  if (is.finite(x$t)) cat("  first-stage t = ", format(round(x$t, 2)), "\n", sep = "")
  invisible(x)
}

#' Anderson-Rubin confidence set for an instrumental-variable coefficient
#'
#' Inverts the Anderson-Rubin test with a heteroskedasticity- or
#' cluster-robust variance: for each candidate `theta` the statistic
#' `C(theta) = n M(theta)' Omega(theta)^-1 M(theta)` with
#' `M(theta) = mean[(Y - theta D) Z]` (after partialling out the controls) is
#' compared with a chi-square critical value with as many degrees of freedom
#' as instruments. The set of non-rejected values is valid whatever the
#' strength of the instrument, and may be unbounded or empty. With one
#' instrument the set is computed analytically; with several, on a grid.
#' The same statistic is Neyman's C(alpha) test, which is why it also applies
#' to the cross-fitted residuals of [est_dml()] (`weak_iv = TRUE`).
#'
#' @inheritParams iv_first_stage
#' @param y Outcome column name.
#' @param conf_level Confidence level.
#' @param theta_grid Grid of candidate values (several instruments); by
#'   default a wide interval around the two-stage least squares estimate.
#'
#' @return A list of class `cm_ar_set` with `intervals` (a data frame of
#'   lower and upper limits; possibly several rows or infinite limits),
#'   `type` (`"bounded"`, `"unbounded"`, `"whole line"`, `"empty"`), the
#'   two-stage least squares `estimate` and `std.error` with their Wald
#'   interval for comparison, `statistic` (the grid, when used), and
#'   `crit_val`. [plot_ar_set()] draws the statistic against `theta`.
#' @references
#' Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a
#' single equation in a complete system of stochastic equations. *Annals of
#' Mathematical Statistics*, 20, 46-63.
#'
#' Andrews, I., Stock, J. H., and Sun, L. (2019). Weak instruments in
#' instrumental variables regression: theory and practice. *Annual Review of
#' Economics*, 11, 727-753.
#' @examples
#' dat <- sim_iv(500, dgp = "weak", seed = 1)
#' iv_ar_confidence_set(dat, y = "y", d = "d", z = "z", x = paste0("x", 1:3))
#' @export
iv_ar_confidence_set <- function(data, y, d, z, x = NULL, cluster = NULL, weights = NULL,
                                 conf_level = 0.95, theta_grid = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, z, x, cluster, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, z, x, cluster, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  w <- if (is.null(weights)) NULL else data[[weights]]
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  ts <- .cm_tsls(as.numeric(data[[y]]), as.numeric(data[[d]]), as.matrix(data[, z, drop = FALSE]), X, w, cl)
  k <- length(z)
  crit <- stats::qchisq(conf_level, df = k)
  if (is.null(theta_grid)) theta_grid <- seq(ts$estimate - 20 * ts$std.error, ts$estimate + 20 * ts$std.error, length.out = 4001L)
  ar <- .cm_ar_statistic(ts$residuals$y, ts$residuals$d, ts$residuals$z, theta_grid, cl, w)
  if (k == 1L) {
    set <- .cm_ar_set_analytic(ar, crit)
  } else {
    keep_g <- ar$statistic <= crit
    iv <- .cm_grid_intervals(theta_grid, keep_g)
    type <- if (nrow(iv) == 0L) "empty" else if (all(keep_g)) "whole line (within the grid)" else if (keep_g[1] || keep_g[length(keep_g)]) "unbounded (reaches the grid edge)" else "bounded"
    set <- list(type = type, intervals = iv)
  }
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  structure(list(intervals = set$intervals, type = set$type, estimate = ts$estimate, std.error = ts$std.error,
                 wald = c(lower = ts$estimate - zq * ts$std.error, upper = ts$estimate + zq * ts$std.error),
                 statistic = data.frame(theta = theta_grid, statistic = ar$statistic), crit_val = crit,
                 conf_level = conf_level, k = k, n = ts$n, n_clusters = ts$n_clusters,
                 first_stage = .cm_first_stage_strength(ts$residuals$d, ts$residuals$z, w, cl)[c("F", "F_robust", "F_effective", "t")]),
            class = "cm_ar_set")
}

#' @export
print.cm_ar_set <- function(x, ...) {
  cat("Anderson-Rubin ", 100 * x$conf_level, "% confidence set (", x$k, " instrument(s), n = ", x$n, ")\n", sep = "")
  cat("  2SLS estimate = ", format(round(x$estimate, 4)), " (SE ", format(round(x$std.error, 4)), "); Wald interval [",
      format(round(x$wald[1], 4)), ", ", format(round(x$wald[2], 4)), "]\n", sep = "")
  cat("  first-stage F (robust) = ", format(round(x$first_stage$F_robust, 2)), "; effective F = ",
      format(round(x$first_stage$F_effective, 2)), "\n", sep = "")
  cat("  Anderson-Rubin set (", x$type, "):\n", sep = "")
  if (nrow(x$intervals) == 0L) cat("    empty\n") else
    for (i in seq_len(nrow(x$intervals))) cat("    [", format(round(x$intervals$lower[i], 4)), ", ", format(round(x$intervals$upper[i], 4)), "]\n", sep = "")
  invisible(x)
}

#' Plot the Anderson-Rubin statistic against the candidate coefficient
#'
#' @param x A `cm_ar_set` object.
#' @return A ggplot object; the horizontal line is the critical value and the
#'   shaded region the confidence set.
#' @export
plot_ar_set <- function(x) {
  df <- x$statistic
  df$inside <- df$statistic <= x$crit_val
  ggplot2::ggplot(df, ggplot2::aes(x = .data$theta, y = .data$statistic)) +
    ggplot2::geom_area(data = df[df$inside, , drop = FALSE], fill = "grey85") +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = x$crit_val, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = x$estimate, colour = "grey40", linetype = "dotted") +
    ggplot2::coord_cartesian(ylim = c(0, min(max(df$statistic), 5 * x$crit_val))) +
    ggplot2::labs(x = "Candidate coefficient", y = "Anderson-Rubin statistic",
                  subtitle = "Dashed: critical value; shaded: confidence set; dotted: 2SLS estimate") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plausibly exogenous instruments: union of intervals
#'
#' Conley, Hansen, and Rossi (2012): allow the instrument a direct effect
#' `gamma` on the outcome, `Y = theta D + gamma Z + ...`. For each `gamma` in
#' `gamma_grid`, subtract `gamma Z` from `Y` and re-estimate by two-stage least
#' squares; the union of the resulting confidence intervals is a confidence
#' set that is valid if the true `gamma` lies in the grid.
#'
#' @inheritParams iv_ar_confidence_set
#' @param gamma_grid Numeric vector of direct effects to entertain.
#'
#' @return A data frame with one row per `gamma` (estimate, standard error,
#'   interval) and attributes `union` (the overall interval) and
#'   `gamma_zero` (the interval at `gamma = 0`).
#' @references
#' Conley, T. G., Hansen, C. B., and Rossi, P. E. (2012). Plausibly exogenous.
#' *Review of Economics and Statistics*, 94(1), 260-272.
#' @examples
#' dat <- sim_iv(1000, dgp = "linear", seed = 1)
#' iv_plausibly_exogenous(dat, "y", "d", "z", paste0("x", 1:5), gamma_grid = seq(-0.2, 0.2, by = 0.1))
#' @export
iv_plausibly_exogenous <- function(data, y, d, z, x = NULL, gamma_grid = seq(-0.5, 0.5, by = 0.1),
                                   cluster = NULL, weights = NULL, conf_level = 0.95) {
  data <- as.data.frame(data)
  if (length(z) != 1L) stop("`z` must name one instrument.", call. = FALSE)
  for (v in c(y, d, z, x, cluster, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, z, x, cluster, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  w <- if (is.null(weights)) NULL else data[[weights]]
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  zv <- as.numeric(data[[z]])
  rows <- lapply(gamma_grid, function(g) {
    ts <- .cm_tsls(as.numeric(data[[y]]) - g * zv, as.numeric(data[[d]]), matrix(zv, ncol = 1), X, w, cl)
    data.frame(gamma = g, estimate = ts$estimate, std.error = ts$std.error,
               conf.low = ts$estimate - zq * ts$std.error, conf.high = ts$estimate + zq * ts$std.error)
  })
  out <- do.call(rbind, rows)
  attr(out, "union") <- c(lower = min(out$conf.low), upper = max(out$conf.high))
  attr(out, "gamma_zero") <- if (any(abs(gamma_grid) < 1e-12)) unlist(out[which.min(abs(gamma_grid)), c("conf.low", "conf.high")]) else NULL
  out
}
