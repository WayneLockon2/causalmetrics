# R/rd-checks.R
#
# The regression discontinuity diagnostic battery: covariate balance with a
# joint test, placebo cutoffs, bandwidth sensitivity, and donut-hole
# estimates, each a loop of rdrobust calls returning a tidy frame.

#' Covariate balance at the cutoff
#'
#' Runs the RD estimator with each covariate as the outcome, which tests
#' whether predetermined characteristics jump at the cutoff, and adds a
#' joint Wald test that all jumps are zero. Individual rows use
#' `rdrobust::rdrobust()` with robust bias-corrected inference; the joint
#' test uses the package's own local linear fit at one common bandwidth
#' with influence functions stacked across covariates, because `rdrobust`
#' does not expose the covariance across outcomes.
#'
#' @param data A data frame.
#' @param covariates Character vector of covariate column names.
#' @param x Running variable column name.
#' @param cutoff Cutoff value.
#' @param h Common bandwidth for the joint test and, if given, for every
#'   covariate; `NULL` lets `rdrobust` select a bandwidth per covariate and
#'   uses their median for the joint test.
#' @param p Local polynomial order.
#' @param kernel Kernel.
#' @param cluster Optional cluster column name.
#' @param conf_level Confidence level.
#'
#' @return A list of class `cm_rd_balance` with `table` (covariate,
#'   estimate, std.error, p.value, conf.low, conf.high, h, n_left,
#'   n_right, all from the robust rows), `joint` (statistic, df, p.value,
#'   the common bandwidth), and `cutoff`.
#' @examples
#' dat <- sim_rd(2000, "covariates", seed = 1)
#' rd_balance(dat, paste0("z", 1:4), "x")
#' @export
rd_balance <- function(data, covariates, x, cutoff = 0, h = NULL, p = 1, kernel = "triangular",
                       cluster = NULL, conf_level = 0.95) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, NULL, x, cutoff, covariates = covariates, cluster = cluster)
  rows <- lapply(covariates, function(v) {
    f <- .cm_rd_fit(as.numeric(pr$data[[v]]), pr$x, cutoff, cluster = pr$cluster, h = h, p = p,
                    kernel = kernel, conf_level = conf_level)
    t <- .cm_rd_tidy_rdrobust(f, term = v)
    t <- t[t$method == "robust", ]
    data.frame(covariate = v, estimate = t$estimate, std.error = t$std.error, statistic = t$statistic,
               p.value = t$p.value, conf.low = t$conf.low, conf.high = t$conf.high,
               h = t$h_left, n_left = t$n_left, n_right = t$n_right, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  h_joint <- if (is.null(h)) stats::median(tab$h) else h
  Y <- as.matrix(pr$data[, covariates, drop = FALSE])
  Y <- scale(Y, center = TRUE, scale = TRUE)
  lj <- .cm_rd_local_jump(Y, pr$x, cutoff, h_joint, kernel = kernel, p = p, cluster = pr$cluster)
  joint <- .cm_if_wald(lj$jumps, lj$inffunc, nrow(Y), pr$cluster)
  joint$h <- h_joint
  structure(list(table = tab, joint = joint, cutoff = cutoff, conf_level = conf_level, x = x),
            class = "cm_rd_balance")
}

#' @export
print.cm_rd_balance <- function(x, ...) {
  cat("Covariate balance at the cutoff (robust bias-corrected rows)\n")
  print(x$table[, c("covariate", "estimate", "std.error", "p.value", "h", "n_left", "n_right")], digits = 4, row.names = FALSE)
  cat("  joint test of no jumps at h = ", format(round(x$joint$h, 4)), ": chi2(", x$joint$df, ") = ",
      format(round(x$joint$statistic, 2)), ", p = ", format(signif(x$joint$p.value, 3)), "\n", sep = "")
  invisible(x)
}

#' Placebo cutoffs
#'
#' Estimates the RD effect at cutoffs where no treatment changes, using
#' only observations on the same side of the true cutoff as the placebo
#' (Imbens and Lemieux 2008), so that the true discontinuity never enters a
#' placebo estimate.
#'
#' @inheritParams rd_balance
#' @param y Outcome column name.
#' @param cutoffs Placebo cutoffs; `NULL` uses the quartiles and medians of
#'   the running variable on each side.
#' @param d Optional treatment column (fuzzy design).
#' @param ... Passed to `rdrobust::rdrobust()`.
#' @return A data frame of class `cm_rd_placebo` with one row per cutoff
#'   (robust rows), including the true cutoff.
#' @export
rd_placebo_cutoffs <- function(data, y, x, cutoff = 0, cutoffs = NULL, d = NULL, cluster = NULL,
                               conf_level = 0.95, ...) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, cluster = cluster)
  if (is.null(cutoffs)) {
    left <- pr$x[pr$x < cutoff]; right <- pr$x[pr$x >= cutoff]
    cutoffs <- c(stats::quantile(left, c(0.25, 0.5, 0.75), names = FALSE),
                 stats::quantile(right, c(0.25, 0.5, 0.75), names = FALSE))
  }
  all_c <- sort(unique(c(cutoffs, cutoff)))
  rows <- lapply(all_c, function(cc) {
    if (cc == cutoff) {
      keep <- rep(TRUE, length(pr$x)); side <- "true cutoff"
    } else if (cc < cutoff) {
      keep <- pr$x < cutoff; side <- "left"
    } else {
      keep <- pr$x >= cutoff; side <- "right"
    }
    f <- tryCatch(.cm_rd_fit(pr$y[keep], pr$x[keep], cc, d = pr$d[keep], cluster = pr$cluster[keep],
                             conf_level = conf_level, ...), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    t <- .cm_rd_tidy_rdrobust(f)
    t <- t[t$method == "robust", ]
    data.frame(cutoff = cc, side = side, estimate = t$estimate, std.error = t$std.error, p.value = t$p.value,
               conf.low = t$conf.low, conf.high = t$conf.high, h = t$h_left, n_left = t$n_left, n_right = t$n_right,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  class(out) <- c("cm_rd_placebo", "data.frame")
  attr(out, "true_cutoff") <- cutoff
  out
}

#' Bandwidth sensitivity
#'
#' Re-estimates the RD effect over a grid of bandwidths, keeping the ratio
#' of the bias bandwidth to the main bandwidth at its data-driven value.
#'
#' @inheritParams rd_placebo_cutoffs
#' @param h_grid Bandwidths; `NULL` uses 0.5 to 2 times the MSE-optimal one.
#' @param rho Ratio `h / b`; `NULL` uses the ratio from the MSE-optimal
#'   selection.
#' @return A data frame of class `cm_rd_sensitivity` (robust rows per
#'   bandwidth) with attribute `"h_mse"`.
#' @export
rd_sensitivity <- function(data, y, x, cutoff = 0, h_grid = NULL, rho = NULL, d = NULL, cluster = NULL,
                           conf_level = 0.95, ...) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, cluster = cluster)
  base <- .cm_rd_fit(pr$y, pr$x, cutoff, d = pr$d, cluster = pr$cluster, conf_level = conf_level, ...)
  h_mse <- base$bws[1, 1]; b_mse <- base$bws[2, 1]
  if (is.null(rho)) rho <- h_mse / b_mse
  if (is.null(h_grid)) h_grid <- h_mse * seq(0.5, 2, by = 0.25)
  rows <- lapply(h_grid, function(hh) {
    f <- tryCatch(.cm_rd_fit(pr$y, pr$x, cutoff, d = pr$d, cluster = pr$cluster, h = hh, b = hh / rho,
                             conf_level = conf_level, ...), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    t <- .cm_rd_tidy_rdrobust(f)
    t <- t[t$method == "robust", ]
    data.frame(h = hh, relative_h = hh / h_mse, estimate = t$estimate, std.error = t$std.error, p.value = t$p.value,
               conf.low = t$conf.low, conf.high = t$conf.high, n_left = t$n_left, n_right = t$n_right)
  })
  out <- do.call(rbind, rows)
  class(out) <- c("cm_rd_sensitivity", "data.frame")
  attr(out, "h_mse") <- h_mse
  out
}

#' Donut-hole estimates
#'
#' Drops observations within a radius of the cutoff before estimating, the
#' check against heaping and manipulation right at the threshold (Barreca
#' et al. 2011).
#'
#' @inheritParams rd_placebo_cutoffs
#' @param radius Radii to drop; `NULL` uses 0 and the 1, 2.5, 5, and 10
#'   percent quantiles of the distance to the cutoff.
#' @return A data frame of class `cm_rd_donut` (robust rows per radius).
#' @export
rd_donut <- function(data, y, x, cutoff = 0, radius = NULL, d = NULL, cluster = NULL, conf_level = 0.95, ...) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, cluster = cluster)
  dist <- abs(pr$x - cutoff)
  if (is.null(radius)) radius <- c(0, stats::quantile(dist, c(0.01, 0.025, 0.05, 0.10), names = FALSE))
  rows <- lapply(radius, function(r) {
    keep <- dist >= r
    f <- tryCatch(.cm_rd_fit(pr$y[keep], pr$x[keep], cutoff, d = pr$d[keep], cluster = pr$cluster[keep],
                             conf_level = conf_level, ...), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    t <- .cm_rd_tidy_rdrobust(f)
    t <- t[t$method == "robust", ]
    data.frame(radius = r, n_dropped = sum(!keep), estimate = t$estimate, std.error = t$std.error, p.value = t$p.value,
               conf.low = t$conf.low, conf.high = t$conf.high, h = t$h_left, n_left = t$n_left, n_right = t$n_right)
  })
  out <- do.call(rbind, rows)
  class(out) <- c("cm_rd_donut", "data.frame")
  out
}

#' The regression discontinuity diagnostic battery
#'
#' Runs the main estimate, covariate balance with a joint test, placebo
#' cutoffs, bandwidth sensitivity, and donut-hole estimates in one call and
#' returns them as one object with a print method and
#' [plot_rd_checks()].
#'
#' @inheritParams rd_placebo_cutoffs
#' @param covariates Optional covariates for the balance table.
#' @param h_grid,radius,cutoffs Grids passed to the individual checks.
#' @return A list of class `cm_rd_checks` with `main` (tidy rows of the
#'   main `rdrobust` fit), `balance`, `placebo`, `sensitivity`, `donut`.
#' @examples
#' dat <- sim_rd(2000, "covariates", seed = 1)
#' ck <- rd_checks(dat, "y", "x", covariates = paste0("z", 1:4))
#' ck
#' plot_rd_checks(ck, "sensitivity")
#' @export
rd_checks <- function(data, y, x, cutoff = 0, covariates = NULL, d = NULL, cluster = NULL,
                      h_grid = NULL, radius = NULL, cutoffs = NULL, conf_level = 0.95, ...) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, cluster = cluster, covariates = covariates)
  main_fit <- .cm_rd_fit(pr$y, pr$x, cutoff, d = pr$d, cluster = pr$cluster, conf_level = conf_level, ...)
  out <- list(
    main = .cm_rd_tidy_rdrobust(main_fit),
    main_fit = main_fit,
    balance = if (!is.null(covariates)) rd_balance(pr$data, covariates, x, cutoff, cluster = cluster, conf_level = conf_level) else NULL,
    placebo = rd_placebo_cutoffs(pr$data, y, x, cutoff, cutoffs = cutoffs, d = d, cluster = cluster, conf_level = conf_level, ...),
    sensitivity = rd_sensitivity(pr$data, y, x, cutoff, h_grid = h_grid, d = d, cluster = cluster, conf_level = conf_level, ...),
    donut = rd_donut(pr$data, y, x, cutoff, radius = radius, d = d, cluster = cluster, conf_level = conf_level, ...),
    cutoff = cutoff, y = y, x = x, conf_level = conf_level
  )
  class(out) <- "cm_rd_checks"
  out
}

#' @export
print.cm_rd_checks <- function(x, ...) {
  cat("RD diagnostics at cutoff ", x$cutoff, "\n", sep = "")
  cat("\nMain estimate:\n")
  print(x$main[, c("method", "estimate", "std.error", "conf.low", "conf.high", "h_left", "n_left", "n_right")], digits = 4, row.names = FALSE)
  if (!is.null(x$balance)) { cat("\n"); print(x$balance) }
  cat("\nPlacebo cutoffs (robust):\n")
  print(as.data.frame(x$placebo)[, c("cutoff", "side", "estimate", "std.error", "p.value", "n_left", "n_right")], digits = 4, row.names = FALSE)
  cat("\nBandwidth sensitivity (robust):\n")
  print(as.data.frame(x$sensitivity)[, c("h", "relative_h", "estimate", "std.error", "conf.low", "conf.high")], digits = 4, row.names = FALSE)
  cat("\nDonut-hole estimates (robust):\n")
  print(as.data.frame(x$donut)[, c("radius", "n_dropped", "estimate", "std.error", "conf.low", "conf.high")], digits = 4, row.names = FALSE)
  invisible(x)
}

#' Plot the RD diagnostic battery
#'
#' @param x A `cm_rd_checks` object, or one of the frames returned by
#'   [rd_balance()], [rd_placebo_cutoffs()], [rd_sensitivity()],
#'   [rd_donut()].
#' @param what Which panel: `"balance"`, `"placebo"`, `"sensitivity"`, or
#'   `"donut"` (ignored when `x` is one of the frames).
#' @return A ggplot object.
#' @export
plot_rd_checks <- function(x, what = c("sensitivity", "placebo", "donut", "balance")) {
  what <- match.arg(what)
  if (inherits(x, "cm_rd_checks")) {
    obj <- switch(what, balance = x$balance, placebo = x$placebo, sensitivity = x$sensitivity, donut = x$donut)
    if (is.null(obj)) stop("The `", what, "` check is not in this object.", call. = FALSE)
    main_est <- x$main$estimate[x$main$method == "robust"]
  } else {
    obj <- x
    main_est <- NULL
    what <- if (inherits(x, "cm_rd_balance")) "balance" else if (inherits(x, "cm_rd_placebo")) "placebo" else
      if (inherits(x, "cm_rd_sensitivity")) "sensitivity" else "donut"
  }
  base <- ggplot2::theme_minimal(base_size = 11)
  if (what == "balance") {
    tab <- obj$table
    tab$t <- tab$estimate / tab$std.error
    tab$covariate <- factor(tab$covariate, levels = rev(tab$covariate))
    return(ggplot2::ggplot(tab, ggplot2::aes(x = .data$t, y = .data$covariate)) +
             ggplot2::geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
             ggplot2::geom_vline(xintercept = 0, colour = "grey30") +
             ggplot2::geom_point(size = 2.5) +
             ggplot2::labs(x = "Jump at the cutoff, robust t-statistic", y = NULL,
                           subtitle = paste0("Joint test p = ", signif(obj$joint$p.value, 3))) + base)
  }
  df <- as.data.frame(obj)
  if (what == "placebo") {
    df$true <- df$side == "true cutoff"
    return(ggplot2::ggplot(df, ggplot2::aes(x = .data$cutoff, y = .data$estimate)) +
             ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
             ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0) +
             ggplot2::geom_point(ggplot2::aes(colour = .data$true), size = 2.5) +
             ggplot2::scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
             ggplot2::labs(x = "Cutoff", y = "Estimate with robust 95% CI") + base)
  }
  if (what == "sensitivity") {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$h, y = .data$estimate)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), fill = "grey80", alpha = 0.6) +
      ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Bandwidth h", y = "Estimate with robust 95% CI") + base
    hm <- attr(obj, "h_mse")
    if (!is.null(hm)) g <- g + ggplot2::geom_vline(xintercept = hm, linetype = "dotted", colour = "grey40")
    return(g)
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$radius, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = "Donut radius (observations within it dropped)", y = "Estimate with robust 95% CI") + base
}
