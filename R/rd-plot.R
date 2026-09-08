# R/rd-plot.R
#
# Regression discontinuity plots as ggplot objects: binned means on each
# side of the cutoff with polynomial or local linear fits.

#' Binned means of an outcome on each side of a cutoff
#'
#' @param data A data frame.
#' @param y,x Outcome and running variable column names.
#' @param cutoff Cutoff value.
#' @param bins `"qs"` (quantile-spaced, default) or `"es"` (evenly spaced).
#' @param n_bins Number of bins per side: one number, a vector `c(left,
#'   right)`, or `NULL` to use the integrated-mean-squared-error optimal
#'   number of `rdrobust::rdplot()` (mimicking-variance variant) when that
#'   package is installed, and 20 per side otherwise.
#' @param weights Optional weights column name.
#' @return A data frame with `side`, `bin`, `bin_x`, `bin_y`, `std.error`,
#'   `n`, `x_min`, `x_max`, and attribute `"n_bins"`.
#' @export
rd_bins <- function(data, y, x, cutoff = 0, bins = c("qs", "es"), n_bins = NULL, weights = NULL) {
  bins <- match.arg(bins)
  pr <- .cm_rd_prepare(data, y, x, cutoff, weights = weights)
  J <- .cm_rd_n_bins(pr$y, pr$x, cutoff, bins, n_bins)
  out <- .cm_rd_bins(pr$y, pr$x, cutoff, J[1], J[2], bins = bins, weights = pr$weights)
  attr(out, "n_bins") <- J
  out
}

.cm_rd_n_bins <- function(y, x, cutoff, bins, n_bins) {
  if (!is.null(n_bins)) {
    n_bins <- as.integer(n_bins)
    return(if (length(n_bins) == 1L) c(n_bins, n_bins) else n_bins[1:2])
  }
  if (requireNamespace("rdrobust", quietly = TRUE)) {
    J <- tryCatch({
      rp <- suppressWarnings(rdrobust::rdplot(y = y, x = x, c = cutoff, hide = TRUE,
                                              binselect = paste0(bins, "mv"), masspoints = "off"))
      as.integer(rp$J)
    }, error = function(e) NULL)
    if (!is.null(J) && all(is.finite(J)) && all(J > 0)) return(J)
  }
  c(20L, 20L)
}

#' Regression discontinuity plot
#'
#' Binned means of the outcome on each side of the cutoff with a fitted
#' curve per side, as a ggplot object that can be modified with the usual
#' `+` layers. Bin counts default to the integrated-mean-squared-error
#' optimal choice of Calonico, Cattaneo, and Titiunik (2015) through
#' `rdrobust::rdplot()` when installed. The fitted curve is a global
#' polynomial of order `p` on each side (the `rdplot` convention) or a
#' local linear fit within bandwidth `h` (the estimator's own view of the
#' data).
#'
#' @inheritParams rd_bins
#' @param p Polynomial order of the global fit (default 4).
#' @param fit `"polynomial"` (default), `"local"` (local linear within `h`),
#'   or `"none"`.
#' @param h Bandwidth for `fit = "local"`; `NULL` uses the MSE-optimal
#'   bandwidth from `rdrobust::rdbwselect()`.
#' @param kernel Kernel for the local fit.
#' @param ci Draw pointwise 95 percent bars around the bin means.
#' @param x_range Optional range of `x` to display.
#' @param colours Two colours for the left and right sides.
#' @param point_size Size of the bin-mean points.
#' @param x_lab,y_lab Axis labels.
#' @return A ggplot object with attributes `"rd_bins"` (the binned data)
#'   and `"rd_fit"` (the fitted curves).
#' @references
#' Calonico, S., Cattaneo, M. D., and Titiunik, R. (2015). Optimal
#' data-driven regression discontinuity plots. *Journal of the American
#' Statistical Association*, 110(512), 1753-1769.
#' @examples
#' dat <- sim_rd(2000, "lee", seed = 1)
#' rd_plot(dat, "y", "x")
#' rd_plot(dat, "y", "x", fit = "local", ci = TRUE) + ggplot2::labs(title = "Local linear fit")
#' @export
rd_plot <- function(data, y, x, cutoff = 0, bins = c("qs", "es"), n_bins = NULL, p = 4,
                    fit = c("polynomial", "local", "none"), h = NULL, kernel = "triangular",
                    ci = FALSE, weights = NULL, x_range = NULL,
                    colours = c("#0072B2", "#D55E00"), point_size = 1.6,
                    x_lab = x, y_lab = y) {
  bins <- match.arg(bins)
  fit <- match.arg(fit)
  pr <- .cm_rd_prepare(data, y, x, cutoff, weights = weights)
  yv <- pr$y; xv <- pr$x
  if (!is.null(x_range)) {
    keep <- xv >= x_range[1] & xv <= x_range[2]
    yv <- yv[keep]; xv <- xv[keep]
    if (!is.null(pr$weights)) pr$weights <- pr$weights[keep]
  }
  J <- .cm_rd_n_bins(yv, xv, cutoff, bins, n_bins)
  bdat <- .cm_rd_bins(yv, xv, cutoff, J[1], J[2], bins = bins, weights = pr$weights)
  bdat$side <- factor(bdat$side, levels = c("left", "right"))

  fdat <- NULL
  if (fit != "none") {
    if (fit == "local" && is.null(h)) {
      h <- .cm_rd_bandwidth(yv, xv, cutoff, kernel = kernel)$h
    }
    fdat <- do.call(rbind, lapply(c("left", "right"), function(s) {
      idx <- if (s == "left") xv < cutoff else xv >= cutoff
      if (fit == "local") idx <- idx & abs(xv - cutoff) <= h
      xs <- xv[idx]; ys <- yv[idx]
      if (length(xs) < p + 2L) return(NULL)
      grid <- seq(min(xs), max(xs), length.out = 100L)
      if (s == "left") grid[length(grid)] <- cutoff else grid[1] <- cutoff
      if (fit == "polynomial") {
        m <- stats::lm(ys ~ stats::poly(xs, degree = p, raw = TRUE))
        pred <- stats::predict(m, newdata = data.frame(xs = grid))
      } else {
        w <- .cm_rd_kernel((xs - cutoff) / h, kernel)
        m <- stats::lm(ys ~ xs, weights = w)
        pred <- stats::predict(m, newdata = data.frame(xs = grid))
      }
      data.frame(side = s, x_grid = grid, fit = as.numeric(pred))
    }))
    fdat$side <- factor(fdat$side, levels = c("left", "right"))
  }

  g <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey40")
  if (ci) {
    g <- g + ggplot2::geom_errorbar(
      data = bdat, ggplot2::aes(x = .data$bin_x, ymin = .data$bin_y - 1.96 * .data$std.error,
                                ymax = .data$bin_y + 1.96 * .data$std.error, colour = .data$side),
      width = 0, alpha = 0.5)
  }
  g <- g + ggplot2::geom_point(data = bdat, ggplot2::aes(x = .data$bin_x, y = .data$bin_y, colour = .data$side),
                               size = point_size)
  if (!is.null(fdat)) {
    g <- g + ggplot2::geom_line(data = fdat, ggplot2::aes(x = .data$x_grid, y = .data$fit, colour = .data$side),
                                linewidth = 0.9)
  }
  g <- g +
    ggplot2::scale_colour_manual(values = c(left = colours[1], right = colours[2]), guide = "none") +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme_minimal(base_size = 11)
  attr(g, "rd_bins") <- bdat
  attr(g, "rd_fit") <- fdat
  attr(g, "n_bins") <- J
  attr(g, "h") <- h
  g
}
