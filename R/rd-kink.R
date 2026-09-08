# R/rd-kink.R
#
# Regression kink designs: slope changes at a kink point, the ratio of the
# outcome kink to the treatment kink, elasticities, and the four standard
# graphs.

#' Regression kink design
#'
#' Estimates the change in the slope of the outcome at the kink (the reduced
#' form) and, when a treatment variable is given, the change in the slope of
#' the treatment (the first stage) and their ratio (the fuzzy kink
#' estimator of Card, Lee, Pei, and Weber 2015), all through
#' `rdrobust::rdrobust(deriv = 1)`. The elasticity of the outcome with
#' respect to the treatment is the ratio scaled by the levels of the
#' treatment and the outcome at the kink.
#'
#' @param data A data frame.
#' @param y,x Outcome and running variable column names.
#' @param cutoff Kink point.
#' @param d Optional treatment column name (fuzzy kink); when `NULL` the
#'   reduced-form slope change is reported and `slope_change` (a known
#'   change in the treatment's slope, for sharp kinks) scales it.
#' @param slope_change Known change in the slope of the treatment at the
#'   kink for sharp designs (default 1: the reduced form is reported).
#' @param elasticity Report the elasticity `kink * D(c) / Y(c)`.
#' @param covariates Optional covariates passed to `rdrobust`.
#' @param p Local polynomial order (default 1; Card et al. use local
#'   linear and local quadratic).
#' @param ... Passed to `rdrobust::rdrobust()` (`h`, `kernel`, `cluster`,
#'   `bwselect`, ...).
#'
#' @return A list of class `cm_rd_kink` with `reduced_form`, `first_stage`,
#'   `kink` (tidy rows of the ratio or the scaled reduced form),
#'   `elasticity`, `levels` (outcome and treatment at the kink), and the
#'   `rdrobust` objects in `fits`.
#' @references
#' Card, D., Lee, D. S., Pei, Z., and Weber, A. (2015). Inference on causal
#' effects in a generalized regression kink design. *Econometrica*, 83(6),
#' 2453-2483.
#' @examples
#' dat <- sim_rd(3000, "kink", seed = 1)
#' rd_kink(dat, "y", "x", d = "d")
#' @export
rd_kink <- function(data, y, x, cutoff = 0, d = NULL, slope_change = 1, elasticity = FALSE,
                    covariates = NULL, p = 1, ...) {
  .cm_check_package("rdrobust")
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, covariates = covariates)
  dots <- list(...)
  rf <- do.call(.cm_rd_fit, c(list(y = pr$y, x = pr$x, cutoff = cutoff, covs = pr$covs, p = p, deriv = 1), dots))
  reduced <- .cm_rd_tidy_rdrobust(rf, term = "slope change of y")
  fits <- list(reduced_form = rf)
  first <- NULL; kink <- NULL
  if (!is.null(d)) {
    fs <- do.call(.cm_rd_fit, c(list(y = pr$d, x = pr$x, cutoff = cutoff, covs = pr$covs, p = p, deriv = 1), dots))
    first <- .cm_rd_tidy_rdrobust(fs, term = "slope change of d")
    fz <- do.call(.cm_rd_fit, c(list(y = pr$y, x = pr$x, cutoff = cutoff, d = pr$d, covs = pr$covs, p = p, deriv = 1), dots))
    kink <- .cm_rd_tidy_rdrobust(fz, term = "kink ratio")
    fits$first_stage <- fs; fits$fuzzy <- fz
  } else {
    kink <- reduced
    kink$term <- "kink (reduced form / slope_change)"
    for (v in c("estimate", "std.error", "conf.low", "conf.high")) kink[[v]] <- kink[[v]] / slope_change
  }
  # levels at the kink from local linear intercepts
  h <- rf$bws[1, 1]
  lev_y <- .cm_rd_local_jump(pr$y, pr$x, cutoff, h)$levels
  y_c <- mean(lev_y)
  d_c <- if (!is.null(d)) mean(.cm_rd_local_jump(pr$d, pr$x, cutoff, h)$levels) else NA_real_
  elas <- NULL
  if (elasticity) {
    if (is.null(d)) stop("`elasticity = TRUE` needs a treatment column `d`.", call. = FALSE)
    scale <- d_c / y_c
    elas <- kink
    elas$term <- "elasticity"
    for (v in c("estimate", "std.error", "conf.low", "conf.high")) elas[[v]] <- elas[[v]] * scale
  }
  structure(list(reduced_form = reduced, first_stage = first, kink = kink, elasticity = elas,
                 levels = c(y = y_c, d = d_c), fits = fits, cutoff = cutoff, y = y, x = x, d = d,
                 data = pr$data, call = match.call()), class = "cm_rd_kink")
}

#' @export
print.cm_rd_kink <- function(x, ...) {
  cat("Regression kink design at ", x$cutoff, "\n", sep = "")
  show <- function(t, lab) {
    cat("\n", lab, ":\n", sep = "")
    print(t[, c("method", "estimate", "std.error", "conf.low", "conf.high", "h_left", "n_left", "n_right")], digits = 4, row.names = FALSE)
  }
  show(x$reduced_form, "Reduced form (change in the slope of y)")
  if (!is.null(x$first_stage)) show(x$first_stage, "First stage (change in the slope of d)")
  show(x$kink, "Kink estimate")
  if (!is.null(x$elasticity)) show(x$elasticity, "Elasticity")
  invisible(x)
}

#' The four regression kink graphs
#'
#' First stage, reduced form, density of the running variable, and a
#' covariate index, each against the running variable with binned means.
#'
#' @param x A `cm_rd_kink` object.
#' @param covariates Optional covariates for the index panel (fitted values
#'   of a regression of the outcome on the covariates).
#' @param n_bins Bins per side for the binned panels.
#' @param p Polynomial order of the fitted lines.
#' @return A patchwork object when `patchwork` is installed, otherwise a
#'   list of ggplot objects.
#' @export
plot_rd_kink <- function(x, covariates = NULL, n_bins = 20, p = 1) {
  dat <- x$data
  panels <- list()
  if (!is.null(x$d)) {
    panels$first_stage <- rd_plot(dat, x$d, x$x, x$cutoff, n_bins = n_bins, p = p, y_lab = x$d) + ggplot2::labs(title = "First stage")
  }
  panels$reduced_form <- rd_plot(dat, x$y, x$x, x$cutoff, n_bins = n_bins, p = p, y_lab = x$y) + ggplot2::labs(title = "Reduced form")
  panels$density <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x$x]])) +
    ggplot2::geom_histogram(bins = 60, fill = "grey70", colour = "white") +
    ggplot2::geom_vline(xintercept = x$cutoff, linetype = "dashed") +
    ggplot2::labs(title = "Density of the running variable", x = x$x, y = "Count") +
    ggplot2::theme_minimal(base_size = 11)
  if (!is.null(covariates)) {
    m <- stats::lm(stats::reformulate(covariates, x$y), data = dat)
    dat$.cm_index <- stats::fitted(m)
    panels$covariate_index <- rd_plot(dat, ".cm_index", x$x, x$cutoff, n_bins = n_bins, p = p, y_lab = "Covariate index") +
      ggplot2::labs(title = "Covariate index")
  }
  if (requireNamespace("patchwork", quietly = TRUE)) return(patchwork::wrap_plots(panels, ncol = 2))
  panels
}
