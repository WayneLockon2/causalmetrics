# R/rd-weak-iv.R
#
# Weak-instrument-robust inference for fuzzy regression discontinuity
# designs: Anderson-Rubin confidence sets on the kernel-weighted local
# sample (Feir, Lemieux, and Marmer 2016).

#' Weak-identification-robust inference for fuzzy RD
#'
#' A fuzzy RD is a local instrumental-variables problem: the indicator of
#' being above the cutoff instruments the treatment, with the running
#' variable (and its interaction with the indicator) as controls, on the
#' observations inside the bandwidth weighted by the kernel. When the
#' first-stage jump is small the Wald interval of the ratio is unreliable
#' (Feir, Lemieux, and Marmer 2016). This function forms the local
#' two-stage least squares estimate, the local first-stage F statistics,
#' and the Anderson-Rubin confidence set, which stays valid however weak
#' the first stage is, by calling [iv_ar_confidence_set()] on the local
#' sample with kernel weights.
#'
#' @param data A data frame.
#' @param y,d,x Outcome, treatment, and running variable column names.
#' @param cutoff Cutoff value.
#' @param h Bandwidth; `NULL` uses the fuzzy MSE-optimal bandwidth of
#'   `rdrobust::rdbwselect()`.
#' @param kernel Kernel for the local weights.
#' @param covariates Optional covariate column names added as controls.
#' @param cluster Optional cluster column name.
#' @param conf_level Confidence level.
#' @param compare Also run `rdrobust::rdrobust(fuzzy = )` at the same
#'   bandwidth and report its rows.
#'
#' @return A list of class `cm_rd_weak_iv` with `ar_set` (the
#'   [iv_ar_confidence_set()] object), `first_stage` (jump, F, robust F,
#'   effective F), `table` (rows: local 2SLS with the Wald interval, the
#'   Anderson-Rubin set, and the `rdrobust` rows when `compare = TRUE`),
#'   `h`, `n_eff`.
#' @references
#' Feir, D., Lemieux, T., and Marmer, V. (2016). Weak identification in
#' fuzzy regression discontinuity designs. *Journal of Business and Economic
#' Statistics*, 34(2), 185-196.
#' @examples
#' dat <- sim_rd(3000, "fuzzy", compliance = 0.7, seed = 1)
#' rd_weak_iv(dat, "y", "d", "x")
#' @export
rd_weak_iv <- function(data, y, d, x, cutoff = 0, h = NULL, kernel = "triangular", covariates = NULL,
                       cluster = NULL, conf_level = 0.95, compare = TRUE) {
  pr <- .cm_rd_prepare(data, y, x, cutoff, d = d, covariates = covariates, cluster = cluster)
  if (is.null(h)) h <- .cm_rd_bandwidth(pr$y, pr$x, cutoff, d = pr$d, kernel = kernel)$h
  w <- .cm_rd_kernel((pr$x - cutoff) / h, kernel)
  keep <- w > 0
  loc <- pr$data[keep, , drop = FALSE]
  loc$.cm_y <- pr$y[keep]; loc$.cm_d <- pr$d[keep]
  loc$.cm_z <- as.integer(pr$x[keep] >= cutoff)
  loc$.cm_xc <- pr$x[keep] - cutoff
  loc$.cm_zxc <- loc$.cm_z * loc$.cm_xc
  loc$.cm_w <- w[keep]
  controls <- c(".cm_xc", ".cm_zxc", covariates)
  ar <- iv_ar_confidence_set(loc, ".cm_y", ".cm_d", ".cm_z", x = controls, cluster = cluster,
                             weights = ".cm_w", conf_level = conf_level)
  fs <- iv_first_stage(loc, ".cm_d", ".cm_z", x = controls, cluster = cluster, weights = ".cm_w")
  rows <- data.frame(
    term = c("fuzzy RD effect", "fuzzy RD effect"),
    method = c("local 2SLS (Wald)", "Anderson-Rubin"),
    estimate = c(ar$estimate, NA_real_), std.error = c(ar$std.error, NA_real_),
    conf.low = c(ar$wald[["lower"]], if (nrow(ar$intervals)) min(ar$intervals$lower) else NA_real_),
    conf.high = c(ar$wald[["upper"]], if (nrow(ar$intervals)) max(ar$intervals$upper) else NA_real_),
    p.value = c(2 * stats::pnorm(-abs(ar$estimate / ar$std.error)), NA_real_),
    h = h, n_eff = sum(keep), stringsAsFactors = FALSE)
  rob <- NULL
  if (compare && requireNamespace("rdrobust", quietly = TRUE)) {
    f <- .cm_rd_fit(pr$y, pr$x, cutoff, d = pr$d, cluster = pr$cluster, h = h, kernel = kernel, conf_level = conf_level)
    t <- .cm_rd_tidy_rdrobust(f, term = "fuzzy RD effect")
    rob <- data.frame(term = t$term, method = paste0("rdrobust ", t$method), estimate = t$estimate, std.error = t$std.error,
                      conf.low = t$conf.low, conf.high = t$conf.high, p.value = t$p.value, h = h, n_eff = sum(keep),
                      stringsAsFactors = FALSE)
    rows <- rbind(rows, rob)
  }
  jump <- if (!is.null(fs$coefficients) && ".cm_z" %in% names(fs$coefficients)) unname(fs$coefficients[[".cm_z"]]) else NA_real_
  structure(list(ar_set = ar, first_stage = list(jump = jump,
                                                  F = fs$F, F_robust = fs$F_robust, F_effective = fs$F_effective),
                 table = rows, h = h, n_eff = sum(keep), kernel = kernel, cutoff = cutoff, conf_level = conf_level,
                 call = match.call()), class = "cm_rd_weak_iv")
}

#' @export
print.cm_rd_weak_iv <- function(x, ...) {
  cat("Fuzzy RD with weak-identification-robust inference (h = ", format(round(x$h, 4)), ", n in window = ", x$n_eff, ")\n", sep = "")
  cat("  local first stage: F = ", format(round(x$first_stage$F, 2)), ", robust F = ", format(round(x$first_stage$F_robust, 2)),
      ", effective F = ", format(round(x$first_stage$F_effective, 2)), "\n", sep = "")
  cat("  Anderson-Rubin set: ", x$ar_set$type, "\n", sep = "")
  print(x$table[, c("method", "estimate", "std.error", "conf.low", "conf.high")], digits = 4, row.names = FALSE)
  invisible(x)
}
