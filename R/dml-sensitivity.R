# R/dml-sensitivity.R
#
# Sensitivity of a partially linear DML estimate to an unobserved confounder:
# the omitted-confounder bias bound in terms of partial R-squareds
# (Chernozhukov et al. 2026, Theorem 12.2.1; Cinelli and Hazlett 2020).

#' Sensitivity of a DML estimate to unobserved confounding
#'
#' In the partially linear model with an unobserved confounder `A`, the
#' bias of the partialling-out estimate is bounded by
#' `|bias| <= sqrt(R2_y R2_d / (1 - R2_d)) * sqrt(E[(Y~ - theta D~)^2] / E[D~^2])`,
#' where `R2_y` is the partial R-squared of `A` with the residualized outcome
#' (after the residualized treatment) and `R2_d` the partial R-squared of
#' `A` with the residualized treatment. The residual moments are available
#' from the `cm_dml` fit. The function returns the bound and the implied
#' interval for `theta` at user-supplied `(R2_y, R2_d)`, the robustness value
#' (the common partial R-squared at which the bound reaches the estimate, or
#' the edge of its confidence interval), and a contour grid for plotting.
#'
#' @param fit A `cm_dml` object with `model = "plr"` or `"pliv"`.
#' @param r2_y,r2_d Partial R-squareds of the confounder with the outcome
#'   and the treatment for the benchmark scenario.
#' @param conf_level Confidence level for the robustness value with
#'   sampling uncertainty.
#' @param grid Values used for the contour grid.
#' @return A list of class `cm_sensitivity` with `bias`, `bounds`,
#'   `robustness_value` (point and with the confidence interval), the
#'   `scale` factor, and `contour` (data frame of `r2_d`, `r2_y`, `bias`).
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., and Syrgkanis, V.
#' (2022). Long story short: omitted variable bias in causal machine
#' learning. NBER Working Paper 30302.
#'
#' Cinelli, C. and Hazlett, C. (2020). Making sense of sensitivity: extending
#' omitted variable bias. *Journal of the Royal Statistical Society B*,
#' 82(1), 39-67.
#' @examples
#' set.seed(1)
#' n <- 800; x <- rnorm(n); d <- 0.5 * x + rnorm(n); y <- d + x^2 + rnorm(n)
#' fit <- est_dml(data.frame(y, d, x), "y", "d", "x", seed = 1)
#' dml_sensitivity(fit, r2_y = 0.05, r2_d = 0.02)
#' @export
dml_sensitivity <- function(fit, r2_y = 0.05, r2_d = 0.05, conf_level = 0.95, grid = seq(0, 0.5, by = 0.01)) {
  if (!inherits(fit, "cm_dml") || !fit$model %in% c("plr", "pliv")) stop("`fit` must be an est_dml() object with model \"plr\" or \"pliv\".", call. = FALSE)
  res <- fit$residuals
  theta <- fit$estimate
  e2 <- mean((res$y_tilde - theta * res$d_tilde)^2)
  d2 <- mean(res$d_tilde^2)
  scale <- sqrt(e2 / d2)
  bias_fun <- function(ry, rd) sqrt(ry * rd / (1 - rd)) * scale
  bias <- bias_fun(r2_y, r2_d)
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  rv <- function(target) {
    q <- (target / scale)^2
    (-q + sqrt(q^2 + 4 * q)) / 2
  }
  rv_point <- rv(abs(theta))
  rv_ci <- if (abs(theta) > zq * fit$std.error) rv(abs(theta) - zq * fit$std.error) else 0
  contour <- expand.grid(r2_d = grid[grid < 1], r2_y = grid)
  contour$bias <- bias_fun(contour$r2_y, contour$r2_d)
  structure(list(estimate = theta, std.error = fit$std.error, bias = bias, r2_y = r2_y, r2_d = r2_d,
                 bounds = c(lower = theta - bias, upper = theta + bias),
                 bounds_ci = c(lower = theta - bias - zq * fit$std.error, upper = theta + bias + zq * fit$std.error),
                 robustness_value = c(point = rv_point, with_ci = rv_ci), scale = scale,
                 residual_variance = e2, treatment_residual_variance = d2, contour = contour,
                 conf_level = conf_level), class = "cm_sensitivity")
}

#' @export
print.cm_sensitivity <- function(x, ...) {
  cat("Sensitivity of the DML estimate to an unobserved confounder\n")
  cat("  estimate = ", format(round(x$estimate, 4)), " (SE ", format(round(x$std.error, 4)), ")\n", sep = "")
  cat("  benchmark R2_y = ", x$r2_y, ", R2_d = ", x$r2_d, ": |bias| <= ", format(round(x$bias, 4)),
      ", bounds [", format(round(x$bounds[1], 4)), ", ", format(round(x$bounds[2], 4)), "]\n", sep = "")
  cat("  robustness value (equal R2 that moves the estimate to zero): ", format(round(x$robustness_value[1], 3)),
      "; with the ", 100 * x$conf_level, "% interval: ", format(round(x$robustness_value[2], 3)), "\n", sep = "")
  invisible(x)
}

#' Contour plot of the omitted-confounder bias bound
#'
#' @param x A `cm_sensitivity` object.
#' @param benchmarks Optional data frame with columns `r2_d`, `r2_y`, `label`
#'   to mark benchmark scenarios.
#' @return A ggplot object with contours of the bias bound; the contour at
#'   the estimate's absolute value marks where the sign is no longer known.
#' @export
plot_dml_sensitivity <- function(x, benchmarks = NULL) {
  df <- x$contour
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$r2_d, y = .data$r2_y, z = .data$bias)) +
    ggplot2::geom_contour(colour = "grey50", bins = 12) +
    ggplot2::geom_contour(breaks = abs(x$estimate), colour = "firebrick", linewidth = 1) +
    ggplot2::labs(x = "Partial R2 of the confounder with the treatment",
                  y = "Partial R2 of the confounder with the outcome",
                  subtitle = "Red contour: bias equal to the estimate (sign no longer identified)") +
    ggplot2::theme_minimal(base_size = 11)
  if (!is.null(benchmarks)) {
    p <- p + ggplot2::geom_point(data = benchmarks, ggplot2::aes(x = .data$r2_d, y = .data$r2_y), inherit.aes = FALSE) +
      ggplot2::geom_text(data = benchmarks, ggplot2::aes(x = .data$r2_d, y = .data$r2_y, label = .data$label), inherit.aes = FALSE, vjust = -0.7, size = 3)
  }
  p
}
