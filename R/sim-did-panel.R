#' Simulate a staggered-adoption panel
#'
#' A data-generating process for teaching and testing difference-in-differences
#' estimators: units are assigned to treatment cohorts (first treatment period
#' `g`) or never treated, potential outcomes have unit and time effects,
#' optional covariate-specific trends, and treatment effects that grow with
#' exposure and differ across cohorts, as in the stylized example of
#' Sant'Anna's lectures.
#'
#' @param n_units Number of units.
#' @param n_periods Number of periods `1, ..., n_periods`.
#' @param groups First treatment periods of the cohorts.
#' @param never_share Share of units never treated.
#' @param mu Cohort-specific effect sizes, recycled over `groups`. The
#'   effect at event time `e >= 0` is `mu[g] * (e + 1)` by default.
#' @param effect Optional function `effect(e, g)` returning the treatment
#'   effect at event time `e` for cohort `g`; overrides `mu`.
#' @param x_trend Slope of the covariate-specific trend: untreated outcomes
#'   grow by `x_trend * x1` per period, so parallel trends holds only
#'   conditional on `x1` when `x_trend != 0`.
#' @param x_select Strength of selection into earlier cohorts on `x1`.
#' @param anticipation Number of periods before `g` in which half of the
#'   eventual first-period effect already appears.
#' @param sd_e Standard deviation of the idiosyncratic error.
#' @param seed Optional seed.
#'
#' @return A data frame in long format with `id`, `time`, `g` (0 for never
#'   treated), `treated` (1 once treated), `event_time`, `x1`, `x2`, `y`, and
#'   `tau` (the realized treatment effect, 0 when untreated). The attribute
#'   `"att_gt"` holds the population `ATT(g, t)` implied by the design.
#'
#' @examples
#' dat <- sim_did_panel(n_units = 200, n_periods = 6, groups = c(3, 5), seed = 1)
#' head(dat)
#' attr(dat, "att_gt")
#' @export
sim_did_panel <- function(n_units = 500, n_periods = 10, groups = c(4, 7),
                          never_share = 0.4, mu = c(3, 2, 1), effect = NULL,
                          x_trend = 0, x_select = 0.5, anticipation = 0, sd_e = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  groups <- sort(unique(as.numeric(groups)))
  if (any(groups <= 1) || any(groups > n_periods)) stop("`groups` must lie in 2..n_periods.", call. = FALSE)
  mu <- rep_len(mu, length(groups))
  if (is.null(effect)) {
    effect <- function(e, g) mu[match(g, groups)] * (e + 1)
  }
  x1 <- stats::rnorm(n_units)
  x2 <- stats::rbinom(n_units, 1, 0.5)
  # cohort assignment: earlier cohorts more likely for high x1
  score <- x_select * x1
  probs <- cbind(never_share, sapply(seq_along(groups), function(k) {
    (1 - never_share) * exp(score * (length(groups) - k + 1) / length(groups))
  }))
  probs[, -1] <- probs[, -1] / rowSums(probs[, -1, drop = FALSE]) * (1 - never_share)
  g <- vapply(seq_len(n_units), function(i) {
    k <- sample.int(length(groups) + 1L, 1L, prob = probs[i, ])
    if (k == 1L) 0 else groups[k - 1L]
  }, numeric(1))
  alpha <- stats::rnorm(n_units)
  lambda <- 0.3 * seq_len(n_periods)
  out <- data.frame(
    id = rep(seq_len(n_units), each = n_periods),
    time = rep(seq_len(n_periods), times = n_units),
    g = rep(g, each = n_periods),
    x1 = rep(x1, each = n_periods),
    x2 = rep(x2, each = n_periods)
  )
  out$event_time <- ifelse(out$g > 0, out$time - out$g, NA_real_)
  out$treated <- as.integer(out$g > 0 & out$time >= out$g)
  tau <- rep(0, nrow(out))
  post <- which(out$treated == 1L)
  tau[post] <- effect(out$event_time[post], out$g[post])
  if (anticipation > 0) {
    ant <- which(out$g > 0 & out$event_time < 0 & out$event_time >= -anticipation)
    tau[ant] <- 0.5 * effect(0, out$g[ant])
  }
  y0 <- rep(alpha, each = n_periods) + lambda[out$time] + x_trend * out$x1 * out$time +
    0.5 * out$x2 + stats::rnorm(nrow(out), sd = sd_e)
  out$tau <- tau
  out$y <- y0 + tau
  cells <- expand.grid(group = groups, time = seq_len(n_periods))
  cells <- cells[cells$time >= cells$group, ]
  cells$att <- effect(cells$time - cells$group, cells$group)
  attr(out, "att_gt") <- cells[order(cells$group, cells$time), ]
  out
}
