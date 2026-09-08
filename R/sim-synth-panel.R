# R/sim-synth-panel.R

#' Simulate a panel from a linear factor model for synthetic control
#'
#' Potential outcomes follow the linear factor model of Ferman and Pinto
#' (2021): `y_jt = c_j + delta_t + lambda_t' mu_j + e_jt`, with unit effects
#' `c_j`, a common trend `delta_t`, common factors `lambda_t` (stationary
#' AR(1) processes), unit loadings `mu_j`, and idiosyncratic shocks. Options
#' place the treated unit in the tail of the loading or of the unit-effect
#' distribution (selection on unobservables), shift the first factor after
#' treatment, and stagger adoption across several treated units.
#'
#' @param n_donors,n_treated Numbers of donor and treated units.
#' @param t_pre,t_post Numbers of pre- and post-treatment periods for a unit
#'   adopting at `t_pre + 1`.
#' @param n_factors Number of common factors.
#' @param sd_e Standard deviation of the idiosyncratic shock; small values
#'   give a near-perfect pre-treatment fit.
#' @param sd_factor,rho Standard deviation and AR(1) coefficient of the
#'   factors.
#' @param trend Slope of the common trend `delta_t = trend * t`.
#' @param select_loading Loading of the treated unit(s) on the first factor,
#'   in standard deviations of the donor loadings (0 keeps it random).
#' @param select_level Unit effect of the treated unit(s), in standard
#'   deviations of the donor effects (0 keeps it random).
#' @param factor_shift Mean shift of the first factor in the post-treatment
#'   periods (a "break", Ferman and Pinto Table 1 Panels C-D); with a
#'   nonzero shift and selection on the loading, DiD and the demeaned
#'   synthetic control are both biased.
#' @param effect Treatment effect: a scalar or a vector over event times
#'   `0, 1, ...` (recycled).
#' @param adoption Optional vector of adoption periods for the treated units
#'   (staggered adoption); default all adopt at `t_pre + 1`.
#' @param seed Optional seed.
#' @return A data frame with `id`, `time`, `y`, `d`, `y0` (untreated
#'   potential outcome), and `tau` (the effect in treated cells), with
#'   attributes `factors`, `loadings`, and `unit_effects`.
#' @examples
#' dat <- sim_synth_panel(n_donors = 20, t_pre = 20, t_post = 5, select_loading = 1.5,
#'                        factor_shift = 1, effect = 2, seed = 1)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
#' fit$estimate
#' @export
sim_synth_panel <- function(n_donors = 20, n_treated = 1, t_pre = 20, t_post = 5, n_factors = 2,
                            sd_e = 1, sd_factor = 1, rho = 0.5, trend = 0.2, select_loading = 0,
                            select_level = 0, factor_shift = 0, effect = 0, adoption = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- n_donors + n_treated
  Tn <- t_pre + t_post
  if (is.null(adoption)) adoption <- rep(t_pre + 1L, n_treated)
  if (length(adoption) != n_treated) stop("`adoption` must have one entry per treated unit.", call. = FALSE)
  lam <- matrix(0, Tn, n_factors)
  for (k in seq_len(n_factors)) {
    lam[1L, k] <- stats::rnorm(1, sd = sd_factor)
    for (t in seq_len(Tn)[-1L]) lam[t, k] <- rho * lam[t - 1L, k] + stats::rnorm(1, sd = sd_factor * sqrt(1 - rho^2))
  }
  lam[(t_pre + 1L):Tn, 1L] <- lam[(t_pre + 1L):Tn, 1L] + factor_shift * sd_factor
  mu <- matrix(stats::rnorm(N * n_factors), N, n_factors)
  cj <- stats::rnorm(N)
  treated_idx <- n_donors + seq_len(n_treated)
  if (select_loading != 0) mu[treated_idx, 1L] <- select_loading
  if (select_level != 0) cj[treated_idx] <- select_level
  delta <- trend * seq_len(Tn)
  y0 <- outer(cj, rep(1, Tn)) + outer(rep(1, N), delta) + mu %*% t(lam) + matrix(stats::rnorm(N * Tn, sd = sd_e), N, Tn)
  d <- matrix(0L, N, Tn); tau <- matrix(0, N, Tn)
  for (i in seq_len(n_treated)) {
    g <- adoption[i]
    post <- seq_len(Tn) >= g
    d[treated_idx[i], post] <- 1L
    ev <- seq_len(Tn)[post] - g
    tau[treated_idx[i], post] <- rep_len(effect, length(ev))[seq_along(ev)]
  }
  y <- y0 + tau
  out <- data.frame(id = rep(seq_len(N), each = Tn), time = rep(seq_len(Tn), N),
                    y = as.vector(t(y)), d = as.vector(t(d)), y0 = as.vector(t(y0)), tau = as.vector(t(tau)))
  attr(out, "factors") <- lam
  attr(out, "loadings") <- mu
  attr(out, "unit_effects") <- cj
  out
}
