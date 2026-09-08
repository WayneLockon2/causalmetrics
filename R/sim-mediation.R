# R/sim-mediation.R
#
# Simulated designs with known direct, indirect, and controlled direct effects,
# for the mediation lecture, the usage vignette, and the tests.

#' Simulate data for mediation analysis
#'
#' Designs with a binary treatment `d`, one or more mediators, and an outcome,
#' with the true total, natural direct, natural indirect, and controlled direct
#' effects attached as `attr(, "truth")` (finite-sample averages of the
#' simulated potential outcomes).
#'
#' @param n Sample size.
#' @param dgp Design:
#'   * `"linear"`: two covariates, confounded treatment, linear mediator and
#'     outcome models without treatment-mediator interaction; product of
#'     coefficients is exact.
#'   * `"interaction"`: as `"linear"` with a treatment-mediator interaction, so
#'     pure and total direct effects differ.
#'   * `"nonlinear"`: nonlinear covariate effects in the treatment, mediator,
#'     and outcome models (the case for machine-learning nuisances).
#'   * `"binary"`: binary mediator and binary outcome (logit models).
#'   * `"correlated_errors"`: `"linear"` with correlated mediator and outcome
#'     errors of correlation `rho` (an unmeasured mediator-outcome confounder);
#'     sequential ignorability fails and the truth records `rho`.
#'   * `"post_treatment_confounder"`: a covariate `z` affected by the treatment
#'     confounds the mediator and the outcome; natural effects are not
#'     identified, the controlled direct effect is.
#'   * `"iv_mediator"`: an unobserved `u` drives both the mediator and the
#'     outcome; an instrument `z` shifts the mediator; treatment randomized.
#'   * `"front_door"`: an unobserved `u` confounds treatment and outcome; the
#'     mediator is unconfounded given the treatment and carries the whole
#'     effect.
#'   * `"parallel"`: two mediators `m1`, `m2` with different treatment effects.
#' @param rho Error correlation for `"correlated_errors"`.
#' @param seed Optional seed.
#'
#' @return A data frame with `y`, `d`, mediator column(s) `m` (or `m1`, `m2`),
#'   covariates `x1`, `x2` (plus `z` or `u` where the design has them), and the
#'   attribute `"truth"`, a list with `total`, `nde`, `nie`, `nde_total`,
#'   `nie_pure`, `cde` (at the mediator reference stored in `m_ref`), and
#'   design-specific entries.
#' @examples
#' dat <- sim_mediation(500, dgp = "linear", seed = 1)
#' attr(dat, "truth")[c("total", "nde", "nie")]
#' @export
sim_mediation <- function(n = 1000L, dgp = c("linear", "interaction", "nonlinear", "binary",
                                             "correlated_errors", "post_treatment_confounder",
                                             "iv_mediator", "front_door", "parallel"),
                          rho = 0.5, seed = NULL) {
  dgp <- match.arg(dgp)
  n <- .cm_check_count(n, "n", min = 20L)
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  out <- switch(dgp,
    linear = , interaction = , correlated_errors = {
      b3 <- if (dgp == "interaction") 0.8 else 0
      p <- stats::plogis(0.4 * x1 - 0.3 * x2)
      d <- stats::rbinom(n, 1, p)
      if (dgp == "correlated_errors") {
        e <- MASS_mvrnorm(n, rho)
        e_m <- e[, 1]; e_y <- e[, 2]
      } else {
        e_m <- stats::rnorm(n); e_y <- stats::rnorm(n)
      }
      m_fun <- function(dd) 0.5 + 1.2 * dd + 0.6 * x1 - 0.4 * x2 + e_m
      y_fun <- function(dd, mm) 1 + 1.0 * dd + 0.9 * mm + b3 * dd * mm + 0.5 * x1 + 0.3 * x2 + e_y
      m0 <- m_fun(0); m1 <- m_fun(1)
      m <- ifelse(d == 1, m1, m0)
      y <- y_fun(d, m)
      m_ref <- 0
      truth <- list(total = mean(y_fun(1, m1) - y_fun(0, m0)),
                    nde = mean(y_fun(1, m0) - y_fun(0, m0)), nie = mean(y_fun(1, m1) - y_fun(1, m0)),
                    nde_total = mean(y_fun(1, m1) - y_fun(0, m1)), nie_pure = mean(y_fun(0, m1) - y_fun(0, m0)),
                    cde = mean(y_fun(1, m_ref) - y_fun(0, m_ref)), m_ref = m_ref,
                    coef = list(a1 = 1.2, b1 = 1.0, b2 = 0.9, b3 = b3), rho = if (dgp == "correlated_errors") rho else 0)
      list(data = data.frame(y = y, d = d, m = m, x1 = x1, x2 = x2), truth = truth)
    },
    nonlinear = {
      p <- stats::plogis(0.8 * sin(x1) - 0.5 * x2^2 + 0.3)
      d <- stats::rbinom(n, 1, p)
      e_m <- stats::rnorm(n); e_y <- stats::rnorm(n)
      m_nl <- function(dd) 0.5 + 1.2 * dd + exp(0.5 * x1) - 0.8 * x2^2 + e_m
      y_nl <- function(dd, mm) 1 + 1.0 * dd + 0.9 * mm + 0.4 * dd * mm + 1.5 * cos(x1) + x1 * x2 + e_y
      m0 <- m_nl(0); m1 <- m_nl(1)
      m <- ifelse(d == 1, m1, m0)
      y <- y_nl(d, m)
      m_ref <- 0
      truth <- list(total = mean(y_nl(1, m1) - y_nl(0, m0)),
                    nde = mean(y_nl(1, m0) - y_nl(0, m0)), nie = mean(y_nl(1, m1) - y_nl(1, m0)),
                    nde_total = mean(y_nl(1, m1) - y_nl(0, m1)), nie_pure = mean(y_nl(0, m1) - y_nl(0, m0)),
                    cde = mean(y_nl(1, m_ref) - y_nl(0, m_ref)), m_ref = m_ref)
      list(data = data.frame(y = y, d = d, m = m, x1 = x1, x2 = x2), truth = truth)
    },
    binary = {
      p <- stats::plogis(0.4 * x1 - 0.3 * x2)
      d <- stats::rbinom(n, 1, p)
      u_m <- stats::runif(n); u_y <- stats::runif(n)
      m_bin <- function(dd) as.numeric(u_m < stats::plogis(-0.5 + 1.5 * dd + 0.5 * x1))
      y_bin <- function(dd, mm) as.numeric(u_y < stats::plogis(-1 + 0.8 * dd + 1.2 * mm + 0.4 * x1 - 0.3 * x2))
      m0 <- m_bin(0); m1 <- m_bin(1)
      m <- ifelse(d == 1, m1, m0)
      y <- y_bin(d, m)
      m_ref <- 0
      truth <- list(total = mean(y_bin(1, m1) - y_bin(0, m0)),
                    nde = mean(y_bin(1, m0) - y_bin(0, m0)), nie = mean(y_bin(1, m1) - y_bin(1, m0)),
                    nde_total = mean(y_bin(1, m1) - y_bin(0, m1)), nie_pure = mean(y_bin(0, m1) - y_bin(0, m0)),
                    cde = mean(y_bin(1, m_ref) - y_bin(0, m_ref)), m_ref = m_ref)
      list(data = data.frame(y = y, d = d, m = m, x1 = x1, x2 = x2), truth = truth)
    },
    post_treatment_confounder = {
      d <- stats::rbinom(n, 1, 0.5)
      e_z <- stats::rnorm(n); e_m <- stats::rnorm(n); e_y <- stats::rnorm(n)
      z_fun <- function(dd) 0.7 * dd + 0.4 * x1 + e_z
      m_ptc <- function(dd, zz) 0.5 + 1.0 * dd + 0.8 * zz + 0.3 * x1 + e_m
      y_ptc <- function(dd, mm, zz) 1 + 0.6 * dd + 0.9 * mm + 0.7 * zz + 0.5 * x1 + 0.3 * x2 + e_y
      z0 <- z_fun(0); z1 <- z_fun(1)
      z <- ifelse(d == 1, z1, z0)
      m0 <- m_ptc(0, z0); m1 <- m_ptc(1, z1)
      m <- ifelse(d == 1, m1, m0)
      y <- y_ptc(d, m, z)
      m_ref <- 0
      truth <- list(total = mean(y_ptc(1, m1, z1) - y_ptc(0, m0, z0)),
                    cde = mean(y_ptc(1, m_ref, z1) - y_ptc(0, m_ref, z0)), m_ref = m_ref,
                    nde = NA_real_, nie = NA_real_, note = "natural effects are not identified: z is a treatment-induced mediator-outcome confounder")
      list(data = data.frame(y = y, d = d, m = m, z = z, x1 = x1, x2 = x2), truth = truth)
    },
    iv_mediator = {
      d <- stats::rbinom(n, 1, 0.5)
      u <- stats::rnorm(n)
      z <- stats::rbinom(n, 1, 0.5)
      e_m <- stats::rnorm(n); e_y <- stats::rnorm(n)
      m_iv <- function(dd, zz) 0.5 - 0.6 * dd + 1.0 * zz + 0.5 * x1 + 1.0 * u + e_m
      y_iv <- function(dd, mm) 1 + 0.8 * dd + 0.7 * mm + 0.4 * x1 + 0.3 * x2 + 1.0 * u + e_y
      m0 <- m_iv(0, z); m1 <- m_iv(1, z)
      m <- ifelse(d == 1, m1, m0)
      y <- y_iv(d, m)
      truth <- list(total = mean(y_iv(1, m1) - y_iv(0, m0)), direct = 0.8, indirect = -0.6 * 0.7,
                    first_stage = -0.6, mediator_effect = 0.7, ols_mediator_bias = "positive (u)",
                    nde = 0.8, nie = -0.6 * 0.7, cde = 0.8, m_ref = 0)
      list(data = data.frame(y = y, d = d, m = m, z = z, x1 = x1, x2 = x2), truth = truth)
    },
    front_door = {
      u <- stats::rnorm(n)
      p <- stats::plogis(0.3 * x1 + 1.2 * u)
      d <- stats::rbinom(n, 1, p)
      u_m <- stats::runif(n)
      m_fd <- function(dd) as.numeric(u_m < stats::plogis(-0.8 + 1.6 * dd + 0.3 * x1))
      e_y <- stats::rnorm(n)
      y_fd <- function(mm) 1 + 1.5 * mm + 0.5 * x1 + 0.3 * x2 + 1.2 * u + e_y
      m0 <- m_fd(0); m1 <- m_fd(1)
      m <- ifelse(d == 1, m1, m0)
      y <- y_fd(m)
      truth <- list(total = mean(y_fd(m1) - y_fd(m0)), nde = 0, nie = mean(y_fd(m1) - y_fd(m0)),
                    cde = 0, m_ref = 0, backdoor_bias = "positive (u raises both d and y)")
      list(data = data.frame(y = y, d = d, m = m, x1 = x1, x2 = x2), truth = truth)
    },
    parallel = {
      p <- stats::plogis(0.4 * x1 - 0.3 * x2)
      d <- stats::rbinom(n, 1, p)
      e1 <- stats::rnorm(n); e2 <- stats::rnorm(n); e_y <- stats::rnorm(n)
      m1_fun <- function(dd) 0.3 + 1.0 * dd + 0.5 * x1 + e1
      m2_fun <- function(dd) -0.2 + 0.4 * dd - 0.3 * x2 + e2
      y_par <- function(dd, a, b) 1 + 0.5 * dd + 0.8 * a + 1.5 * b + 0.4 * x1 + 0.2 * x2 + e_y
      a0 <- m1_fun(0); a1 <- m1_fun(1); b0 <- m2_fun(0); b1 <- m2_fun(1)
      m1 <- ifelse(d == 1, a1, a0); m2 <- ifelse(d == 1, b1, b0)
      y <- y_par(d, m1, m2)
      truth <- list(total = mean(y_par(1, a1, b1) - y_par(0, a0, b0)),
                    nde = 0.5, nie = mean(y_par(1, a1, b1) - y_par(1, a0, b0)),
                    nie_m1 = 1.0 * 0.8, nie_m2 = 0.4 * 1.5, cde = 0.5, m_ref = c(0, 0))
      list(data = data.frame(y = y, d = d, m1 = m1, m2 = m2, x1 = x1, x2 = x2), truth = truth)
    }
  )
  dat <- out$data
  attr(dat, "truth") <- out$truth
  attr(dat, "dgp") <- dgp
  dat
}

# Bivariate normal draws with correlation rho (no MASS dependency).
MASS_mvrnorm <- function(n, rho) {
  z1 <- stats::rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
  cbind(z1, z2)
}
