# R/sim-hte.R
#
# Simulated data with known conditional average treatment effects, for the
# usage vignette, the lecture, and the tests.

#' Simulate data with heterogeneous treatment effects
#'
#' Designs used in the lecture on heterogeneous effects and policy learning.
#' Every design returns the true conditional average treatment effect
#' (`tau_true`), propensity score (`p_true`), and conditional means
#' (`mu0_true`, `mu1_true`) next to the observed data.
#'
#' @param n Sample size.
#' @param dgp Design:
#'   * `"smooth"`: five standard normal covariates, propensity
#'     `plogis(0.5 x1 - 0.5 x2)`, baseline `x1 + x2^2 + 0.5 x3`, and a CATE
#'     `1 + x1 + 0.5 (x2^2 - 1)` that depends on `x1` and `x2` only.
#'   * `"simple_cate"`: DGP 1 of Chernozhukov et al. (2026, ch. 15). One
#'     covariate `x1 ~ U(0, 1)`, rare treatment (`p_treat = 0.05`), constant
#'     effect 0.5, baseline `0.3 * 1{x1 in [0.6, 0.8]}`, noise sd 0.05.
#'   * `"complex_cate"`: DGP 2. Effect `0.5 * 1{x1 in [0.6, 0.8]}`, baseline
#'     0.1, `p_treat = 0.05`.
#'   * `"unbalanced"`: DGP 3, the same as DGP 2 with `p_treat = 0.95`.
#'   * `"binary_outcome"`: the conversion example of Facure (2023, ch. 23).
#'     A randomized nudge, `age ~ Gamma(10, 4)`, `income ~ 100 Gamma(20, 2)`,
#'     latent index `-4.5 + 0.001 income + nudge (1 + 0.01 age) + N(0, 1)`,
#'     outcome `1{latent > threshold}`. Income moves the baseline only; age
#'     moves the effect only.
#'   * `"policy"`: five uniform covariates on `[-1, 1]`, a CATE that is
#'     piecewise constant on a depth-two tree in `x1`, `x2`, `x3` plus a
#'     small smooth term, confounded assignment `plogis(0.5 x1 + 0.5 x4)`.
#'     The optimal policy is the depth-two tree stored in
#'     `attr(, "optimal_tree")`.
#' @param p_treat Treatment probability for the three Chernozhukov et al.
#'   designs (ignored otherwise).
#' @param threshold Latent-index threshold for `"binary_outcome"` (0 gives
#'   about 50 percent conversion, 2 about 12 percent, -2 about 93 percent).
#' @param sd_e Noise standard deviation (designs `"smooth"` and `"policy"`).
#' @param seed Optional seed.
#'
#' @return A data frame with columns `y`, `d`, the covariates, `tau_true`,
#'   `p_true`, `mu0_true`, `mu1_true`, and attribute `"dgp"`.
#'
#' @references
#' Chernozhukov, V., Hansen, C., Kallus, N., Spindler, M., and Syrgkanis, V.
#' (2026). *Applied Causal Inference Powered by ML and AI*, chapter 15.
#'
#' Facure, M. (2023). *Causal Inference for the Brave and True*, chapter 23.
#'
#' @examples
#' dat <- sim_hte(500, dgp = "smooth", seed = 1)
#' head(dat)
#' mean(dat$tau_true)
#' @export
sim_hte <- function(n = 1000L, dgp = c("smooth", "simple_cate", "complex_cate", "unbalanced",
                                       "binary_outcome", "policy"),
                    p_treat = NULL, threshold = 0, sd_e = 1, seed = NULL) {
  dgp <- match.arg(dgp)
  n <- .cm_check_count(n, "n", min = 10L)
  if (!is.null(seed)) set.seed(seed)
  out <- switch(dgp,
    smooth = {
      X <- matrix(stats::rnorm(n * 5), n, 5)
      colnames(X) <- paste0("x", 1:5)
      p <- stats::plogis(0.5 * X[, 1] - 0.5 * X[, 2])
      d <- stats::rbinom(n, 1, p)
      mu0 <- X[, 1] + X[, 2]^2 + 0.5 * X[, 3]
      tau <- 1 + X[, 1] + 0.5 * (X[, 2]^2 - 1)
      y <- mu0 + d * tau + stats::rnorm(n, sd = sd_e)
      data.frame(y = y, d = d, X, tau_true = tau, p_true = p, mu0_true = mu0, mu1_true = mu0 + tau)
    },
    simple_cate = , complex_cate = , unbalanced = {
      if (is.null(p_treat)) p_treat <- if (dgp == "unbalanced") 0.95 else 0.05
      x1 <- stats::runif(n)
      bump <- as.numeric(x1 >= 0.6 & x1 <= 0.8)
      if (dgp == "simple_cate") {
        tau <- rep(0.5, n)
        mu0 <- 0.3 * bump
      } else {
        tau <- 0.5 * bump
        mu0 <- rep(0.1, n)
      }
      p <- rep(p_treat, n)
      d <- stats::rbinom(n, 1, p)
      y <- mu0 + d * tau + stats::rnorm(n, sd = 0.05)
      data.frame(y = y, d = d, x1 = x1, tau_true = tau, p_true = p, mu0_true = mu0, mu1_true = mu0 + tau)
    },
    binary_outcome = {
      d <- stats::rbinom(n, 1, 0.5)
      age <- stats::rgamma(n, shape = 10, scale = 4)
      income <- stats::rgamma(n, shape = 20, scale = 2) * 100
      mean0 <- -4.5 + 0.001 * income
      mean1 <- mean0 + 1 + 0.01 * age
      latent <- ifelse(d == 1, mean1, mean0) + stats::rnorm(n)
      y <- as.integer(latent > threshold)
      mu0 <- 1 - stats::pnorm(threshold - mean0)
      mu1 <- 1 - stats::pnorm(threshold - mean1)
      data.frame(y = y, d = d, age = age, income = income, latent = latent,
                 tau_true = mu1 - mu0, p_true = rep(0.5, n), mu0_true = mu0, mu1_true = mu1,
                 tau_latent = 1 + 0.01 * age)
    },
    policy = {
      X <- matrix(stats::runif(n * 5, -1, 1), n, 5)
      colnames(X) <- paste0("x", 1:5)
      p <- stats::plogis(0.5 * X[, 1] + 0.5 * X[, 4])
      d <- stats::rbinom(n, 1, p)
      mu0 <- 0.5 * X[, 1] + X[, 2]^2 + X[, 4]
      tau <- ifelse(X[, 1] <= 0,
                    ifelse(X[, 3] <= -0.25, 1.5, -0.5),
                    ifelse(X[, 2] <= 0.25, -1, 1)) + 0.3 * X[, 4]
      y <- mu0 + d * tau + stats::rnorm(n, sd = sd_e)
      data.frame(y = y, d = d, X, tau_true = tau, p_true = p, mu0_true = mu0, mu1_true = mu0 + tau)
    }
  )
  attr(out, "dgp") <- dgp
  if (dgp == "policy") {
    attr(out, "optimal_tree") <- list(
      root = "x1 <= 0",
      left = "x3 <= -0.25: treat (tau = 1.5 + 0.3 x4); else do not treat (tau = -0.5 + 0.3 x4)",
      right = "x2 <= 0.25: do not treat (tau = -1 + 0.3 x4); else treat (tau = 1 + 0.3 x4)"
    )
  }
  out
}
