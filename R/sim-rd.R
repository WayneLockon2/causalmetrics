# R/sim-rd.R
#
# Simulated regression discontinuity designs with known effects at the
# cutoff, for the usage vignette, the lecture, and the tests.

#' Simulate regression discontinuity designs
#'
#' All designs use the running variable `x ~ 2 Beta(2, 4) - 1` with the
#' cutoff at 0 and the smooth conditional mean of the Lee (2008) design as
#' calibrated by Calonico, Cattaneo, and Titiunik (2014), whose jump at the
#' cutoff is 0.04. Each design returns the true effect at the cutoff in the
#' attribute `"tau_true"` and the two conditional means `mu0` and `mu1`.
#'
#' @param n Sample size.
#' @param dgp Design:
#'   * `"lee"`: sharp RD, `d = 1{x >= 0}`, effect `tau`.
#'   * `"fuzzy"`: the cutoff moves the treatment probability from
#'     `1 - compliance` to `compliance`; the effect of `d` on `y` is `tau`.
#'   * `"covariates"`: sharp RD with four covariates `z1`-`z4` that shift the
#'     outcome (`z1`, `z2` strongly, nonlinearly in `z2`) and are smooth in
#'     `x`; the covariate adjustment of [rd_adjust()] removes their noise.
#'   * `"cia"`: sharp RD where the outcome depends on `x` only through a
#'     covariate `z` correlated with `x`, so the conditional independence
#'     assumption of Angrist and Rokkanen (2015) holds and effects away from
#'     the cutoff equal `tau`.
#'   * `"kink"`: a sharp regression kink, `d = 1 + slope * max(x, 0)`, effect
#'     `tau` per unit of `d` (the slope of `y` changes by `tau * slope`; the
#'     level 1 keeps the treatment positive at the kink so an elasticity is
#'     defined).
#'   * `"discrete"`: the `"lee"` design with `x` rounded to a grid of width
#'     `grid`, as with grade point averages.
#'   * `"manipulated"`: the `"lee"` design in which a share `manip` of the
#'     units just below the cutoff move just above it, with their untreated
#'     outcomes raised, creating a density jump and a biased naive contrast.
#'   * `"multi_cutoff"`: three cutoffs at -0.5, 0, 0.5 (column `cutoff`)
#'     with effects `tau`, `tau + 0.03`, `tau + 0.06`.
#' @param tau Effect at the cutoff (default 0.04 in outcome units, 0.5 for
#'   `"cia"`, 2 for `"kink"`).
#' @param compliance Treatment probability above the cutoff in `"fuzzy"`
#'   (below it is `1 - compliance`).
#' @param sd_e Noise standard deviation (0.1295 in the CCT calibration).
#' @param grid Rounding width for `"discrete"`.
#' @param manip Share of units in `(-0.1, 0)` moved above the cutoff in
#'   `"manipulated"`.
#' @param slope Slope of the treatment above the kink in `"kink"`.
#' @param seed Optional seed.
#'
#' @return A data frame with `x`, `d`, `y`, design-specific columns, `mu0`,
#'   `mu1`, and attributes `"tau_true"`, `"dgp"`, `"cutoff"`.
#' @references
#' Calonico, S., Cattaneo, M. D., and Titiunik, R. (2014). Robust
#' nonparametric confidence intervals for regression-discontinuity designs.
#' *Econometrica*, 82(6), 2295-2326.
#' @examples
#' dat <- sim_rd(1000, "lee", seed = 1)
#' attr(dat, "tau_true")
#' @export
sim_rd <- function(n = 1000L, dgp = c("lee", "fuzzy", "covariates", "cia", "kink", "discrete",
                                      "manipulated", "multi_cutoff"),
                   tau = NULL, compliance = 0.8, sd_e = 0.1295, grid = 0.1, manip = 0.5,
                   slope = 1, seed = NULL) {
  dgp <- match.arg(dgp)
  n <- .cm_check_count(n, "n", min = 20L)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(tau)) tau <- switch(dgp, cia = 0.5, kink = 2, 0.04)
  x <- 2 * stats::rbeta(n, 2, 4) - 1
  mu_lee <- function(x) {
    ifelse(x < 0,
           0.48 + 1.27 * x + 7.18 * x^2 + 20.21 * x^3 + 21.54 * x^4 + 7.33 * x^5,
           0.52 + 0.84 * x - 3.00 * x^2 + 7.99 * x^3 - 9.01 * x^4 + 3.56 * x^5)
  }
  # untreated mean: continuous at 0 (remove the 0.04 jump of the calibration)
  mu0 <- ifelse(x < 0, mu_lee(x), mu_lee(x) - 0.04)
  e <- stats::rnorm(n, sd = sd_e)
  cutoff_attr <- 0
  out <- switch(dgp,
    lee = {
      d <- as.integer(x >= 0)
      data.frame(x = x, d = d, y = mu0 + tau * d + e, mu0 = mu0, mu1 = mu0 + tau)
    },
    discrete = {
      xg <- round(x / grid) * grid
      d <- as.integer(xg >= 0)
      mu0g <- ifelse(xg < 0, mu_lee(xg), mu_lee(xg) - 0.04)
      data.frame(x = xg, d = d, y = mu0g + tau * d + e, mu0 = mu0g, mu1 = mu0g + tau)
    },
    fuzzy = {
      above <- x >= 0
      d <- stats::rbinom(n, 1, ifelse(above, compliance, 1 - compliance))
      data.frame(x = x, d = d, z = as.integer(above), y = mu0 + tau * d + e, mu0 = mu0, mu1 = mu0 + tau)
    },
    covariates = {
      z1 <- x + stats::rnorm(n); z2 <- stats::rnorm(n); z3 <- stats::rnorm(n); z4 <- stats::rnorm(n)
      d <- as.integer(x >= 0)
      g <- 0.3 * z1 + 0.3 * (z2^2 - 1) + 0.05 * z3
      data.frame(x = x, d = d, z1 = z1, z2 = z2, z3 = z3, z4 = z4,
                 y = mu0 + g + tau * d + e, mu0 = mu0 + g, mu1 = mu0 + g + tau)
    },
    cia = {
      z <- x + stats::rnorm(n, sd = 0.5)
      d <- as.integer(x >= 0)
      m0 <- 1 + 2 * z
      data.frame(x = x, d = d, z = z, y = m0 + tau * d + e, mu0 = m0, mu1 = m0 + tau)
    },
    kink = {
      d <- 1 + slope * pmax(x, 0)
      data.frame(x = x, d = d, y = mu0 + tau * d + e, mu0 = mu0 + tau, mu1 = mu0 + tau * d)
    },
    manipulated = {
      d <- as.integer(x >= 0)
      y0 <- mu0
      movers <- which(x > -0.1 & x < 0 & stats::runif(n) < manip)
      x[movers] <- -x[movers]
      y0[movers] <- y0[movers] + 0.1
      d <- as.integer(x >= 0)
      data.frame(x = x, d = d, y = y0 + tau * d + e, mu0 = y0, mu1 = y0 + tau, mover = seq_len(n) %in% movers)
    },
    multi_cutoff = {
      cuts <- c(-0.5, 0, 0.5)
      cutoff <- sample(cuts, n, replace = TRUE)
      xx <- cutoff + 0.6 * x
      d <- as.integer(xx >= cutoff)
      eff <- tau + 0.03 * match(cutoff, cuts) - 0.03
      m0 <- 0.5 + 0.8 * (xx - cutoff) - 0.6 * (xx - cutoff)^2
      cutoff_attr <- cuts
      data.frame(x = xx, cutoff = cutoff, d = d, y = m0 + eff * d + e, mu0 = m0, mu1 = m0 + eff, tau_cutoff = eff)
    }
  )
  attr(out, "tau_true") <- if (dgp == "multi_cutoff") tau + c(0, 0.03, 0.06) else tau
  attr(out, "dgp") <- dgp
  attr(out, "cutoff") <- cutoff_attr
  out
}
