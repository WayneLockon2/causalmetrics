# R/sim-iv.R
#
# Simulated designs for the instrumental-variables lecture, each with the
# true parameter attached.

#' Simulate instrumental-variable designs
#'
#' @param n Sample size.
#' @param dgp Design:
#'   * `"linear"`: partially linear IV. Five standard normal covariates, an
#'     unobserved confounder `a`, `h(x) = 0.5 x1^2 + 0.5 sin(2 x2)`, instrument
#'     `z = h(x) + e_z`, treatment `d = pi z + 0.8 a + 0.5 x1^2 + e_d`, outcome
#'     `y = theta d + a + h(x) + 0.5 x3 + e_y`. The instrument is
#'     excluded only conditional on nonlinear functions of `x1` and `x2`, so
#'     linear controls leave bias. Attribute `theta`.
#'   * `"weak"`: three covariates entering linearly (`h(x) = 0.5 x1 + 0.5 x2`),
#'     the same confounder with small idiosyncratic noise (so the endogeneity
#'     correlation is about 0.9), and a first-stage coefficient chosen so that
#'     the concentration parameter is about `concentration`
#'     (`pi = sqrt(concentration / n)`).
#'   * `"late"`: binary instrument assigned with probability `plogis(0.5 x1)`,
#'     compliance types drawn with covariate-dependent probabilities
#'     (always-takers, never-takers, compliers), complier effect
#'     `2 + x1`, never-taker and always-taker effects that differ from it.
#'     Attributes `late` (the population LATE) and `ate`.
#'   * `"mte"`: a Roy model with a continuous instrument. Unobserved
#'     resistance `v ~ U(0, 1)`, propensity `p(z, x) = plogis(-0.5 + z + 0.3 x1)`,
#'     `d = 1{p > v}`, and a marginal treatment effect `mte(u) = 1.5 - 2 (u - 0.5) - 3 (u - 0.5)^2`
#'     plus a covariate term `0.5 x1`. Attributes `mte` (the function), `ate`,
#'     and `att`.
#'   * `"probit_cf"`: binary outcome with an endogenous continuous regressor.
#'     `d = 0.8 z + 0.5 x1 + v`, `y* = -0.3 + 0.7 d + 0.5 x1 + u`,
#'     `(u, v)` bivariate normal with correlation `rho`, `y = 1{y* > 0}`.
#'     Attribute `ape` (the true average partial effect of `d`).
#'   * `"shift_share"`: `n` regions with Dirichlet shares over `n_industries`
#'     industries, national shocks `g_k`, regional treatment
#'     `d_i = sum_k s_ik (g_k + eta_ik) + u_i`, outcome
#'     `y_i = beta d_i + 0.5 x1 + e_i` with `u_i` correlated with `e_i`.
#'     Shares are columns `s1, ..., sK`, shocks are in attribute `shocks`.
#'   * `"judge"`: `n_judges` judges with leniencies drawn from `U(0.1, 0.9)`,
#'     random assignment, `d = 1{leniency + 0.3 x1 + 0.5 a > u}` and
#'     `y = theta d + a + 0.5 x1 + e`. Attribute `theta`.
#' @param pi First-stage coefficient (designs `"linear"`, `"probit_cf"`).
#' @param concentration Concentration parameter for `"weak"`.
#' @param rho Error correlation (`"probit_cf"`).
#' @param n_industries,n_judges Sizes for the shift-share and judge designs.
#' @param seed Optional seed.
#'
#' @return A data frame with attribute `"dgp"` and the design's truth as
#'   attributes.
#' @examples
#' dat <- sim_iv(500, dgp = "late", seed = 1)
#' attr(dat, "late")
#' @export
sim_iv <- function(n = 1000L, dgp = c("linear", "weak", "late", "mte", "probit_cf", "shift_share", "judge"),
                   pi = 1, concentration = 10, rho = 0.5, n_industries = 20L, n_judges = 30L, seed = NULL) {
  dgp <- match.arg(dgp)
  n <- .cm_check_count(n, "n", min = 20L)
  if (!is.null(seed)) set.seed(seed)
  out <- switch(dgp,
    linear = , weak = {
      p_x <- if (dgp == "linear") 5L else 3L
      X <- matrix(stats::rnorm(n * p_x), n, p_x, dimnames = list(NULL, paste0("x", seq_len(p_x))))
      a <- stats::rnorm(n)
      h <- if (dgp == "linear") 0.5 * X[, 1]^2 + 0.5 * sin(2 * X[, 2]) else 0.5 * X[, 1] + 0.5 * X[, 2]
      z <- h + stats::rnorm(n)
      pi_use <- if (dgp == "weak") sqrt(concentration / n) else pi
      sd_e <- if (dgp == "linear") 1 else 0.3
      d <- pi_use * z + 0.8 * a + (if (dgp == "linear") 0.5 * X[, 1]^2 else 0.5 * X[, 1]) + stats::rnorm(n, sd = sd_e)
      y <- 1 * d + a + h + 0.5 * X[, 3] + stats::rnorm(n, sd = sd_e)
      df <- data.frame(y = y, d = d, z = z, X)
      attr(df, "theta") <- 1
      attr(df, "pi") <- pi_use
      df
    },
    late = {
      x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
      pz <- stats::plogis(0.5 * x1)
      z <- stats::rbinom(n, 1, pz)
      # compliance types
      p_always <- 0.15 + 0.1 * stats::plogis(x2)
      p_never <- 0.25 - 0.1 * stats::plogis(x1)
      u <- stats::runif(n)
      type <- ifelse(u < p_always, "always", ifelse(u < p_always + p_never, "never", "complier"))
      d0 <- as.integer(type == "always")
      d1 <- as.integer(type != "never")
      d <- ifelse(z == 1, d1, d0)
      tau <- ifelse(type == "complier", 2 + x1, ifelse(type == "always", 3 + x1, 0.5 + x1))
      y0 <- x1 + x2^2 + stats::rnorm(n)
      y <- y0 + d * tau
      df <- data.frame(y = y, d = d, z = z, x1 = x1, x2 = x2, type = type, tau_true = tau)
      attr(df, "late") <- mean(tau[type == "complier"])
      attr(df, "ate") <- mean(tau)
      attr(df, "complier_share") <- mean(type == "complier")
      df
    },
    mte = {
      x1 <- stats::rnorm(n)
      z <- stats::rnorm(n)
      p <- stats::plogis(-0.5 + z + 0.3 * x1)
      v <- stats::runif(n)
      d <- as.integer(p > v)
      mte_fun <- function(u) 1.5 - 2 * (u - 0.5) - 3 * (u - 0.5)^2
      y0 <- 0.5 * x1 + stats::rnorm(n, sd = 0.5)
      y1 <- y0 + mte_fun(v) + 0.5 * x1 + stats::rnorm(n, sd = 0.3)
      y <- ifelse(d == 1, y1, y0)
      df <- data.frame(y = y, d = d, z = z, x1 = x1, p_true = p, v = v)
      attr(df, "mte") <- function(u) mte_fun(u) + 0.5 * mean(x1)
      attr(df, "ate") <- mean(mte_fun(v)) + 0.5 * mean(x1)
      attr(df, "att") <- mean((mte_fun(v) + 0.5 * x1)[d == 1])
      df
    },
    probit_cf = {
      x1 <- stats::rnorm(n)
      z <- stats::rnorm(n)
      uv <- matrix(stats::rnorm(2 * n), n, 2)
      v <- uv[, 2]
      u <- rho * v + sqrt(1 - rho^2) * uv[, 1]
      d <- 0.8 * z + 0.5 * x1 + v
      ystar <- -0.3 + 0.7 * d + 0.5 * x1 + u
      y <- as.integer(ystar > 0)
      df <- data.frame(y = y, d = d, z = z, x1 = x1, v_true = v)
      # true APE of d: E[phi(-0.3 + 0.7 d + 0.5 x1 + rho v) * 0.7] averaged over (d, x1, v)
      # with the structural error u given v having sd sqrt(1 - rho^2)
      s <- sqrt(1 - rho^2)
      attr(df, "ape") <- mean(stats::dnorm((-0.3 + 0.7 * d + 0.5 * x1 + rho * v) / s) * 0.7 / s)
      df
    },
    shift_share = {
      K <- .cm_check_count(n_industries, "n_industries", min = 2L)
      S <- matrix(stats::rgamma(n * K, shape = 0.5), n, K)
      S <- S / rowSums(S)
      colnames(S) <- paste0("s", seq_len(K))
      g <- stats::rnorm(K, sd = 2)
      x1 <- stats::rnorm(n)
      eta <- matrix(stats::rnorm(n * K, sd = 0.5), n, K)
      e <- stats::rnorm(n)
      u <- 0.7 * e + stats::rnorm(n, sd = 0.5)
      d <- as.numeric(rowSums(S * (matrix(g, n, K, byrow = TRUE) + eta))) + 0.3 * x1 + u
      y <- 1 * d + 0.5 * x1 + e
      df <- data.frame(y = y, d = d, x1 = x1, S)
      df$b <- as.numeric(S %*% g)
      attr(df, "beta") <- 1
      attr(df, "shocks") <- g
      df
    },
    judge = {
      J <- .cm_check_count(n_judges, "n_judges", min = 2L)
      len <- stats::runif(J, 0.1, 0.9)
      judge <- sample.int(J, n, replace = TRUE)
      x1 <- stats::rnorm(n)
      a <- stats::rnorm(n)
      d <- as.integer(len[judge] + 0.3 * stats::plogis(x1) + 0.4 * a > stats::runif(n) + 0.35)
      y <- 1 * d + a + 0.5 * x1 + stats::rnorm(n)
      df <- data.frame(y = y, d = d, judge = judge, x1 = x1, leniency_true = len[judge])
      attr(df, "theta") <- 1
      df
    }
  )
  attr(out, "dgp") <- dgp
  out
}
