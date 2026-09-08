# Helpers for lecture 06 (regression discontinuity). Sourced by the lecture.
# All figure text is ASCII.

library(ggplot2)
library(dplyr)
library(tibble)

rd_colors <- c(
  "Global polynomial, order 4" = "#B22222",
  "Global polynomial, order 2" = "#E69F00",
  "Local linear, MSE bandwidth" = "#0072B2",
  "Local linear, half bandwidth" = "#56B4E9",
  "Local quadratic" = "#009E73",
  "Truth" = "#000000"
)

rd_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Global polynomial jump estimate at the cutoff (Gelman-Imbens comparison).
global_poly_jump <- function(y, x, cutoff = 0, order = 4) {
  d <- as.integer(x >= cutoff)
  xc <- x - cutoff
  m <- lm(y ~ d * poly(xc, order, raw = TRUE))
  V <- sandwich::vcovHC(m, type = "HC1")
  c(estimate = unname(coef(m)["d"]), std.error = sqrt(V["d", "d"]))
}

# Local linear jump at a fixed bandwidth with the package's own local fit
# (conventional standard errors), for the coverage experiment.
local_linear_jump <- function(y, x, cutoff = 0, h, kernel = "triangular", p = 1) {
  lj <- causalmetrics:::.cm_rd_local_jump(y, x, cutoff, h, kernel = kernel, p = p)
  c(estimate = lj$jumps, std.error = lj$std.error)
}

# Monte Carlo of Section 1 and 3: estimators of the jump on the Lee design.
run_estimator_experiment <- function(n_sims = 200, n = 1000, tau = 0.04, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_rd(n, "lee", tau = tau)
    bw <- rdrobust::rdbwselect(dat$y, dat$x, c = 0)$bws
    h <- bw[1, 1]
    g4 <- global_poly_jump(dat$y, dat$x, order = 4)
    g2 <- global_poly_jump(dat$y, dat$x, order = 2)
    l1 <- local_linear_jump(dat$y, dat$x, h = h)
    lh <- local_linear_jump(dat$y, dat$x, h = h / 2)
    l2 <- local_linear_jump(dat$y, dat$x, h = h, p = 2)
    rows[[s]] <- tibble(
      sim = s,
      estimator = names(rd_colors)[1:5],
      estimate = c(g4["estimate"], g2["estimate"], l1["estimate"], lh["estimate"], l2["estimate"]),
      std.error = c(g4["std.error"], g2["std.error"], l1["std.error"], lh["std.error"], l2["std.error"]),
      h = h
    )
  }
  bind_rows(rows) %>% mutate(truth = tau)
}

summarise_estimator_experiment <- function(res) {
  res %>%
    group_by(estimator) %>%
    summarise(bias = mean(estimate - truth), sd = sd(estimate), rmse = sqrt(mean((estimate - truth)^2)),
              mean_se = mean(std.error),
              coverage = mean(abs(estimate - truth) <= 1.96 * std.error), .groups = "drop") %>%
    mutate(estimator = factor(estimator, levels = names(rd_colors))) %>%
    arrange(estimator)
}

# Monte Carlo of Section 4: coverage of conventional, robust bias-corrected,
# and honest intervals at the MSE-optimal bandwidth.
run_coverage_experiment <- function(n_sims = 200, n = 1000, tau = 0.04, seed = 2) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_rd(n, "lee", tau = tau)
    f <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
    hn <- RDHonest::RDHonest(y ~ x, data = dat, cutoff = 0)
    co <- hn$coefficients
    rows[[s]] <- tibble(
      sim = s,
      method = c("Conventional", "Bias-corrected", "Robust bias-corrected", "Honest (RDHonest)"),
      estimate = c(f$coef[1, 1], f$coef[2, 1], f$coef[3, 1], co$estimate),
      conf.low = c(f$ci[1, 1], f$ci[2, 1], f$ci[3, 1], co$conf.low),
      conf.high = c(f$ci[1, 2], f$ci[2, 2], f$ci[3, 2], co$conf.high),
      h = c(rep(f$bws[1, 1], 3), co$bandwidth)
    )
  }
  bind_rows(rows) %>% mutate(truth = tau, covered = conf.low <= truth & conf.high >= truth,
                             length = conf.high - conf.low)
}

summarise_coverage_experiment <- function(res) {
  res %>%
    group_by(method) %>%
    summarise(bias = mean(estimate - truth), coverage = mean(covered), mean_length = mean(length),
              mean_h = mean(h), .groups = "drop") %>%
    mutate(method = factor(method, levels = c("Conventional", "Bias-corrected", "Robust bias-corrected", "Honest (RDHonest)"))) %>%
    arrange(method)
}

# Monte Carlo of Section 6: precision gains from covariate adjustment.
run_adjustment_experiment <- function(n_sims = 100, n = 2000, tau = 0.2, seed = 3) {
  set.seed(seed)
  rf <- mlr3::lrn("regr.ranger", num.trees = 200, min.node.size = 20)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_rd(n, "covariates", tau = tau)
    z <- paste0("z", 1:4)
    raw <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
    lin <- rdrobust::rdrobust(dat$y, dat$x, c = 0, covs = as.matrix(dat[, z]))
    a_lin <- rd_adjust(dat, "y", "x", covariates = z, fit = TRUE)
    a_rf <- rd_adjust(dat, "y", "x", covariates = z, learner = rf, pooled = FALSE, fit = TRUE)
    get <- function(f) c(f$coef[3, 1], f$se[3, 1])
    est <- rbind(get(raw), get(lin), get(a_lin$fits$adjusted), get(a_rf$fits$adjusted))
    rows[[s]] <- tibble(sim = s, method = c("No covariates", "Linear (rdrobust covs)", "Linear adjustment (rd_adjust)",
                                            "Forest adjustment (rd_adjust)"),
                        estimate = est[, 1], std.error = est[, 2])
  }
  bind_rows(rows) %>% mutate(truth = tau)
}

summarise_adjustment_experiment <- function(res) {
  res %>%
    group_by(method) %>%
    summarise(bias = mean(estimate - truth), sd = sd(estimate), mean_se = mean(std.error),
              coverage = mean(abs(estimate - truth) <= 1.96 * std.error), .groups = "drop") %>%
    mutate(method = factor(method, levels = c("No covariates", "Linear (rdrobust covs)", "Linear adjustment (rd_adjust)",
                                              "Forest adjustment (rd_adjust)"))) %>%
    arrange(method)
}

# Section 9: kink estimates across bandwidths and polynomial orders (the
# Card et al. Figure 4 pattern).
run_kink_experiment <- function(n = 4000, tau = 2, seed = 4, h_grid = seq(0.1, 0.6, by = 0.05)) {
  set.seed(seed)
  dat <- sim_rd(n, "kink", tau = tau)
  rows <- list()
  for (p in 1:2) for (h in h_grid) {
    f <- tryCatch(rdrobust::rdrobust(dat$y, dat$x, c = 0, fuzzy = dat$d, deriv = 1, p = p, h = h), error = function(e) NULL)
    if (is.null(f)) next
    rows[[length(rows) + 1L]] <- tibble(order = paste0("Local polynomial order ", p), h = h,
                                        estimate = f$coef[3, 1], std.error = f$se[3, 1],
                                        conf.low = f$ci[3, 1], conf.high = f$ci[3, 2])
  }
  bind_rows(rows) %>% mutate(truth = tau)
}
