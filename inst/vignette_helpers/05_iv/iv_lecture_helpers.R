# Helpers for lecture 05 (instrumental variables and unobserved confounding).
# Sourced by the lecture. All figure text is ASCII.

library(ggplot2)
library(dplyr)
library(tibble)
library(mlr3)
library(mlr3learners)
lgr::get_logger("mlr3")$set_threshold("warn")

iv_colors <- c(
  "OLS" = "#B22222",
  "2SLS" = "#0072B2",
  "2SLS, linear controls" = "#0072B2",
  "DML-IV, forests" = "#009E73",
  "DML-LATE" = "#009E73",
  "Control function" = "#009E73",
  "Forbidden regression" = "#E69F00",
  "Naive probit" = "#B22222",
  "Wald" = "#0072B2",
  "Anderson-Rubin" = "#009E73",
  "Truth" = "#000000"
)

iv_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

iv_forest <- function(min_node = 10, trees = 300) lrn("regr.ranger", num.trees = trees, min.node.size = min_node)

# Monte Carlo of Section 5.2: the partially linear IV design has a nonlinear
# confounder path. Compare OLS, 2SLS with linear controls, and DML-IV with
# forests.
run_pliv_experiment <- function(n_sims = 100, n = 2000, seed = 1) {
  set.seed(seed)
  rows <- list()
  x <- paste0("x", 1:5)
  for (s in seq_len(n_sims)) {
    dat <- sim_iv(n, dgp = "linear")
    ols <- estimatr::lm_robust(y ~ d + x1 + x2 + x3 + x4 + x5, data = dat)
    tsls <- estimatr::iv_robust(y ~ d + x1 + x2 + x3 + x4 + x5 | z + x1 + x2 + x3 + x4 + x5, data = dat)
    dml <- est_dml(dat, "y", "d", x, z = "z", model = "pliv",
                   learner_l = iv_forest(), learner_m = iv_forest(), learner_z = iv_forest())
    rows[[s]] <- tibble(
      sim = s,
      estimator = c("OLS", "2SLS, linear controls", "DML-IV, forests"),
      estimate = c(coef(ols)[["d"]], coef(tsls)[["d"]], dml$estimate),
      std.error = c(ols$std.error[["d"]], tsls$std.error[["d"]], dml$std.error)
    )
  }
  bind_rows(rows)
}

# Monte Carlo of Section 5.3: compliance varies with the covariates. 2SLS
# with linear controls versus the doubly robust LATE.
run_late_experiment <- function(n_sims = 100, n = 2000, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_iv(n, dgp = "late")
    truth <- attr(dat, "late")
    tsls <- estimatr::iv_robust(y ~ d + x1 + x2 | z + x1 + x2, data = dat)
    dml <- est_dml(dat, "y", "d", c("x1", "x2"), z = "z", model = "iivm",
                   learner_mu0 = iv_forest(), learner_mu1 = iv_forest())
    rows[[s]] <- tibble(
      sim = s, truth = truth,
      estimator = c("2SLS, linear controls", "DML-LATE"),
      estimate = c(coef(tsls)[["d"]], dml$estimate),
      std.error = c(tsls$std.error[["d"]], dml$std.error)
    )
  }
  bind_rows(rows)
}

summarise_experiment <- function(res, truth = 1) {
  if (!"truth" %in% names(res)) res$truth <- truth
  res %>%
    group_by(estimator) %>%
    summarise(bias = mean(estimate - truth), sd = sd(estimate), rmse = sqrt(mean((estimate - truth)^2)),
              mean_se = mean(std.error),
              coverage = mean(abs(estimate - truth) <= 1.96 * std.error), .groups = "drop")
}

# Monte Carlo of Section 6: the weak-instrument problem. For several
# concentration parameters, the sampling distribution of 2SLS and the
# coverage of the Wald interval versus the Anderson-Rubin set.
run_weak_iv_experiment <- function(n_sims = 500, n = 500, concentrations = c(2, 10, 50), seed = 1) {
  set.seed(seed)
  rows <- list()
  for (mu2 in concentrations) {
    for (s in seq_len(n_sims)) {
      dat <- sim_iv(n, dgp = "weak", concentration = mu2)
      ar <- iv_ar_confidence_set(dat, "y", "d", "z", x = paste0("x", 1:3))
      ols <- coef(lm(y ~ d + x1 + x2 + x3, data = dat))[["d"]]
      inside_ar <- any(ar$intervals$lower <= 1 & ar$intervals$upper >= 1)
      rows[[length(rows) + 1L]] <- tibble(
        concentration = mu2, sim = s, tsls = ar$estimate, se = ar$std.error, ols = ols,
        first_stage_F = ar$first_stage$F_robust,
        wald_covers = abs(ar$estimate - 1) <= 1.96 * ar$std.error,
        ar_covers = inside_ar, ar_bounded = ar$type == "bounded",
        ar_length = if (ar$type == "bounded") ar$intervals$upper - ar$intervals$lower else Inf
      )
    }
  }
  bind_rows(rows)
}

# Monte Carlo of Section 4.3: probit outcome with an endogenous regressor.
# Naive probit, the forbidden regression, and the control function, judged by
# the average partial effect.
run_cf_experiment <- function(n_sims = 200, n = 2000, rho = 0.6, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_iv(n, dgp = "probit_cf", rho = rho)
    truth <- attr(dat, "ape")
    naive <- suppressWarnings(glm(y ~ d + x1, data = dat, family = binomial("probit")))
    dat$d_hat <- fitted(lm(d ~ z + x1, data = dat))
    forb <- suppressWarnings(glm(y ~ d_hat + x1, data = dat, family = binomial("probit")))
    aug <- cf_residuals(dat, "d", "z", "x1")
    cf <- suppressWarnings(glm(y ~ d + x1 + v_hat, data = aug, family = binomial("probit")))
    rows[[s]] <- tibble(
      sim = s, truth = truth,
      estimator = c("Naive probit", "Forbidden regression", "Control function"),
      ape = c(cf_ape(naive, dat, "d"), cf_ape(forb, dat, "d_hat"), cf_ape(cf, aug, "d"))
    )
  }
  bind_rows(rows)
}

plot_sampling <- function(res, truth = 1, colors = iv_colors, ncol = 3) {
  ggplot(res, aes(x = estimate, fill = estimator)) +
    geom_histogram(bins = 40, alpha = 0.8, colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = truth, linetype = "dashed") +
    facet_wrap(~ estimator, scales = "free_x", ncol = ncol) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(x = "Estimate", y = "Simulations") + iv_theme()
}
