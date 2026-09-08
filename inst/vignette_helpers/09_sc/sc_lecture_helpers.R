# Helpers for lecture 09 (synthetic control). Sourced by the lecture.
# All figure text is ASCII so the pdf device under --as-cran is happy.

library(ggplot2)
library(dplyr)
library(tibble)

sc_colors <- c(
  "DiD" = "#B22222",
  "Synthetic control" = "#0072B2",
  "Demeaned synthetic control" = "#009E73",
  "Augmented synthetic control" = "#56B4E9",
  "Synthetic DiD" = "#E69F00",
  "Unconstrained regression" = "#CC79A7"
)

sc_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# The two designs of the opening experiment (Ferman and Pinto 2021, Table 1).
# A: the treated unit's fixed effect sits two standard deviations above the
#    donors', with no post-treatment break in the factors.
# B: the treated unit loads 1.5 standard deviations above the donors on the
#    first factor, and that factor shifts up by one standard deviation after
#    treatment.
sc_designs <- list(
  "A: level in the tail" = list(select_level = 2, select_loading = 0, factor_shift = 0),
  "B: loading in the tail, factor break" = list(select_level = 0, select_loading = 1.5, factor_shift = 1)
)

# Six estimators on one simulated panel. Returns one row per estimator.
fit_sc_estimators <- function(dat) {
  sc <- synth_control(dat, "id", "time", "y", "d")
  dm <- synth_control(dat, "id", "time", "y", "d", demean = TRUE)
  aug <- synth_control(dat, "id", "time", "y", "d", demean = TRUE, augment = "ridge")
  ols <- synth_control(dat, "id", "time", "y", "d", demean = TRUE, constraints = "none")
  did <- sdid_weights(dat, "id", "time", "y", "d", estimator = "did")
  sdid <- sdid_weights(dat, "id", "time", "y", "d", estimator = "sdid")
  tibble(
    estimator = names(sc_colors),
    estimate = c(did$estimate, sc$estimate, dm$estimate, aug$estimate, sdid$estimate, ols$estimate),
    pre_rmspe = c(NA, sc$pre_rmspe, dm$pre_rmspe, aug$pre_rmspe, NA, ols$pre_rmspe),
    l2 = c(1 / did$N0, sc$concentration[["l2"]], dm$concentration[["l2"]], aug$concentration[["l2"]],
           sum(sdid$omega$weight^2), sum(ols$weights$weight^2))
  )
}

run_sc_experiment <- function(n_sims = 200, designs = sc_designs, n_donors = 20, t_pre = 20, t_post = 5,
                              effect = 1, sd_e = 1, seed = 1) {
  set.seed(seed)
  bind_rows(lapply(names(designs), function(nm) {
    d <- designs[[nm]]
    bind_rows(lapply(seq_len(n_sims), function(s) {
      dat <- sim_synth_panel(n_donors = n_donors, t_pre = t_pre, t_post = t_post, sd_e = sd_e,
                             select_level = d$select_level, select_loading = d$select_loading,
                             factor_shift = d$factor_shift, effect = effect)
      fit_sc_estimators(dat) %>% mutate(design = nm, sim = s)
    }))
  }))
}

summarise_sc_experiment <- function(results, effect = 1) {
  results %>%
    group_by(design, estimator) %>%
    summarise(
      bias = mean(estimate - effect),
      sd = sd(estimate),
      rmse = sqrt(mean((estimate - effect)^2)),
      pre_rmspe = mean(pre_rmspe),
      l2 = mean(l2),
      .groups = "drop"
    ) %>%
    mutate(estimator = factor(estimator, levels = names(sc_colors))) %>%
    arrange(design, estimator)
}

plot_sc_experiment <- function(results, effect = 1, colors = sc_colors, xlim = 4) {
  d <- results %>%
    mutate(error = estimate - effect,
           estimator = factor(estimator, levels = names(colors)))
  ggplot(d %>% filter(abs(error) <= xlim), aes(x = error, fill = estimator)) +
    geom_histogram(binwidth = 0.25, boundary = 0, colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    facet_grid(design ~ estimator, labeller = labeller(design = label_wrap_gen(14), estimator = label_wrap_gen(14))) +
    scale_fill_manual(values = colors, guide = "none") +
    coord_cartesian(xlim = c(-xlim, xlim)) +
    labs(x = "Estimate minus true effect", y = "Simulations") +
    sc_theme() +
    theme(strip.text = element_text(size = 8), strip.text.y = element_text(angle = 0))
}

# A treated unit whose latent structure (unit effect, trend, factor loadings)
# is an exact convex combination of the donors', plus its own noise: weights
# w_true reproduce the treated unit's loadings. Used to show that the fitted
# weights do not converge to w_true when the noise variance is positive
# (Ferman and Pinto 2021, Proposition 1).
simulate_convex_treated <- function(n_donors = 8, t_pre = 50, t_post = 5, sd_e = 1, w_true, trend = 0.2, seed = NULL) {
  dat <- sim_synth_panel(n_donors = n_donors, t_pre = t_pre, t_post = t_post, sd_e = sd_e, trend = trend, seed = seed)
  Tn <- t_pre + t_post
  lam <- attr(dat, "factors"); mu <- attr(dat, "loadings"); cj <- attr(dat, "unit_effects")
  latent <- outer(cj[seq_len(n_donors)], rep(1, Tn)) + outer(rep(1, n_donors), trend * seq_len(Tn)) +
    mu[seq_len(n_donors), , drop = FALSE] %*% t(lam)
  dat$y[dat$id == n_donors + 1] <- as.numeric(crossprod(latent, w_true)) + rnorm(Tn, sd = sd_e)
  attr(dat, "donor_loadings") <- mu[seq_len(n_donors), , drop = FALSE]
  dat
}

# Overfitting: pre-treatment fit improves mechanically with the number of
# donors, while the post-treatment error does not. No treatment effect.
run_overfit_experiment <- function(n_sims = 100, n_donors = c(5, 10, 20, 40), t_pre = 20, t_post = 5,
                                   sd_e = 1, seed = 1) {
  set.seed(seed)
  bind_rows(lapply(n_donors, function(J) {
    bind_rows(lapply(seq_len(n_sims), function(s) {
      dat <- sim_synth_panel(n_donors = J, t_pre = t_pre, t_post = t_post, sd_e = sd_e)
      f <- synth_control(dat, "id", "time", "y", "d", demean = TRUE)
      tibble(n_donors = J, sim = s, estimate = f$estimate, pre_rmspe = f$pre_rmspe,
             post_rmspe = f$post_rmspe, l2 = f$concentration[["l2"]])
    }))
  }))
}

summarise_overfit <- function(results) {
  results %>%
    group_by(n_donors) %>%
    summarise(pre_rmspe = mean(pre_rmspe), post_rmspe = mean(post_rmspe), bias = mean(estimate),
              sd = sd(estimate), rmse = sqrt(mean(estimate^2)), l2 = mean(l2), .groups = "drop")
}

# Size and power of the inference tools: in-space placebo rank test,
# moving-block conformal test, and the Ferman-Pinto specification test.
run_inference_experiment <- function(n_sims = 200, effect = 0, n_donors = 20, t_pre = 20, t_post = 5,
                                     sd_e = 1, seed = 1) {
  set.seed(seed)
  bind_rows(lapply(seq_len(n_sims), function(s) {
    dat <- sim_synth_panel(n_donors = n_donors, t_pre = t_pre, t_post = t_post, sd_e = sd_e, effect = effect)
    fit <- synth_control(dat, "id", "time", "y", "d", demean = TRUE)
    pl <- synth_placebo(fit)
    ci <- synth_conformal(fit, grid = fit$estimate)
    st <- synth_spec_test(fit)
    tibble(sim = s, effect = effect, estimate = fit$estimate,
           p_placebo = pl$p_value, p_conformal = ci$p_value, p_spec = st$p_value)
  }))
}

summarise_inference <- function(results) {
  results %>%
    group_by(effect) %>%
    summarise(
      `In-space placebo` = mean(p_placebo <= 0.10),
      `Conformal` = mean(p_conformal <= 0.10),
      `Specification test` = mean(p_spec <= 0.10),
      .groups = "drop"
    )
}
