# Helpers for lecture 07 (difference-in-differences). Sourced by the vignette.
# All figure text is ASCII so the pdf device under --as-cran is happy.

library(ggplot2)
library(dplyr)
library(tibble)

did_colors <- c(
  "Static TWFE" = "#B22222",
  "Dynamic TWFE, binned" = "#E69F00",
  "TWFE + covariates" = "#CC79A7",
  "Sun-Abraham" = "#009E73",
  "Callaway-Sant'Anna, not yet treated" = "#0072B2",
  "Imputation (BJS)" = "#56B4E9"
)

did_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Sant'Anna's stylized staggered design: four cohorts, cohort-specific
# intercepts, unit and time effects, effects mu_g * (e + 1). The panel ends
# the year before the last cohort adopts, so that cohort is the never-treated
# comparison within the sample (the "last-treated cohort" comparison of the
# slides).
simulate_staggered_dgp <- function(n_units = 1000, n_states = 40, years = 1980:2003,
                                   groups = c(1986, 1992, 1998, 2004), mu = c(3, 2, 1, 3),
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  state <- sample(seq_len(n_states), n_units, replace = TRUE)
  state_group <- sample(groups, n_states, replace = TRUE)
  g_true <- state_group[state]
  g <- ifelse(g_true == max(groups), 0, g_true)
  alpha <- rnorm(n_units, mean = state / 5, sd = 1)
  nT <- length(years)
  dat <- tibble(
    id = rep(seq_len(n_units), each = nT),
    state = rep(state, each = nT),
    year = rep(years, times = n_units),
    g = rep(g, each = nT),
    alpha = rep(alpha, each = nT)
  ) %>%
    mutate(
      event_time = ifelse(g > 0, year - g, NA_real_),
      treated = as.integer(g > 0 & year >= g),
      lambda = (year - 1980) / 10 + rnorm(n()),
      tau = ifelse(treated == 1, mu[match(g, groups)] * (event_time + 1), 0),
      y = (2010 - rep(g_true, each = nT)) + alpha + lambda + tau + rnorm(n(), sd = 0.5)
    )
  treated_groups <- groups[groups < max(groups)]
  w <- table(factor(g[g > 0], levels = treated_groups))
  attr(dat, "truth") <- tibble(event_time = 0:5) %>%
    mutate(att = sapply(event_time, function(e) {
      keep <- treated_groups + e <= max(years)
      sum(mu[match(treated_groups[keep], groups)] * (e + 1) * w[keep]) / sum(w[keep])
    }))
  dat
}

# One replication of the Section 1 experiment: five estimators of the
# event-study parameters at event times 0..5 (truth attached to the data).
run_staggered_experiment <- function(n_sims = 50, n_units = 1000, seed = 1, max_e = 5) {
  set.seed(seed)
  out <- vector("list", n_sims)
  for (s in seq_len(n_sims)) {
    dat <- simulate_staggered_dgp(n_units = n_units)
    truth <- attr(dat, "truth")
    dat$rel_bin <- ifelse(dat$g == 0, -1000, pmin(pmax(dat$event_time, -5), max_e + 1))
    static <- fixest::feols(y ~ treated | id + year, data = dat)
    dyn <- fixest::feols(y ~ i(rel_bin, ref = c(-1, -1000)) | id + year, data = dat)
    sa <- fixest::feols(y ~ sunab(g, year) | id + year, data = dat)
    cs <- causalmetrics::att_gt(dat, id = "id", time = "year", group = "g", y = "y",
                                control_group = "notyet", n_boot = 0)
    cs_dyn <- causalmetrics::aggregate_att(cs, type = "dynamic", min_e = -5, max_e = max_e, n_boot = 0)
    imp <- causalmetrics::did_imputation(dat, id = "id", time = "year", group = "g", y = "y",
                                         horizons = 0:max_e, pre_window = 0)
    take <- function(fr, nm) {
      fr <- fr[fr$event_time %in% 0:max_e, ]
      tibble(sim = s, estimator = nm, event_time = fr$event_time, estimate = fr$estimate)
    }
    fr_dyn <- causalmetrics::event_study_frame(a = dyn)
    fr_sa <- causalmetrics::event_study_frame(a = sa)
    out[[s]] <- bind_rows(
      tibble(sim = s, estimator = "Static TWFE", event_time = 0:max_e, estimate = unname(coef(static)[["treated"]])),
      take(fr_dyn, "Dynamic TWFE, binned"),
      take(fr_sa, "Sun-Abraham"),
      take(tibble(event_time = cs_dyn$by$event_time, estimate = cs_dyn$by$estimate), "Callaway-Sant'Anna, not yet treated"),
      take(tibble(event_time = imp$by_event$event_time, estimate = imp$by_event$estimate), "Imputation (BJS)")
    ) %>% left_join(truth, by = "event_time")
  }
  bind_rows(out)
}

summarise_staggered_experiment <- function(results) {
  results %>%
    group_by(estimator, event_time) %>%
    summarise(truth = mean(att), mean = mean(estimate), bias = mean(estimate - att),
              sd = sd(estimate), .groups = "drop") %>%
    mutate(estimator = factor(estimator, levels = names(did_colors)))
}

plot_staggered_experiment <- function(results) {
  s <- summarise_staggered_experiment(results)
  ggplot(s, aes(x = event_time, y = mean, colour = estimator, group = estimator)) +
    geom_line(aes(y = truth), colour = "black", linewidth = 1, linetype = "longdash") +
    geom_pointrange(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = did_colors, drop = TRUE) +
    labs(x = "Event time", y = "Estimate (mean +/- 2 SD across simulations)",
         colour = NULL, subtitle = "Black dashed line: true event-study parameters") +
    did_theme() + guides(colour = guide_legend(nrow = 3))
}

# Kang-Schafer style two-period design of Sant'Anna's lecture 5: parallel
# trends holds conditional on X, TWFE with linear X is badly biased.
simulate_covariate_dgp <- function(n = 1000, seed = NULL, misspecify = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(4 * n), n, 4)
  Z <- cbind(exp(X[, 1] / 2), X[, 2] / (1 + exp(X[, 1])) + 10, (X[, 1] * X[, 3] / 25 + 0.6)^3, (X[, 2] + X[, 4] + 20)^2)
  Z <- scale(Z)
  W <- if (misspecify) Z else X
  f_reg <- 210 + 27.4 * W[, 1] + 13.7 * (W[, 2] + W[, 3] + W[, 4])
  f_ps <- 0.75 * (-W[, 1] + 0.5 * W[, 2] - 0.25 * W[, 3] - 0.1 * W[, 4])
  d <- as.integer(plogis(f_ps) >= runif(n))
  v <- rnorm(n, mean = d * f_reg, sd = 1)
  y1 <- f_reg + v + rnorm(n)
  y2 <- 2 * f_reg + v + rnorm(n)
  tibble(id = rep(seq_len(n), each = 2), time = rep(1:2, n), g = rep(ifelse(d == 1, 2, 0), each = 2),
         y = as.vector(rbind(y1, y2)),
         x1 = rep(X[, 1], each = 2), x2 = rep(X[, 2], each = 2), x3 = rep(X[, 3], each = 2), x4 = rep(X[, 4], each = 2))
}

run_covariate_experiment <- function(n_sims = 200, n = 1000, seed = 1, misspecify = FALSE) {
  set.seed(seed)
  rows <- vector("list", n_sims)
  for (s in seq_len(n_sims)) {
    dat <- simulate_covariate_dgp(n = n, misspecify = misspecify)
    dat$post <- as.integer(dat$time == 2)
    dat$treat <- as.integer(dat$g == 2)
    tw <- fixest::feols(y ~ treat:post + x1 + x2 + x3 + x4 | id + time, data = dat, cluster = ~id)
    twx <- fixest::feols(y ~ treat:post + post:(x1 + x2 + x3 + x4) | id + time, data = dat, cluster = ~id)
    est <- sapply(c("reg", "ipw", "dr"), function(m) {
      f <- causalmetrics::att_gt(dat, id = "id", time = "time", group = "g", y = "y",
                                 x = c("x1", "x2", "x3", "x4"), method = m, n_boot = 0)
      c(f$att_gt$att, f$att_gt$std.error)
    })
    rows[[s]] <- tibble(
      sim = s,
      estimator = c("TWFE + X", "TWFE + post:X", "Outcome regression", "IPW", "Doubly robust"),
      estimate = c(coef(tw)[["treat:post"]], coef(twx)[["treat:post"]], est[1, ]),
      std.error = c(se(tw)[["treat:post"]], se(twx)[["treat:post"]], est[2, ])
    )
  }
  bind_rows(rows)
}

summarise_covariate_experiment <- function(results, truth = 0) {
  results %>%
    group_by(estimator) %>%
    summarise(bias = mean(estimate - truth), sd = sd(estimate), rmse = sqrt(mean((estimate - truth)^2)),
              mean_se = mean(std.error),
              coverage = mean(abs(estimate - truth) <= 1.96 * std.error), .groups = "drop") %>%
    mutate(estimator = factor(estimator, levels = c("TWFE + X", "TWFE + post:X", "Outcome regression", "IPW", "Doubly robust"))) %>%
    arrange(estimator)
}
