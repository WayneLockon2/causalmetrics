# vignette_helpers_soo.R
# ------------------------------------------------------------
# Additional visualization helpers for the
# "Selection on Observables" vignette.
#
# This file complements causal_estimators_ggplot_template.R.
# Together they provide:
#   - the five-estimator panel (template)
#   - the opening scatter / naive-gap figure
#   - the propensity-score overlap histogram
#   - the IPW weight histogram
#   - the AIPW Series R (Robinson view, 4 plots)
#   - the AIPW Series A (correction view, 4 plots)
#   - a simple balance plot (standardized mean differences)
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# Theme + palette (kept consistent with the template)
# ------------------------------------------------------------

soo_colors <- c(Control = "#993C1D", Treated = "#185FA5")
soo_fills  <- c(Control = "#F5C4B3", Treated = "#B5D4F4")

soo_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey35")
    )
}

# ------------------------------------------------------------
# 1. Opening scatter with naive-gap annotation
# ------------------------------------------------------------

plot_opening_scatter <- function(df) {
  group_means <- df %>%
    mutate(group = if_else(treat == 1, "Treated", "Control")) %>%
    group_by(group) %>%
    summarise(y_mean = mean(y), x_mean = mean(x), .groups = "drop")

  naive_diff <- with(df, mean(y[treat == 1]) - mean(y[treat == 0]))

  df_plot <- df %>%
    mutate(group = if_else(treat == 1, "Treated", "Control"))

  ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(aes(fill = group, color = group),
               shape = 21, size = 2.6, stroke = 0.6, alpha = 0.85) +
    geom_hline(data = group_means,
               aes(yintercept = y_mean, color = group),
               linetype = "dashed", linewidth = 0.6) +
    geom_segment(
      data = tibble(
        x = max(df$x) + 1,
        xend = max(df$x) + 1,
        y = group_means$y_mean[group_means$group == "Control"],
        yend = group_means$y_mean[group_means$group == "Treated"]
      ),
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.18, "cm"), ends = "both"),
      inherit.aes = FALSE,
      color = "grey25", linewidth = 0.6
    ) +
    annotate(
      "label",
      x = max(df$x) + 3,
      y = mean(group_means$y_mean),
      label = paste0("Naive gap\n= ", round(naive_diff, 2)),
      hjust = 0, size = 3.3,
      border.colour = NA, fill = alpha("white", 0.8)
    ) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "The selection problem",
      subtitle = "Treated and control groups differ in X; the raw Y gap is not the causal effect.",
      x = "Pre-treatment covariate X",
      y = "Outcome Y",
      fill = NULL, color = NULL
    ) +
    coord_cartesian(xlim = c(min(df$x) - 2, max(df$x) + 18)) +
    soo_theme()
}

# ------------------------------------------------------------
# 2. Propensity-score overlap histogram
# ------------------------------------------------------------

plot_ps_overlap <- function(df) {
  ps_model <- glm(treat ~ x, data = df, family = binomial())
  df_plot <- df %>%
    mutate(
      ps = predict(ps_model, type = "response"),
      group = if_else(treat == 1, "Treated", "Control")
    )

  ggplot(df_plot, aes(x = ps, fill = group, color = group)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 25, alpha = 0.55, position = "identity"
    ) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Overlap check: propensity-score distributions",
      subtitle = "Comparable units exist only where both histograms have mass.",
      x = "Estimated propensity score ê(X)",
      y = "Density",
      fill = NULL, color = NULL
    ) +
    soo_theme()
}

# ------------------------------------------------------------
# 3. IPW weight histogram
# ------------------------------------------------------------

plot_ipw_weights <- function(df) {
  ps_model <- glm(treat ~ x, data = df, family = binomial())
  df_plot <- df %>%
    mutate(
      ps = pmin(pmax(predict(ps_model, type = "response"), 0.02), 0.98),
      ipw = if_else(treat == 1, 1 / ps, 1 / (1 - ps)),
      group = if_else(treat == 1, "Treated", "Control")
    )

  ggplot(df_plot, aes(x = ipw, fill = group, color = group)) +
    geom_histogram(bins = 25, alpha = 0.55, position = "identity") +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "IPW weight distribution",
      subtitle = "Right tail: a few units carry a large share of the estimate.",
      x = "Inverse-propensity weight",
      y = "Count",
      fill = NULL, color = NULL
    ) +
    soo_theme()
}

# ------------------------------------------------------------
# 4. AIPW Series R — Robinson view (4 plots)
# ------------------------------------------------------------
# All four plots are built from the same residualized data so they
# tell a single story. We use a flexible polynomial as the "ML"
# learner; in production one would use mlr3/grf/etc.

.soo_residualize <- function(df, poly_deg = 3) {
  fit_y <- lm(y ~ poly(x, poly_deg), data = df)
  fit_d <- glm(treat ~ poly(x, poly_deg), data = df, family = binomial())

  df %>%
    mutate(
      y_hat   = predict(fit_y),
      d_hat   = predict(fit_d, type = "response"),
      y_resid = y - y_hat,
      d_resid = treat - d_hat,
      group   = if_else(treat == 1, "Treated", "Control")
    )
}

plot_robinson_step1 <- function(df) {
  df_plot <- df %>% mutate(group = if_else(treat == 1, "Treated", "Control"))
  group_means <- df_plot %>% group_by(group) %>%
    summarise(y_mean = mean(y), .groups = "drop")
  naive_diff <- with(df_plot, mean(y[treat == 1]) - mean(y[treat == 0]))

  ggplot(df_plot, aes(x = x, y = y, color = group, fill = group)) +
    geom_point(shape = 21, size = 2.4, alpha = 0.75) +
    geom_hline(data = group_means,
               aes(yintercept = y_mean, color = group),
               linetype = "dashed", linewidth = 0.6) +
    annotate("label",
             x = min(df$x), y = max(df$y),
             label = paste0("Naive gap = ", round(naive_diff, 2),
                            "\n(biased: X differs by group)"),
             hjust = 0, vjust = 1, size = 3.2, border.colour = NA,
             fill = alpha("white", 0.85)) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Step 1: The naive comparison",
      subtitle = "Group-mean difference, ignoring X.",
      x = "X", y = "Y", color = NULL, fill = NULL
    ) +
    soo_theme()
}

plot_robinson_step2 <- function(df, poly_deg = 3) {
  df_r <- .soo_residualize(df, poly_deg)
  x_grid <- seq(min(df_r$x), max(df_r$x), length.out = 200)
  smooth_y <- tibble(
    x = x_grid,
    y_hat = predict(lm(y ~ poly(x, poly_deg), data = df_r),
                    newdata = data.frame(x = x_grid))
  )

  ggplot(df_r, aes(x = x, y = y, color = group, fill = group)) +
    geom_segment(aes(xend = x, yend = y_hat), alpha = 0.35, linewidth = 0.3) +
    geom_point(shape = 21, size = 2.4, alpha = 0.75) +
    geom_line(data = smooth_y, aes(x = x, y = y_hat),
              inherit.aes = FALSE, color = "grey20", linewidth = 0.9) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Step 2: Predict Y from X",
      subtitle = expression(paste(
        "Fit ", hat(E) * "[Y | X]", " on all units; vertical lines are residuals ",
        tilde(Y) == Y - hat(E) * "[Y | X]", "."
      )),
      x = "X", y = "Y", color = NULL, fill = NULL
    ) +
    soo_theme()
}

plot_robinson_step3 <- function(df, poly_deg = 3) {
  df_r <- .soo_residualize(df, poly_deg)
  x_grid <- seq(min(df_r$x), max(df_r$x), length.out = 200)
  smooth_d <- tibble(
    x = x_grid,
    d_hat = predict(glm(treat ~ poly(x, poly_deg),
                        data = df_r, family = binomial()),
                    newdata = data.frame(x = x_grid),
                    type = "response")
  )

  ggplot(df_r, aes(x = x, y = treat, color = group, fill = group)) +
    geom_segment(aes(xend = x, yend = d_hat), alpha = 0.35, linewidth = 0.3) +
    geom_point(shape = 21, size = 2.4, alpha = 0.75) +
    geom_line(data = smooth_d, aes(x = x, y = d_hat),
              inherit.aes = FALSE, color = "grey20", linewidth = 0.9) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Step 3: Predict D from X (propensity)",
      subtitle = expression(paste(
        "Fit ", hat(e)(X) == hat(E) * "[D | X]", "; vertical lines are residuals ",
        tilde(D) == D - hat(e)(X), "."
      )),
      x = "X", y = "D (0/1)", color = NULL, fill = NULL
    ) +
    soo_theme()
}

plot_robinson_step4 <- function(df, poly_deg = 3) {
  df_r <- .soo_residualize(df, poly_deg)
  fit_resid <- lm(y_resid ~ d_resid, data = df_r)
  slope <- coef(fit_resid)["d_resid"]

  ggplot(df_r, aes(x = d_resid, y = y_resid, color = group, fill = group)) +
    geom_point(shape = 21, size = 2.4, alpha = 0.75) +
    geom_smooth(aes(x = d_resid, y = y_resid),
                inherit.aes = FALSE,
                method = "lm", formula = y ~ x,
                color = "grey20", se = FALSE, linewidth = 0.9) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    annotate("label",
             x = min(df_r$d_resid), y = max(df_r$y_resid),
             label = paste0("Residualized slope = ", round(slope, 2),
                            "\n(estimate of the ATE)"),
             hjust = 0, vjust = 1, size = 3.2, border.colour = NA,
             fill = alpha("white", 0.85)) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Step 4: Regress residuals on residuals",
      subtitle = expression(paste(
        "Slope of ", tilde(Y), " on ", tilde(D), " is the partialled-out treatment effect."
      )),
      x = expression(tilde(D) == D - hat(e)(X)),
      y = expression(tilde(Y) == Y - hat(E) * "[Y | X]"),
      color = NULL, fill = NULL
    ) +
    soo_theme()
}

# ------------------------------------------------------------
# 5. AIPW Series A — correction view (4 plots)
# ------------------------------------------------------------
# Built from arm-specific outcome models and a propensity model.

.soo_aipw_pieces <- function(df, poly_deg = 2) {
  fit_c <- lm(y ~ poly(x, poly_deg), data = df %>% filter(treat == 0))
  fit_t <- lm(y ~ poly(x, poly_deg), data = df %>% filter(treat == 1))
  ps_mod <- glm(treat ~ poly(x, poly_deg), data = df, family = binomial())

  df %>%
    mutate(
      mu0 = predict(fit_c, newdata = .),
      mu1 = predict(fit_t, newdata = .),
      ps  = pmin(pmax(predict(ps_mod, newdata = ., type = "response"), 0.02), 0.98),
      ipw = if_else(treat == 1, 1 / ps, 1 / (1 - ps)),
      mu_own = if_else(treat == 1, mu1, mu0),
      resid  = y - mu_own,
      group  = if_else(treat == 1, "Treated", "Control")
    )
}

plot_aipw_step1 <- function(df, poly_deg = 2) {
  d <- .soo_aipw_pieces(df, poly_deg)
  x_grid <- seq(min(d$x), max(d$x), length.out = 200)
  fit_c <- lm(y ~ poly(x, poly_deg), data = d %>% filter(treat == 0))
  fit_t <- lm(y ~ poly(x, poly_deg), data = d %>% filter(treat == 1))
  curves <- bind_rows(
    tibble(x = x_grid, y = predict(fit_c, newdata = data.frame(x = x_grid)),
           model = "control"),
    tibble(x = x_grid, y = predict(fit_t, newdata = data.frame(x = x_grid)),
           model = "treated")
  )
  reg_est <- mean(d$mu1 - d$mu0)

  ggplot(d, aes(x = x, y = y)) +
    geom_point(aes(fill = group, color = group),
               shape = 21, size = 2.4, alpha = 0.6) +
    geom_line(data = curves, aes(x = x, y = y, linetype = model),
              color = "grey20", linewidth = 0.9, inherit.aes = FALSE) +
    annotate("label",
             x = min(d$x), y = max(d$y),
             label = paste0("Regression-imputation\nestimate = ", round(reg_est, 2)),
             hjust = 0, vjust = 1, size = 3.2, border.colour = NA,
             fill = alpha("white", 0.85)) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    scale_linetype_manual(
      values = c(control = "dashed", treated = "solid"),
      labels = c(control = expression(hat(mu)[0](x) * ": control model"),
                 treated = expression(hat(mu)[1](x) * ": treated model"))
    ) +
    labs(
      title = "Step 1: Outcome-regression imputation",
      subtitle = "Two arm-specific fits; the gap between curves is the regression estimate.",
      x = "X", y = "Y", color = NULL, fill = NULL, linetype = NULL
    ) +
    soo_theme()
}

plot_aipw_step2 <- function(df, poly_deg = 2) {
  d <- .soo_aipw_pieces(df, poly_deg)

  ggplot(d, aes(x = x, y = resid, color = group, fill = group)) +
    geom_segment(aes(xend = x, yend = 0), alpha = 0.35, linewidth = 0.3) +
    geom_point(shape = 21, size = 2.4, alpha = 0.75) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    labs(
      title = "Step 2: Residuals from the own-arm regression",
      subtitle = expression(paste(
        hat(epsilon)[i] == Y[i] - hat(mu)[D[i]](X[i]), ". These are what AIPW corrects for."
      )),
      x = "X",
      y = expression("Residual " * hat(epsilon)),
      color = NULL, fill = NULL
    ) +
    soo_theme()
}

plot_aipw_step3 <- function(df, poly_deg = 2) {
  d <- .soo_aipw_pieces(df, poly_deg)

  ggplot(d, aes(x = x, y = resid, color = group, fill = group)) +
    geom_segment(aes(xend = x, yend = 0), alpha = 0.30, linewidth = 0.3) +
    geom_point(aes(size = ipw), shape = 21, alpha = 0.7) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = soo_fills) +
    scale_color_manual(values = soo_colors) +
    scale_size_area(max_size = 8, name = "IPW weight") +
    labs(
      title = "Step 3: Up-weight residuals where overlap is thin",
      subtitle = expression(paste(
        "Point size " %prop% " 1/", hat(e)(X), " for treated, 1/(1 - ", hat(e)(X), ") for controls."
      )),
      x = "X",
      y = expression("Residual " * hat(epsilon)),
      color = NULL, fill = NULL
    ) +
    soo_theme()
}

plot_aipw_step4 <- function(df, poly_deg = 2) {
  d <- .soo_aipw_pieces(df, poly_deg)

  reg_term <- mean(d$mu1 - d$mu0)
  correction <- mean(
    d$treat * (d$y - d$mu1) / d$ps -
      (1 - d$treat) * (d$y - d$mu0) / (1 - d$ps)
  )
  aipw <- reg_term + correction

  bars <- tibble(
    component = factor(
      c("Regression\nimputation",
        "Propensity-weighted\nresidual correction",
        "AIPW\nestimate"),
      levels = c("Regression\nimputation",
                 "Propensity-weighted\nresidual correction",
                 "AIPW\nestimate")
    ),
    value = c(reg_term, correction, aipw),
    role  = c("term", "term", "total")
  )

  ggplot(bars, aes(x = component, y = value, fill = role)) +
    geom_col(width = 0.55, color = "grey20") +
    geom_text(aes(label = round(value, 2)),
              vjust = -0.4, size = 3.5) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    scale_fill_manual(values = c(term = "#B5D4F4", total = "#185FA5"),
                      guide = "none") +
    labs(
      title = "Step 4: AIPW = regression term + propensity-weighted correction",
      subtitle = "If the regression model is wrong, the correction term debiases the estimate.",
      x = NULL, y = "Component value"
    ) +
    coord_cartesian(ylim = c(min(bars$value) - 0.2, max(bars$value) + 0.4)) +
    soo_theme()
}

# ------------------------------------------------------------
# 6. Balance plot (standardized mean differences before/after)
# ------------------------------------------------------------
# Inputs: a named list of data frames, each with `treat` and the
# covariates to check. Useful for comparing balance before/after
# matching or weighting.

compute_smd <- function(d, treat_col = "treat", vars = NULL) {
  if (is.null(vars)) vars <- setdiff(names(d), treat_col)

  bind_rows(lapply(vars, function(v) {
    x_t <- d[[v]][d[[treat_col]] == 1]
    x_c <- d[[v]][d[[treat_col]] == 0]
    pooled_sd <- sqrt((var(x_t) + var(x_c)) / 2)
    tibble(
      variable = v,
      smd = (mean(x_t) - mean(x_c)) / pooled_sd
    )
  }))
}

plot_balance <- function(data_list, treat_col = "treat", vars = NULL) {
  smd_long <- bind_rows(lapply(names(data_list), function(label) {
    compute_smd(data_list[[label]], treat_col, vars) %>%
      mutate(sample = label)
  }))

  ggplot(smd_long, aes(x = abs(smd), y = variable, color = sample, shape = sample)) +
    geom_vline(xintercept = c(0.1, 0.25), linetype = "dashed", color = "grey60") +
    geom_point(size = 3.2) +
    labs(
      title = "Covariate balance: |standardized mean difference|",
      subtitle = "Dashed lines mark common 0.10 and 0.25 thresholds.",
      x = "|SMD|", y = NULL, color = NULL, shape = NULL
    ) +
    soo_theme()
}
