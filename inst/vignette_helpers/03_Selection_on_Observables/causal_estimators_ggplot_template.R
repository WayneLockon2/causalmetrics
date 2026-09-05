# causal_estimators_ggplot_template.R
# ------------------------------------------------------------
# Reusable ggplot2 template for visualizing causal estimators.
#
# Methods included:
#   1. Outcome regression / linear regression
#   2. Propensity score matching
#   3. Inverse propensity weighting (IPW)
#   4. Coarsened exact matching (CEM)
#   5. Augmented inverse propensity weighting (AIPW)
#
# The goal is to mimic the logic of an interactive JavaScript/SVG figure,
# but in an RMarkdown-friendly ggplot2 workflow.
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# 1. Simulate a treatment/control dataset
# ------------------------------------------------------------

simulate_causal_data <- function(
  n_control = 80,
  n_treated = 60,
  true_ate = 1.5,
  seed = 123
) {
  set.seed(seed)

  x_control <- pmin(pmax(rnorm(n_control, mean = 42, sd = 18), 0), 100)
  x_treated <- pmin(pmax(rnorm(n_treated, mean = 65, sd = 16), 0), 100)

  bind_rows(
    tibble(id = paste0("C", seq_along(x_control)), treat = 0, x = x_control),
    tibble(id = paste0("T", seq_along(x_treated)), treat = 1, x = x_treated)
  ) %>%
    mutate(
      group = if_else(treat == 1, "Treated", "Control"),
      y0 = 2 + 0.045 * x + 0.00025 * (x - 50)^2 + rnorm(n(), 0, 0.28),
      y = y0 + true_ate * treat
    )
}

# ------------------------------------------------------------
# 2. Propensity score helper
# ------------------------------------------------------------

add_propensity <- function(df) {
  ps_model <- glm(treat ~ x, data = df, family = binomial())

  df %>%
    mutate(
      ps = predict(ps_model, type = "response"),
      # Trim extreme propensities for numerical stability and readable plots.
      ps = pmin(pmax(ps, 0.02), 0.98)
    )
}

# ------------------------------------------------------------
# 3. Fitted line helper
# ------------------------------------------------------------
# If weight_col is NULL, this fits ordinary least squares lines separately
# for treated and control groups.
# If weight_col is supplied, this fits weighted least squares lines.
# This is useful for IPW and AIPW visualization because the weighted fitted
# line can differ from the unweighted regression line.

fit_group_lines <- function(df, weight_col = NULL, x_grid_length = 100) {
  x_all <- range(df$x, na.rm = TRUE)

  make_one_group <- function(d, group_name) {
    w <- if (is.null(weight_col)) NULL else d[[weight_col]]

    fit <- if (is.null(w)) {
      lm(y ~ x, data = d)
    } else {
      lm(y ~ x, data = d, weights = w)
    }

    support <- range(d$x, na.rm = TRUE)

    make_piece <- function(xmin, xmax, line_type, piece_name) {
      if (xmax <= xmin) return(tibble())

      x_grid <- seq(xmin, xmax, length.out = x_grid_length)

      tibble(
        group = group_name,
        x = x_grid,
        y = predict(fit, newdata = data.frame(x = x_grid)),
        line_type = line_type,
        line_group = paste(group_name, piece_name)
      )
    }

    bind_rows(
      make_piece(x_all[1], support[1], "extrapolated", "left"),
      make_piece(support[1], support[2], "observed support", "middle"),
      make_piece(support[2], x_all[2], "extrapolated", "right")
    )
  }

  bind_rows(
    make_one_group(df %>% filter(treat == 0), "Control"),
    make_one_group(df %>% filter(treat == 1), "Treated")
  )
}

# ------------------------------------------------------------
# 4. Main ggplot template
# ------------------------------------------------------------
# The plotting function takes two objects:
#   df   = unit-level data with id, treat, x, y
#   spec = estimator-specific visual information
#
# A spec can contain:
#   title         : plot title
#   ate           : estimated effect
#   lines         : fitted lines
#   pairs         : matched treatment-control links
#   bins          : CEM bins / common support rectangles
#   bin_effects   : within-bin treatment-control mean differences
#   point_weights : unit weights for point size
#   effect        : one vertical effect annotation
#   dim_ids       : units to fade out

plot_causal_template <- function(
  df,
  spec,
  xlab = "Covariate X",
  ylab = "Outcome Y"
) {
  df_plot <- df %>%
    mutate(
      group = if_else(treat == 1, "Treated", "Control"),
      point_alpha = ifelse(
        !is.null(spec$dim_ids) & id %in% spec$dim_ids,
        0.18,
        1
      )
    )

  if (!is.null(spec$point_weights)) {
    df_plot <- df_plot %>%
      left_join(spec$point_weights, by = "id") %>%
      mutate(plot_weight = ifelse(is.na(plot_weight), 1, plot_weight))
  }

  p <- ggplot()

  # CEM-style bins: active bins are shaded blue, inactive bins are grey/dashed.
  if (!is.null(spec$bins)) {
    active_bins <- spec$bins %>% filter(active)
    dropped_bins <- spec$bins %>% filter(!active)

    p <- p +
      geom_rect(
        data = active_bins,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "#E6F1FB",
        color = "#378ADD",
        alpha = 0.25,
        linewidth = 0.3
      ) +
      geom_rect(
        data = dropped_bins,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "grey85",
        color = "grey60",
        alpha = 0.35,
        linetype = "dashed",
        linewidth = 0.3
      )
  }

  # CEM within-bin contrasts.
  # These are not one-to-one matches. The bin is the matched stratum.
  if (!is.null(spec$bin_effects)) {
    p <- p +
      geom_segment(
        data = spec$bin_effects,
        aes(x = x_mid, xend = x_mid, y = y_control, yend = y_treated),
        inherit.aes = FALSE,
        linewidth = 0.8,
        linetype = "dotted",
        color = "grey25"
      ) +
      geom_point(
        data = spec$bin_effects,
        aes(x = x_mid, y = y_control),
        inherit.aes = FALSE,
        shape = 21,
        size = 3,
        fill = "#F5C4B3",
        color = "#993C1D"
      ) +
      geom_point(
        data = spec$bin_effects,
        aes(x = x_mid, y = y_treated),
        inherit.aes = FALSE,
        shape = 21,
        size = 3,
        fill = "#B5D4F4",
        color = "#185FA5"
      ) +
      geom_label(
        data = spec$bin_effects,
        aes(x = x_mid, y = (y_treated + y_control) / 2, label = label),
        parse = TRUE,
        inherit.aes = FALSE,
        size = 3,
        border.colour = NA,
        alpha = 0.9
      )
  }

  # PSM links.
  if (!is.null(spec$pairs)) {
    p <- p +
      geom_segment(
        data = spec$pairs,
        aes(x = x_t, y = y_t, xend = x_c, yend = y_c),
        inherit.aes = FALSE,
        color = "grey45",
        alpha = 0.5,
        linewidth = 0.35
      )
  }

  # Main unit-level points.
  # If point weights are supplied, point size represents weight.
  if (!is.null(spec$point_weights)) {
    p <- p +
      geom_point(
        data = df_plot,
        aes(
          x = x,
          y = y,
          fill = group,
          color = group,
          alpha = point_alpha,
          size = plot_weight
        ),
        shape = 21,
        stroke = 0.6
      )
  } else {
    p <- p +
      geom_point(
        data = df_plot,
        aes(x = x, y = y, fill = group, color = group, alpha = point_alpha),
        shape = 21,
        size = 2.6,
        stroke = 0.6
      )
  }

  # Fitted lines.
  if (!is.null(spec$lines)) {
    p <- p +
      geom_line(
        data = spec$lines,
        aes(
          x = x,
          y = y,
          color = group,
          linetype = line_type,
          group = line_group
        ),
        linewidth = 0.9
      )
  }

  # Optional effect annotation.
  if (!is.null(spec$effect)) {
    p <- p +
      geom_segment(
        data = spec$effect,
        aes(x = x, xend = x, y = y0, yend = y1),
        inherit.aes = FALSE,
        linetype = "dotted",
        linewidth = 0.5
      ) +
      geom_label(
        data = spec$effect,
        aes(x = x, y = (y0 + y1) / 2, label = label),
        inherit.aes = FALSE,
        nudge_x = 4,
        size = 3,
        border.colour = NA
      )
  }

  p +
    scale_fill_manual(values = c(Control = "#F5C4B3", Treated = "#B5D4F4")) +
    scale_color_manual(values = c(Control = "#993C1D", Treated = "#185FA5")) +
    scale_alpha_identity(guide = "none") +
    scale_size_area(max_size = 10, name = "Weight") +
    # Only panels with fitted lines map linetype; adding the scale elsewhere
    # triggers a "no shared levels" warning in ggplot2 >= 4.0.
    (if (!is.null(spec$lines)) {
      scale_linetype_manual(
        values = c(
          "observed support" = "solid",
          "extrapolated" = "dashed"
        ),
        drop = FALSE
      )
    }) +
    labs(
      title = spec$title,
      subtitle = paste0("Estimated effect = ", round(spec$ate, 2)),
      x = xlab,
      y = ylab,
      fill = NULL,
      color = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey35")
    )
}

# ------------------------------------------------------------
# 5. Estimator specs
# ------------------------------------------------------------

# 5.1 Outcome regression / linear regression

build_lm_spec <- function(df, x_ref = median(df$x, na.rm = TRUE)) {
  fit_c <- lm(y ~ x, data = df %>% filter(treat == 0))
  fit_t <- lm(y ~ x, data = df %>% filter(treat == 1))

  y0_ref <- predict(fit_c, newdata = data.frame(x = x_ref))
  y1_ref <- predict(fit_t, newdata = data.frame(x = x_ref))

  ate <- mean(
    predict(fit_t, newdata = df) -
      predict(fit_c, newdata = df)
  )

  list(
    title = "Outcome regression / linear model",
    ate = ate,
    lines = fit_group_lines(df),
    effect = tibble(
      x = x_ref,
      y0 = as.numeric(y0_ref),
      y1 = as.numeric(y1_ref),
      label = paste0("gap = ", round(y1_ref - y0_ref, 2))
    )
  )
}

# 5.2 Propensity score matching
# This implements nearest-neighbor matching with replacement for ATT.

build_psm_spec <- function(df) {
  d <- add_propensity(df)

  treated <- d %>% filter(treat == 1)
  control <- d %>% filter(treat == 0)

  nearest_control <- sapply(
    treated$ps,
    function(p) which.min(abs(control$ps - p))
  )

  pairs <- tibble(
    id_t = treated$id,
    id_c = control$id[nearest_control],
    x_t = treated$x,
    y_t = treated$y,
    x_c = control$x[nearest_control],
    y_c = control$y[nearest_control]
  )

  matched_controls <- unique(pairs$id_c)
  dim_ids <- setdiff(control$id, matched_controls)

  list(
    title = "Propensity score matching: ATT",
    ate = mean(pairs$y_t - pairs$y_c),
    pairs = pairs,
    dim_ids = dim_ids
  )
}

# 5.3 Inverse propensity weighting
# Point size = IPW weight.
# Lines = weighted least-squares lines, not ordinary regression lines.

build_ipw_spec <- function(df) {
  d <- add_propensity(df) %>%
    mutate(
      ipw = if_else(treat == 1, 1 / ps, 1 / (1 - ps))
    )

  # Cap plotted weights for readability. The estimator still uses uncapped ipw.
  cap <- quantile(d$ipw, 0.95, na.rm = TRUE)
  d <- d %>%
    mutate(ipw_plot = pmin(ipw, cap))

  # Hájek-style weighted means.
  mu1 <- with(d, sum((treat == 1) * ipw * y) / sum((treat == 1) * ipw))
  mu0 <- with(d, sum((treat == 0) * ipw * y) / sum((treat == 0) * ipw))

  list(
    title = "Inverse propensity weighting",
    ate = mu1 - mu0,
    point_weights = d %>% transmute(id, plot_weight = ipw_plot),
    lines = fit_group_lines(d, weight_col = "ipw")
  )
}

# 5.4 Coarsened exact matching
# In standard CEM, the bin/stratum is the match.
# There is no required one-to-one matching inside the bin.
# The effect is the weighted average of within-bin treated-control mean differences.

assign_cem_bin <- function(x, breaks) {
  idx <- findInterval(x, breaks, rightmost.closed = TRUE)
  pmax(1, pmin(idx, length(breaks) - 1))
}

build_cem_spec <- function(
  df,
  breaks = seq(0, 100, by = 20),
  estimand = c("ATT", "ATE", "ATC"),
  use_regression = FALSE
) {
  estimand <- match.arg(estimand)

  d <- df %>%
    mutate(
      bin_id = cut(
        x,
        breaks = breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
      )
    )

  bins <- tibble(
    bin_id = seq_len(length(breaks) - 1),
    xmin = breaks[-length(breaks)],
    xmax = breaks[-1]
  )

  bin_counts <- d %>%
    group_by(bin_id) %>%
    summarise(
      n_treated = sum(treat == 1),
      n_control = sum(treat == 0),
      active = n_treated > 0 & n_control > 0,
      .groups = "drop"
    )

  bins <- bins %>%
    left_join(bin_counts, by = "bin_id") %>%
    mutate(
      n_treated = ifelse(is.na(n_treated), 0, n_treated),
      n_control = ifelse(is.na(n_control), 0, n_control),
      active = ifelse(is.na(active), FALSE, active)
    )

  kept_bins <- bins %>%
    filter(active) %>%
    pull(bin_id)

  kept <- d %>%
    filter(bin_id %in% kept_bins)

  bin_effects <- kept %>%
    group_by(bin_id) %>%
    group_modify(~ {
      bin_data <- .x

      if (!use_regression) {
        # Pure CEM: mean difference.
        y0 <- mean(bin_data$y[bin_data$treat == 0])
        y1 <- mean(bin_data$y[bin_data$treat == 1])
        tau <- y1 - y0
      } else {
        # Equivalent to mean difference when formula is y ~ treat.
        fit <- lm(y ~ treat, data = bin_data)
        y0 <- coef(fit)["(Intercept)"]
        tau <- coef(fit)["treat"]
        y1 <- y0 + tau
      }

      tibble(
        y_control = as.numeric(y0),
        y_treated = as.numeric(y1),
        tau = as.numeric(tau),
        n_treated = sum(bin_data$treat == 1),
        n_control = sum(bin_data$treat == 0)
      )
    }) %>%
    ungroup() %>%
    left_join(bins %>% select(bin_id, xmin, xmax), by = "bin_id") %>%
    mutate(
      x_mid = (xmin + xmax) / 2,
      bin_weight = switch(
        estimand,
        ATT = n_treated,
        ATE = n_treated + n_control,
        ATC = n_control,
        stop("`estimand` must be one of 'ATT', 'ATE', or 'ATC'.")
      ),
      label = paste0("tau == ", round(tau, 2))
    )

  ate <- weighted.mean(bin_effects$tau, bin_effects$bin_weight)

  list(
    title = paste0("Coarsened exact matching: ", estimand),
    ate = ate,
    bins = bins,
    bin_effects = bin_effects,
    dim_ids = setdiff(df$id, kept$id)
  )
}

# 5.5 Augmented inverse propensity weighting
# ATE = AIPW estimating equation.
# Point size = IPW weight.
# Lines = weighted least-squares lines for visualization.

build_aipw_spec <- function(df) {
  d <- add_propensity(df) %>%
    mutate(
      ipw = if_else(treat == 1, 1 / ps, 1 / (1 - ps))
    )

  # Outcome models for the AIPW estimator.
  fit_c <- lm(y ~ x, data = d %>% filter(treat == 0))
  fit_t <- lm(y ~ x, data = d %>% filter(treat == 1))

  m0_hat <- predict(fit_c, newdata = d)
  m1_hat <- predict(fit_t, newdata = d)

  d <- d %>%
    mutate(
      m0 = m0_hat,
      m1 = m1_hat,
      aipw_score = (m1 - m0) +
        treat * (y - m1) / ps -
        (1 - treat) * (y - m0) / (1 - ps)
    )

  # Cap plotted weights for readability. The estimator still uses uncapped ipw.
  cap <- quantile(d$ipw, 0.95, na.rm = TRUE)
  d <- d %>%
    mutate(ipw_plot = pmin(ipw, cap))

  # Weighted lines for display.
  weighted_lines <- fit_group_lines(d, weight_col = "ipw")

  # Optional visual gap at median x using weighted fits.
  x_ref <- median(d$x, na.rm = TRUE)

  fit_c_w <- lm(y ~ x, data = d %>% filter(treat == 0), weights = ipw)
  fit_t_w <- lm(y ~ x, data = d %>% filter(treat == 1), weights = ipw)

  y0_ref <- predict(fit_c_w, newdata = data.frame(x = x_ref))
  y1_ref <- predict(fit_t_w, newdata = data.frame(x = x_ref))

  list(
    title = "Augmented IPW / doubly robust",
    ate = mean(d$aipw_score),
    point_weights = d %>% transmute(id, plot_weight = ipw_plot),
    lines = weighted_lines,
    effect = tibble(
      x = x_ref,
      y0 = as.numeric(y0_ref),
      y1 = as.numeric(y1_ref),
      label = paste0("weighted fit gap = ", round(y1_ref - y0_ref, 2))
    )
  )
}

# ------------------------------------------------------------
# 6. Convenience wrappers
# ------------------------------------------------------------

build_all_specs <- function(df, cem_estimand = "ATT") {
  list(
    lm = build_lm_spec(df),
    psm = build_psm_spec(df),
    ipw = build_ipw_spec(df),
    cem = build_cem_spec(df, estimand = cem_estimand),
    aipw = build_aipw_spec(df)
  )
}

plot_all_estimators <- function(df, cem_estimand = "ATT") {
  specs <- build_all_specs(df, cem_estimand = cem_estimand)
  lapply(specs, function(spec) plot_causal_template(df, spec))
}

# ------------------------------------------------------------
# 7. Example usage
# ------------------------------------------------------------
# Uncomment this section in an interactive R session or RMarkdown document.
#
# df <- simulate_causal_data(seed = 123)
#
# spec_lm   <- build_lm_spec(df)
# spec_psm  <- build_psm_spec(df)
# spec_ipw  <- build_ipw_spec(df)
# spec_cem  <- build_cem_spec(df, estimand = "ATT")
# spec_aipw <- build_aipw_spec(df)
#
# plot_causal_template(df, spec_lm)
# plot_causal_template(df, spec_psm)
# plot_causal_template(df, spec_ipw)
# plot_causal_template(df, spec_cem)
# plot_causal_template(df, spec_aipw)
#
# If you use patchwork, you can combine the plots:
#
# library(patchwork)
# plots <- plot_all_estimators(df, cem_estimand = "ATT")
# wrap_plots(plots, ncol = 2)
#
# In RMarkdown, you can source this file:
#
# source("causal_estimators_ggplot_template.R")
# df <- simulate_causal_data()
# plot_causal_template(df, build_ipw_spec(df))
# ------------------------------------------------------------
