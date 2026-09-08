# Helpers for lecture 10 (mediation and mechanisms). Sourced by the lecture.
# All figure text is ASCII.

library(ggplot2)
library(dplyr)
library(tibble)
library(dagitty)
library(ggdag)
library(mlr3)
library(mlr3learners)

med_colors <- c(
  "Total" = "#000000",
  "Direct" = "#0072B2",
  "Indirect" = "#D55E00",
  "Naive" = "#999999",
  "Truth" = "#000000",
  "Regression" = "#E69F00",
  "DML linear" = "#56B4E9",
  "DML forest" = "#009E73"
)

med_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# ---------------------------------------------------------------------------
# DAGs. `draw_dag()` takes a dagitty object with coordinates, draws observed
# nodes as filled circles and the nodes in `latent` as hollow dashed circles,
# and labels the edges in `edge_labels` (a named vector "from->to" = label).
# ---------------------------------------------------------------------------
draw_dag <- function(dag, latent = NULL, edge_labels = NULL, title = NULL, node_size = 11,
                     text_size = 3.6, expand = 0.35) {
  td <- tidy_dagitty(dag)
  df <- as.data.frame(td)
  df$latent <- df$name %in% latent
  nodes <- df[!duplicated(df$name), c("name", "x", "y", "latent")]
  edges <- df[!is.na(df$to), ]
  p <- ggplot(edges, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges(edge_width = 0.6, arrow_directed = grid::arrow(length = grid::unit(7, "pt"), type = "closed")) +
    geom_point(data = nodes[!nodes$latent, ], aes(x = x, y = y), inherit.aes = FALSE, size = node_size, shape = 21,
               fill = "white", colour = "black", stroke = 0.8) +
    geom_point(data = nodes[nodes$latent, ], aes(x = x, y = y), inherit.aes = FALSE, size = node_size, shape = 21,
               fill = "grey92", colour = "grey40", stroke = 0.8) +
    geom_text(data = nodes, aes(x = x, y = y, label = name), inherit.aes = FALSE, size = text_size) +
    coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(add = expand)) +
    scale_y_continuous(expand = expansion(add = expand)) +
    theme_dag() +
    theme(plot.title = element_text(size = 11, hjust = 0.5))
  if (!is.null(edge_labels)) {
    lab <- edges %>%
      mutate(key = paste0(name, "->", to)) %>%
      filter(key %in% names(edge_labels)) %>%
      mutate(label = edge_labels[key], mx = (x + xend) / 2, my = (y + yend) / 2)
    p <- p + geom_label(data = lab, aes(x = mx, y = my, label = label), inherit.aes = FALSE, size = 3.2,
                        linewidth = 0, fill = "white", label.padding = unit(0.1, "lines"))
  }
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

# The graphs used in the lecture.
dag_mediation <- function() {
  dagify(Y ~ D + M + X, M ~ D + X, D ~ X,
         coords = list(x = c(D = 0, M = 1, Y = 2, X = 1), y = c(D = 0, M = 1, Y = 0, X = -1)))
}
dag_mediation_u <- function() {
  dagify(Y ~ D + M + X + U, M ~ D + X + U, D ~ X,
         coords = list(x = c(D = 0, M = 1, Y = 2, X = 1, U = 2.2), y = c(D = 0, M = 1, Y = 0, X = -1, U = 1.1)))
}
dag_post_treatment <- function() {
  dagify(Y ~ D + M + L + X, M ~ D + L + X, L ~ D + X, D ~ X,
         coords = list(x = c(D = 0, L = 0.8, M = 1.4, Y = 2.2, X = 1.1), y = c(D = 0, L = 1, M = 0.6, Y = 0, X = -1)))
}
dag_iv_mediator <- function() {
  dagify(Y ~ D + M + U, M ~ D + Z + U,
         coords = list(x = c(D = 0, M = 1, Y = 2, Z = 1, U = 2.2), y = c(D = 0, M = 1, Y = 0, Z = 2, U = 1.1)))
}
dag_front_door <- function() {
  dagify(Y ~ M + U, M ~ D, D ~ U,
         coords = list(x = c(D = 0, M = 1, Y = 2, U = 1), y = c(D = 0, M = 0, Y = 0, U = 1)))
}
dag_parallel <- function() {
  dagify(Y ~ D + M1 + M2, M1 ~ D, M2 ~ D,
         coords = list(x = c(D = 0, M1 = 1, M2 = 1, Y = 2), y = c(D = 0, M1 = 0.8, M2 = -0.8, Y = 0)))
}
dag_collider <- function() {
  dagify(Y ~ U, M ~ D + U,
         coords = list(x = c(D = 0, M = 1, Y = 2, U = 1.6), y = c(D = 0, M = 0, Y = 0, U = 1)))
}

# ---------------------------------------------------------------------------
# Section 1: the post-treatment trap. Repeated samples from the linear
# design with an unmeasured M-Y confounder U of strength `u_strength`;
# returns the coefficient on D in Y ~ D (total) and in Y ~ D + M (naive
# direct) together with the truth.
# ---------------------------------------------------------------------------
sim_bad_control <- function(n = 2000, u_strength = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  u <- rnorm(n)
  m <- 0.5 + 1.2 * d + 0.5 * x + u_strength * u + rnorm(n)
  y <- 1 + 1.0 * d + 0.9 * m + 0.5 * x + u_strength * u + rnorm(n)
  data.frame(y = y, d = d, m = m, x = x, u = u)
}

mc_bad_control <- function(n_sims = 200, n = 2000, u_strength = c(0, 0.5, 1, 1.5), seed = 1) {
  set.seed(seed)
  out <- list()
  for (us in u_strength) {
    for (s in seq_len(n_sims)) {
      dat <- sim_bad_control(n, us)
      total <- coef(lm(y ~ d + x, data = dat))[["d"]]
      naive <- coef(lm(y ~ d + m + x, data = dat))[["d"]]
      oracle <- coef(lm(y ~ d + m + x + u, data = dat))[["d"]]
      out[[length(out) + 1]] <- tibble(u_strength = us, sim = s, total = total, naive = naive, oracle = oracle)
    }
  }
  bind_rows(out) %>%
    tidyr::pivot_longer(c(total, naive, oracle), names_to = "coefficient", values_to = "estimate") %>%
    mutate(coefficient = factor(coefficient, levels = c("total", "naive", "oracle"),
                                labels = c("Y ~ D + X (total)", "Y ~ D + M + X (naive direct)", "Y ~ D + M + X + U (oracle)")))
}

# ---------------------------------------------------------------------------
# Section 7: regression versus double machine learning in the nonlinear
# design. Bias, standard deviation, RMSE, and coverage of the NDE and NIE
# for mediate_reg (linear, with interaction) and mediate_dml with linear
# learners and with random forests.
# ---------------------------------------------------------------------------
mc_reg_vs_dml <- function(n_sims = 100, n = 2000, seed = 1) {
  set.seed(seed)
  rf <- mlr3::lrn("regr.ranger", num.trees = 200, min.node.size = 10)
  rfc <- mlr3::lrn("classif.ranger", num.trees = 200, min.node.size = 10, predict_type = "prob")
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- causalmetrics::sim_mediation(n, dgp = "nonlinear", seed = seed * 1000 + s)
    tr <- attr(dat, "truth")
    fits <- list(
      "Regression" = causalmetrics::mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), interaction = TRUE,
                                                method = "simulation", n_sim = 200, seed = s),
      "DML linear" = causalmetrics::mediate_dml(dat, "y", "d", "m", x = c("x1", "x2"), seed = s),
      "DML forest" = causalmetrics::mediate_dml(dat, "y", "d", "m", x = c("x1", "x2"), learner_y = rf, learner_d = rfc, seed = s)
    )
    for (nm in names(fits)) {
      e <- fits[[nm]]$effects
      for (term in c("nde", "nie")) {
        r <- e[e$term == term, ]
        rows[[length(rows) + 1]] <- tibble(sim = s, estimator = nm, term = term, estimate = r$estimate,
                                           std.error = r$std.error, truth = tr[[term]],
                                           covered = r$conf.low <= tr[[term]] & tr[[term]] <= r$conf.high)
      }
    }
  }
  bind_rows(rows)
}

summarise_mc <- function(mc) {
  mc %>%
    group_by(estimator, term) %>%
    summarise(truth = mean(truth), bias = mean(estimate - truth), sd = sd(estimate),
              rmse = sqrt(mean((estimate - truth)^2)), coverage = mean(covered), .groups = "drop")
}

# ---------------------------------------------------------------------------
# Section 9: back door versus front door across repeated samples.
# ---------------------------------------------------------------------------
mc_front_door <- function(n_sims = 200, n = 2000, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- causalmetrics::sim_mediation(n, dgp = "front_door", seed = seed * 1000 + s)
    fd <- causalmetrics::front_door(dat, "y", "d", "m", x = c("x1", "x2"), n_boot = 0)
    rows[[length(rows) + 1]] <- tibble(sim = s, truth = attr(dat, "truth")$total,
                                       front_door = fd$effects$estimate[fd$effects$term == "total"],
                                       backdoor = fd$effects$estimate[fd$effects$term == "backdoor"])
  }
  bind_rows(rows)
}


# Section 8: two mediators whose errors are correlated (shared component
# `share`), so that sequential addition attributes the coefficient change
# differently in different orders. Indirect effects are 1.0 * 0.8 = 0.8
# through m1 and 0.4 * 1.5 = 0.6 through m2; the direct effect is 0.5.
sim_correlated_mediators <- function(n = 3000, share = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  d <- rbinom(n, 1, plogis(0.4 * x1 - 0.3 * x2))
  e1 <- rnorm(n); e2 <- rnorm(n)
  m1 <- 0.3 + 1.0 * d + 0.5 * x1 + e1
  m2 <- -0.2 + 0.4 * d - 0.3 * x2 + share * e1 + e2
  y <- 1 + 0.5 * d + 0.8 * m1 + 1.5 * m2 + 0.4 * x1 + 0.2 * x2 + rnorm(n)
  data.frame(y = y, d = d, m1 = m1, m2 = m2, x1 = x1, x2 = x2)
}

# Sequential addition of mediators: the treatment coefficient after each
# step, in the order given, for the "add and watch" figure of Section 8.
sequential_coefficients <- function(data, y, d, x_base, order) {
  rhs <- x_base
  out <- tibble(step = 0L, added = "(base)", coefficient = coef(lm(reformulate(c(d, rhs), y), data = data))[[d]])
  for (k in seq_along(order)) {
    rhs <- c(rhs, order[k])
    out <- bind_rows(out, tibble(step = k, added = order[k],
                                 coefficient = coef(lm(reformulate(c(d, rhs), y), data = data))[[d]]))
  }
  out
}

# Format helpers for inline numbers.
fmt <- function(x, digits = 2) formatC(x, digits = digits, format = "f")
fmt_pct <- function(x, digits = 0) paste0(formatC(100 * x, digits = digits, format = "f"), "%")
