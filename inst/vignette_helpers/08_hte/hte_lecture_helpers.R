# Helpers for lecture 08 (heterogeneous treatment effects and policy
# learning). Sourced by the lecture. All figure text is ASCII.

library(ggplot2)
library(dplyr)
library(tibble)
library(mlr3)
library(mlr3learners)

hte_colors <- c(
  "S-learner" = "#B22222",
  "T-learner" = "#E69F00",
  "X-learner" = "#CC79A7",
  "DA X-learner" = "#7B3F99",
  "DR-learner" = "#0072B2",
  "R-learner" = "#009E73",
  "Q-aggregation" = "#000000",
  "Best single" = "#999999",
  "Truth" = "#000000"
)

hte_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Learners used throughout: random forests for regressions, logistic
# regression for the propensity score (the designs of Section 3 have a
# constant propensity, which a logit fits exactly).
hte_forest <- function(min_node = 10, trees = 200) {
  lrn("regr.ranger", num.trees = trees, min.node.size = min_node)
}

# Fit the five meta-learners plus the domain-adapted X-learner on `train`
# and return their predictions on `newdata`.
fit_all_learners <- function(train, x, x_het = x, newdata, learner = hte_forest(),
                             learner_final = learner, seed = NULL) {
  sc <- dr_scores(train, "y", "d", x, learner_mu = learner, seed = seed)
  spec <- list(
    "S-learner" = list(method = "s"),
    "T-learner" = list(method = "t"),
    "X-learner" = list(method = "x"),
    "DA X-learner" = list(method = "x", adapt = TRUE),
    "DR-learner" = list(method = "dr"),
    "R-learner" = list(method = "r")
  )
  fits <- lapply(spec, function(s) {
    do.call(cate_learner, c(list(scores = sc, x = x, x_het = x_het, learner = learner,
                                 learner_final = learner_final), s))
  })
  preds <- sapply(fits, function(f) predict(f, newdata))
  list(fits = fits, preds = preds, scores = sc)
}

# Monte Carlo of Section 3: RMSE of each learner and of two ensembles
# across the three designs of Chernozhukov et al. (2026, Example 15.1.1).
# Learners are fitted on 70% of a sample of size n, ensembles are stacked on
# the remaining 30%, and the RMSE against the true CATE is computed on a
# fresh sample of size n_test from the same design (so no learner is judged
# on its own training rows). First stages use small leaves so that the
# discontinuity in the baseline can be learned; final stages use larger
# leaves because their labels are noisy.
run_metalearner_experiment <- function(n_sims = 100, n_grid = c(500, 2000), n_test = 2000,
                                       dgps = c("simple_cate", "complex_cate", "unbalanced"),
                                       final_min_node = 100, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (g in dgps) {
    test <- sim_hte(n_test, dgp = g)
    for (n in n_grid) for (s in seq_len(n_sims)) {
      dat <- sim_hte(n, dgp = g)
      n_train <- round(0.7 * n)
      train <- dat[seq_len(n_train), ]
      score_set <- dat[-seq_len(n_train), ]
      fl <- fit_all_learners(train, "x1", newdata = test, learner = hte_forest(min_node = 5, trees = 100),
                             learner_final = hte_forest(min_node = final_min_node, trees = 100))
      sc_score <- dr_scores(score_set, "y", "d", "x1", learner_mu = hte_forest(min_node = 5, trees = 100))
      base <- fl$fits
      names(base) <- c("s_learner", "t_learner", "x_learner", "dax_learner", "dr_learner", "r_learner")
      ens <- list(
        "Q-aggregation" = do.call(cate_ensemble, c(list(sc_score, method = "q"), base)),
        "Best single" = do.call(cate_ensemble, c(list(sc_score, method = "best"), base))
      )
      preds <- cbind(fl$preds, sapply(ens, function(e) predict(e, test)))
      rmse <- sqrt(colMeans((preds - test$tau_true)^2))
      rows[[length(rows) + 1L]] <- tibble(dgp = g, n = n, sim = s, learner = names(rmse), rmse = unname(rmse))
    }
  }
  bind_rows(rows)
}

summarise_metalearner_experiment <- function(res) {
  res %>%
    group_by(dgp, n, learner) %>%
    summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse), median_rmse = median(rmse), .groups = "drop") %>%
    mutate(learner = factor(learner, levels = names(hte_colors)),
           dgp = factor(dgp, levels = c("simple_cate", "complex_cate", "unbalanced"),
                        labels = c("DGP 1: simple CATE, hard baseline, 5% treated",
                                   "DGP 2: hard CATE, simple baseline, 5% treated",
                                   "DGP 3: hard CATE, simple baseline, 95% treated"))) %>%
    arrange(dgp, n, learner)
}

# One draw of a design with the fitted CATE of each learner along x1.
plot_learner_fits <- function(dat, preds, which = colnames(preds), title = NULL) {
  df <- bind_rows(lapply(which, function(w) tibble(x1 = dat$x1, learner = w, cate = preds[, w])))
  truth <- tibble(x1 = dat$x1, cate = dat$tau_true)
  ggplot() +
    geom_line(data = truth %>% arrange(x1), aes(x1, cate), colour = "black", linewidth = 1) +
    geom_line(data = df %>% arrange(x1), aes(x1, cate, colour = learner), linewidth = 0.7, alpha = 0.9) +
    scale_colour_manual(values = hte_colors, name = NULL) +
    labs(x = "x1", y = "CATE", title = title, subtitle = "Black: true CATE") +
    hte_theme()
}

# Cumulative effect curve of Facure (2023, ch. 23): order units by a
# prediction, take the difference in means between treated and controls among
# the top share, and trace it as the share grows.
cumulative_effect_curve <- function(df, treatment, outcome, prediction, min_rows = 500, steps = 100) {
  ord <- df[order(-df[[prediction]]), ]
  n <- nrow(ord)
  rows <- unique(c(seq(min_rows, n, length.out = steps), n))
  tibble(share = rows / n,
         effect = sapply(rows, function(r) {
           top <- ord[seq_len(floor(r)), ]
           mean(top[[outcome]][top[[treatment]] == 1]) - mean(top[[outcome]][top[[treatment]] == 0])
         }))
}

# Monte Carlo of Section 7: regret of policies learned by empirical welfare
# maximization (depth-1 and depth-2 trees) and by thresholding CATE models
# (T-learner and DR-learner signs), evaluated with the true CATE on a large
# test sample.
run_policy_experiment <- function(n_sims = 50, n_grid = c(500, 1000, 2000, 4000), n_test = 20000, seed = 1) {
  set.seed(seed)
  test <- sim_hte(n_test, dgp = "policy")
  best_value <- mean(pmax(test$tau_true, 0))
  x <- paste0("x", 1:5)
  rows <- list()
  for (n in n_grid) {
    for (s in seq_len(n_sims)) {
      dat <- sim_hte(n, dgp = "policy")
      sc <- dr_scores(dat, "y", "d", x, learner_mu = hte_forest(min_node = 20, trees = 100))
      pols <- list(
        "Tree, depth 1" = policy_learn(sc, method = "tree", depth = 1, holdout = 0),
        "Tree, depth 2" = policy_learn(sc, method = "tree", depth = 2, holdout = 0),
        "Sign of T-learner CATE" = cate_learner(scores = sc, method = "t", learner = hte_forest(min_node = 20, trees = 100)),
        "Sign of DR-learner CATE" = cate_learner(scores = sc, method = "dr", learner = hte_forest(min_node = 20, trees = 100))
      )
      for (nm in names(pols)) {
        assign_test <- if (inherits(pols[[nm]], "cm_policy")) predict(pols[[nm]], test) else as.numeric(predict(pols[[nm]], test) > 0)
        value <- mean(assign_test * test$tau_true)
        rows[[length(rows) + 1L]] <- tibble(n = n, sim = s, policy = nm, value = value, regret = best_value - value)
      }
    }
  }
  bind_rows(rows)
}

# Monte Carlo of Section 3.6: the "smooth" design with linear oracles. The
# outcome regressions are misspecified (the baseline is quadratic in x2 and
# the learners are linear in x1, ..., x5) while the logistic propensity model
# is correctly specified. The final stage is linear in x1 and x2^2, which
# can represent the true CATE exactly. Learners that rely on outcome
# modelling alone (S, T, X) inherit the misspecification; the DR- and
# R-learners are protected by the correct propensity.
run_misspecification_experiment <- function(n_sims = 200, n = 2000, n_test = 5000, seed = 1) {
  set.seed(seed)
  test <- sim_hte(n_test, dgp = "smooth")
  test$x2sq <- test$x2^2
  x <- paste0("x", 1:5)
  sum_err <- NULL
  sum_sq <- NULL
  rmse_rows <- list()
  for (s in seq_len(n_sims)) {
    dat <- sim_hte(n, dgp = "smooth")
    dat$x2sq <- dat$x2^2
    fl <- fit_all_learners(dat, x, x_het = c("x1", "x2sq"), newdata = test,
                           learner = lrn("regr.lm"), learner_final = lrn("regr.lm"))
    err <- fl$preds - test$tau_true
    if (is.null(sum_err)) { sum_err <- err; sum_sq <- err^2 } else { sum_err <- sum_err + err; sum_sq <- sum_sq + err^2 }
    rmse_rows[[s]] <- tibble(sim = s, learner = colnames(err), rmse = sqrt(colMeans(err^2)))
  }
  mean_err <- sum_err / n_sims
  mean_sq <- sum_sq / n_sims
  list(
    rmse = bind_rows(rmse_rows),
    pointwise = tibble(learner = colnames(err),
                       bias_rms = sqrt(colMeans(mean_err^2)),
                       sd_rms = sqrt(colMeans(mean_sq - mean_err^2)),
                       rmse = sqrt(colMeans(mean_sq)))
  )
}

summarise_misspecification_experiment <- function(res) {
  res$pointwise %>%
    mutate(learner = factor(learner, levels = names(hte_colors))) %>%
    arrange(learner)
}
