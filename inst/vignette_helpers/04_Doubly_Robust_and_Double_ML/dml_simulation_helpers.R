# dml_simulation_helpers.R
# ------------------------------------------------------------
# Simulation helpers for the "Doubly Robust Estimation and Double
# Machine Learning" lecture (Section 04).
#
# Provides:
#   - simulate_dml_dgp():        a confounded design with nonlinear
#                                 nuisance functions and a constant effect
#   - naive_plugin_forest():     the plug-in estimator that DML replaces
#   - run_dml_experiment():      Monte Carlo over the competing estimators
#   - summarise_dml_experiment(): bias, dispersion, RMSE, coverage table
#   - plot_dml_experiment():     centred sampling distributions with the
#                                 normal approximation each estimator reports
#
# All figure text is ASCII or plotmath: the lecture PDF is rendered with
# the pdf() device, which cannot encode most Unicode glyphs.
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)

dml_colors <- c(
  "OLS with linear controls"              = "grey45",
  "Naive plug-in forest"                  = "#993C1D",
  "Naive plug-in forest, cross-fitted"    = "#C97B5A",
  "Orthogonal score, full-sample forests" = "#B8860B",
  "Orthogonal score, cross-fitted"        = "#185FA5",
  "AIPW, cross-fitted"                    = "#2E7D32"
)

dml_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey35")
    )
}

# ------------------------------------------------------------
# 1. Data-generating process
# ------------------------------------------------------------
# Binary treatment, constant effect `theta`, so the ATE, the ATT, and the
# partially linear coefficient coincide. The confounders x1 and x2 enter
# both the propensity score and the outcome nonlinearly, so a linear
# adjustment is misspecified while a forest can learn both surfaces.

simulate_dml_dgp <- function(n = 500, p = 10, theta = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]

  m0 <- plogis(0.8 * x1 - 0.4 * (x2^2 - 1))
  g0 <- 2 * plogis(2 * x1) + 0.5 * x2^2 + 0.25 * x3

  d <- rbinom(n, 1, m0)
  y <- theta * d + g0 + rnorm(n)

  data.frame(y = y, d = d, X, m_true = m0, g_true = g0)
}

# ------------------------------------------------------------
# 2. The naive plug-in estimator
# ------------------------------------------------------------
# One random forest of Y on (D, X). The ATE is estimated by predicting
# every unit under D = 1 and D = 0 and averaging the difference, i.e.
# regression imputation with a flexible learner. The reported standard
# error is the one a naive analyst would use: the dispersion of the
# imputed differences divided by sqrt(n).

naive_plugin_forest <- function(dat, x, num.trees = 100) {
  train <- dat[, c("y", "d", x)]
  rf <- ranger::ranger(y ~ ., data = train, num.trees = num.trees)
  d1 <- train
  d1$d <- 1
  d0 <- train
  d0$d <- 0
  contrast <- predict(rf, d1)$predictions - predict(rf, d0)$predictions
  c(estimate = mean(contrast), std.error = stats::sd(contrast) / sqrt(nrow(train)))
}

# The same plug-in with cross-fitting: the forest is fitted on K - 1 folds
# and the contrast is predicted on the held-out fold. Sample splitting
# alone does not repair a non-orthogonal moment.
naive_plugin_forest_crossfit <- function(dat, x, num.trees = 100, folds = 5) {
  train <- dat[, c("y", "d", x)]
  n <- nrow(train)
  fold_id <- sample(rep(seq_len(folds), length.out = n))
  contrast <- numeric(n)
  for (k in seq_len(folds)) {
    rf <- ranger::ranger(y ~ ., data = train[fold_id != k, ], num.trees = num.trees)
    test <- train[fold_id == k, ]
    d1 <- test
    d1$d <- 1
    d0 <- test
    d0$d <- 0
    contrast[fold_id == k] <- predict(rf, d1)$predictions - predict(rf, d0)$predictions
  }
  c(estimate = mean(contrast), std.error = stats::sd(contrast) / sqrt(n))
}

# ------------------------------------------------------------
# 3. Monte Carlo experiment
# ------------------------------------------------------------
# Six estimators on every simulated sample:
#   - OLS with linear controls (Section 03's regression adjustment)
#   - naive plug-in forest, and the same plug-in cross-fitted
#   - orthogonal (partialling-out) score, forests fit on the full sample
#   - orthogonal (partialling-out) score, cross-fitted forests (est_dml)
#   - AIPW score, cross-fitted forests (est_dml, interactive model)

run_dml_experiment <- function(n_sims = 100, n = 500, p = 10, theta = 1,
                               num.trees = 100, folds = 5, seed = 1,
                               verbose = FALSE) {
  x <- paste0("x", seq_len(p))
  learner_regr <- function() mlr3::lrn("regr.ranger", num.trees = num.trees)
  learner_prob <- function() mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = num.trees)

  set.seed(seed)
  rows <- vector("list", n_sims)
  for (s in seq_len(n_sims)) {
    if (verbose && s %% 10 == 0) message("simulation ", s, " of ", n_sims)
    dat <- simulate_dml_dgp(n = n, p = p, theta = theta)

    ols <- estimatr::lm_robust(
      stats::reformulate(c("d", x), response = "y"), data = dat, se_type = "HC2"
    )
    naive <- naive_plugin_forest(dat, x, num.trees = num.trees)
    naive_cf <- naive_plugin_forest_crossfit(dat, x, num.trees = num.trees, folds = folds)
    full <- est_dml(
      dat, y = "y", d = "d", x = x, model = "plr",
      learner_l = learner_regr(), learner_m = learner_prob(),
      cross_fit = FALSE
    )
    cross <- est_dml(
      dat, y = "y", d = "d", x = x, model = "plr",
      learner_l = learner_regr(), learner_m = learner_prob(),
      folds = folds
    )
    aipw <- est_dml(
      dat, y = "y", d = "d", x = x, model = "irm", estimand = "ATE",
      learner_p = learner_prob(), learner_mu0 = learner_regr(), learner_mu1 = learner_regr(),
      folds = folds
    )

    rows[[s]] <- tibble(
      sim = s,
      estimator = names(dml_colors),
      estimate = c(ols$coefficients[["d"]], naive[["estimate"]], naive_cf[["estimate"]],
                   full$estimate, cross$estimate, aipw$estimate),
      std.error = c(ols$std.error[["d"]], naive[["std.error"]], naive_cf[["std.error"]],
                    full$std.error, cross$std.error, aipw$std.error)
    )
  }
  out <- bind_rows(rows)
  out$estimator <- factor(out$estimator, levels = names(dml_colors))
  attr(out, "theta") <- theta
  out
}

# ------------------------------------------------------------
# 4. Summary table and figure
# ------------------------------------------------------------

summarise_dml_experiment <- function(results, theta = attr(results, "theta"), conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  results %>%
    group_by(estimator) %>%
    summarise(
      bias = mean(estimate - theta),
      sd = stats::sd(estimate),
      rmse = sqrt(mean((estimate - theta)^2)),
      mean_se = mean(std.error),
      coverage = mean(abs(estimate - theta) <= z * std.error),
      .groups = "drop"
    )
}

plot_dml_experiment <- function(results, theta = attr(results, "theta"),
                                estimators = levels(results$estimator), bins = 25,
                                colors = dml_colors, ncol = 2) {
  res <- results %>%
    filter(estimator %in% estimators) %>%
    mutate(centered = estimate - theta, estimator = factor(estimator, levels = estimators))

  # The normal approximation each estimator itself reports: centred at zero
  # with the average reported standard error. A histogram that is shifted
  # or wider than its curve is a failure of bias or of the standard error.
  curves <- res %>%
    group_by(estimator) %>%
    summarise(se = mean(std.error), .groups = "drop") %>%
    rowwise() %>%
    do({
      grid <- seq(min(res$centered), max(res$centered), length.out = 200)
      tibble(estimator = .$estimator, centered = grid, density = stats::dnorm(grid, 0, .$se))
    }) %>%
    ungroup() %>%
    mutate(estimator = factor(estimator, levels = estimators))

  ggplot(res, aes(x = centered)) +
    geom_histogram(aes(y = after_stat(density), fill = estimator),
                   bins = bins, color = "white", alpha = 0.85) +
    geom_line(data = curves, aes(y = density), color = "grey20", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    facet_wrap(~ estimator, ncol = ncol, scales = "free_y") +
    scale_fill_manual(values = colors, guide = "none") +
    labs(
      x = "Estimate minus true effect",
      y = "Density",
      subtitle = "Curve: the normal approximation implied by the estimator's own reported standard error"
    ) +
    dml_theme()
}


# ------------------------------------------------------------
# 5. Heterogeneous effects: what the partially linear coefficient targets
# ------------------------------------------------------------
# Binary treatment with a covariate-dependent effect tau(X) = 1 + x1^2,
# largest exactly where the propensity score is extreme and the overlap
# weights m(1 - m) are small. The true nuisances are returned so estimators
# can be compared without any learning error: the PLM coefficient is the
# overlap-weighted average of tau(X); the interactive model targets the
# plain average.

simulate_hte_dgp <- function(n = 5000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  m0 <- plogis(1.5 * x1 - 0.5 * x2)
  d <- rbinom(n, 1, m0)
  tau <- 1 + x1^2
  g0 <- 0.5 * x1 + x2^2
  y <- g0 + d * tau + rnorm(n)
  data.frame(
    y = y, d = d, x1 = x1, x2 = x2,
    tau_true = tau, m_true = m0,
    g0_true = g0, g1_true = g0 + tau,
    l_true = g0 + m0 * tau
  )
}

# ------------------------------------------------------------
# 6. High-dimensional linear model: single selection vs double lasso
# ------------------------------------------------------------
# Example 4.3.1 of Chernozhukov et al.: p = n = 100, coefficients 1/j^2 on
# the controls in both equations, small residual treatment noise sd_v.

simulate_hd_linear_dgp <- function(n = 100, p = 100, theta = 1, sd_v = 0.25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(W) <- paste0("w", seq_len(p))
  coefs <- 1 / seq_len(p)^2
  m_true <- as.numeric(W %*% coefs)
  d <- m_true + rnorm(n) * sd_v
  y <- theta * d + as.numeric(W %*% coefs) + rnorm(n)
  data.frame(y = y, d = d, W, m_true = m_true)
}

# Indices of the controls selected by a cross-validated lasso of `y` on `X`.
.lasso_selected <- function(y, X, penalty.factor = NULL, s = "lambda.min") {
  if (is.null(penalty.factor)) penalty.factor <- rep(1, ncol(X))
  cv <- glmnet::cv.glmnet(X, y, penalty.factor = penalty.factor)
  which(as.numeric(stats::coef(cv, s = s))[-1] != 0)
}

# Naive single selection: lasso of Y on (D, W) with D unpenalised, keep the
# selected controls, refit by least squares, report the usual robust SE.
single_selection_lasso <- function(dat, w) {
  X <- as.matrix(dat[, c("d", w)])
  sel <- .lasso_selected(dat$y, X, penalty.factor = c(0, rep(1, length(w))))
  selected <- w[setdiff(sel, 1) - 1]
  fit <- estimatr::lm_robust(
    stats::reformulate(c("d", selected), response = "y"), data = dat, se_type = "HC1"
  )
  c(estimate = fit$coefficients[["d"]], std.error = fit$std.error[["d"]],
    n_selected = length(selected))
}

# Post-double-selection (Belloni, Chernozhukov, and Hansen 2014): lasso Y on W
# and D on W, take the union of the selected controls, refit by least squares.
double_selection_lasso <- function(dat, w) {
  X <- as.matrix(dat[, w])
  sel <- union(.lasso_selected(dat$y, X), .lasso_selected(dat$d, X))
  fit <- estimatr::lm_robust(
    stats::reformulate(c("d", w[sel]), response = "y"), data = dat, se_type = "HC1"
  )
  c(estimate = fit$coefficients[["d"]], std.error = fit$std.error[["d"]],
    n_selected = length(sel))
}

lasso_colors <- c(
  "Single selection"        = "#993C1D",
  "Post-double-selection"   = "#2E7D32",
  "CV lasso, cross-fitted"  = "#185FA5"
)

run_double_lasso_experiment <- function(n_sims = 60, n = 100, p = 100, theta = 1,
                                        sd_v = 0.25, folds = 5, seed = 1) {
  w <- paste0("w", seq_len(p))
  set.seed(seed)
  rows <- vector("list", n_sims)
  for (s in seq_len(n_sims)) {
    dat <- simulate_hd_linear_dgp(n = n, p = p, theta = theta, sd_v = sd_v)
    single <- single_selection_lasso(dat, w)
    double <- double_selection_lasso(dat, w)
    dl <- est_dml(
      dat, y = "y", d = "d", x = w, model = "plr",
      learner_l = mlr3::lrn("regr.cv_glmnet", nfolds = 5),
      learner_m = mlr3::lrn("regr.cv_glmnet", nfolds = 5),
      folds = folds
    )
    rows[[s]] <- tibble(
      sim = s,
      estimator = names(lasso_colors),
      estimate = c(single[["estimate"]], double[["estimate"]], dl$estimate),
      std.error = c(single[["std.error"]], double[["std.error"]], dl$std.error),
      n_selected = c(single[["n_selected"]], double[["n_selected"]], NA_real_),
      msq_treatment_residual = c(NA_real_, NA_real_,
                                 dl$diagnostics$identification$mean_sq_treatment_residual),
      msq_m_error = c(NA_real_, NA_real_, mean((dl$nuisance$m_hat - dat$m_true)^2))
    )
  }
  out <- bind_rows(rows)
  out$estimator <- factor(out$estimator, levels = names(lasso_colors))
  attr(out, "theta") <- theta
  attr(out, "sd_v") <- sd_v
  out
}
