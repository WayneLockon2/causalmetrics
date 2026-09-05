sim_plr_data <- function(n = 800, theta = 1.5, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  d <- 0.5 * x1 + sin(x2) + rnorm(n)
  y <- theta * d + x1^2 + exp(x2 / 2) + rnorm(n)
  data.frame(
    y = y, d = d, x1 = x1, x2 = x2,
    l_true = x1^2 + exp(x2 / 2) + theta * (0.5 * x1 + sin(x2)),
    m_true = 0.5 * x1 + sin(x2)
  )
}

sim_irm_data <- function(n = 800, tau = 2, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  p <- plogis(-0.2 + 0.6 * x1 - 0.4 * x2)
  d <- rbinom(n, 1, p)
  mu0 <- 1 + x1 + x2
  y <- mu0 + tau * d + rnorm(n)
  data.frame(y = y, d = d, x1 = x1, x2 = x2, p = p, mu0 = mu0, mu1 = mu0 + tau)
}

sim_att_data <- function(n = 3000, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  p <- plogis(-0.2 + 0.8 * x1 - 0.4 * x2)
  d <- rbinom(n, 1, p)
  mu0 <- 1 + x1 + x2
  tau <- 1 + x1
  y <- mu0 + d * tau + rnorm(n)
  data.frame(y = y, d = d, x1 = x1, x2 = x2, p = p, mu0 = mu0, mu1 = mu0 + tau, tau = tau)
}

test_that("partially linear model with supplied true nuisances recovers theta", {
  dat <- sim_plr_data(n = 2000)
  fit <- est_dml(dat, "y", "d", l_hat = "l_true", m_hat = "m_true")

  expect_s3_class(fit, "cm_dml")
  expect_equal(fit$model, "plr")
  expect_equal(fit$estimand, "theta")
  expect_lt(abs(fit$estimate - 1.5), 4 * fit$std.error)

  # Closed form: residual-on-residual slope and HC1-type standard error.
  yt <- dat$y - dat$l_true
  dt_ <- dat$d - dat$m_true
  n <- nrow(dat)
  expect_equal(fit$estimate, sum(yt * dt_) / sum(dt_^2))
  e <- yt - fit$estimate * dt_
  expect_equal(fit$std.error, sqrt(sum(dt_^2 * e^2) / (n - 1)) / (mean(dt_^2) * sqrt(n)))
  expect_equal(fit$residuals$y_tilde, yt)
  expect_equal(fit$residuals$d_tilde, dt_)
  expect_equal(fit$diagnostics$prediction_mode, "supplied")
  expect_equal(fit$n_rep, 1L)
  expect_true(is.na(fit$n_treated))
  expect_output(print(fit), "Double ML estimate \\(partially linear regression")

  # Vectors and column names are interchangeable.
  fit_vec <- est_dml(dat[, c("y", "d")], "y", "d", l_hat = dat$l_true, m_hat = dat$m_true)
  expect_equal(fit_vec$estimate, fit$estimate)
})

test_that("full-sample OLS nuisances reproduce OLS exactly (FWL) with HC1 standard errors", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  set.seed(2)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  d <- 0.3 * x1 - 0.2 * x2 + rnorm(n)
  y <- 2 * d + x1 + 0.5 * x2 + rnorm(n)
  dat <- data.frame(y = y, d = d, x1 = x1, x2 = x2)

  fit <- est_dml(dat, "y", "d", x = c("x1", "x2"), cross_fit = FALSE)
  ols <- lm(y ~ d + x1 + x2, data = dat)
  expect_equal(fit$estimate, unname(coef(ols)[["d"]]))

  yt <- unname(resid(lm(y ~ x1 + x2, data = dat)))
  dt_ <- unname(resid(lm(d ~ x1 + x2, data = dat)))
  e <- yt - fit$estimate * dt_
  expect_equal(fit$std.error, sqrt(sum(dt_^2 * e^2) / (n - 1)) / (mean(dt_^2) * sqrt(n)))
  expect_equal(fit$diagnostics$prediction_mode, "full_sample")
  expect_false(fit$diagnostics$fold_summary$training_excludes_test)
  expect_equal(fit$learners$l_hat, "regr.lm")
})

test_that("cross-fitted flexible learners fit the nonlinear nuisances better than linear ones", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  dat <- sim_plr_data(n = 1500)
  fit_lin <- est_dml(dat, "y", "d", x = c("x1", "x2"), seed = 1)
  fit_rf <- est_dml(
    dat, "y", "d", x = c("x1", "x2"),
    learner_l = mlr3::lrn("regr.ranger", num.trees = 200),
    learner_m = mlr3::lrn("regr.ranger", num.trees = 200),
    seed = 1
  )
  q_lin <- fit_lin$diagnostics$nuisance
  q_rf <- fit_rf$diagnostics$nuisance
  expect_lt(q_rf$rmse[q_rf$nuisance == "l_hat"], q_lin$rmse[q_lin$nuisance == "l_hat"])
  expect_lt(abs(fit_rf$estimate - 1.5), 4 * fit_rf$std.error)
  expect_equal(fit_rf$learners$l_hat, "regr.ranger")

  fs <- fit_rf$diagnostics$fold_summary
  expect_equal(nrow(fs), 5L)
  expect_true(all(fs$training_excludes_test))
  expect_equal(sum(fs$test_n), nrow(dat))
  expect_equal(nrow(fit_rf$fold_estimates), 5L)
  expect_equal(fit_rf$diagnostics$prediction_mode, "out_of_fold")

  # Deterministic learners plus a seed give reproducible results.
  fit_lin2 <- est_dml(dat, "y", "d", x = c("x1", "x2"), seed = 1)
  expect_equal(fit_lin$estimate, fit_lin2$estimate)
  expect_identical(fit_lin$fold_id, fit_lin2$fold_id)
})

test_that("the fold-average solution averages the fold estimates that the pooled solution stacks", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_plr_data(n = 600)
  f2 <- est_dml(dat, "y", "d", x = c("x1", "x2"), seed = 3, solve = "pooled")
  f1 <- est_dml(dat, "y", "d", x = c("x1", "x2"), seed = 3, solve = "fold_average")

  expect_equal(f1$fold_estimates$estimate, f2$fold_estimates$estimate)
  expect_equal(f1$estimate, mean(f1$fold_estimates$estimate))
  expect_equal(f2$estimate, with(f2$residuals, sum(y_tilde * d_tilde) / sum(d_tilde^2)))
  expect_false(isTRUE(all.equal(f1$estimate, f2$estimate)))
  expect_lt(abs(f1$estimate - f2$estimate), 0.5)
  expect_equal(f1$solve, "fold_average")
  expect_equal(f2$solve, "pooled")
  expect_output(print(f1), "average of 5 fold estimate")
  expect_output(print(f2), "pooled across 5 fold")
})

test_that("repeated cross-fitting aggregates with the median rule", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_plr_data(n = 600)
  fit <- est_dml(dat, "y", "d", x = c("x1", "x2"), n_rep = 3, seed = 4)

  expect_equal(nrow(fit$repetitions), 3L)
  expect_equal(fit$n_rep, 3L)
  est <- fit$repetitions$estimate
  se <- fit$repetitions$std.error
  expect_equal(fit$estimate, median(est))
  expect_equal(fit$std.error, sqrt(median(se^2 + (est - median(est))^2)))
  expect_equal(fit$rep_selected, which.min(abs(est - median(est))))
  expect_equal(nrow(fit$fold_estimates), 15L)
  expect_false(identical(est[1], est[2]))
  expect_output(print(fit), "3 repetition\\(s\\) \\(median aggregation\\)")

  fit2 <- est_dml(dat, "y", "d", x = c("x1", "x2"), n_rep = 3, seed = 4)
  expect_equal(fit$repetitions, fit2$repetitions)

  expect_error(est_dml(dat, "y", "d", l_hat = "l_true", m_hat = "m_true", n_rep = 2), "n_rep > 1")
  expect_error(est_dml(dat, "y", "d", x = c("x1", "x2"), n_rep = 2, cross_fit = FALSE), "cross_fit = TRUE")
})

test_that("binary treatment in the partially linear model uses a classifier for E[D | X]", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_irm_data(n = 800)
  fit <- est_dml(dat, "y", "d", x = c("x1", "x2"), model = "plr", seed = 1)
  expect_equal(fit$learners$m_hat, "classif.log_reg")
  expect_true(all(fit$nuisance$m_hat > 0 & fit$nuisance$m_hat < 1))
  expect_equal(fit$n_treated + fit$n_control, fit$n)
  expect_lt(abs(fit$estimate - 2), 4 * fit$std.error)
  expect_true(all(fit$diagnostics$fold_summary$training_excludes_test))
})

test_that("interactive model reproduces est_aipw for the ATE", {
  dat <- sim_irm_data(n = 1000)
  a <- est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  i <- est_dml(dat, "y", "d", model = "irm", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  expect_equal(i$estimate, a$estimate)
  expect_equal(i$std.error, a$std.error)
  expect_equal(i$score + i$estimate, a$score)
  expect_equal(i$estimand, "ATE")
  expect_equal(i$score_type, "AIPW_ATE")
  expect_equal(i$n_treated, a$n_treated)

  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  a2 <- est_aipw(dat, "y", "d", x = c("x1", "x2"), seed = 7)
  i2 <- est_dml(dat, "y", "d", model = "irm", x = c("x1", "x2"), seed = 7)
  expect_equal(i2$estimate, a2$estimate)
  expect_equal(i2$std.error, a2$std.error)
  expect_identical(i2$fold_id, a2$fold_id)
  expect_equal(i2$learners$p_hat, "classif.log_reg")
  expect_equal(nrow(i2$diagnostics$nuisance), 2L)
})

test_that("ATT score targets the treated population", {
  dat <- sim_att_data()
  att_true <- mean(dat$tau[dat$d == 1])
  ate_true <- mean(dat$tau)
  expect_gt(att_true - ate_true, 0.15)

  fit_att <- est_dml(dat, "y", "d", model = "irm", estimand = "ATT", p_hat = "p", mu0_hat = "mu0")
  fit_ate <- est_dml(dat, "y", "d", model = "irm", estimand = "ATE",
                     p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  expect_lt(abs(fit_att$estimate - att_true), 4 * fit_att$std.error)
  expect_lt(abs(fit_ate$estimate - ate_true), 4 * fit_ate$std.error)

  w <- with(dat, d - (1 - d) * p / (1 - p))
  expect_equal(fit_att$estimate, sum(w * (dat$y - dat$mu0)) / sum(dat$d))
  expect_null(fit_att$nuisance$mu1_hat)
  expect_equal(fit_att$learners$mu1_hat, "not needed for ATT")
  expect_equal(fit_att$score_type, "AIPW_ATT")

  a <- est_aipw(dat, "y", "d", estimand = "ATT", p_hat = "p", mu0_hat = "mu0")
  expect_equal(a$estimate, fit_att$estimate)
  expect_equal(a$std.error, fit_att$std.error)
  expect_null(a$mu1_hat)
  expect_equal(a$diagnostics$nuisance$source$mu1_hat, "not needed for ATT")
  expect_output(print(a), "AIPW estimate \\(ATT\\)")

  expect_error(est_dml(dat, "y", "d", model = "plr", estimand = "ATT", x = c("x1", "x2")),
               "only available")
})

test_that("est_dml validates its inputs", {
  dat <- sim_plr_data(n = 200)
  expect_error(est_dml(dat, "y", "d", x = c("x1", "x2"), model = "irm"), "binary")
  expect_error(est_dml(dat, "y", "d"), "`x` must be supplied")
  expect_error(est_dml(dat, "y", "d", l_hat = "l_true", m_hat = "m_true", p_hat = "l_true"),
               "belong to")
  expect_error(est_dml(dat, "y", "d", x = c("x1", "x2"), n_rep = 0), "n_rep")
  expect_error(est_dml(dat, "y", "d", x = c("x1", "x2"), solve = "dml2"), "should be one of")
  expect_error(est_dml(dat, "y", "nope", x = "x1"), "not found")
  expect_error(est_dml(dat, "y", "d", x = c("x1", "y"), l_hat = "l_true"), "must not include")

  dat$y[1] <- NA
  expect_error(est_dml(dat, "y", "d", l_hat = "l_true", m_hat = "m_true"), "na_action")
  expect_warning(
    fit <- est_dml(dat, "y", "d", l_hat = "l_true", m_hat = "m_true", na_action = "omit"),
    "omitted"
  )
  expect_equal(fit$n, 199L)
  expect_equal(fit$diagnostics$sample$omitted_missing, 1L)
})
