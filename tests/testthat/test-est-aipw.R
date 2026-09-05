sim_aipw_data <- function(n = 600, tau = 2, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  p <- plogis(-0.2 + 0.6 * x1 - 0.4 * x2)
  d <- rbinom(n, 1, p)
  mu0 <- 1 + x1 + x2
  mu1 <- mu0 + tau
  y <- mu0 + tau * d + rnorm(n)
  data.frame(y = y, d = d, x1 = x1, x2 = x2, p = p, mu0 = mu0, mu1 = mu1)
}

test_that("est_aipw recovers the ATE with true nuisances and reports IF standard errors", {
  dat <- sim_aipw_data(n = 2000, tau = 2)
  fit <- est_aipw(dat, y = "y", d = "d", p_hat = "p", mu0_hat = "mu0",
                  mu1_hat = "mu1", p_clip = NULL)

  expect_s3_class(fit, "cm_aipw")
  expect_equal(fit$estimand, "ATE")
  expect_equal(fit$n, nrow(dat))
  expect_equal(fit$n_treated + fit$n_control, nrow(dat))
  expect_lt(abs(fit$estimate - 2), 4 * fit$std.error)

  # The estimate is the mean AIPW score and the SE is its IF standard error.
  score <- with(dat, (mu1 - mu0) + d * (y - mu1) / p - (1 - d) * (y - mu0) / (1 - p))
  expect_equal(fit$estimate, mean(score))
  expect_equal(fit$score, score)
  expect_equal(fit$std.error, sd(score) / sqrt(nrow(dat)))
  expect_equal(fit$conf.low, fit$estimate - qnorm(0.975) * fit$std.error)
  expect_equal(fit$conf.high, fit$estimate + qnorm(0.975) * fit$std.error)

  expect_named(fit$weights, c("treated", "control", "ipw"))
  expect_equal(fit$diagnostics$nuisance$prediction_mode, "supplied")
  expect_equal(fit$diagnostics$sample$omitted_missing, 0L)
  expect_output(print(fit), "AIPW estimate \\(ATE\\)")
})

test_that("nuisances can be supplied as vectors or column names with identical results", {
  dat <- sim_aipw_data(n = 500)
  fit_cols <- est_aipw(dat, y = "y", d = "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  fit_vecs <- est_aipw(dat[, c("y", "d")], y = "y", d = "d",
                       p_hat = dat$p, mu0_hat = dat$mu0, mu1_hat = dat$mu1)
  expect_equal(fit_cols$estimate, fit_vecs$estimate)
  expect_equal(fit_cols$std.error, fit_vecs$std.error)

  # data.table input behaves like data.frame input.
  fit_dt <- est_aipw(data.table::as.data.table(dat), y = "y", d = "d",
                     p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  expect_equal(fit_dt$estimate, fit_cols$estimate)
})

test_that("treatment coding is validated and coerced", {
  dat <- sim_aipw_data(n = 300)
  args <- list(y = "y", d = "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  base <- do.call(est_aipw, c(list(dat), args))$estimate

  expect_equal(do.call(est_aipw, c(list(transform(dat, d = as.logical(d))), args))$estimate, base)
  expect_equal(do.call(est_aipw, c(list(transform(dat, d = as.character(d))), args))$estimate, base)
  expect_equal(do.call(est_aipw, c(list(transform(dat, d = factor(d))), args))$estimate, base)
  expect_error(do.call(est_aipw, c(list(transform(dat, d = d + 1)), args)), "binary")
  expect_error(do.call(est_aipw, c(list(dat[dat$d == 1, ]), args)), "Both treated and control")
})

test_that("input validation is informative", {
  dat <- sim_aipw_data(n = 200)
  full <- list(p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")

  expect_error(est_aipw(1:10, y = "y", d = "d"), "data frame")
  expect_error(do.call(est_aipw, c(list(dat, y = "nope", d = "d"), full)), "not found")
  expect_error(do.call(est_aipw, c(list(dat, y = "y", d = "d", estimand = "overlap"), full)),
               "must be \\\"ATE\\\" or \\\"ATT\\\"")
  expect_error(est_aipw(dat, y = "y", d = "d", p_hat = "p", mu0_hat = "mu0"),
               "`x` must be supplied")
  expect_error(do.call(est_aipw, c(list(dat, y = "y", d = "d", conf_level = 1.2), full)),
               "conf_level")
  expect_error(est_aipw(dat, y = "y", d = "d", p_hat = dat$p[-1], mu0_hat = "mu0", mu1_hat = "mu1"),
               "length nrow\\(data\\)")
  expect_error(do.call(est_aipw, c(list(transform(dat, p = p * 2), y = "y", d = "d"), full)),
               "between 0 and 1")
  expect_error(do.call(est_aipw, c(list(dat, y = "y", d = "d", p_clip = c(0, 0.99)), full)),
               "strictly inside")
  expect_error(est_aipw(dat, y = "y", d = "d", x = c("x1", "y"), p_hat = "p"),
               "must not include")
  expect_error(est_aipw(dat, y = "y", d = "d", x = c("x1", "x1"), p_hat = "p"),
               "duplicated")
  expect_error(est_aipw(dat, y = "y", d = "d", x = ".cm_y", p_hat = "p"),
               "reserved")
})

test_that("na_action controls missing-value handling", {
  dat <- sim_aipw_data(n = 300)
  dat$y[c(3, 10)] <- NA
  expect_error(est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1"),
               "na_action = 'omit'")
  expect_warning(
    fit <- est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1",
                    na_action = "omit"),
    "2 row\\(s\\) omitted"
  )
  expect_equal(fit$n, 298L)
  expect_equal(fit$diagnostics$sample$omitted_missing, 2L)
  expect_equal(fit$diagnostics$sample$n_before_missing, 300L)
})

test_that("trimming and clipping are applied and recorded", {
  dat <- sim_aipw_data(n = 500)
  common <- list(y = "y", d = "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")

  fit <- do.call(est_aipw, c(list(dat, trim = c(0.2, 0.8)), common))
  expect_equal(fit$diagnostics$trimming$trim, c(0.2, 0.8))
  expect_equal(fit$diagnostics$trimming$n_trimmed, sum(dat$p < 0.2 | dat$p > 0.8))
  expect_equal(fit$n, sum(dat$p >= 0.2 & dat$p <= 0.8))
  expect_true(all(fit$p_hat >= 0.2 & fit$p_hat <= 0.8))

  fit_clip <- do.call(est_aipw, c(list(dat, p_clip = c(0.3, 0.7)), common))
  expect_equal(fit_clip$diagnostics$propensity$n_clipped_low, sum(dat$p < 0.3))
  expect_equal(fit_clip$diagnostics$propensity$n_clipped_high, sum(dat$p > 0.7))
  expect_true(all(fit_clip$p_hat >= 0.3 & fit_clip$p_hat <= 0.7))
  expect_equal(fit_clip$n, nrow(dat))

  expect_error(do.call(est_aipw, c(list(dat, trim = c(0.49999, 0.5)), common)),
               "Trimming removed too many")
})

test_that(".cm_make_folds stratifies by treatment and is reproducible", {
  d <- rep(c(0L, 1L), times = c(70, 30))
  f1 <- .cm_make_folds(d, folds = 5, seed = 42)
  f2 <- .cm_make_folds(d, folds = 5, seed = 42)
  expect_identical(f1, f2)
  expect_setequal(unique(f1), 1:5)

  tab <- table(f1, d)
  expect_true(all(tab > 0))
  expect_true(all(abs(tab[, "1"] - 6) <= 1))
  expect_true(all(abs(tab[, "0"] - 14) <= 1))

  expect_error(.cm_make_folds(d, folds = 1, seed = 1), "at least 2")
  expect_error(.cm_make_folds(d, folds = 101, seed = 1), "cannot exceed")
})

test_that(".cm_make_folds leaves the global RNG stream untouched when seeded", {
  set.seed(123)
  before <- runif(1)
  set.seed(123)
  invisible(.cm_make_folds(rep(0:1, 20), folds = 2, seed = 7))
  after <- runif(1)
  expect_equal(before, after)
})

test_that("weight and support diagnostics behave", {
  expect_equal(.cm_ess(c(1, 1, 1, 1)), 4)
  expect_equal(.cm_ess(c(1, 0, 0, 0)), 1)
  expect_true(is.na(.cm_ess(numeric(0))))

  p <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  d <- c(0, 0, 0, 1, 0, 1, 1, 1)
  cs <- .cm_common_support(p, d)
  expect_equal(cs$low, 0.4)
  expect_equal(cs$high, 0.5)
  expect_equal(cs$share_outside, 6 / 8)
})

test_that("internal mlr3 nuisance estimation cross-fits by default", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_aipw_data(n = 600)
  fit <- est_aipw(dat, "y", "d", x = c("x1", "x2"), folds = 4, seed = 11)

  expect_equal(fit$diagnostics$nuisance$prediction_mode, "out_of_fold")
  expect_equal(unname(unlist(fit$diagnostics$nuisance$source)), rep("mlr3_out_of_fold", 3))
  fs <- fit$diagnostics$nuisance$fold_summary
  expect_equal(nrow(fs), 4L)
  expect_true(all(fs$training_excludes_test))
  expect_equal(sum(fs$test_n), nrow(dat))
  expect_lt(abs(fit$estimate - 2), 4 * fit$std.error)

  fit2 <- est_aipw(dat, "y", "d", x = c("x1", "x2"), folds = 4, seed = 11)
  expect_equal(fit$estimate, fit2$estimate)
  expect_identical(fit$fold_id, fit2$fold_id)
})

test_that("cross_fit = FALSE fits on the full sample and supplied fold_id is honoured", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_aipw_data(n = 400)
  fit_full <- est_aipw(dat, "y", "d", x = c("x1", "x2"), cross_fit = FALSE)
  expect_equal(fit_full$diagnostics$nuisance$prediction_mode, "full_sample")
  expect_false(fit_full$diagnostics$nuisance$fold_summary$training_excludes_test)

  dat$fold <- rep(1:3, length.out = nrow(dat))
  fit_fid <- est_aipw(dat, "y", "d", x = c("x1", "x2"), fold_id = "fold")
  expect_equal(fit_fid$diagnostics$call$folds, 3L)
  expect_identical(fit_fid$fold_id, as.integer(as.factor(dat$fold)))
})

test_that("partially supplied nuisances are completed by learners", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  dat <- sim_aipw_data(n = 400)
  fit <- est_aipw(dat, "y", "d", x = c("x1", "x2"), p_hat = "p", folds = 3, seed = 2)
  expect_equal(fit$diagnostics$nuisance$source$p_hat, "supplied")
  expect_equal(fit$diagnostics$nuisance$source$mu0_hat, "mlr3_out_of_fold")
  expect_equal(fit$diagnostics$nuisance$source$mu1_hat, "mlr3_out_of_fold")
})

test_that("custom learners and binary outcomes work", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  dat <- sim_aipw_data(n = 400)
  fit <- est_aipw(
    dat, "y", "d", x = c("x1", "x2"),
    learner_p   = mlr3::lrn("classif.ranger", predict_type = "prob", num.trees = 50),
    learner_mu0 = mlr3::lrn("regr.ranger", num.trees = 50),
    learner_mu1 = mlr3::lrn("regr.ranger", num.trees = 50),
    folds = 3, seed = 3
  )
  expect_s3_class(fit, "cm_aipw")
  expect_lt(abs(fit$estimate - 2), 6 * fit$std.error)

  dat$yb <- as.integer(dat$y > median(dat$y))
  fitb <- est_aipw(
    dat, "yb", "d", x = c("x1", "x2"),
    learner_mu0 = mlr3::lrn("classif.log_reg", predict_type = "prob"),
    learner_mu1 = mlr3::lrn("classif.log_reg", predict_type = "prob"),
    folds = 3, seed = 4
  )
  expect_equal(fitb$diagnostics$nuisance$outcome_type, "binary")
  expect_true(all(fitb$mu0_hat >= 0 & fitb$mu0_hat <= 1))

  expect_error(
    est_aipw(dat, "y", "d", x = c("x1", "x2"),
             learner_mu0 = mlr3::lrn("classif.log_reg", predict_type = "prob"),
             folds = 3),
    "binary `outcome_type`"
  )
})
