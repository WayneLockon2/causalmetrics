sim_tidy_data <- function(n = 600, seed = 11) {
  set.seed(seed)
  x1 <- rnorm(n)
  p <- plogis(0.5 * x1)
  d <- rbinom(n, 1, p)
  y <- 1 + x1 + 2 * d + rnorm(n)
  data.frame(y = y, d = d, x1 = x1, p = p, mu0 = 1 + x1, mu1 = 3 + x1)
}

test_that("tidy and glance follow broom conventions for cm_aipw", {
  dat <- sim_tidy_data()
  a <- est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")

  td <- tidy(a)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 1L)
  expect_named(td, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_equal(td$term, "d")
  expect_equal(td$estimate, a$estimate)
  expect_equal(td$statistic, a$estimate / a$std.error)
  expect_true(td$p.value >= 0 && td$p.value <= 1)
  expect_equal(td$conf.low, a$conf.low)

  td90 <- tidy(a, conf.level = 0.9)
  expect_gt(td90$conf.low, td$conf.low)
  expect_named(tidy(a, conf.int = FALSE), c("term", "estimate", "std.error", "statistic", "p.value"))

  gl <- glance(a)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$nobs, a$n)
  expect_equal(gl$estimand, "ATE")
  expect_equal(gl$method, "AIPW")
})

test_that("tidy and glance work for cm_dml and feed modelsummary", {
  dat <- sim_tidy_data()
  a <- est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
  f <- est_dml(dat, "y", "d", model = "irm", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")

  td <- tidy(f)
  expect_equal(td$term, "d")
  expect_equal(td$estimate, f$estimate)
  gl <- glance(f)
  expect_equal(gl$model, "irm")
  expect_equal(gl$nobs, f$n)
  expect_true("rmse_p_hat" %in% names(gl))

  withr_opts <- options(modelsummary_get = "broom")
  on.exit(options(withr_opts), add = TRUE)
  ms <- modelsummary::modelsummary(list("AIPW" = a, "DML" = f), output = "data.frame")
  expect_true(any(ms$term == "d"))
  expect_true(any(grepl("Num.Obs", ms$term)))
})
