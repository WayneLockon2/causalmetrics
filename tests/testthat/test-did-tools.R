test_that("Bacon decomposition reproduces the static TWFE coefficient", {
  dat <- sim_did_panel(n_units = 200, n_periods = 8, groups = c(3, 6), seed = 1)
  b <- bacon_decomp(dat, id = "id", time = "time", y = "y", d = "treated")
  expect_equal(b$sum_weights, 1, tolerance = 1e-8)
  expect_equal(b$check, b$twfe, tolerance = 1e-8)
  expect_true(all(b$decomposition$weight >= 0))
  expect_setequal(b$by_type$type, c("treated vs never treated", "earlier vs later treated", "later vs earlier treated"))
  expect_output(print(b), "Goodman-Bacon")
  # without never-treated units the decomposition still holds
  dat2 <- dat[dat$g > 0, ]
  b2 <- bacon_decomp(dat2, id = "id", time = "time", y = "y", d = "treated")
  expect_equal(b2$check, b2$twfe, tolerance = 1e-8)
})

test_that("dCDH weights reproduce the TWFE coefficient and flag negative weights", {
  dat <- sim_did_panel(n_units = 200, n_periods = 10, groups = c(3, 5, 8), never_share = 0.1, seed = 2)
  w <- twfe_weights(dat, id = "id", time = "time", d = "treated", y = "y")
  expect_equal(sum(w$cells$weight), 1, tolerance = 1e-10)
  expect_equal(w$check, w$twfe, tolerance = 1e-8)
  expect_gt(w$share_negative, 0)
  expect_output(print(w), "negative weight")
})

test_that("imputation estimator recovers the dynamic effects and matches a manual fixest imputation", {
  dat <- sim_did_panel(n_units = 400, n_periods = 8, groups = c(3, 6), seed = 3)
  fit <- did_imputation(dat, id = "id", time = "time", group = "g", y = "y", pre_window = 2)
  truth <- attr(dat, "att_gt")
  # manual: first stage on untreated cells, average tau by event time
  untreated <- dat[dat$treated == 0, ]
  s1 <- fixest::feols(y ~ 1 | id + time, data = untreated)
  tr <- dat[dat$treated == 1, ]
  tr$tau <- tr$y - predict(s1, newdata = tr)
  manual <- tapply(tr$tau, tr$event_time, mean)
  expect_equal(fit$by_event$estimate, as.numeric(manual), tolerance = 1e-8)
  expect_equal(fit$overall$estimate, mean(tr$tau), tolerance = 1e-8)
  # close to the truth (effects mu * (e + 1) averaged over cohorts present)
  truth_e <- tapply(truth$att, truth$time - truth$group, mean)
  expect_lt(max(abs(fit$by_event$estimate - unname(truth_e))), 0.5)
  expect_true(all(is.finite(fit$by_event$std.error)))
  expect_equal(nrow(fit$pre), 2)
  expect_output(print(fit), "Imputation")
  expect_s3_class(tidy(fit), "data.frame")
})

test_that("imputation standard errors are close to the sampling variability", {
  ests <- replicate(40, {
    d <- sim_did_panel(n_units = 300, n_periods = 6, groups = 4, never_share = 0.5)
    f <- did_imputation(d, id = "id", time = "time", group = "g", y = "y", pre_window = 0)
    c(f$overall$estimate, f$overall$std.error)
  })
  expect_lt(abs(sd(ests[1, ]) / mean(ests[2, ]) - 1), 0.4)
})

test_that("event_study_frame combines estimators and plots", {
  dat <- sim_did_panel(n_units = 300, n_periods = 8, groups = c(3, 6), seed = 4)
  fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", n_boot = 49, seed = 1)
  cs <- aggregate_att(fit, type = "dynamic", n_boot = 49, seed = 1)
  dat$rel <- ifelse(dat$g > 0, dat$time - dat$g, -1000)
  tw <- fixest::feols(y ~ i(rel, ref = c(-1, -1000)) | id + time, data = dat, cluster = ~id)
  imp <- did_imputation(dat, id = "id", time = "time", group = "g", y = "y", pre_window = 2)
  fr <- event_study_frame("CS" = cs, "TWFE" = tw, "Imputation" = imp, ref = -1)
  expect_equal(levels(fr$method), c("CS", "TWFE", "Imputation"))
  expect_true(all(c("event_time", "estimate", "std.error", "conf.low", "band.low") %in% names(fr)))
  expect_equal(fr$estimate[fr$method == "TWFE" & fr$event_time == -1], 0)
  expect_true(any(is.finite(fr$band.low[fr$method == "CS"])))
  p <- plot_event_study(fr)
  expect_s3_class(p, "ggplot")
  # sunab coefficients are parsed too
  sa <- fixest::feols(y ~ sunab(g, time) | id + time, data = dat[dat$g > 0 | TRUE, ], cluster = ~id)
  fr2 <- event_study_frame("Sun-Abraham" = sa)
  expect_gt(nrow(fr2), 3)
})

test_that("pretrend_power behaves sensibly", {
  dat <- sim_did_panel(n_units = 300, n_periods = 8, groups = c(3, 6), seed = 5)
  cs <- aggregate_att(att_gt(dat, id = "id", time = "time", group = "g", y = "y", n_boot = 0), type = "dynamic", n_boot = 0)
  pp <- pretrend_power(cs, slope = c(0, 0.2), n_sim = 2000, seed = 1)
  expect_lt(pp$slopes[[1]]$power, 0.25)
  expect_gt(pp$slopes[[2]]$power, pp$slopes[[1]]$power)
  expect_true(all(is.finite(pp$detectable$slope)))
  expect_gt(pp$detectable$slope[2], pp$detectable$slope[1])
  expect_equal(pp$slopes[[2]]$bias$unconditional_bias, 0.2 * (pp$slopes[[2]]$bias$event_time + 1))
  expect_output(print(pp), "Pre-trends test power")
})

test_that("permutation test rejects a large effect and not a null one", {
  dat <- sim_did_panel(n_units = 40, n_periods = 6, groups = 4, mu = 5, seed = 6)
  pt <- did_permutation_test(dat, id = "id", time = "time", group = "g", y = "y", n_perm = 199, seed = 1)
  expect_lt(pt$p.value, 0.05)
  null <- sim_did_panel(n_units = 40, n_periods = 6, groups = 4, mu = 0, seed = 7)
  pt0 <- did_permutation_test(null, id = "id", time = "time", group = "g", y = "y", n_perm = 199, seed = 1)
  expect_gt(pt0$p.value, 0.05)
  expect_output(print(pt), "Permutation test")
})
