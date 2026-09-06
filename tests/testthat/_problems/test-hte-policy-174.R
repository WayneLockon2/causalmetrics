# Extracted from test-hte-policy.R:174

# prequel ----------------------------------------------------------------------
skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")
x5 <- paste0("x", 1:5)

# test -------------------------------------------------------------------------
dat <- sim_hte(1200, dgp = "smooth", seed = 7)
sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
sc$x <- x5
v <- cate_validate(sc, "tau_true", n_boot = 99, seed = 1)
m <- lm(sc$score ~ I(dat$tau_true - mean(dat$tau_true)))
expect_equal(v$blp$estimate, unname(coef(m)), tolerance = 1e-10)
expect_lt(abs(v$blp$estimate[2] - 1), 3 * v$blp$std.error[2])
cal <- v$calibration$table
expect_equal(nrow(cal), 4)
expect_equal(sum(cal$n), 1200)
expect_true(all(diff(cal$mean_tau) > 0))
expect_lt(v$calibration$test_equal$p.value, 0.01)
cv <- v$curves
expect_equal(cv$q[nrow(cv)], 1)
expect_equal(cv$toc[nrow(cv)], 0, tolerance = 1e-10)
expect_equal(cv$qini, cv$toc * cv$share, tolerance = 1e-10)
top <- dat$tau_true >= quantile(dat$tau_true, 0.8, type = 1)
expect_equal(cv$toc[cv$q == 0.2], mean(sc$score[top]) - mean(sc$score), tolerance = 1e-10)
