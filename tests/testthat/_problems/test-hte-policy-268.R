# Extracted from test-hte-policy.R:268

# prequel ----------------------------------------------------------------------
skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")
x5 <- paste0("x", 1:5)

# test -------------------------------------------------------------------------
dat <- sim_hte(1500, dgp = "smooth", seed = 9)
sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
fr <- policy_frontier(sc, "tau_true", y0_hat = "mu0_true", budget = 0.3, curvature = c(0, 1, 5, 50))
expect_equal(fr$table$overlap_impact[1], 1)
expect_true(all(diff(fr$table$overlap_impact) <= 1e-12))
expect_true(all(diff(fr$table$overlap_deprivation) >= -1e-12))
expect_equal(colMeans(fr$selected), rep(0.3, 4), tolerance = 1e-3, ignore_attr = TRUE)
expect_s3_class(plot_policy_frontier(fr, 5), "ggplot")
expect_output(print(fr), "Targeting frontier")
dat$pos <- exp(dat$mu0_true / 4)
sc2 <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
fr2 <- policy_frontier(sc2, "tau_true", y0_hat = "pos", utility = "crra", curvature = c(0.5, 1, 2))
