# Extracted from test-hte-policy.R:103

# prequel ----------------------------------------------------------------------
skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")
x5 <- paste0("x", 1:5)

# test -------------------------------------------------------------------------
dat <- sim_hte(1000, dgp = "smooth", seed = 5)
sc <- dr_scores(dat, "y", "d", x5, seed = 1)
g <- cut(dat$x1, quantile(dat$x1, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
ga <- cate_gate(sc, g, n_boot = 199, seed = 1)
expect_equal(ga$table$estimate, unname(tapply(sc$score, g, mean)), tolerance = 1e-10)
