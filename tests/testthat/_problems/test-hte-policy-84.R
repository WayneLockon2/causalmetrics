# Extracted from test-hte-policy.R:84

# prequel ----------------------------------------------------------------------
skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")
x5 <- paste0("x", 1:5)

# test -------------------------------------------------------------------------
dat <- sim_hte(1000, dgp = "smooth", seed = 4)
sc <- dr_scores(dat, "y", "d", x5, seed = 1)
blp <- cate_blp(sc, ~ x1 + I(x2^2), n_boot = 199, seed = 1)
m <- lm(sc$score ~ x1 + I(x2^2), data = dat)
expect_equal(blp$beta, unname(coef(m)), tolerance = 1e-10)
X <- model.matrix(m)
e <- resid(m)
hc1 <- solve(crossprod(X)) %*% crossprod(X * e) %*% solve(crossprod(X)) * nrow(X) / (nrow(X) - ncol(X))
expect_equal(blp$coefficients$std.error, sqrt(diag(hc1)), tolerance = 1e-10)
