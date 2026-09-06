# Extracted from test-hte-policy.R:123

# prequel ----------------------------------------------------------------------
skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")
x5 <- paste0("x", 1:5)

# test -------------------------------------------------------------------------
dat <- sim_hte(1800, dgp = "smooth", seed = 6)
train <- dat[1:1000, ]
score_set <- dat[1001:1800, ]
m_dr <- cate_learner(train, "y", "d", x5, x_het = c("x1", "x2"), method = "dr", seed = 1)
m_t <- cate_learner(train, "y", "d", x5, x_het = c("x1", "x2"), method = "t", seed = 1)
m_const <- cate_learner(train, "y", "d", x5, x_het = "x5", method = "t", seed = 1)
sc <- dr_scores(score_set, "y", "d", x5, seed = 2)
cs <- cate_score(sc, dr_learner = m_dr, t_learner = m_t, weak = m_const)
expect_equal(cs$table$model, c("dr_learner", "t_learner", "weak"))
expect_equal(cs$table$loss, colMeans((sc$score - cs$predictions)^2))
