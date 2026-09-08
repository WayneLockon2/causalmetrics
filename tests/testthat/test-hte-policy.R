skip_if_not_installed("mlr3")
skip_if_not_installed("mlr3learners")

x5 <- paste0("x", 1:5)

test_that("dr_scores reproduces the AIPW score and the ATE with known nuisances", {
  dat <- sim_hte(800, dgp = "smooth", seed = 1)
  sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  H <- dat$d / dat$p_true - (1 - dat$d) / (1 - dat$p_true)
  manual <- dat$mu1_true - dat$mu0_true + H * (dat$y - ifelse(dat$d == 1, dat$mu1_true, dat$mu0_true))
  expect_equal(sc$score, manual, tolerance = 1e-12)
  aipw <- est_aipw(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  expect_equal(sc$ate$estimate, aipw$estimate, tolerance = 1e-10)
  expect_equal(tidy(sc)$estimate, sc$ate$estimate)
  # the conditional mean of the score is the CATE: slope near one
  b <- coef(lm(sc$score ~ dat$tau_true))
  expect_lt(abs(b[2] - 1), 0.25)
  expect_output(print(sc), "Doubly robust pseudo-outcomes")
  expect_equal(glance(sc)$nobs, 800)
  # residuals for the R-learner
  expect_equal(sc$residuals$d_tilde, dat$d - dat$p_true)
  expect_equal(sc$residuals$y_tilde, dat$y - (dat$p_true * dat$mu1_true + (1 - dat$p_true) * dat$mu0_true))
})

test_that("dr_scores cross-fits internally and handles the alternative signals", {
  dat <- sim_hte(600, dgp = "smooth", seed = 2)
  sc <- dr_scores(dat, "y", "d", x5, folds = 3, seed = 1)
  expect_s3_class(sc, "cm_scores")
  expect_equal(length(unique(sc$fold_id)), 3L)
  expect_equal(nrow(sc$nuisance), 600)
  expect_true(all(sc$nuisance$p >= 0.01 & sc$nuisance$p <= 0.99))
  sc_ipw <- dr_scores(dat, "y", "d", x5, type = "ipw", seed = 1)
  sc_reg <- dr_scores(dat, "y", "d", x5, type = "reg", seed = 1)
  expect_equal(sc_reg$score, sc_reg$nuisance$mu1 - sc_reg$nuisance$mu0)
  expect_gt(sd(sc_ipw$score), sd(sc$score))
  expect_error(dr_scores(dat, "y", "d"), "`x` is required")
  dat$y[1] <- NA
  expect_error(dr_scores(dat, "y", "d", x5), "Missing values")
  sc_na <- dr_scores(dat, "y", "d", x5, na_action = "omit", seed = 1)
  expect_equal(sc_na$n, 599)
})

test_that("cate_learner runs every method, predicts on new data, and the DR-learner recovers a linear CATE", {
  dat <- sim_hte(1500, dgp = "smooth", seed = 3)
  train <- dat[1:1000, ]
  test <- dat[1001:1500, ]
  sc <- dr_scores(train, "y", "d", x5, seed = 1)
  for (m in c("dr", "r", "x", "t", "s")) {
    fit <- cate_learner(scores = sc, x_het = c("x1", "x2"), method = m)
    expect_s3_class(fit, "cm_cate")
    expect_equal(length(fit$tau_hat), 1000)
    pr <- predict(fit, test)
    expect_equal(length(pr), 500)
    expect_true(all(is.finite(pr)))
    expect_output(print(fit), "CATE model")
    expect_equal(nrow(tidy(fit)), 1000)
  }
  # linear final stage on x1 with the true nuisances: slope on x1 near one
  sc0 <- dr_scores(train, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  sc0$x <- x5
  fit <- cate_learner(scores = sc0, x_het = "x1", method = "dr")
  expect_lt(abs(coef(fit$final$learner$model)[["x1"]] - 1), 0.3)
  # R-learner final stage is the weighted regression of the ratio label
  fit_r <- cate_learner(scores = sc0, x_het = "x1", method = "r")
  w <- sc0$residuals$d_tilde^2
  lab <- sc0$residuals$y_tilde / sc0$residuals$d_tilde
  manual <- lm(lab ~ train$x1, weights = w)
  expect_equal(unname(coef(fit_r$final$learner$model)), unname(coef(manual)), tolerance = 1e-8)
  # the X-learner with adaptation and the S-learner with interactions run
  expect_s3_class(cate_learner(scores = sc, x_het = "x1", method = "x", adapt = TRUE), "cm_cate")
  expect_s3_class(cate_learner(scores = sc, method = "s", interactions = TRUE), "cm_cate")
  expect_error(cate_learner(train, "y", "d"), "Supply")
})

test_that("cate_blp equals OLS with HC1 and its bands and tests behave", {
  dat <- sim_hte(1000, dgp = "smooth", seed = 4)
  sc <- dr_scores(dat, "y", "d", x5, seed = 1)
  blp <- cate_blp(sc, ~ x1 + I(x2^2), n_boot = 199, seed = 1)
  m <- lm(sc$score ~ x1 + I(x2^2), data = dat)
  expect_equal(blp$beta, unname(coef(m)), tolerance = 1e-10)
  X <- model.matrix(m)
  e <- resid(m)
  hc1 <- solve(crossprod(X)) %*% crossprod(X * e) %*% solve(crossprod(X)) * nrow(X) / (nrow(X) - ncol(X))
  expect_equal(blp$coefficients$std.error, unname(sqrt(diag(hc1))), tolerance = 1e-10)
  expect_gte(blp$crit_val, blp$crit_pointwise)
  expect_lt(blp$test_heterogeneity$p.value, 0.01)
  grid <- data.frame(x1 = seq(-2, 2, by = 0.5), x2 = 0)
  pr <- predict(blp, grid, seed = 1)
  expect_equal(pr$estimate, as.numeric(cbind(1, grid$x1, 0) %*% blp$beta))
  expect_true(all(pr$band.low <= pr$conf.low))
  expect_s3_class(plot_cate_blp(blp, grid, seed = 1), "ggplot")
  # intercept only returns the ATE
  b0 <- cate_blp(sc, NULL)
  expect_equal(b0$beta, sc$ate$estimate)
  expect_equal(tidy(blp)$term, c("(Intercept)", "x1", "I(x2^2)"))
})

test_that("cate_gate returns group means of the score with a joint test", {
  dat <- sim_hte(1000, dgp = "smooth", seed = 5)
  sc <- dr_scores(dat, "y", "d", x5, seed = 1)
  g <- cut(dat$x1, quantile(dat$x1, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
  ga <- cate_gate(sc, g, n_boot = 199, seed = 1)
  expect_equal(ga$table$estimate, as.numeric(tapply(sc$score, g, mean)), tolerance = 1e-10)
  expect_equal(ga$table$n, as.integer(table(g)))
  expect_lt(ga$test_equal$p.value, 0.01)
  expect_s3_class(plot_cate_gate(ga), "ggplot")
  expect_output(print(ga), "Group average")
  dat$grp <- g
  sc2 <- dr_scores(dat, "y", "d", x5, seed = 1)
  expect_equal(cate_gate(sc2, "grp", n_boot = 99)$table$estimate, ga$table$estimate)
})

test_that("cate_score and cate_ensemble select and stack models", {
  dat <- sim_hte(1800, dgp = "smooth", seed = 6)
  train <- dat[1:1000, ]
  score_set <- dat[1001:1800, ]
  m_dr <- cate_learner(train, "y", "d", x5, x_het = c("x1", "x2"), method = "dr", seed = 1)
  m_t <- cate_learner(train, "y", "d", x5, x_het = c("x1", "x2"), method = "t", seed = 1)
  m_const <- cate_learner(train, "y", "d", x5, x_het = "x5", method = "t", seed = 1)
  sc <- dr_scores(score_set, "y", "d", x5, seed = 2)
  cs <- cate_score(sc, dr_learner = m_dr, t_learner = m_t, weak = m_const)
  expect_equal(cs$table$model, c("dr_learner", "t_learner", "weak"))
  expect_equal(cs$table$loss, unname(colMeans((sc$score - cs$predictions)^2)))
  expect_equal(cs$table$diff, cs$loss_constant - cs$table$loss, tolerance = 1e-10)
  expect_equal(cs$table$score_norm, 1 - cs$table$loss / cs$loss_constant)
  # the informative models beat the weak one out of sample
  expect_gt(cs$table$diff[1], cs$table$diff[3])
  cs2 <- cate_score(sc, dr_learner = m_dr, t_learner = m_t, baseline = "t_learner")
  expect_equal(cs2$table$diff[2], 0)
  expect_error(cate_score(sc, s = m_dr), "partially")
  expect_output(print(cs), "Doubly robust loss")
  # ensembles
  for (rule in c("q", "convex", "best", "ols")) {
    e <- cate_ensemble(sc, dr_learner = m_dr, t_learner = m_t, weak = m_const, method = rule)
    expect_s3_class(e, "cm_cate")
    expect_equal(length(predict(e, train)), 1000)
    if (rule != "ols") {
      expect_equal(sum(e$ensemble$weights), 1, tolerance = 1e-8)
      expect_true(all(e$ensemble$weights >= 0))
    }
  }
  e_best <- cate_ensemble(sc, dr_learner = m_dr, t_learner = m_t, method = "best")
  expect_equal(max(e_best$ensemble$weights), 1)
  # convex weights solve the same quadratic program as quadprog
  skip_if_not_installed("quadprog")
  e_cv <- cate_ensemble(sc, dr_learner = m_dr, t_learner = m_t, weak = m_const, method = "convex")
  P <- sweep(cs$predictions, 2, colMeans(cs$predictions))
  tg <- sc$score - mean(sc$score)
  M <- ncol(P)
  qp <- quadprog::solve.QP(2 * crossprod(P) / nrow(P) + diag(1e-10, M), 2 * as.numeric(crossprod(P, tg)) / nrow(P),
                           cbind(rep(1, M), diag(M)), c(1, rep(0, M)), meq = 1)
  expect_equal(unname(e_cv$ensemble$weights), qp$solution, tolerance = 1e-4)
  expect_equal(e_cv$ensemble$loss, mean((tg - P %*% qp$solution)^2), tolerance = 1e-8)
})

test_that("cate_validate reproduces the score regression and matches grf's rank-weighted effects", {
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
  # the top-q set treats units strictly above the type-1 quantile threshold
  top <- dat$tau_true > quantile(dat$tau_true, 0.8, type = 1)
  expect_equal(cv$toc[cv$q == 0.2], mean(sc$score[top]) - mean(sc$score), tolerance = 1e-10)
  expect_true(all(cv$toc_band_low <= cv$toc_low + 1e-12))
  expect_gt(v$areas$conf.low[1], 0)
  expect_equal(nrow(v$group_diff), 5)
  expect_output(print(v), "Validation of a CATE model")
  for (w in c("calibration", "toc", "qini")) expect_s3_class(plot_cate_validation(v, w), "ggplot")
  # constant predictions: ties handled, curves zero
  v0 <- cate_validate(sc, rep(1, 1200), n_boot = 49)
  expect_true(all(abs(v0$curves$toc) < 1e-10))

  skip_if_not_installed("grf")
  X <- as.matrix(dat[, x5])
  cf <- grf::causal_forest(X, dat$y, dat$d, num.trees = 300, seed = 1)
  tau_cf <- predict(cf)$predictions
  sc_cf <- dr_scores(dat, "y", "d", p_hat = cf$W.hat, mu0_hat = cf$Y.hat - cf$W.hat * tau_cf,
                     mu1_hat = cf$Y.hat + (1 - cf$W.hat) * tau_cf, p_clip = c(0, 1))
  expect_equal(sc_cf$score, as.numeric(grf::get_scores(cf)), tolerance = 1e-10)
  prio <- round(dat$x1, 1)
  vg <- cate_validate(sc_cf, prio, quantile_grid = seq(0.1, 1, by = 0.1), n_boot = 49, seed = 1)
  r_autoc <- grf::rank_average_treatment_effect(cf, prio, target = "AUTOC", q = seq(0.1, 1, by = 0.1), R = 50)
  r_qini <- grf::rank_average_treatment_effect(cf, prio, target = "QINI", q = seq(0.1, 1, by = 0.1), R = 50)
  expect_equal(vg$curves$toc, r_autoc$TOC$estimate, tolerance = 1e-8)
  expect_equal(vg$areas$estimate, c(r_autoc$estimate, r_qini$estimate), tolerance = 1e-8)
})

test_that("policy_value and policy_learn recover the optimal tree and agree with policytree", {
  dat <- sim_hte(3000, dgp = "policy", seed = 8)
  sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  sc$x <- x5
  oracle <- as.numeric(dat$tau_true > 0)
  pv <- policy_value(sc, list(oracle = oracle, all = rep(1, 3000), none = rep(0, 3000)))
  expect_equal(pv$estimate[1], mean(oracle * sc$score))
  expect_equal(pv$estimate[3], 0)
  expect_equal(pv$diff_vs_first[2], pv$estimate[2] - pv$estimate[1])
  pv_boot <- policy_value(sc, list(oracle = oracle, none = rep(0, 3000)), n_boot = 200, seed = 1)
  B <- attr(pv_boot, "boot")
  expect_equal(dim(B), c(200, 2))
  expect_equal(pv_boot$boot.low[1], unname(quantile(B[, 1], 0.025)))
  expect_lt(abs(pv_boot$boot.se[1] - pv_boot$std.error[1]) / pv_boot$std.error[1], 0.35)
  expect_equal(unname(pv_boot$boot.se[2]), 0)
  expect_equal(dim(attr(pv_boot, "boot_diff")), c(200, 2))
  # weighted values: the adaptively weighted estimator and its standard error
  w <- 1 / seq_len(3000)^0.5
  pv_w <- policy_value(sc, list(oracle = oracle, none = rep(0, 3000)), weights = w, n_boot = 50, seed = 1)
  c_i <- oracle * sc$score
  expect_equal(pv_w$estimate[1], sum(w * c_i) / sum(w))
  expect_equal(pv_w$std.error[1], sqrt(sum(w^2 * (c_i - pv_w$estimate[1])^2)) / sum(w))
  expect_equal(pv_w$diff_vs_first[2], -pv_w$estimate[1])
  expect_equal(dim(attr(pv_w, "boot")), c(50, 2))
  expect_error(policy_value(sc, oracle, weights = -w), "weights")
  pv_all <- policy_value(sc, oracle, baseline = "all")
  expect_equal(pv_all$estimate, mean((oracle - 1) * sc$score))
  expect_equal(policy_value(sc, oracle, cost = 0.5)$estimate, mean(oracle * (sc$score - 0.5)))

  pol <- policy_learn(sc, method = "tree", depth = 2, seed = 1)
  expect_s3_class(pol, "cm_policy")
  expect_equal(pol$rule$var, "x1")
  expect_lt(abs(pol$rule$threshold), 0.1)
  expect_gt(mean(pol$assign == oracle), 0.95)
  expect_equal(nrow(pol$value), 2)
  expect_equal(predict(pol, dat), pol$assign)
  expect_output(print(pol), "Treatment policy")
  expect_s3_class(plot_policy_tree(pol), "ggplot")
  expect_equal(tidy(pol), pol$value)
  # depth 1 exact search equals brute force over all splits of x1
  p1 <- policy_learn(sc, x = "x1", depth = 1, holdout = 0)
  g <- sc$score
  brute <- max(sapply(sort(unique(dat$x1)), function(t) {
    l <- dat$x1 <= t
    max(sum(g[l]), 0) + max(sum(g[!l]), 0)
  }), max(sum(g), 0))
  expect_equal(p1$rule$value, brute, tolerance = 1e-10)
  expect_equal(p1$value$value_vs_none, brute / 3000, tolerance = 1e-10)
  # other methods run
  for (m in c("linear", "classifier")) {
    p <- policy_learn(sc, method = m, seed = 1)
    expect_true(all(p$assign %in% c(0, 1)))
    expect_gt(p$value$value_vs_none[2], 0)
  }
  expect_s3_class(policy_learn(sc, method = "linear", seed = 1)$model, "glm")
  expect_true(inherits(policy_learn(sc, method = "classifier", seed = 1)$model, "Learner"))
  expect_equal(pol$model$var, "x1")
  m_dr <- cate_learner(scores = sc, x_het = x5, method = "dr")
  pb <- policy_learn(sc, method = "budget", budget = 0.3, tau_hat = m_dr, seed = 1)
  expect_lt(abs(mean(pb$assign) - 0.3), 0.05)
  expect_error(policy_learn(sc, method = "budget"), "budget")
  expect_error(policy_learn(sc, depth = 3), "depth 1 or 2")

  # weighted tree search equals the unweighted search on duplicated rows
  wi <- c(rep(2, 400), rep(1, 800))
  sc_w <- dr_scores(dat[1:1200, ], "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  sc_dup <- dr_scores(dat[c(1:1200, 1:400), ], "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  tw <- policy_learn(sc_w, x = c("x1", "x2", "x3"), depth = 2, holdout = 0, weights = wi, max_root_splits = 1e6)
  td <- policy_learn(sc_dup, x = c("x1", "x2", "x3"), depth = 2, holdout = 0, max_root_splits = 1e6)
  expect_equal(tw$rule$threshold, td$rule$threshold)
  expect_equal(tw$rule$value * mean(wi), td$rule$value, tolerance = 1e-8)
  expect_equal(tw$value$value_vs_none, policy_value(sc_w, tw$assign, weights = wi)$estimate)
  expect_equal(length(tw$weights), 1200)

  skip_if_not_installed("policytree")
  small <- dat[1:300, ]
  scs <- dr_scores(small, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  scs$x <- x5
  a <- policy_learn(scs, x = c("x1", "x2", "x3"), depth = 2, holdout = 0, max_root_splits = 1e6)
  b <- policy_learn(scs, x = c("x1", "x2", "x3"), depth = 2, holdout = 0, engine = "policytree")
  expect_equal(a$value$value_vs_none, b$value$value_vs_none, tolerance = 1e-10)
  a1 <- policy_learn(scs, x = c("x1", "x2", "x3"), depth = 1, holdout = 0)
  b1 <- policy_learn(scs, x = c("x1", "x2", "x3"), depth = 1, holdout = 0, engine = "policytree")
  expect_equal(a1$value$value_vs_none, b1$value$value_vs_none, tolerance = 1e-10)
  expect_equal(length(predict(b, small)), 300)
  b_step <- policy_learn(scs, x = c("x1", "x2", "x3"), depth = 2, holdout = 0, engine = "policytree",
                         split_step = 5, weights = rep(c(1, 2), 150))
  expect_true(all(b_step$assign %in% c(0, 1)))
})

test_that("policy_frontier moves from impact-only to deprivation-only targeting", {
  dat <- sim_hte(1500, dgp = "smooth", seed = 9)
  sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  fr <- policy_frontier(sc, "tau_true", y0_hat = "mu0_true", budget = 0.3, curvature = c(0, 1, 5, 50))
  expect_equal(fr$table$overlap_impact[1], 1)
  expect_true(all(diff(fr$table$overlap_impact) <= 1e-12))
  expect_true(all(diff(fr$table$overlap_deprivation) >= -1e-12))
  expect_equal(colMeans(fr$selected), rep(0.3, 4), tolerance = 1e-3, ignore_attr = TRUE)
  expect_s3_class(plot_policy_frontier(fr, 5), "ggplot")
  expect_output(print(fr), "Targeting frontier")
  dat$pos <- exp(dat$mu0_true / 4) + 5
  sc2 <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true")
  fr2 <- policy_frontier(sc2, "tau_true", y0_hat = "pos", utility = "crra", curvature = c(0.5, 1, 2))
  expect_equal(nrow(fr2$table), 3)
})

test_that("sim_hte designs have the documented properties", {
  d1 <- sim_hte(400, "simple_cate", seed = 1)
  expect_true(all(d1$tau_true == 0.5))
  d2 <- sim_hte(400, "complex_cate", seed = 1)
  expect_true(all(d2$tau_true %in% c(0, 0.5)))
  d3 <- sim_hte(400, "unbalanced", seed = 1)
  expect_gt(mean(d3$d), 0.85)
  d4 <- sim_hte(400, "binary_outcome", seed = 1)
  expect_true(all(d4$y %in% c(0, 1)))
  d5 <- sim_hte(400, "policy", seed = 1)
  expect_type(attr(d5, "optimal_tree"), "list")
})

test_that("bind_scores pools score objects and glance.cm_blp feeds modelsummary", {
  dat <- sim_hte(900, dgp = "smooth", seed = 11)
  a <- dr_scores(dat[1:300, ], "y", "d", x5, seed = 1)
  b <- dr_scores(dat[301:600, ], "y", "d", x5, seed = 2)
  cc <- dr_scores(dat[601:900, ], "y", "d", x5, seed = 3)
  pooled <- bind_scores(first = a, second = b, cc)
  expect_s3_class(pooled, "cm_scores")
  expect_equal(pooled$n, 900)
  expect_equal(pooled$score, c(a$score, b$score, cc$score))
  expect_equal(as.character(unique(pooled$data$set)), c("first", "second", "set3"))
  expect_equal(length(unique(pooled$fold_id)), 15L)
  expect_equal(pooled$ate$estimate, mean(pooled$score))
  expect_equal(bind_scores(list(a, b))$n, 600)
  expect_error(bind_scores(a, dr_scores(dat[1:300, ], "y", "d", x5, type = "reg", seed = 1)), "same `type`")
  # projection on a set-level dictionary
  comp <- data.frame(set = c("first", "second", "set3"), z = c(-1, 0, 1))
  extra <- comp[match(pooled$data$set, comp$set), "z", drop = FALSE]
  blp <- cate_blp(pooled, ~ z, data = extra, uniform = FALSE)
  m <- lm(pooled$score ~ extra$z)
  expect_equal(blp$beta, unname(coef(m)), tolerance = 1e-10)
  expect_equal(blp$r.squared, summary(m)$r.squared, tolerance = 1e-6)
  expect_equal(blp$adj.r.squared, summary(m)$adj.r.squared, tolerance = 1e-6)
  g <- glance(blp)
  expect_equal(g$nobs, 900)
  skip_if_not_installed("modelsummary")
  tab <- modelsummary::modelsummary(list(blp), output = "data.frame")
  expect_true(any(grepl("z", tab$term)))
  expect_true(any(grepl("R2", tab$term)))
})

test_that("dr_scores trims on the propensity and records the share dropped", {
  dat <- sim_hte(800, dgp = "smooth", seed = 12)
  dat$p_true[1:10] <- 0.005
  sc <- dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true", trim = c(0.01, 0.99))
  expect_equal(sc$n, 790)
  expect_equal(sc$diagnostics$trimmed_share, 10 / 800)
  expect_equal(nrow(sc$data), 790)
  expect_equal(length(sc$fold_id), 790)
  expect_output(print(sc), "trimmed share")
  expect_error(dr_scores(dat, "y", "d", p_hat = "p_true", mu0_hat = "mu0_true", mu1_hat = "mu1_true", trim = c(0.5, 0.4)), "increasing")
})

test_that("multi-arm scores reproduce the IPS reward, agree with policytree, and contrast to pairwise scores", {
  set.seed(21)
  n <- 3000
  x1 <- rnorm(n); x2 <- rnorm(n)
  grp <- sample(c("a", "b", "c"), n, TRUE, c(0.6, 0.25, 0.15))
  arm <- sample(c("7", "14", "30"), n, TRUE, c(0.15, 0.15, 0.7))
  mu <- cbind(`7` = 0.15 + 0.05 * (grp == "a") + 0.02 * x1, `14` = 0.15 + 0.05 * (grp == "b") + 0.02 * x1,
              `30` = 0.15 + 0.05 * (grp == "c") + 0.02 * x1)
  y <- rbinom(n, 1, mu[cbind(seq_len(n), match(arm, colnames(mu)))])
  dat <- data.frame(y = y, arm = arm, x1 = x1, x2 = x2, ga = as.integer(grp == "a"), gb = as.integer(grp == "b"))
  xs <- c("x1", "x2", "ga", "gb")
  e <- matrix(rep(c(0.15, 0.15, 0.7), each = n), n, dimnames = list(NULL, c("7", "14", "30")))
  sc <- dr_scores(dat, "y", "arm", xs, arms = c("7", "14", "30"), p_hat = e, seed = 1)
  expect_s3_class(sc, "cm_scores_multi")
  expect_equal(dim(sc$gamma), c(n, 3))
  I <- matrix(0, n, 3); I[cbind(seq_len(n), match(arm, sc$arms))] <- 1
  expect_equal(sc$gamma, sc$nuisance$mu + I * (y - sc$nuisance$mu) / sc$nuisance$e, ignore_attr = TRUE)
  expect_equal(nrow(tidy(sc)), 3)
  expect_output(print(sc), "3 arms")
  # IPS reward of the off-policy literature
  sci <- dr_scores(dat, "y", "arm", xs, arms = c("7", "14", "30"), p_hat = e, type = "ipw")
  pol <- ifelse(dat$ga == 1, "7", ifelse(dat$gb == 1, "14", "30"))
  ips <- sum((arm == pol) * y / e[cbind(seq_len(n), match(pol, colnames(e)))]) / n
  pv <- policy_value(sci, list(rule = pol, all30 = rep("30", n)), baseline = "30")
  expect_equal(pv$estimate[1] + pv$baseline_value[1], ips, tolerance = 1e-12)
  expect_equal(pv$estimate[2], 0)
  expect_equal(pv$gain_pct[1], 100 * pv$estimate[1] / pv$baseline_value[1])
  expect_equal(unname(unlist(pv[1, c("share_7", "share_14", "share_30")])), as.numeric(prop.table(table(factor(pol, levels = c("7", "14", "30"))))))
  expect_error(policy_value(sc, pol, baseline = "99"), "one of the arms")
  # trees over three arms
  p1 <- policy_learn(sc, x = xs, depth = 1, holdout = 0)
  expect_true(all(p1$assign %in% sc$arms))
  small <- dat[1:400, ]
  scs <- dr_scores(small, "y", "arm", xs, arms = c("7", "14", "30"), p_hat = e[1:400, ], seed = 1)
  a <- policy_learn(scs, x = xs, depth = 2, holdout = 0, max_root_splits = 1e6)
  expect_output(print(a), "3 arms")
  skip_if_not_installed("policytree")
  b <- policy_learn(scs, x = xs, depth = 2, holdout = 0, engine = "policytree")
  expect_equal(a$value$value, b$value$value, tolerance = 1e-10)
  expect_true(all(predict(b, small) %in% sc$arms))
  # a pairwise contrast equals the pairwise binary object built from the same nuisances
  cs <- contrast_scores(sc, "7", "30")
  expect_s3_class(cs, "cm_scores")
  rows <- arm %in% c("7", "30")
  pair <- dat[rows, ]; pair$W <- as.integer(pair$arm == "7")
  ref <- dr_scores(pair, "y", "W", p_hat = rep(0.15 / 0.85, nrow(pair)),
                   mu0_hat = sc$nuisance$mu[rows, "30"], mu1_hat = sc$nuisance$mu[rows, "7"])
  expect_equal(cs$score, ref$score, tolerance = 1e-12)
  expect_equal(cs$d, ref$d)
  expect_s3_class(cate_learner(scores = cs, x_het = xs, method = "dr"), "cm_cate")
  expect_error(contrast_scores(sc, "7", "7"), "two different arms")
})
