truth_of <- function(dat, what) unlist(attr(dat, "truth")[what])

test_that("mediate_reg reproduces the product of coefficients and matches the linear truth", {
  dat <- sim_mediation(3000, dgp = "linear", seed = 1)
  fit <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), method = "delta")
  a1 <- coef(lm(m ~ d + x1 + x2, data = dat))[["d"]]
  b <- coef(lm(y ~ d + m + x1 + x2, data = dat))
  expect_equal(fit$effects$estimate[fit$effects$term == "nie"], a1 * b[["m"]], tolerance = 1e-10)
  expect_equal(fit$effects$estimate[fit$effects$term == "nde"], b[["d"]], tolerance = 1e-10)
  expect_equal(fit$effects$estimate[fit$effects$term == "total"], coef(lm(y ~ d + x1 + x2, data = dat))[["d"]], tolerance = 1e-10)
  tr <- truth_of(dat, c("total", "nde", "nie"))
  est <- fit$effects$estimate[match(c("total", "nde", "nie"), fit$effects$term)]
  se <- fit$effects$std.error[match(c("total", "nde", "nie"), fit$effects$term)]
  expect_true(all(abs(est - tr) < 3 * se))
  # Sobel variance
  va1 <- causalmetrics:::.cm_vcov_robust(lm(m ~ d + x1 + x2, data = dat))["d", "d"]
  vb2 <- causalmetrics:::.cm_vcov_robust(lm(y ~ d + m + x1 + x2, data = dat))["m", "m"]
  expect_equal(se[3], sqrt(a1^2 * vb2 + b[["m"]]^2 * va1), tolerance = 1e-10)
  # shares: proportion mediated equals S1 in the linear case; delta-method SEs are finite
  expect_equal(fit$shares$estimate[1], fit$shares$estimate[2], tolerance = 1e-10)
  expect_true(all(is.finite(fit$shares$std.error)))
  # factor covariates are expanded as the regressions expand them
  dat$g <- factor(sample(letters[1:3], nrow(dat), replace = TRUE))
  ff <- mediate_reg(dat, "y", "d", "m", x = c("x1", "g"), method = "delta")
  expect_equal(ff$effects$estimate[3], coef(lm(m ~ d + x1 + g, data = dat))[["d"]] * coef(lm(y ~ d + m + x1 + g, data = dat))[["m"]], tolerance = 1e-10)
  # simulation and bootstrap inference agree with the delta method within reason
  fs <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), method = "simulation", n_sim = 500, seed = 1)
  fb <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), method = "bootstrap", n_boot = 99, seed = 1)
  expect_lt(abs(fs$effects$std.error[3] / se[3] - 1), 0.25)
  expect_lt(abs(fb$effects$std.error[3] / se[3] - 1), 0.35)
  expect_output(print(fit), "Causal mediation by regression")
  expect_s3_class(plot_mediation(fit), "ggplot")
  expect_equal(nrow(tidy(fit)), nrow(fit$effects))
  expect_error(mediate_reg(dat, "y", "d", "m", interaction = TRUE, method = "delta"), "delta")
})

test_that("mediate_reg handles interactions, binary models, and parallel mediators", {
  di <- sim_mediation(3000, dgp = "interaction", seed = 2)
  fi <- mediate_reg(di, "y", "d", "m", x = c("x1", "x2"), interaction = TRUE, n_sim = 300, seed = 1)
  tr <- truth_of(di, c("total", "nde", "nie", "nde_total", "nie_pure", "cde"))
  est <- fi$effects$estimate[match(names(tr), fi$effects$term)]
  se <- fi$effects$std.error[match(names(tr), fi$effects$term)]
  expect_true(all(abs(est - tr) < 3.5 * se))
  expect_false(isTRUE(all.equal(est[2], est[4])))
  db <- sim_mediation(3000, dgp = "binary", seed = 3)
  fb <- mediate_reg(db, "y", "d", "m", x = c("x1", "x2"), outcome = "logit", mediator = "logit", n_sim = 300, seed = 1)
  tr <- truth_of(db, c("total", "nde", "nie"))
  est <- fb$effects$estimate[match(names(tr), fb$effects$term)]
  se <- fb$effects$std.error[match(names(tr), fb$effects$term)]
  expect_true(all(abs(est - tr) < 3.5 * se))
  # continuous mediator with a binary outcome uses Monte Carlo integration
  db$mc <- db$m + rnorm(nrow(db), sd = 0.1)
  fmc <- mediate_reg(db, "y", "d", "mc", x = c("x1", "x2"), outcome = "logit", n_sim = 100, seed = 1)
  expect_true(all(is.finite(fmc$effects$estimate)))
  dp <- sim_mediation(3000, dgp = "parallel", seed = 4)
  fp <- mediate_reg(dp, "y", "d", c("m1", "m2"), x = c("x1", "x2"), n_sim = 300, seed = 1)
  expect_true(all(c("nie_m1", "nie_m2") %in% fp$effects$term))
  e <- fp$effects$estimate
  expect_equal(e[fp$effects$term == "nie"], e[fp$effects$term == "nie_m1"] + e[fp$effects$term == "nie_m2"], tolerance = 1e-8)
  tr <- truth_of(dp, c("nie_m1", "nie_m2"))
  expect_lt(abs(e[fp$effects$term == "nie_m1"] - tr[1]), 0.15)
  expect_lt(abs(e[fp$effects$term == "nie_m2"] - tr[2]), 0.2)
})

test_that("mediate_sensitivity recovers the indirect effect at the true error correlation", {
  dc <- sim_mediation(6000, dgp = "correlated_errors", rho = 0.5, seed = 5)
  fit <- mediate_reg(dc, "y", "d", "m", x = c("x1", "x2"), method = "delta")
  sens <- mediate_sensitivity(fit, rho = c(0, 0.5))
  expect_lt(abs(sens$curve$nie[2] - truth_of(dc, "nie")), 0.1)
  expect_gt(abs(sens$curve$nie[1] - truth_of(dc, "nie")), 0.3)
  expect_true(is.finite(sens$rho_zero))
  expect_s3_class(plot_mediate_sensitivity(sens), "ggplot")
  expect_output(print(sens), "Sensitivity")
  dp <- sim_mediation(500, dgp = "parallel", seed = 4)
  expect_error(mediate_sensitivity(mediate_reg(dp, "y", "d", c("m1", "m2"), n_sim = 50)), "one mediator")
})

test_that("mediate_dml recovers the truth with linear learners and with forests", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  dat <- sim_mediation(3000, dgp = "linear", seed = 1)
  fit <- mediate_dml(dat, "y", "d", "m", x = c("x1", "x2"), seed = 1)
  expect_s3_class(fit, "cm_med_dml")
  expect_equal(fit$effects$term, c("total", "nde", "nie", "nde_total", "nie_pure"))
  tr <- truth_of(dat, c("total", "nde", "nie", "nde_total", "nie_pure"))
  expect_true(all(abs(fit$effects$estimate - tr) < 3.5 * fit$effects$std.error))
  # the four potential-outcome means are consistent with the effects
  po <- fit$potential$estimate
  expect_equal(fit$effects$estimate[1], po[1] - po[2], tolerance = 1e-10)
  expect_equal(nrow(fit$psi), fit$n_used)
  expect_output(print(fit), "Double machine learning mediation")
  # binary mediator with a controlled direct effect
  db <- sim_mediation(3000, dgp = "binary", seed = 3)
  fb <- mediate_dml(db, "y", "d", "m", x = c("x1", "x2"), m_ref = 0, seed = 1)
  expect_true("cde" %in% fb$effects$term)
  tr <- truth_of(db, c("total", "nde", "nie"))
  expect_true(all(abs(fb$effects$estimate[1:3] - tr) < 3.5 * fb$effects$std.error[1:3]))
  expect_lt(abs(fb$effects$estimate[fb$effects$term == "cde"] - truth_of(db, "cde")), 0.08)
  # repeated cross-fitting
  fr <- mediate_dml(dat, "y", "d", "m", x = c("x1", "x2"), n_rep = 2, folds = 3, seed = 1)
  expect_equal(fr$n_rep, 2)
  skip_if_not_installed("ranger")
  dn <- sim_mediation(3000, dgp = "nonlinear", seed = 2)
  rf <- mlr3::lrn("regr.ranger", num.trees = 100, min.node.size = 10)
  rfc <- mlr3::lrn("classif.ranger", num.trees = 100, min.node.size = 10, predict_type = "prob")
  fn <- mediate_dml(dn, "y", "d", "m", x = c("x1", "x2"), learner_y = rf, learner_d = rfc, seed = 1)
  fl <- mediate_reg(dn, "y", "d", "m", x = c("x1", "x2"), interaction = TRUE, n_sim = 50, seed = 1)
  tr <- truth_of(dn, c("nde", "nie"))
  err_dml <- abs(fn$effects$estimate[2:3] - tr)
  err_reg <- abs(fl$effects$estimate[match(c("nde", "nie"), fl$effects$term)] - tr)
  expect_true(all(err_dml < err_reg))
})

test_that("gelbach_decomp is exact and order invariant", {
  dp <- sim_mediation(1500, dgp = "parallel", seed = 4)
  g <- gelbach_decomp(dp, "y", "d", x_base = c("x1", "x2"), x_add = c("m1", "m2"), n_boot = 0)
  expect_lt(abs(g$identity_gap), 1e-10)
  expect_equal(sum(g$contributions$contribution), g$coefficients$estimate[g$coefficients$term == "change"], tolerance = 1e-10)
  g2 <- gelbach_decomp(dp, "y", "d", x_base = c("x1", "x2"), x_add = c("m2", "m1"), n_boot = 0)
  expect_equal(g2$contributions$contribution[g2$contributions$covariate == "m1"],
               g$contributions$contribution[g$contributions$covariate == "m1"], tolerance = 1e-10)
  # contributions equal the sequential-ignorability indirect effects in this linear design
  tr <- truth_of(dp, c("nie_m1", "nie_m2"))
  expect_lt(abs(g$contributions$contribution[1] - tr[1]), 0.15)
  gb <- gelbach_decomp(dp, "y", "d", x_base = c("x1", "x2"), x_add = c("m1", "m2"),
                       groups = list(both = c("m1", "m2")), n_boot = 49, seed = 1)
  expect_equal(gb$groups$contribution, sum(gb$contributions$contribution), tolerance = 1e-10)
  expect_true(all(is.finite(gb$contributions$std.error)))
  expect_s3_class(plot_gelbach(gb), "ggplot")
  expect_s3_class(plot_gelbach(gb, groups = TRUE), "ggplot")
  expect_output(print(gb), "Gelbach")
  dp$cl <- rep(1:50, length.out = nrow(dp))
  gf <- gelbach_decomp(dp, "y", "d", x_add = c("m1", "m2"), fe = "cl", cluster = "cl", n_boot = 19, seed = 1)
  expect_lt(abs(gf$identity_gap), 1e-8)
})

test_that("mediate_cde recovers the controlled direct effect with a treatment-induced confounder", {
  dz <- sim_mediation(4000, dgp = "post_treatment_confounder", seed = 5)
  fit <- mediate_cde(dz, "y", "d", "m", x_pre = c("x1", "x2"), x_post = "z", n_boot = 49, seed = 1)
  cde <- fit$effects[fit$effects$term == "cde", ]
  expect_lt(abs(cde$estimate - truth_of(dz, "cde")), 3 * cde$std.error + 0.02)
  naive <- fit$effects$estimate[fit$effects$term == "naive_direct"]
  expect_gt(abs(naive - truth_of(dz, "cde")), abs(cde$estimate - truth_of(dz, "cde")))
  expect_output(print(fit), "sequential g-estimation")
  # without post-treatment confounders it equals the regression CDE
  dl <- sim_mediation(1500, dgp = "linear", seed = 1)
  f1 <- mediate_cde(dl, "y", "d", "m", x_pre = c("x1", "x2"), n_boot = 0 + 9)
  f2 <- mediate_reg(dl, "y", "d", "m", x = c("x1", "x2"), method = "delta")
  expect_equal(f1$effects$estimate[f1$effects$term == "cde"], f2$effects$estimate[f2$effects$term == "cde"], tolerance = 1e-8)
})

test_that("mediate_iv decomposes the effect through an instrumented mediator", {
  di <- sim_mediation(4000, dgp = "iv_mediator", seed = 6)
  fit <- mediate_iv(di, "y", "d", "m", z = "z", x = c("x1", "x2"), n_boot = 49, seed = 1)
  e <- fit$effects
  get <- function(t) e$estimate[e$term == t]
  se <- function(t) e$std.error[e$term == t]
  expect_lt(abs(get("direct") - truth_of(di, "direct")), 3 * se("direct") + 0.02)
  expect_lt(abs(get("indirect") - truth_of(di, "indirect")), 3 * se("indirect") + 0.02)
  expect_equal(get("total_z"), get("direct") + get("indirect"), tolerance = 1e-8)
  expect_gt(abs(get("naive_mediator_effect") - 0.7), abs(get("mediator_effect") - 0.7))
  iv <- fixest::feols(y ~ d + x1 + x2 | m ~ z, data = di)
  expect_equal(get("direct"), unname(coef(iv)["d"]), tolerance = 1e-10)
  expect_gt(fit$first_stage_F, 50)
  di$g <- as.integer(di$x2 > 0)
  fh <- mediate_iv(di, "y", "d", "m", z = "z", x = c("x1", "x2"), homogeneity_by = "g", n_boot = 9)
  expect_equal(nrow(fh$homogeneity), 2)
  expect_output(print(fh), "homogeneity")
  expect_s3_class(plot_mediation(fit, terms = c("total", "direct", "indirect")), "ggplot")
})

test_that("front_door recovers the effect under unobserved treatment-outcome confounding", {
  skip_if_not_installed("mlr3")
  df <- sim_mediation(4000, dgp = "front_door", seed = 7)
  fr <- front_door(df, "y", "d", "m", x = c("x1", "x2"), n_boot = 49, seed = 1)
  tot <- fr$effects[fr$effects$term == "total", ]
  expect_lt(abs(tot$estimate - truth_of(df, "total")), 3 * tot$std.error + 0.03)
  bd <- fr$effects$estimate[fr$effects$term == "backdoor"]
  expect_gt(abs(bd - truth_of(df, "total")), 0.3)
  fa <- front_door(df, "y", "d", "m", x = c("x1", "x2"), method = "aipw", seed = 1)
  tot_a <- fa$effects[fa$effects$term == "total", ]
  expect_lt(abs(tot_a$estimate - truth_of(df, "total")), 3 * tot_a$std.error + 0.03)
  expect_equal(nrow(fa$psi), nrow(df))
  expect_output(print(fa), "Front-door")
  df$mc <- df$m + rnorm(nrow(df), sd = 0.01)
  expect_error(front_door(df, "y", "d", "mc", method = "aipw"), "binary")
  fc <- front_door(df, "y", "d", "mc", x = c("x1", "x2"), n_boot = 9)
  expect_true(is.finite(fc$effects$estimate[1]))
})

test_that("sim_mediation designs carry consistent truths", {
  for (g in c("linear", "interaction", "nonlinear", "binary", "correlated_errors", "post_treatment_confounder",
              "iv_mediator", "front_door", "parallel")) {
    d <- sim_mediation(300, dgp = g, seed = 1)
    expect_true(is.list(attr(d, "truth")))
    expect_true(all(d$d %in% c(0, 1)))
  }
  dl <- sim_mediation(300, dgp = "linear", seed = 1)
  tr <- attr(dl, "truth")
  expect_equal(tr$total, tr$nde + tr$nie, tolerance = 1e-10)
})

test_that("weights are carried through mediate_reg, mediate_iv, and gelbach_decomp", {
  dat <- sim_mediation(2000, dgp = "linear", seed = 11)
  dat$w <- runif(nrow(dat), 0.5, 2)
  fw <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), weights = "w", method = "delta")
  a1 <- coef(lm(m ~ d + x1 + x2, data = dat, weights = w))[["d"]]
  b2 <- coef(lm(y ~ d + m + x1 + x2, data = dat, weights = w))[["m"]]
  expect_equal(fw$effects$estimate[fw$effects$term == "nie"], a1 * b2, tolerance = 1e-10)
  f1 <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), method = "delta")
  expect_false(isTRUE(all.equal(fw$effects$estimate[3], f1$effects$estimate[3])))
  # unit weights reproduce the unweighted fit, including the robust covariance
  dat$one <- 1
  f_one <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), weights = "one", method = "delta")
  expect_equal(f_one$effects, f1$effects, tolerance = 1e-10)
  expect_equal(f_one$shares$std.error, f1$shares$std.error, tolerance = 1e-10)
  # weighted robust covariance equals sandwich's HC1 for a weighted lm
  fit_w <- lm(y ~ d + m + x1 + x2, data = dat, weights = w)
  V <- causalmetrics:::.cm_vcov_robust(fit_w)
  X <- model.matrix(fit_w); e <- residuals(fit_w); ww <- dat$w
  bread <- solve(crossprod(X * sqrt(ww))); meat <- crossprod(X * (ww * e))
  expect_equal(V, bread %*% meat %*% bread * nrow(X) / (nrow(X) - ncol(X)), tolerance = 1e-10, ignore_attr = TRUE)
  # logit with weights
  db <- sim_mediation(1500, dgp = "binary", seed = 3); db$w <- runif(nrow(db), 0.5, 2)
  fb <- suppressWarnings(mediate_reg(db, "y", "d", "m", x = c("x1", "x2"), outcome = "logit", mediator = "logit", weights = "w", n_sim = 50, seed = 1))
  expect_true(all(is.finite(fb$effects$estimate)))
  # mediate_iv
  di <- sim_mediation(2000, dgp = "iv_mediator", seed = 6); di$w <- runif(nrow(di), 0.5, 2)
  fi <- mediate_iv(di, "y", "d", "m", z = "z", x = c("x1", "x2"), weights = "w", n_boot = 9)
  iv <- fixest::feols(y ~ d + x1 + x2 | m ~ z, data = di, weights = ~w)
  expect_equal(fi$effects$estimate[fi$effects$term == "direct"], unname(coef(iv)["d"]), tolerance = 1e-10)
  # gelbach
  dp <- sim_mediation(1500, dgp = "parallel", seed = 4); dp$w <- runif(nrow(dp), 0.5, 2)
  g <- gelbach_decomp(dp, "y", "d", x_base = c("x1", "x2"), x_add = c("m1", "m2"), weights = "w", n_boot = 0)
  expect_lt(abs(g$identity_gap), 1e-10)
  expect_equal(g$coefficients$estimate[1], coef(lm(y ~ d + x1 + x2, data = dp, weights = w))[["d"]], tolerance = 1e-10)
})
