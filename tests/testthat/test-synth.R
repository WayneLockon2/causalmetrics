# Synthetic control: identities with sdid_weights and quadprog, exact recovery,
# placebo, conformal, specification test, leave-one-out, separate mode.

ferman_sc <- function(y_before, y_after, demean = FALSE) {
  # Ferman and Pinto's replication code (_aux.R), quadprog on raw outcomes
  y <- y_before[, 1]; X <- y_before[, -1]
  if (demean) X <- cbind(1, X)
  Dmat <- t(X) %*% X; dvec <- t(X) %*% y
  if (demean) {
    Amat <- t(rbind(c(0, rep(1, ncol(X) - 1)), cbind(0, diag(ncol(X) - 1)))); bvec <- c(1, rep(0, ncol(X) - 1))
  } else {
    Amat <- t(rbind(rep(1, ncol(X)), diag(ncol(X)))); bvec <- c(1, rep(0, ncol(X)))
  }
  m <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  if (demean) list(w = m$solution[-1], effects = -m$solution[1] + y_after %*% c(1, -m$solution[-1]))
  else list(w = m$solution, effects = y_after %*% c(1, -m$solution))
}

test_that("all-lag synthetic control matches sdid_weights and the quadprog solution", {
  skip_if_not_installed("quadprog")
  dat <- sim_synth_panel(n_donors = 12, t_pre = 20, t_post = 5, effect = 2, seed = 11)
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
  expect_equal(sum(fit$weights$weight), 1, tolerance = 1e-8)
  expect_true(all(fit$weights$weight >= -1e-10))
  Y <- fit$block$Y
  pre <- t(Y[, 1:fit$T0]); post <- t(Y[, -(1:fit$T0)])
  ref <- ferman_sc(cbind(pre[, ncol(pre)], pre[, -ncol(pre)]), cbind(post[, ncol(post)], post[, -ncol(post)]))
  expect_equal(fit$weights$weight, ref$w, tolerance = 1e-6)
  expect_equal(fit$estimate, mean(ref$effects), tolerance = 1e-6)
  sw <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d", estimator = "sc", sparsify = FALSE)
  expect_equal(fit$estimate, sw$estimate, tolerance = 1e-3)
  # demeaned SC equals the intercept QP and the difp average
  fit_dm <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
  ref_dm <- ferman_sc(cbind(pre[, ncol(pre)], pre[, -ncol(pre)]), cbind(post[, ncol(post)], post[, -ncol(post)]), demean = TRUE)
  expect_equal(fit_dm$weights$weight, ref_dm$w, tolerance = 1e-6)
  expect_equal(fit_dm$estimate, mean(ref_dm$effects), tolerance = 1e-6)
  difp <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d", estimator = "difp", sparsify = FALSE)
  expect_equal(fit_dm$estimate, difp$estimate, tolerance = 1e-3)
  expect_equal(fit_dm$intercept, mean(fit_dm$effects$treated[1:fit$T0]) - sum(fit_dm$weights$weight * rowMeans(Y[1:fit$N0, 1:fit$T0])))
})

test_that("a treated unit built as a convex combination of donors is recovered exactly", {
  dat <- sim_synth_panel(n_donors = 8, t_pre = 12, t_post = 4, sd_e = 0.5, seed = 12)
  w_true <- c(0.5, 0.3, 0.2, rep(0, 5))
  Y0 <- sapply(split(dat$y0, dat$id)[1:8], identity)
  y_tr <- as.numeric(Y0 %*% w_true) + c(rep(0, 12), rep(3, 4))
  dat$y[dat$id == 9] <- y_tr
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
  expect_equal(fit$weights$weight, w_true, tolerance = 1e-5)
  expect_equal(fit$estimate, 3, tolerance = 1e-5)
  expect_lt(fit$pre_rmspe, 1e-5)
  expect_true(all(fit$balance$in_range))
  expect_equal(unname(fit$concentration[["l2"]]), sum(w_true^2), tolerance = 1e-5)
})

test_that("covariates, chosen lags, V weights, constraints, ridge, and augmentation run", {
  dat <- sim_synth_panel(n_donors = 10, t_pre = 15, t_post = 4, effect = 2, seed = 13)
  dat$x1 <- 0.5 * dat$y0 + rnorm(nrow(dat))
  for (v in c("equal", "regression", "mspe")) {
    f <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", x = "x1", lags = c(5, 10, 15),
                       pre_window = 11:15, v = v, v_control = list(maxit = 100))
    expect_equal(sum(f$v$v), 1, tolerance = 1e-8)
    expect_equal(nrow(f$balance), 4L)
    expect_equal(f$balance$type, c("covariate", rep("lag", 3)))
    expect_true(f$settings$standardize)
  }
  f_user <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", x = "x1", lags = c(5, 10, 15), v = c(1, 1, 1, 7))
  expect_equal(f_user$v$v, c(0.1, 0.1, 0.1, 0.7))
  expect_error(synth_control(dat, id = "id", time = "time", y = "y", d = "d", lags = c(5, 10, 15), v = c(1, 1)), "one entry per predictor")
  f_none <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", constraints = "none")
  expect_false(isTRUE(all.equal(sum(f_none$weights$weight), 1)))
  skip_if_not_installed("quadprog")
  f_nn <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", constraints = "nonnegative")
  expect_true(all(f_nn$weights$weight >= 0))
  f_ridge <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", zeta = 1)
  f_plain <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
  expect_lt(sum(f_ridge$weights$weight^2), sum(f_plain$weights$weight^2) + 1e-8)
  f_aug <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", augment = "ridge")
  expect_true(is.numeric(f_aug$lambda) && nrow(f_aug$augmentation) == 19L)
  f_aug_big <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", augment = "ridge", lambda = 1e12)
  expect_equal(f_aug_big$estimate, f_plain$estimate, tolerance = 1e-6)
  # reuse weights on another outcome
  dat$y2 <- dat$y + 1
  f_re <- synth_control(dat, id = "id", time = "time", y = "y2", d = "d", weights = f_plain)
  expect_equal(f_re$weights$weight, f_plain$weights$weight)
  expect_equal(f_re$estimate, f_plain$estimate, tolerance = 1e-8)
})

test_that("placebo, conformal, specification, and leave-one-out tools behave", {
  dat <- sim_synth_panel(n_donors = 10, t_pre = 15, t_post = 4, effect = 40, seed = 14)
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
  pl <- synth_placebo(fit)
  expect_s3_class(pl, "cm_synth_placebo")
  expect_equal(dim(pl$gaps), c(19L, 11L))
  expect_equal(pl$mspe$rank[1L], 1L)
  expect_equal(pl$p_value, 1 / 11)
  expect_equal(nrow(pl$p_value_by_period), 4L)
  pl_cut <- synth_placebo(fit, mspe_limit = 2)
  expect_true(sum(pl_cut$kept) <= 11L)
  pt <- synth_placebo(fit, type = "time", placebo_time = 10)
  expect_s3_class(pt, "cm_synth")
  expect_equal(pt$T0, 9L)
  expect_equal(pt$T1, 6L)
  expect_error(synth_placebo(fit, type = "time", placebo_time = 2), "at least two periods")
  ci <- synth_conformal(fit, grid = seq(30, 50, by = 1))
  expect_lt(ci$p_value, 0.2)
  expect_true(ci$conf.low <= fit$estimate && ci$conf.high >= fit$estimate)
  expect_equal(nrow(ci$grid_p), 21L)
  ci_true <- synth_conformal(fit, null = 40)
  expect_gt(ci_true$p_value, 0.05)
  ci_pp <- synth_conformal(fit, per_period = TRUE, grid = seq(30, 50, by = 2))
  expect_equal(nrow(ci_pp$per_period), 4L)
  ci_iid <- synth_conformal(fit, permutations = "iid", n_perm = 99, seed = 1, grid = seq(30, 50, by = 5))
  expect_true(ci_iid$p_value >= 0 && ci_iid$p_value <= 1)
  st <- synth_spec_test(fit)
  expect_true(st$p_value >= 0 && st$p_value <= 1)
  expect_equal(nrow(st$paths), 19L)
  lo <- synth_loo(fit)
  expect_equal(nrow(lo$table), sum(fit$weights$weight > 0.001))
  expect_equal(ncol(lo$paths), nrow(lo$table))
  for (ty in c("path", "gap", "demeaned_path", "weights", "balance")) expect_s3_class(plot_synth(fit, ty), "ggplot")
  expect_s3_class(plot_synth(pl), "ggplot")
  expect_s3_class(plot_synth(pl, "ratio"), "ggplot")
  expect_s3_class(plot_synth(lo), "ggplot")
  expect_s3_class(plot_synth(st), "ggplot")
  expect_output(print(fit), "Synthetic control")
  expect_output(print(pl), "placebo")
  expect_output(print(ci), "Conformal")
  expect_output(print(st), "Ferman-Pinto")
  expect_output(print(lo), "Leave-one-donor-out")
  td <- tidy(fit)
  expect_equal(nrow(td), 5L)
  expect_equal(td$estimate[5L], fit$estimate)
  expect_equal(nrow(glance(fit)), 1L)
})

test_that("the specification test reproduces Ferman and Pinto's conformal_fp", {
  skip_if_not_installed("quadprog")
  conformal_fp <- function(pret, postt, q = 2) {
    full <- rbind(pret, postt)
    sc_mod <- ferman_sc(full, full, demean = TRUE)
    did_mod <- (full[, 1] - rowMeans(full[, -1])) - (mean(full[, 1]) - mean(full[, -1]))
    res <- sc_mod$effects - did_mod
    T0 <- nrow(pret); T1 <- nrow(postt); Tt <- T0 + T1
    distr <- t(sapply(0:(Tt - 1), function(j) {
      index <- 1:Tt + j; index <- ifelse(index > Tt, index - Tt, index)
      res_post <- res[index][(T0 + 1):Tt]
      c((sqrt(1 / T1) * sum(abs(res_post)^q))^(1 / q), sqrt(1 / T1) * abs(sum(res_post)))
    }))
    c(mean(distr[1, 1] <= distr[, 1]), mean(distr[1, 2] <= distr[, 2]))
  }
  dat <- sim_synth_panel(n_donors = 8, t_pre = 15, t_post = 5, select_loading = 1.5, factor_shift = 1.5, seed = 15)
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
  st <- synth_spec_test(fit, q = 2)
  Y <- fit$block$Y
  pre <- t(Y[, 1:fit$T0]); post <- t(Y[, -(1:fit$T0)])
  ref <- conformal_fp(cbind(pre[, ncol(pre)], pre[, -ncol(pre)]), cbind(post[, ncol(post)], post[, -ncol(post)]))
  expect_equal(st$p_value, ref[2L], tolerance = 1e-8)
  expect_equal(st$p_value_norm, ref[1L], tolerance = 1e-8)
})

test_that("the separate mode handles staggered adoption", {
  dat <- sim_synth_panel(n_donors = 12, n_treated = 3, t_pre = 15, t_post = 6, effect = 2,
                         adoption = c(13, 16, 19), seed = 16)
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", treated_units = "separate", horizon = 3)
  expect_equal(fit$mode, "separate")
  expect_equal(nrow(fit$unit_table), 3L)
  expect_equal(fit$unit_table$adoption, c(13, 16, 19))
  # later adopters serve as donors for earlier ones only while untreated over the window
  expect_equal(fit$unit_table$n_donors, c(14L, 13L, 12L))
  expect_true(all(fit$units[[1]]$T1 == 3L))
  expect_true(is.finite(fit$std.error))
  expect_true(all(0:2 %in% fit$event_gap$event_time))
  expect_s3_class(plot_synth(fit), "ggplot")
  expect_output(print(fit), "separately")
  expect_equal(nrow(tidy(fit)), 1L)
  expect_error(synth_control(dat, id = "id", time = "time", y = "y", d = "d"), "different periods")
  expect_error(synth_placebo(fit), "average mode")
})

test_that("sim_synth_panel returns a consistent panel", {
  dat <- sim_synth_panel(n_donors = 5, n_treated = 2, t_pre = 6, t_post = 3, effect = c(1, 2, 3), seed = 1)
  expect_equal(nrow(dat), 7 * 9)
  expect_equal(unique(dat$tau[dat$d == 1]), c(1, 2, 3))
  expect_true(all(dat$y[dat$d == 0] == dat$y0[dat$d == 0]))
  expect_equal(dim(attr(dat, "factors")), c(9L, 2L))
  expect_error(sim_synth_panel(n_treated = 2, adoption = 5), "one entry per treated")
})

test_that("Andersson (2019) synthetic Sweden weights are reproduced when the data are present", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("haven")
  path <- file.path(testthat::test_path("..", ".."), "data_raw", "09_Synthetic_Control", "andersson2019carbon", "carbontax_data.dta")
  skip_if_not(file.exists(path), "Andersson (2019) data not available")
  d <- as.data.frame(haven::read_dta(path))
  d$treated <- as.integer(d$country == "Sweden" & d$year >= 1990)
  set.seed(1)
  fit <- synth_control(d, id = "country", time = "year", y = "CO2_transport_capita", d = "treated",
                       x = c("GDP_per_capita", "vehicles_capita", "gas_cons_capita", "urban_pop"),
                       lags = c(1989, 1980, 1970), pre_window = 1980:1989, v = "mspe",
                       v_control = list(maxit = 2000, starts = 3))
  published <- c(Denmark = 0.384, Belgium = 0.195, `New Zealand` = 0.177, Greece = 0.090, `United States` = 0.088, Switzerland = 0.061)
  w <- setNames(fit$weights$weight, fit$weights$country)
  expect_true(all(abs(w[names(published)] - published) < 0.05))
  expect_true(all(w[setdiff(names(w), names(published))] < 0.01))
  expect_lt(fit$pre_rmspe, 0.05)
  expect_equal(fit$balance$treated[1L], 20121.479, tolerance = 1e-3)
})

test_that("the default conformal grid widens so the confidence set contains non-rejected nulls", {
  dat <- sim_synth_panel(n_donors = 10, t_pre = 12, t_post = 6, effect = 0.3, seed = 21)
  fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
  ci <- synth_conformal(fit)
  if (ci$p_value > 1 - ci$level) expect_true(ci$conf.low <= 0 && ci$conf.high >= 0)
  expect_true(ci$grid_p$p_value[1L] <= 1 - ci$level || nrow(ci$grid_p) > 61L)
})
