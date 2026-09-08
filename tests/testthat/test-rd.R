skip_if_not_installed("rdrobust")

test_that("sim_rd designs have the documented properties", {
  d <- sim_rd(500, "lee", seed = 1)
  expect_equal(attr(d, "tau_true"), 0.04)
  expect_equal(d$mu1 - d$mu0, rep(0.04, 500))
  expect_true(all(d$d == as.integer(d$x >= 0)))
  f <- sim_rd(500, "fuzzy", compliance = 0.8, seed = 1)
  expect_gt(mean(f$d[f$x >= 0]), mean(f$d[f$x < 0]))
  k <- sim_rd(500, "kink", seed = 1)
  expect_true(all(k$d[k$x < 0] == 1))
  expect_equal(attr(k, "tau_true"), 2)
  m <- sim_rd(500, "multi_cutoff", seed = 1)
  expect_equal(sort(unique(m$cutoff)), c(-0.5, 0, 0.5))
  dd <- sim_rd(500, "discrete", grid = 0.1, seed = 1)
  expect_true(all(abs(dd$x / 0.1 - round(dd$x / 0.1)) < 1e-8))
  mp <- sim_rd(2000, "manipulated", seed = 1)
  expect_gt(sum(mp$x > 0 & mp$x < 0.1), sum(mp$x < 0 & mp$x > -0.1))
})

test_that("rd_bins and rd_plot bin the data on each side and return ggplot objects", {
  dat <- sim_rd(2000, "lee", seed = 1)
  b <- rd_bins(dat, "y", "x", n_bins = 10)
  expect_equal(as.integer(table(b$side)), c(10L, 10L))
  expect_true(all(b$bin_x[b$side == "left"] < 0) && all(b$bin_x[b$side == "right"] >= 0))
  expect_equal(sum(b$n), 2000)
  # evenly spaced bins have equal widths
  be <- rd_bins(dat, "y", "x", bins = "es", n_bins = 8)
  widths <- be$x_max - be$x_min
  expect_lt(diff(range(widths[be$side == "left"])), 1e-10)
  p <- rd_plot(dat, "y", "x")
  expect_s3_class(p, "ggplot")
  expect_true(all(attr(p, "n_bins") > 0))
  # the default bin count comes from rdrobust::rdplot
  rp <- rdrobust::rdplot(dat$y, dat$x, c = 0, hide = TRUE, binselect = "qsmv", masspoints = "off")
  expect_equal(attr(p, "n_bins"), as.integer(rp$J))
  p2 <- rd_plot(dat, "y", "x", fit = "local", ci = TRUE, n_bins = 15)
  expect_s3_class(p2, "ggplot")
  expect_true(is.numeric(attr(p2, "h")))
  fd <- attr(p2, "rd_fit")
  # the local fits meet the cutoff from both sides
  expect_equal(max(fd$x_grid[fd$side == "left"]), 0)
  expect_equal(min(fd$x_grid[fd$side == "right"]), 0)
  p3 <- rd_plot(dat, "y", "x", fit = "none", n_bins = c(5, 8))
  expect_equal(attr(p3, "n_bins"), c(5L, 8L))
})

test_that("tidiers and rd_frame reproduce the printed numbers of the RD packages", {
  dat <- sim_rd(2000, "lee", seed = 2)
  fit <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
  fr <- rd_frame(a = fit)
  expect_equal(fr$estimate, as.numeric(fit$coef))
  expect_equal(fr$conf.low, as.numeric(fit$ci[, 1]))
  skip_if_not_installed("RDHonest")
  hn <- RDHonest::RDHonest(y ~ x, data = dat, cutoff = 0)
  th <- tidy(hn)
  expect_equal(th$estimate, hn$coefficients$estimate)
  fr2 <- rd_frame(rdrobust = fit, honest = hn, methods = "robust")
  expect_equal(nrow(fr2), 2)
  expect_equal(fr2$method, c("robust", "honest"))
})

test_that("rd_balance matches rdrobust covariate by covariate and its joint test behaves", {
  dat <- sim_rd(2000, "covariates", seed = 3)
  bal <- rd_balance(dat, paste0("z", 1:4), "x")
  for (v in paste0("z", 1:4)) {
    f <- rdrobust::rdrobust(dat[[v]], dat$x, c = 0)
    expect_equal(bal$table$estimate[bal$table$covariate == v], f$coef[3, 1], tolerance = 1e-8)
    expect_equal(bal$table$p.value[bal$table$covariate == v], f$pv[3, 1], tolerance = 1e-8)
  }
  expect_equal(bal$joint$df, 4L)
  expect_gt(bal$joint$p.value, 0.01)
  # a covariate that jumps at the cutoff is detected by the joint test
  dat$z5 <- rnorm(2000) + 2 * (dat$x >= 0)
  bal2 <- rd_balance(dat, c("z1", "z5"), "x")
  expect_lt(bal2$joint$p.value, 1e-4)
  expect_equal(bal2$joint$df, 2L)
  expect_lt(bal2$table$p.value[bal2$table$covariate == "z5"], 1e-4)
  expect_s3_class(plot_rd_checks(bal2), "ggplot")
})

test_that("the local jump with influence functions reproduces a weighted least squares jump", {
  dat <- sim_rd(1500, "lee", seed = 4)
  h <- 0.3
  lj <- .cm_rd_local_jump(dat$y, dat$x, 0, h, kernel = "uniform", p = 1)
  w <- as.numeric(abs(dat$x) <= h)
  m <- lm(y ~ I(x >= 0) * x, data = dat, weights = w)
  expect_equal(lj$jumps, unname(coef(m)[["I(x >= 0)TRUE"]]), tolerance = 1e-8)
  # the influence-function SE equals the HC0-type sandwich of that regression
  X <- model.matrix(m)[w > 0, ]
  e <- resid(m)[w > 0]
  V <- solve(crossprod(X)) %*% crossprod(X * e) %*% solve(crossprod(X))
  expect_equal(lj$std.error, sqrt(V[2, 2]), tolerance = 1e-6)
})

test_that("placebo cutoffs, sensitivity, donut, and rd_checks run and are consistent with rdrobust", {
  dat <- sim_rd(2500, "lee", seed = 5)
  pl <- rd_placebo_cutoffs(dat, "y", "x", cutoffs = c(-0.4, 0.3))
  expect_equal(nrow(pl), 3)
  expect_equal(pl$side, c("left", "true cutoff", "right"))
  f_true <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
  expect_equal(pl$estimate[pl$side == "true cutoff"], f_true$coef[3, 1], tolerance = 1e-8)
  f_left <- rdrobust::rdrobust(dat$y[dat$x < 0], dat$x[dat$x < 0], c = -0.4)
  expect_equal(pl$estimate[1], f_left$coef[3, 1], tolerance = 1e-8)
  se <- rd_sensitivity(dat, "y", "x", h_grid = c(0.15, 0.25))
  expect_equal(nrow(se), 2)
  expect_equal(se$h, c(0.15, 0.25))
  rho <- f_true$bws[1, 1] / f_true$bws[2, 1]
  f_h <- rdrobust::rdrobust(dat$y, dat$x, c = 0, h = 0.25, b = 0.25 / rho)
  expect_equal(se$estimate[2], f_h$coef[3, 1], tolerance = 1e-8)
  dn <- rd_donut(dat, "y", "x", radius = c(0, 0.02))
  expect_equal(dn$n_dropped[1], 0)
  expect_equal(dn$n_dropped[2], sum(abs(dat$x) < 0.02))
  ck <- rd_checks(dat, "y", "x", h_grid = c(0.2, 0.3), radius = c(0, 0.01), cutoffs = c(-0.3, 0.3))
  expect_s3_class(ck, "cm_rd_checks")
  expect_null(ck$balance)
  expect_output(print(ck), "Placebo cutoffs")
  for (w in c("placebo", "sensitivity", "donut")) expect_s3_class(plot_rd_checks(ck, w), "ggplot")
  expect_error(plot_rd_checks(ck, "balance"), "not in this object")
})

test_that("rd_adjust removes covariate noise and equals rdrobust when nothing is adjusted", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  dat <- sim_rd(3000, "covariates", seed = 6)
  adj <- rd_adjust(dat, "y", "x", covariates = paste0("z", 1:4), seed = 1)
  expect_s3_class(adj, "cm_rd_adjust")
  expect_equal(nrow(adj$data), 3000)
  expect_true(adj$r2 > 0.2)
  cmp <- adj$comparison[adj$comparison$method == "robust", ]
  expect_lt(cmp$std.error[cmp$outcome == "adjusted"], cmp$std.error[cmp$outcome == "raw"])
  expect_lt(abs(cmp$estimate[cmp$outcome == "adjusted"] - 0.04), 3 * cmp$std.error[cmp$outcome == "adjusted"])
  # the raw fit is rdrobust
  f <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
  expect_equal(cmp$estimate[cmp$outcome == "raw"], f$coef[3, 1], tolerance = 1e-8)
  # adjusting on an irrelevant covariate changes little
  dat$noise <- rnorm(3000)
  adj0 <- rd_adjust(dat, "y", "x", covariates = "noise", seed = 1)
  expect_lt(abs(adj0$r2), 0.02)
  # side-averaged variant and a fuzzy design
  adj2 <- rd_adjust(dat, "y", "x", covariates = paste0("z", 1:4), pooled = FALSE, h_train = 0.5, weighted = TRUE, seed = 1)
  expect_true(all(is.finite(adj2$data$y_adj)))
  fz <- sim_rd(2000, "fuzzy", seed = 1)
  fz$z1 <- fz$x + rnorm(2000)
  adjf <- rd_adjust(fz, "y", "x", covariates = "z1", d = "d", seed = 1)
  expect_true("d_adj" %in% names(adjf$data))
  expect_output(print(adjf), "Flexible covariate adjustment")
  expect_equal(nrow(rd_frame(adj = adjf, methods = "robust")), 2)
})

test_that("rd_extrapolate recovers effects away from the cutoff under conditional independence", {
  dat <- sim_rd(4000, "cia", seed = 7)
  ex <- rd_extrapolate(dat, "y", "x", covariates = "z", window = 0.5)
  expect_s3_class(ex, "cm_rd_extrapolate")
  expect_true(all(ex$test$p.value > 0.001))
  expect_lt(abs(ex$effects$estimate[1] - 0.5), 3 * ex$effects$std.error[1])
  expect_lt(abs(ex$effects$estimate[2] - 0.5), 3 * ex$effects$std.error[2])
  expect_true(all(abs(ex$curve$estimate - 0.5) < 4 * ex$curve$std.error + 0.05))
  expect_s3_class(plot_rd_extrapolate(ex, truth = 0.5), "ggplot")
  # the test rejects when the outcome depends on x directly
  bad <- sim_rd(4000, "lee", seed = 7)
  bad$z <- rnorm(4000)
  ex_bad <- rd_extrapolate(bad, "y", "x", covariates = "z", window = 0.5)
  expect_lt(min(ex_bad$test$p.value), 1e-4)
  skip_if_not_installed("mlr3")
  ex2 <- rd_extrapolate(dat, "y", "x", covariates = "z", window = 0.5, method = "aipw", seed = 1)
  expect_lt(abs(ex2$effects$estimate[1] - 0.5), 3 * ex2$effects$std.error[1])
  expect_output(print(ex2), "aipw")
})

test_that("rd_weak_iv matches the Wald interval when the first stage is strong and widens when weak", {
  strong <- sim_rd(3000, "fuzzy", compliance = 0.9, seed = 8)
  wk <- rd_weak_iv(strong, "y", "d", "x")
  expect_s3_class(wk, "cm_rd_weak_iv")
  expect_gt(wk$first_stage$F_effective, 50)
  tw <- wk$table[wk$table$method == "local 2SLS (Wald)", ]
  ta <- wk$table[wk$table$method == "Anderson-Rubin", ]
  expect_lt(abs((ta$conf.high - ta$conf.low) / (tw$conf.high - tw$conf.low) - 1), 0.25)
  # local 2SLS at the same bandwidth equals rdrobust's conventional fuzzy estimate
  rob <- wk$table[wk$table$method == "rdrobust conventional", ]
  expect_lt(abs(tw$estimate - rob$estimate), 0.02)
  weak <- sim_rd(3000, "fuzzy", compliance = 0.52, seed = 8)
  wk2 <- rd_weak_iv(weak, "y", "d", "x")
  expect_lt(wk2$first_stage$F_effective, 10)
  ta2 <- wk2$table[wk2$table$method == "Anderson-Rubin", ]
  tw2 <- wk2$table[wk2$table$method == "local 2SLS (Wald)", ]
  expect_true(is.na(ta2$conf.low) || (ta2$conf.high - ta2$conf.low) > (tw2$conf.high - tw2$conf.low) * 0.9)
  expect_output(print(wk2), "Anderson-Rubin set")
})

test_that("rd_kink reproduces rdrobust with deriv = 1 and the four graphs draw", {
  dat <- sim_rd(4000, "kink", seed = 9)
  kk <- rd_kink(dat, "y", "x", d = "d", elasticity = TRUE)
  rf <- rdrobust::rdrobust(dat$y, dat$x, c = 0, deriv = 1, p = 1)
  expect_equal(kk$reduced_form$estimate, as.numeric(rf$coef), tolerance = 1e-8)
  fz <- rdrobust::rdrobust(dat$y, dat$x, c = 0, fuzzy = dat$d, deriv = 1, p = 1)
  expect_equal(kk$kink$estimate, as.numeric(fz$coef), tolerance = 1e-8)
  expect_lt(abs(kk$kink$estimate[3] - 2), 4 * kk$kink$std.error[3])
  expect_true(is.finite(kk$levels[["y"]]))
  # sharp version scales the reduced form by a known slope change
  ks <- rd_kink(dat, "y", "x", slope_change = 1)
  expect_equal(ks$kink$estimate, ks$reduced_form$estimate)
  expect_output(print(kk), "Elasticity")
  skip_if_not_installed("patchwork")
  expect_s3_class(plot_rd_kink(kk), "patchwork")
})
