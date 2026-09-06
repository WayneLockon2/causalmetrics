sim_block <- function(n0 = 30, n1 = 5, t0 = 8, t1 = 4, tau = 2, seed = 1) {
  set.seed(seed)
  n <- n0 + n1; Tn <- t0 + t1
  f <- rnorm(Tn)
  load <- rnorm(n)
  alpha <- rnorm(n)
  y <- outer(alpha, rep(1, Tn)) + outer(load, f) + matrix(rnorm(n * Tn, sd = 0.3), n)
  treated <- seq_len(n) > n0
  y[treated, (t0 + 1):Tn] <- y[treated, (t0 + 1):Tn] + tau
  data.frame(id = rep(seq_len(n), each = Tn), time = rep(seq_len(Tn), n),
             y = as.vector(t(y)), d = as.integer(rep(treated, each = Tn) & rep(seq_len(Tn) > t0, n)))
}

test_that("sdid weights lie on the simplex and the estimate equals the weighted fixest regression", {
  dat <- sim_block()
  w <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d")
  expect_equal(sum(w$omega$weight), 1, tolerance = 1e-8)
  expect_equal(sum(w$lambda$weight), 1, tolerance = 1e-8)
  expect_true(all(w$omega$weight >= -1e-12) && all(w$lambda$weight >= -1e-12))
  fit <- fixest::feols(y ~ d | id + time, data = merge(dat, w$cell_weights, by = c("id", "time")), weights = ~weight)
  expect_equal(unname(coef(fit)[["d"]]), w$estimate, tolerance = 1e-8)
  expect_lt(abs(w$estimate - 2), 0.5)
  # DID weights are uniform and the estimate is the plain DiD
  d0 <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d", estimator = "did")
  expect_equal(unique(round(d0$omega$weight, 10)), round(1 / d0$N0, 10))
  plain <- fixest::feols(y ~ d | id + time, data = dat)
  expect_equal(d0$estimate, unname(coef(plain)[["d"]]), tolerance = 1e-8)
  sc <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d", estimator = "sc")
  expect_equal(sum(sc$lambda$weight), 0)
  expect_output(print(w), "Synthetic difference-in-differences")
  expect_s3_class(plot_sdid(w), "ggplot")
  expect_s3_class(plot_sdid(w, type = "weights"), "ggplot")
})

test_that("sdid matches synthdid when available and the standard errors run", {
  dat <- sim_block(n0 = 25, n1 = 4)
  w <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d")
  se_pl <- sdid_se(w, method = "placebo", n_reps = 50, seed = 1)
  se_bs <- sdid_se(w, method = "bootstrap", n_reps = 50, seed = 1)
  se_jk <- sdid_se(w, method = "jackknife")
  expect_true(all(c(se_pl$std.error, se_bs$std.error, se_jk$std.error) > 0))
})

test_that("att_dose recovers level effects and its TWFE identities hold", {
  set.seed(3)
  n <- 2000
  dose <- ifelse(runif(n) < 0.3, 0, runif(n, 0.5, 3))
  dat <- data.frame(id = rep(1:n, each = 2), time = rep(1:2, n), dose = rep(dose, each = 2))
  alpha <- rnorm(n)
  dat$y <- rep(alpha, each = 2) + 0.5 * dat$time + (dat$time == 2) * (2 * dat$dose - 0.3 * dat$dose^2) + rnorm(2 * n, sd = 0.5)
  fit <- att_dose(dat, id = "id", time = "time", y = "y", dose = "dose", n_bins = 4)
  truth <- 2 * fit$by_dose$dose - 0.3 * fit$by_dose$dose^2
  expect_lt(max(abs(fit$by_dose$att - truth) / fit$by_dose$std.error), 4)
  expect_equal(fit$checks[["sum_w_level_x_att"]], fit$checks[["twfe_binned"]], tolerance = 1e-8)
  expect_equal(fit$checks[["sum_w_slope_x_increment"]], fit$checks[["twfe_binned"]], tolerance = 1e-8)
  expect_equal(fit$checks[["sum_w_slope_x_dose_gap"]], 1, tolerance = 1e-8)
  expect_true(all(fit$twfe_weights$w_slope[-1] >= 0))
  expect_true(any(fit$twfe_weights$w_level < 0))
  fs <- att_dose(dat, id = "id", time = "time", y = "y", dose = "dose", dose_type = "spline", n_boot = 49, seed = 1)
  expect_equal(nrow(fs$curve), 50)
  expect_lt(max(abs(fs$curve$acrt - (2 - 0.6 * fs$curve$dose))[5:45]), 0.5)
  expect_s3_class(plot_att_dose(fs), "ggplot")
  expect_s3_class(plot_att_dose(fs, "acr"), "ggplot")
  expect_output(print(fit), "continuous treatment")
})

test_that("Frank-Wolfe unit weights agree with a quadratic-programming solution", {
  skip_if_not_installed("quadprog")
  dat <- sim_block(n0 = 20, n1 = 3, t0 = 8, t1 = 3, seed = 2)
  w <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d", sparsify = FALSE)
  Y <- w$Y; N0 <- w$N0; T0 <- w$T0
  A <- t(Y[seq_len(N0), seq_len(T0)])
  b <- colMeans(Y[(N0 + 1):nrow(Y), seq_len(T0), drop = FALSE])
  Ac <- scale(A, center = TRUE, scale = FALSE); bc <- b - mean(b)
  zeta <- w$zeta[["omega"]]
  D <- 2 * (crossprod(Ac) + T0 * zeta^2 * diag(N0))
  qp <- quadprog::solve.QP(D + diag(1e-10, N0), 2 * crossprod(Ac, bc), cbind(1, diag(N0)), c(1, rep(0, N0)), meq = 1)
  expect_lt(max(abs(qp$solution - w$omega$weight)), 0.02)
  obj <- function(om) sum((bc - Ac %*% om)^2) + T0 * zeta^2 * sum(om^2)
  expect_lt(abs(obj(w$omega$weight) - obj(qp$solution)) / obj(qp$solution), 1e-3)
})
