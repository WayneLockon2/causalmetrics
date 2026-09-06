test_that("att_gt reproduces did::att_gt cells and analytic standard errors", {
  skip_if_not_installed("did")
  dat <- sim_did_panel(n_units = 400, n_periods = 7, groups = c(3, 5), x_trend = 0.2, seed = 11)
  for (ctrl in c("never", "notyet")) {
    for (m in c("dr", "reg", "ipw")) {
      ours <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", x = c("x1", "x2"),
                     method = m, control_group = ctrl, n_boot = 0)
      ref <- suppressWarnings(did::att_gt(
        yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~ x1 + x2,
        data = dat, control_group = if (ctrl == "never") "nevertreated" else "notyettreated",
        est_method = m, bstrap = FALSE, cband = FALSE
      ))
      expect_equal(ours$att_gt$group, ref$group)
      expect_equal(ours$att_gt$time, ref$t)
      expect_equal(ours$att_gt$att, ref$att, tolerance = 1e-6)
      expect_equal(ours$att_gt$std.error, ref$se, tolerance = 1e-5)
    }
  }
})

test_that("universal base period and anticipation match did", {
  skip_if_not_installed("did")
  dat <- sim_did_panel(n_units = 300, n_periods = 7, groups = c(4, 6), anticipation = 1, seed = 3)
  ours <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", x = "x1",
                 base_period = "universal", anticipation = 1, control_group = "notyet", n_boot = 0)
  ref <- suppressWarnings(did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~x1,
                                      data = dat, control_group = "notyettreated", base_period = "universal",
                                      anticipation = 1, bstrap = FALSE, cband = FALSE))
  expect_equal(ours$att_gt$att, ref$att, tolerance = 1e-6)
  expect_equal(ours$att_gt$std.error, ref$se, tolerance = 1e-5)
})

test_that("aggregations match did::aggte", {
  skip_if_not_installed("did")
  dat <- sim_did_panel(n_units = 400, n_periods = 8, groups = c(3, 5, 7), x_trend = 0.1, seed = 5)
  ours <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", x = "x1", control_group = "notyet", n_boot = 0)
  ref <- suppressWarnings(did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~x1,
                                      data = dat, control_group = "notyettreated", bstrap = FALSE, cband = FALSE))
  for (type in c("simple", "group", "calendar", "dynamic")) {
    a <- aggregate_att(ours, type = type, n_boot = 0)
    r <- did::aggte(ref, type = type, bstrap = FALSE, cband = FALSE)
    expect_equal(a$overall$estimate, r$overall.att, tolerance = 1e-6, info = type)
    expect_equal(a$overall$std.error, r$overall.se, tolerance = 1e-5, info = type)
    if (type != "simple") {
      expect_equal(a$by$estimate, r$att.egt, tolerance = 1e-6, info = type)
      expect_equal(a$by$std.error, r$se.egt, tolerance = 1e-5, info = type)
    }
  }
  a <- aggregate_att(ours, type = "dynamic", balance_e = 1, n_boot = 0)
  r <- did::aggte(ref, type = "dynamic", balance_e = 1, bstrap = FALSE, cband = FALSE)
  expect_equal(a$by$event_time, r$egt)
  expect_equal(a$by$estimate, r$att.egt, tolerance = 1e-6)
  expect_equal(a$by$std.error, r$se.egt, tolerance = 1e-5)
})

test_that("panel cells reproduce DRDID and the differenced est_dml ATT score", {
  skip_if_not_installed("DRDID")
  dat <- sim_did_panel(n_units = 500, n_periods = 2, groups = 2, x_trend = 0.3, seed = 7)
  wide <- reshape(dat[, c("id", "time", "g", "x1", "x2", "y")], idvar = c("id", "g", "x1", "x2"),
                  timevar = "time", direction = "wide")
  X <- cbind(1, wide$x1, wide$x2)
  d <- as.numeric(wide$g == 2)
  ref <- DRDID::drdid_panel(wide$y.2, wide$y.1, d, covariates = X, inffunc = TRUE)
  cell <- causalmetrics:::.cm_did_panel_cell(wide$y.2, wide$y.1, d, X, rep(1, nrow(wide)), "dr")
  expect_equal(cell$att, ref$ATT, tolerance = 1e-8)
  expect_equal(cell$inffunc, as.numeric(ref$att.inf.func), tolerance = 1e-6)
  fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", x = c("x1", "x2"), n_boot = 0)
  expect_equal(fit$att_gt$att, ref$ATT, tolerance = 1e-8)
  expect_equal(fit$att_gt$std.error, ref$se, tolerance = 1e-6)
  # Hajek-normalized DR-DiD equals the est_dml ATT score up to the control normalization
  wide$dy <- wide$y.2 - wide$y.1
  wide$d <- d
  ps <- fitted(glm(d ~ x1 + x2, data = wide, family = binomial))
  mu0 <- predict(lm(dy ~ x1 + x2, data = wide[wide$d == 0, ]), newdata = wide)
  dml <- est_dml(wide, y = "dy", d = "d", model = "irm", estimand = "ATT",
                 p_hat = ps, mu0_hat = mu0, mu1_hat = mu0)
  expect_lt(abs(dml$estimate - ref$ATT), 0.05)
})

test_that("repeated cross-sections match did and DRDID", {
  skip_if_not_installed("did")
  dat <- sim_did_panel(n_units = 500, n_periods = 4, groups = c(3), x_trend = 0.2, seed = 9)
  set.seed(1)
  rcs <- dat[sample(nrow(dat), 1200), ]
  rcs$row <- seq_len(nrow(rcs))
  ours <- att_gt(rcs, id = "row", time = "time", group = "g", y = "y", x = "x1",
                 sampling = "rcs", control_group = "never", n_boot = 0)
  ref <- suppressWarnings(did::att_gt(yname = "y", tname = "time", idname = "row", gname = "g", xformla = ~x1,
                                      data = rcs, panel = FALSE, control_group = "nevertreated",
                                      bstrap = FALSE, cband = FALSE))
  expect_equal(ours$att_gt$att, ref$att, tolerance = 1e-6)
  expect_equal(ours$att_gt$std.error, ref$se, tolerance = 1e-5)
})

test_that("bootstrap standard errors and bands are sensible and clustered inference runs", {
  dat <- sim_did_panel(n_units = 300, n_periods = 6, groups = c(3, 5), seed = 2)
  dat$state <- (dat$id %% 30) + 1
  fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", control_group = "notyet",
                n_boot = 199, seed = 1)
  fit0 <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", control_group = "notyet", n_boot = 0)
  expect_equal(fit$att_gt$att, fit0$att_gt$att)
  expect_lt(max(abs(fit$att_gt$std.error / fit0$att_gt$std.error - 1)), 0.35)
  expect_gt(fit$crit_val, qnorm(0.975))
  fitc <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", cluster = "state", n_boot = 99, seed = 1)
  expect_true(all(is.finite(fitc$att_gt$std.error)))
  expect_s3_class(tidy(fit), "data.frame")
  expect_output(print(fit), "Group-time average treatment effects")
  agg <- aggregate_att(fit, type = "dynamic", n_boot = 99, seed = 1)
  expect_output(print(agg), "Overall ATT")
  expect_equal(nrow(tidy(agg)), nrow(agg$by))
  truth <- attr(dat, "att_gt")
  post <- fit$att_gt[fit$att_gt$post == 1, ]
  expect_lt(max(abs(post$att - truth$att)), 4 * max(post$std.error))
})

test_that("learner-based cells recover the effect", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  dat <- sim_did_panel(n_units = 600, n_periods = 3, groups = 3, x_trend = 0.3, seed = 4)
  fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", x = c("x1", "x2"),
                learner_p = mlr3::lrn("classif.log_reg", predict_type = "prob"),
                learner_or = mlr3::lrn("regr.lm"), folds = 3, seed = 1, n_boot = 0)
  expect_equal(fit$method, "dr_learner")
  truth <- attr(dat, "att_gt")
  expect_lt(abs(fit$att_gt$att[fit$att_gt$post == 1] - truth$att), 4 * fit$att_gt$std.error[fit$att_gt$post == 1])
})
