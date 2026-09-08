sim_structured <- function(n = 3000, seed = 1) {
  set.seed(seed)
  x <- matrix(runif(n * 3), n)
  theta <- cbind(x[, 1] - 0.5, 1 + x[, 2], -0.5 + x[, 3], 4)
  t_obs <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.5))
  u <- theta[, 1] + theta[, 2] * t_obs[, 1] + theta[, 3] * t_obs[, 2]
  y <- theta[, 4] * plogis(u) + rnorm(n, sd = 0.1)
  truth <- sapply(list(c(1, 0), c(0, 1), c(1, 1)), function(tt) {
    mean(theta[, 4] * plogis(theta[, 1] + theta[, 2] * tt[1] + theta[, 3] * tt[2]) - theta[, 4] * plogis(theta[, 1]))
  })
  list(x = x, theta = theta, t_obs = t_obs, y = y, truth = truth)
}

test_that("with the true nuisance the correction is centred and the estimator matches the plug-in contrast", {
  s <- sim_structured()
  fit <- est_dml_structured(s$theta, s$y, s$t_obs, targets = c("10", "01", "11"),
                            fold_id = rep(1:3, length.out = length(s$y)))
  pooled <- fit$estimates[fit$estimates$fold == "pooled", ]
  expect_equal(pooled$target, c("(1, 0)", "(0, 1)", "(1, 1)"))
  expect_equal(pooled$plugin, s$truth, tolerance = 1e-10)
  expect_true(all(abs(pooled$estimate - s$truth) < 4 * pooled$std.error + 1e-3))
  expect_true(all(abs(colMeans(fit$psi - fit$plugin)) < 0.02))
  expect_equal(dim(fit$psi), c(length(s$y), 3L))
  expect_equal(fit$best, "(1, 0)")
  # contrasts: the best minus itself is zero, others are differences of estimates
  ctr <- fit$contrasts[fit$contrasts$fold == "pooled", ]
  expect_equal(ctr$estimate[ctr$target == "(1, 0)"], 0)
  expect_equal(ctr$estimate[ctr$target == "(0, 1)"],
               pooled$estimate[1] - pooled$estimate[2], tolerance = 1e-10)
  # pooled se equals the independent-fold combination and the averaged bounds are wider
  folds <- fit$estimates[fit$estimates$fold != "pooled" & fit$estimates$target == "(1, 1)", ]
  expect_equal(pooled$std.error[3], sqrt(sum(folds$std.error^2)) / 3)
  expect_lt(pooled$ci_low_avg[3], pooled$conf.low[3])
  expect_output(print(fit), "Structured double machine learning")
  td <- tidy(fit)
  expect_equal(td$term, pooled$target)
})

test_that("the blockwise solve reproduces a hand-written influence function", {
  s <- sim_structured(n = 300)
  fit <- est_dml_structured(s$theta, s$y, s$t_obs, targets = "11")
  lk <- causalmetrics:::.cm_structured_link("gen_sigmoid", 2L, 4L)
  dist <- fit$t_dist
  hand <- vapply(1:300, function(r) {
    th <- s$theta[r, , drop = FALSE]
    L <- Reduce(`+`, lapply(seq_len(nrow(dist$t)), function(j) {
      g <- as.numeric(lk$grad(th, dist$t[j, ]))
      2 * dist$prob[j] * tcrossprod(g)
    }))
    ell <- 2 * as.numeric(lk$grad(th, s$t_obs[r, ])) * (lk$G(th, s$t_obs[r, ]) - s$y[r])
    H <- lk$G(th, c(1, 1)) - lk$G(th, c(0, 0))
    Hg <- as.numeric(lk$grad(th, c(1, 1)) - lk$grad(th, c(0, 0)))
    H - sum(Hg * solve(L, ell))
  }, numeric(1))
  expect_equal(as.numeric(fit$psi[, 1]), hand, tolerance = 1e-9)
  # ridge changes the solve
  fr <- est_dml_structured(s$theta, s$y, s$t_obs, targets = "11", ridge = 0.1)
  expect_false(isTRUE(all.equal(fr$psi, fit$psi)))
})

test_that("the linear link with a constant OLS theta reproduces the OLS contrast, and links agree", {
  s <- sim_structured(n = 1000)
  b <- coef(lm(s$y ~ s$t_obs))
  th_lin <- matrix(b, 1000, 3, byrow = TRUE)
  fl <- est_dml_structured(th_lin, s$y, s$t_obs, targets = "11", link = "linear")
  expect_equal(fl$estimates$estimate[fl$estimates$fold == "pooled"], unname(b[2] + b[3]), tolerance = 1e-8)
  yb <- as.numeric(s$y > 2)
  tl <- est_dml_structured(s$theta[, 1:3], yb, s$t_obs, targets = "11", link = "logit")
  tc <- est_dml_structured(s$theta[, 1:3], yb, s$t_obs, targets = "11", link = "custom",
                           G = function(th, t) plogis(th[, 1] + rowSums(th[, -1] * t)),
                           G_grad = function(th, t) {
                             p <- plogis(th[, 1] + rowSums(th[, -1] * t))
                             p * (1 - p) * cbind(1, t)
                           })
  expect_equal(tl$psi, tc$psi, tolerance = 1e-12)
  # explicit assignment distribution and input checks
  td <- data.frame(a = c(0, 1, 0, 1), b = c(0, 0, 1, 1), prob = rep(0.25, 4))
  fd <- est_dml_structured(s$theta, s$y, s$t_obs, targets = "11", t_dist = td)
  expect_s3_class(fd, "cm_dml_structured")
  expect_error(est_dml_structured(s$theta[, 1:3], s$y, s$t_obs, targets = "11"), "columns")
  expect_error(est_dml_structured(s$theta, s$y, s$t_obs, targets = "111"), "length")
  expect_warning(est_dml_structured(s$theta, s$y, s$t_obs, targets = "11", t_dist = data.frame(a = c(0, 1, 0, 1), b = c(0, 0, 1, 1), prob = rep(0.3, 4))), "sum to one")
  # a rank-deficient assignment distribution makes Lambda(x) singular
  expect_error(est_dml_structured(s$theta, s$y, s$t_obs, targets = "11", t_dist = data.frame(a = c(0, 1), b = c(0, 1), prob = c(0.5, 0.5))), "singular")
})
