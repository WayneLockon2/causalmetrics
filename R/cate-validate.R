# R/cate-validate.R
#
# Validation of a fixed CATE model on a held-out test sample: the best linear
# predictor heterogeneity test, calibration by CATE quantile groups, and the
# targeting operating characteristic (TOC) and QINI curves with their areas.

#' Validate a CATE model on a held-out sample
#'
#' Runs the diagnostics of Chernozhukov et al. (2026, Section 15.3) for a
#' fixed CATE model `tau_hat` using the doubly robust pseudo-outcomes of a
#' test sample:
#'
#' * **Heterogeneity test.** OLS of the pseudo-outcome on
#'   `(1, tau_hat - mean(tau_hat))`. The slope estimates
#'   `Cov(Y(1) - Y(0), tau_hat(X)) / Var(tau_hat(X))`; a significant slope
#'   means the model is correlated with the true effect, and the ideal value
#'   is 1. The intercept is the ATE.
#' * **Calibration.** Groups by quantiles of `tau_hat`; in each group the
#'   doubly robust group average effect is compared with the mean prediction.
#'   `cal1` and `cal2` are the group-weighted absolute and squared
#'   discrepancies; the Wald test checks whether the groups differ.
#' * **Targeting curves.** For each share `q`, `TOC(q)` is the average effect
#'   among the top-`q` predicted units minus the ATE and `QINI(q) = TOC(q) q`
#'   (both are covariances between the individual effect and the targeting
#'   indicator). Pointwise standard errors come from the influence functions,
#'   simultaneous one-sided bands from the multiplier bootstrap. The areas
#'   `AUTOC` and `AUQC` integrate the curves over all ranks (the rank-weighted
#'   average treatment effects of Yadlowsky et al. 2021, as in
#'   `grf::rank_average_treatment_effect()`) and are reported with one-sided
#'   intervals from their influence functions. Ties in
#'   `tau_hat` are broken at random with the rule of their Remark 15.3.1.
#' * **Group differences.** Covariate means in the top and bottom groups.
#'
#' @param scores A `cm_scores` object on the test sample.
#' @param tau_hat A `cm_cate` object (predicted on `scores$data`), a numeric
#'   vector, or a column name.
#' @param n_groups Number of quantile groups for calibration.
#' @param quantile_grid Shares `q` at which the curves are evaluated.
#' @param breaks Optional quantile thresholds of `tau_hat` computed on
#'   non-test data (a list with `groups` and `curve` elements, or `NULL` to
#'   use the test-sample quantiles).
#' @param covariates Covariates for the group-difference table (default the
#'   `x` of `scores`).
#' @param n_boot Multiplier bootstrap draws.
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_cate_val` with `blp` (the heterogeneity
#'   regression table), `calibration` (group table, `cal1`, `cal2`, and
#'   `test_equal`), `curves` (data frame with `q`, `toc`, `qini`, standard
#'   errors, pointwise and one-sided simultaneous lower bounds), `areas`
#'   (`autoc`, `auqc` with standard errors and lower bounds), `group_diff`,
#'   and `tau_hat`. [plot_cate_validation()] draws the pieces.
#' @examples
#' dat <- sim_hte(2000, dgp = "smooth", seed = 1)
#' x <- paste0("x", 1:5)
#' train <- dat[1:1200, ]; test <- dat[1201:2000, ]
#' m <- cate_learner(train, "y", "d", x, x_het = c("x1", "x2"), method = "dr", seed = 1)
#' sc_test <- dr_scores(test, "y", "d", x, seed = 2)
#' v <- cate_validate(sc_test, m, n_boot = 199, seed = 1)
#' v
#' @references
#' Yadlowsky, S., Fleming, S., Shah, N., Brunskill, E., and Wager, S. (2021).
#' Evaluating treatment prioritization rules via rank-weighted average
#' treatment effects. arXiv:2111.07966.
#' @seealso [cate_score()], [plot_cate_validation()]
#' @export
cate_validate <- function(scores, tau_hat, n_groups = 4L, quantile_grid = seq(0.05, 1, by = 0.05),
                          breaks = NULL, covariates = NULL, n_boot = 999L, conf_level = 0.95,
                          seed = NULL) {
  .cm_check_scores(scores)
  n_groups <- .cm_check_count(n_groups, "n_groups", min = 2L)
  y <- scores$score
  n <- scores$n
  tau <- .cm_cate_predictions(tau_hat, scores, "tau_hat")
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  crit1 <- stats::qnorm(conf_level)

  # 1. heterogeneity test -----------------------------------------------------
  tau_c <- tau - mean(tau)
  blp <- if (stats::var(tau) > 0) {
    fit <- .cm_ols_if(cbind(1, tau_c), y, conf_level = conf_level, uniform = FALSE)
    data.frame(term = c("ATE", "tau_hat (centered)"), estimate = fit$beta, std.error = fit$se,
               statistic = fit$beta / fit$se, p.value = 2 * stats::pnorm(-abs(fit$beta / fit$se)),
               conf.low = fit$beta - crit * fit$se, conf.high = fit$beta + crit * fit$se,
               stringsAsFactors = FALSE)
  } else {
    data.frame(term = c("ATE", "tau_hat (centered)"), estimate = c(mean(y), NA), std.error = c(stats::sd(y) / sqrt(n), NA),
               statistic = NA, p.value = NA, conf.low = NA, conf.high = NA, stringsAsFactors = FALSE)
  }

  # 2. calibration ------------------------------------------------------------
  gb <- breaks$groups
  if (is.null(gb)) gb <- stats::quantile(tau, probs = seq(0, 1, length.out = n_groups + 1L), names = FALSE)
  gb[1] <- -Inf
  gb[length(gb)] <- Inf
  grp <- cut(tau, breaks = unique(gb), include.lowest = TRUE, labels = FALSE)
  grp <- factor(grp, levels = seq_len(max(grp)))
  G <- nlevels(grp)
  cal <- do.call(rbind, lapply(seq_len(G), function(k) {
    i <- grp == levels(grp)[k]
    m <- sum(i)
    gate <- mean(y[i])
    se <- if (m > 1L) stats::sd(y[i]) / sqrt(m) else NA_real_
    data.frame(group = k, n = m, share = m / n, mean_tau = mean(tau[i]), gate = gate, std.error = se,
               conf.low = gate - crit * se, conf.high = gate + crit * se)
  }))
  cal1 <- sum(abs(cal$gate - cal$mean_tau) * cal$share)
  cal2 <- sum((cal$gate - cal$mean_tau)^2 * cal$share)
  if_gate <- vapply(seq_len(G), function(k) {
    i <- grp == levels(grp)[k]
    ifelse(i, (y - cal$gate[k]) / cal$share[k], 0)
  }, numeric(n))
  test_equal <- if (G > 1L) {
    C <- cbind(-1, diag(G - 1L))
    .cm_if_wald(as.numeric(C %*% cal$gate), if_gate %*% t(C), n)
  } else list(statistic = NA_real_, df = 0L, p.value = NA_real_)
  cal1_const <- sum(abs(cal$gate - mean(y)) * cal$share)

  # 3. TOC and QINI curves ----------------------------------------------------
  q <- sort(unique(quantile_grid))
  if (any(q <= 0 | q > 1)) stop("`quantile_grid` must lie in (0, 1].", call. = FALSE)
  thr <- breaks$curve
  if (is.null(thr)) thr <- .cm_quantile_threshold(tau, q)
  theta <- mean(y)
  ind <- vapply(seq_along(q), function(l) {
    above <- tau > thr[l]
    tie <- tau == thr[l]
    share_above <- mean(above)
    share_tie <- mean(tie)
    lambda <- if (share_tie > 0) .cm_clip((q[l] - share_above) / share_tie, 0, 1) else 0
    above + lambda * tie
  }, numeric(n))
  ind <- matrix(ind, nrow = n)
  pi_q <- colMeans(ind)
  toc <- as.numeric(colMeans((y - theta) * sweep(ind, 2, pi_q, "/")))
  qini <- as.numeric(colMeans((y - theta) * ind))
  if_toc <- sweep((y - theta) * sweep(ind, 2, pi_q, "/") - (y - theta), 2, toc)
  if_qini <- sweep((y - theta) * sweep(ind, 2, pi_q) , 2, qini)
  se_toc <- sqrt(colMeans(if_toc^2) / n)
  se_qini <- sqrt(colMeans(if_qini^2) / n)
  sup_toc <- .cm_sup_crit(if_toc, n, n_boot = n_boot, conf_level = conf_level, one_sided = TRUE, seed = seed)
  sup_qini <- .cm_sup_crit(if_qini, n, n_boot = n_boot, conf_level = conf_level, one_sided = TRUE, seed = seed)
  ct <- if (is.finite(sup_toc$crit)) sup_toc$crit else crit1
  cq <- if (is.finite(sup_qini$crit)) sup_qini$crit else crit1
  curves <- data.frame(q = q, share = pi_q, threshold = thr,
                       toc = toc, toc_se = se_toc, toc_low = toc - crit1 * se_toc, toc_band_low = toc - ct * se_toc,
                       qini = qini, qini_se = se_qini, qini_low = qini - crit1 * se_qini, qini_band_low = qini - cq * se_qini)
  # Areas integrated over all ranks (the rank-weighted average treatment
  # effects of Yadlowsky et al. 2021, as in grf::rank_average_treatment_effect):
  # AUTOC = int_0^1 TOC(u) du, AUQC = int_0^1 QINI(u) du. With u = k/n the
  # TOC weight of the unit ranked r is sum_{k >= r} 1/k - 1 and the QINI
  # weight is (n - r + 1)/n - (n + 1)/(2n); tied predictions receive the
  # average weight of their block, the expectation under random tie-breaking.
  rk <- rank(-tau, ties.method = "first")
  h_tail <- rev(cumsum(1 / rev(seq_len(n))))  # sum_{k >= r} 1/k
  w_toc <- h_tail[rk] - 1
  w_qini <- (n - rk + 1) / n - (n + 1) / (2 * n)
  tie_block <- match(tau, tau)
  w_toc <- stats::ave(w_toc, tie_block)
  w_qini <- stats::ave(w_qini, tie_block)
  autoc <- mean((y - theta) * w_toc)
  auqc <- mean((y - theta) * w_qini)
  if_autoc <- (y - theta) * w_toc - autoc
  if_auqc <- (y - theta) * w_qini - auqc
  areas <- data.frame(
    statistic = c("AUTOC", "AUQC"), estimate = c(autoc, auqc),
    std.error = c(sqrt(mean(if_autoc^2) / n), sqrt(mean(if_auqc^2) / n))
  )
  areas$conf.low <- areas$estimate - crit1 * areas$std.error
  areas$p.value <- stats::pnorm(-areas$estimate / areas$std.error)
  het_stat <- data.frame(curve = c("TOC", "QINI"),
                         max_lower_band = c(max(curves$toc_band_low), max(curves$qini_band_low)),
                         at_q = c(q[which.max(curves$toc_band_low)], q[which.max(curves$qini_band_low)]))

  # 4. group differences ------------------------------------------------------
  covs <- covariates %||% scores$x
  group_diff <- NULL
  if (!is.null(covs) && G > 1L) {
    top <- grp == levels(grp)[G]
    bottom <- grp == levels(grp)[1]
    group_diff <- do.call(rbind, lapply(covs, function(v) {
      z <- as.numeric(scores$data[[v]])
      m1 <- mean(z[top]); m0 <- mean(z[bottom])
      s1 <- stats::sd(z[top]) / sqrt(sum(top)); s0 <- stats::sd(z[bottom]) / sqrt(sum(bottom))
      data.frame(covariate = v, mean_top = m1, se_top = s1, mean_bottom = m0, se_bottom = s0,
                 diff = m1 - m0, se_diff = sqrt(s1^2 + s0^2), stringsAsFactors = FALSE)
    }))
  }

  structure(list(
    blp = blp,
    calibration = list(table = cal, cal1 = cal1, cal2 = cal2,
                       cal1_normalized = if (cal1_const > 0) 1 - cal1 / cal1_const else NA_real_,
                       test_equal = test_equal, groups = grp),
    curves = curves, areas = areas, heterogeneity_statistic = het_stat,
    group_diff = group_diff, tau_hat = tau, ate = theta, n = n,
    conf_level = conf_level, crit_band = c(toc = ct, qini = cq), call = match.call()
  ), class = "cm_cate_val")
}

#' @export
print.cm_cate_val <- function(x, ...) {
  cat("Validation of a CATE model on ", x$n, " held-out observations\n", sep = "")
  cat("\nHeterogeneity test (OLS of the score on the centered prediction):\n")
  print(x$blp[, c("term", "estimate", "std.error", "p.value", "conf.low", "conf.high")], digits = 4, row.names = FALSE)
  cat("\nCalibration by prediction quantile groups:\n")
  print(x$calibration$table[, c("group", "n", "mean_tau", "gate", "std.error", "conf.low", "conf.high")], digits = 4, row.names = FALSE)
  cat("  cal1 = ", format(round(x$calibration$cal1, 4)), ", normalized = ",
      format(round(x$calibration$cal1_normalized, 3)), "; equal-effects test p = ",
      format(signif(x$calibration$test_equal$p.value, 3)), "\n", sep = "")
  cat("\nAreas under the targeting curves (one-sided ", 100 * x$conf_level, "% lower bounds):\n", sep = "")
  print(x$areas, digits = 4, row.names = FALSE)
  cat("  largest lower simultaneous band: TOC ", format(round(x$heterogeneity_statistic$max_lower_band[1], 4)),
      " at q = ", x$heterogeneity_statistic$at_q[1], "; QINI ",
      format(round(x$heterogeneity_statistic$max_lower_band[2], 4)), " at q = ",
      x$heterogeneity_statistic$at_q[2], "\n", sep = "")
  invisible(x)
}

#' Plot the validation diagnostics of a CATE model
#'
#' @param x A `cm_cate_val` object.
#' @param what `"calibration"` (group average effect against mean
#'   prediction), `"toc"`, or `"qini"` (curve with pointwise and one-sided
#'   simultaneous lower bounds).
#' @return A ggplot object.
#' @export
plot_cate_validation <- function(x, what = c("calibration", "toc", "qini")) {
  what <- match.arg(what)
  if (what == "calibration") {
    tab <- x$calibration$table
    rng <- range(c(tab$mean_tau, tab$conf.low, tab$conf.high), na.rm = TRUE)
    return(
      ggplot2::ggplot(tab, ggplot2::aes(x = .data$mean_tau, y = .data$gate)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::coord_cartesian(xlim = rng, ylim = rng) +
        ggplot2::labs(x = "Mean predicted CATE in group", y = "Doubly robust group average effect") +
        ggplot2::theme_minimal(base_size = 11)
    )
  }
  cv <- x$curves
  if (what == "toc") {
    df <- data.frame(q = cv$q, estimate = cv$toc, low = cv$toc_low, band = cv$toc_band_low)
    ylab <- "TOC(q): effect among top q minus ATE"
  } else {
    df <- data.frame(q = cv$q, estimate = cv$qini, low = cv$qini_low, band = cv$qini_band_low)
    ylab <- "QINI(q): gain over treating a random share q"
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$q)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$band, ymax = .data$estimate), fill = "grey80", alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = .data$low), linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate), linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = .data$estimate), size = 1.5) +
    ggplot2::labs(x = "Share of the population treated (q)", y = ylab,
                  subtitle = "Dotted: pointwise lower bound; shaded: one-sided simultaneous band") +
    ggplot2::theme_minimal(base_size = 11)
}
