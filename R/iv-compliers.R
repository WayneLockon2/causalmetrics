# R/iv-compliers.R
#
# Who the instrument moves: Abadie's kappa weights and complier
# characteristics, the LATE score components as pseudo-outcomes for
# heterogeneity tools, instrument-specific LATEs and 2SLS weights, bounds on
# the ATE, and leave-out leniency instruments.

#' Complier share and characteristics (Abadie's kappa)
#'
#' Under the LATE assumptions with a binary instrument `Z` and binary
#' treatment `D`, Abadie (2003) shows that the weights
#' `kappa = 1 - D (1 - Z) / (1 - pi(X)) - (1 - D) Z / pi(X)`, with
#' `pi(X) = P(Z = 1 | X)`, average to the complier share and that
#' `E[kappa g(X)] / E[kappa]` is the complier mean of any function of the
#' covariates. Always-takers are the units with `D = 1, Z = 0` (weighted by
#' `1 / (1 - pi)`) and never-takers those with `D = 0, Z = 1` (weighted by
#' `1 / pi`). The function reports the three group shares and the covariate
#' means of each group next to the full-sample means, with bootstrap standard
#' errors that re-estimate the instrument propensity.
#'
#' @param data A data frame.
#' @param d,z Binary treatment and binary instrument column names.
#' @param covariates Character vector of covariates to profile.
#' @param x Optional covariates for the instrument propensity (logistic
#'   regression); `NULL` uses the unconditional share of `Z = 1`.
#' @param p_hat Optional supplied `P(Z = 1 | X)` (vector or column name).
#' @param n_boot Bootstrap replications for standard errors (0 to skip).
#' @param seed Optional seed.
#'
#' @return A list of class `cm_compliers` with `shares` (complier,
#'   always-taker, never-taker shares with standard errors), `table`
#'   (covariate means by group, complier/population ratios, standard errors
#'   of the complier means), `kappa` (the weights), and the propensity.
#' @references
#' Abadie, A. (2003). Semiparametric instrumental variable estimation of
#' treatment response models. *Journal of Econometrics*, 113(2), 231-263.
#' @examples
#' dat <- sim_iv(2000, dgp = "late", seed = 1)
#' complier_profile(dat, d = "d", z = "z", covariates = c("x1", "x2"), x = "x1", n_boot = 49)
#' @export
complier_profile <- function(data, d, z, covariates, x = NULL, p_hat = NULL, n_boot = 199L, seed = NULL) {
  data <- as.data.frame(data)
  for (v in c(d, z, covariates, x)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(d, z, covariates, x), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  dv <- .cm_as_binary(data[[d]], d)
  zv <- .cm_as_binary(data[[z]], z)
  p_in <- .cm_get_optional_numeric(p_hat, data, nrow(data), "p_hat")
  est <- function(idx) {
    dd <- dv[idx]; zz <- zv[idx]
    if (!is.null(p_in)) {
      pi <- p_in[idx]
    } else if (is.null(x)) {
      pi <- rep(mean(zz), length(idx))
    } else {
      df <- data[idx, x, drop = FALSE]; df$.cm_z <- zz
      pi <- stats::fitted(stats::glm(stats::reformulate(x, ".cm_z"), data = df, family = stats::binomial()))
    }
    pi <- pmin(pmax(pi, 0.01), 0.99)
    kappa <- 1 - dd * (1 - zz) / (1 - pi) - (1 - dd) * zz / pi
    w_at <- dd * (1 - zz) / (1 - pi)
    w_nt <- (1 - dd) * zz / pi
    shares <- c(complier = mean(kappa), always_taker = mean(w_at), never_taker = mean(w_nt))
    X <- as.matrix(data[idx, covariates, drop = FALSE])
    means <- cbind(all = colMeans(X),
                   compliers = colSums(kappa * X) / sum(kappa),
                   always_takers = colSums(w_at * X) / sum(w_at),
                   never_takers = colSums(w_nt * X) / sum(w_nt))
    list(shares = shares, means = means, kappa = kappa, pi = pi)
  }
  n <- nrow(data)
  point <- est(seq_len(n))
  se_shares <- rep(NA_real_, 3); se_means <- matrix(NA_real_, length(covariates), 4)
  if (n_boot > 0) {
    boots <- .cm_with_seed(seed, lapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      tryCatch(est(idx), error = function(e) NULL)
    }))
    boots <- boots[!vapply(boots, is.null, logical(1))]
    if (length(boots) > 1L) {
      se_shares <- apply(sapply(boots, function(b) b$shares), 1, stats::sd)
      arr <- simplify2array(lapply(boots, function(b) b$means))
      se_means <- apply(arr, c(1, 2), stats::sd)
    }
  }
  shares <- data.frame(group = names(point$shares), share = unname(point$shares), std.error = unname(se_shares))
  tab <- data.frame(covariate = covariates, point$means, check.names = FALSE)
  tab$ratio_compliers <- tab$compliers / tab$all
  tab$se_compliers <- se_means[, 2]
  tab$se_always <- se_means[, 3]
  tab$se_never <- se_means[, 4]
  rownames(tab) <- NULL
  structure(list(shares = shares, table = tab, kappa = point$kappa, propensity = point$pi,
                 n = n, n_boot = n_boot), class = "cm_compliers")
}

#' @export
print.cm_compliers <- function(x, ...) {
  cat("Compliance groups (Abadie's kappa), n = ", x$n, "\n", sep = "")
  print(x$shares, digits = 3, row.names = FALSE)
  cat("Covariate means by group:\n")
  print(x$table[, c("covariate", "all", "compliers", "always_takers", "never_takers", "ratio_compliers", "se_compliers")],
        digits = 3, row.names = FALSE)
  invisible(x)
}

#' LATE score components as pseudo-outcomes
#'
#' Extracts from an interactive-IV `est_dml()` fit the two cross-fitted
#' scores whose ratio is the LATE: `psi_b`, the doubly robust effect of the
#' instrument on the outcome, and `psi_a`, the doubly robust effect of the
#' instrument on the treatment (complier indicator in expectation). Both are
#' returned as `cm_scores` objects so that the heterogeneity tools of the
#' package ([cate_blp()], [cate_gate()], [cate_learner()], [policy_learn()])
#' apply to compliers: a rule with value `mean(pi(X) (psi_b - c psi_a))`
#' measures gains among compliers net of cost `c`.
#'
#' @param fit A `cm_dml` object with `model = "iivm"`.
#' @return A list with `outcome` and `treatment` (both `cm_scores`), the
#'   LATE `estimate`, and the covariate names.
#' @seealso [late_blp()]
#' @export
late_scores <- function(fit) {
  if (!inherits(fit, "cm_dml") || !identical(fit$model, "iivm")) stop("`fit` must be an est_dml(model = \"iivm\") object.", call. = FALSE)
  sc <- fit$score_components
  data <- fit$data
  n <- nrow(data)
  make <- function(score, type) {
    structure(list(score = as.numeric(score), nuisance = fit$nuisance, residuals = NULL,
                   y = data[[fit$outcome]], d = as.integer(.cm_as_binary(data[[fit$treatment]], "d")),
                   x = fit$covariates, y_name = fit$outcome, d_name = fit$treatment,
                   data = data, fold_id = fit$fold_id, n = n, type = type,
                   ate = list(estimate = mean(score), std.error = stats::sd(score) / sqrt(n)),
                   learners = fit$learners, outcome_type = "continuous", p_clip = NULL,
                   diagnostics = list(score = .cm_summary(score))), class = "cm_scores")
  }
  list(outcome = make(sc$psi_b, "late_outcome"), treatment = make(sc$psi_a, "late_treatment"),
       estimate = fit$estimate, x = fit$covariates)
}

#' Best linear predictor of the complier CATE
#'
#' Projects the two LATE score components of [late_scores()] on a dictionary
#' `p(X)` and forms the ratio `p(x)' b_outcome / p(x)' b_treatment` on
#' `newdata`, the conditional LATE `E[Y(1) - Y(0) | X = x, complier]` when the
#' complier effect is linear in the dictionary (Chernozhukov et al. 2026,
#' eq. 12.3.1). Standard errors come from the joint influence functions of
#' the two regressions by the delta method.
#'
#' @param scores The list returned by [late_scores()].
#' @param formula One-sided formula in the covariates.
#' @param newdata Data frame on which to evaluate the ratio.
#' @param conf_level Confidence level.
#' @return A data frame: `newdata` with `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, plus the numerator and denominator
#'   predictions; the two `cm_blp` fits are attached as attributes.
#' @export
late_blp <- function(scores, formula, newdata, conf_level = 0.95) {
  num <- cate_blp(scores$outcome, formula, uniform = FALSE, conf_level = conf_level)
  den <- cate_blp(scores$treatment, formula, uniform = FALSE, conf_level = conf_level)
  newdata <- as.data.frame(newdata)
  Xn <- stats::model.matrix(stats::delete.response(num$terms), newdata)
  a <- as.numeric(Xn %*% num$beta)
  b <- as.numeric(Xn %*% den$beta)
  tau <- a / b
  n <- num$n
  if_a <- num$inffunc %*% t(Xn)
  if_b <- den$inffunc %*% t(Xn)
  if_tau <- sweep(if_a - sweep(if_b, 2, tau, "*"), 2, b, "/")
  se <- sqrt(colMeans(if_tau^2) / n)
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  out <- cbind(newdata, data.frame(estimate = tau, std.error = se, conf.low = tau - zq * se, conf.high = tau + zq * se,
                                   outcome_effect = a, compliance = b))
  attr(out, "blp_outcome") <- num
  attr(out, "blp_treatment") <- den
  out
}

#' Instrument-specific LATEs and the weights two-stage least squares puts on them
#'
#' With several instruments, the two-stage least squares coefficient equals a
#' weighted average of the just-identified Wald estimates that use one
#' instrument at a time: `beta_2SLS = sum_k omega_k beta_k` with
#' `beta_k = Cov(Z_k, Y) / Cov(Z_k, D)` and
#' `omega_k = pi_k Cov(Z_k, D) / sum_j pi_j Cov(Z_j, D)`, where `pi` is the
#' first-stage coefficient vector (all instruments residualized on the
#' controls). The weights sum to one and are negative when an instrument's
#' first-stage coefficient and its covariance with the treatment disagree in
#' sign, which is what makes the 2SLS estimand hard to interpret with
#' heterogeneous effects (Mogstad, Torgovitsky, and Walters 2021).
#'
#' @param data A data frame.
#' @param y,d Outcome and treatment column names.
#' @param z Character vector of at least two instrument column names.
#' @param x Optional exogenous control column names.
#' @param weights Optional weights column name.
#' @return A list of class `cm_late_weights` with `table` (instrument,
#'   first-stage coefficient and t, reduced-form coefficient, `beta_k`,
#'   `omega_k`), `estimate` (2SLS), and `check` (`sum(omega_k beta_k)`).
#' @references
#' Mogstad, M., Torgovitsky, A., and Walters, C. R. (2021). The causal
#' interpretation of two-stage least squares with multiple instrumental
#' variables. *American Economic Review*, 111(11), 3663-3698.
#' @export
iv_late_weights <- function(data, y, d, z, x = NULL, weights = NULL) {
  data <- as.data.frame(data)
  if (length(z) < 2L) stop("`z` must name at least two instruments.", call. = FALSE)
  for (v in c(y, d, z, x, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, z, x, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  n <- nrow(data)
  w <- if (is.null(weights)) rep(1, n) else data[[weights]]
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  ts <- .cm_tsls(as.numeric(data[[y]]), as.numeric(data[[d]]), as.matrix(data[, z, drop = FALSE]), X, w)
  yt <- ts$residuals$y; dt <- ts$residuals$d; Zt <- ts$residuals$z
  pi_hat <- ts$first_stage
  cov_zd <- colSums(w * Zt * dt) / n
  cov_zy <- colSums(w * Zt * yt) / n
  beta_k <- cov_zy / cov_zd
  omega <- pi_hat * cov_zd / sum(pi_hat * cov_zd)
  fs <- .cm_first_stage_strength(dt, Zt, w)
  t_k <- pi_hat / sqrt(diag(fs$vcov))
  rf <- vapply(seq_along(z), function(k) sum(w * Zt[, k] * yt) / sum(w * Zt[, k]^2), numeric(1))
  tab <- data.frame(instrument = z, first_stage = unname(pi_hat), first_stage_t = unname(t_k),
                    reduced_form = rf, cov_zd = unname(cov_zd), beta_k = unname(beta_k), omega_k = unname(omega))
  structure(list(table = tab, estimate = ts$estimate, std.error = ts$std.error,
                 check = sum(omega * beta_k), n = n, first_stage_F = fs$F_robust), class = "cm_late_weights")
}

#' @export
print.cm_late_weights <- function(x, ...) {
  cat("Two-stage least squares with ", nrow(x$table), " instruments: estimate = ", format(round(x$estimate, 4)),
      " (SE ", format(round(x$std.error, 4)), "), first-stage F = ", format(round(x$first_stage_F, 2)), "\n", sep = "")
  cat("Instrument-specific IV estimates and 2SLS weights (sum(omega_k beta_k) = ", format(round(x$check, 4)), "):\n", sep = "")
  print(x$table, digits = 4, row.names = FALSE)
  if (any(x$table$omega_k < 0)) cat("  Negative weights: ", paste(x$table$instrument[x$table$omega_k < 0], collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Bounds on the ATE with a binary instrument
#'
#' Manski's (1990) instrumental-variable bounds for a binary outcome, binary
#' treatment, and binary instrument: for each instrument value the
#' no-assumption bounds on `E[Y(1)]` and `E[Y(0)]` are computed and, under
#' instrument independence, intersected across values of `Z`. The ATE bounds
#' are the difference of the intersected bounds. These are valid but not the
#' sharp bounds of Balke and Pearl (1997), which additionally use the joint
#' restrictions implied by exclusion; the LATE assumptions point-identify only
#' the complier effect, so these bounds show what the data say about the
#' whole population.
#'
#' @param data A data frame.
#' @param y,d,z Binary outcome, treatment, and instrument column names.
#' @param n_boot Bootstrap replications for the bounds' sampling variability.
#' @param seed Optional seed.
#' @return A list of class `cm_iv_bounds` with `bounds` (lower and upper on
#'   `E[Y(1)]`, `E[Y(0)]`, and the ATE), `no_instrument` (the bounds without
#'   using `Z`), the `late` Wald estimate, and bootstrap standard errors.
#' @references
#' Manski, C. F. (1990). Nonparametric bounds on treatment effects.
#' *American Economic Review*, 80(2), 319-323.
#'
#' Balke, A. and Pearl, J. (1997). Bounds on treatment effects from studies
#' with imperfect compliance. *Journal of the American Statistical
#' Association*, 92(439), 1171-1176.
#' @export
iv_ate_bounds <- function(data, y, d, z, n_boot = 199L, seed = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, z)) .cm_check_column(v, data)
  yv <- .cm_as_binary(data[[y]], y); dv <- .cm_as_binary(data[[d]], d); zv <- .cm_as_binary(data[[z]], z)
  keep <- !is.na(yv) & !is.na(dv) & !is.na(zv)
  yv <- yv[keep]; dv <- dv[keep]; zv <- zv[keep]
  calc <- function(idx) {
    yy <- yv[idx]; dd <- dv[idx]; zz <- zv[idx]
    b <- function(sel) {
      p11 <- mean(yy[sel] == 1 & dd[sel] == 1); p10 <- mean(yy[sel] == 1 & dd[sel] == 0)
      pd1 <- mean(dd[sel] == 1); pd0 <- 1 - pd1
      c(l1 = p11, u1 = p11 + pd0, l0 = p10, u0 = p10 + pd1)
    }
    all <- b(rep(TRUE, length(idx)))
    b0 <- b(zz == 0); b1 <- b(zz == 1)
    l1 <- max(b0["l1"], b1["l1"]); u1 <- min(b0["u1"], b1["u1"])
    l0 <- max(b0["l0"], b1["l0"]); u0 <- min(b0["u0"], b1["u0"])
    wald <- (mean(yy[zz == 1]) - mean(yy[zz == 0])) / (mean(dd[zz == 1]) - mean(dd[zz == 0]))
    c(ey1_lower = l1, ey1_upper = u1, ey0_lower = l0, ey0_upper = u0,
      ate_lower = l1 - u0, ate_upper = u1 - l0,
      noiv_lower = unname(all["l1"] - all["u0"]), noiv_upper = unname(all["u1"] - all["l0"]), late = wald)
  }
  n <- length(yv)
  pt <- calc(seq_len(n))
  se <- rep(NA_real_, length(pt))
  if (n_boot > 0) {
    B <- .cm_with_seed(seed, replicate(n_boot, calc(sample.int(n, n, replace = TRUE))))
    se <- apply(B, 1, stats::sd)
  }
  names(se) <- names(pt)
  bounds <- data.frame(quantity = c("E[Y(1)]", "E[Y(0)]", "ATE"),
                       lower = pt[c("ey1_lower", "ey0_lower", "ate_lower")],
                       upper = pt[c("ey1_upper", "ey0_upper", "ate_upper")],
                       se_lower = se[c("ey1_lower", "ey0_lower", "ate_lower")],
                       se_upper = se[c("ey1_upper", "ey0_upper", "ate_upper")])
  rownames(bounds) <- NULL
  structure(list(bounds = bounds, no_instrument = c(lower = unname(pt["noiv_lower"]), upper = unname(pt["noiv_upper"])),
                 late = unname(pt["late"]), late_se = unname(se["late"]), n = n), class = "cm_iv_bounds")
}

#' @export
print.cm_iv_bounds <- function(x, ...) {
  cat("Manski instrumental-variable bounds (binary Y, D, Z), n = ", x$n, "\n", sep = "")
  print(x$bounds, digits = 3, row.names = FALSE)
  cat("  ATE bounds without the instrument: [", format(round(x$no_instrument[1], 3)), ", ",
      format(round(x$no_instrument[2], 3)), "]; Wald LATE = ", format(round(x$late, 3)), "\n", sep = "")
  invisible(x)
}

#' Leave-out leniency instrument for examiner designs
#'
#' Builds the standard judge-leniency instrument: for each unit, the mean
#' treatment rate of its examiner computed over all *other* cases of that
#' examiner (leave-one-out), optionally after residualizing the treatment on
#' controls and within groups (court by year, for example) so that leniency
#' is compared among examiners who see the same case pool.
#'
#' @param data A data frame.
#' @param examiner Examiner identifier column name.
#' @param d Treatment (decision) column name.
#' @param x Optional control column names; the treatment is residualized on
#'   them by OLS before averaging.
#' @param group Optional grouping column name (randomization strata); the
#'   leave-out mean is computed within examiner-by-group cells and the
#'   residualization includes group fixed effects.
#' @param name Name of the new column.
#' @return `data` with the leniency column added, plus attribute
#'   `"leniency_summary"` (examiner-level standard deviation and number of
#'   examiners).
#' @export
leniency_instrument <- function(data, examiner, d, x = NULL, group = NULL, name = "leniency") {
  data <- as.data.frame(data)
  for (v in c(examiner, d, x, group)) .cm_check_column(v, data)
  dv <- as.numeric(data[[d]])
  rhs <- c(x, if (!is.null(group)) paste0("factor(", group, ")"))
  resid <- if (length(rhs)) {
    stats::residuals(stats::lm(stats::reformulate(rhs, "dv"), data = cbind(data, dv = dv)))
  } else dv
  cell <- if (is.null(group)) as.character(data[[examiner]]) else paste(data[[examiner]], data[[group]], sep = "_")
  sums <- tapply(resid, cell, sum)
  counts <- tapply(resid, cell, length)
  loo <- (sums[cell] - resid) / (counts[cell] - 1)
  loo[counts[cell] <= 1] <- NA_real_
  data[[name]] <- as.numeric(loo)
  attr(data, "leniency_summary") <- list(n_examiners = length(unique(data[[examiner]])),
                                         sd_across_cells = stats::sd(sums / counts),
                                         min_cases = min(counts))
  data
}
