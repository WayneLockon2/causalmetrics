# R/mediate-dml.R
#
# Doubly robust / double machine learning causal mediation for a binary
# treatment: the efficient influence function of Tchetgen Tchetgen and
# Shpitser (2012) in the form used by Farbmacher, Huber, Laffers, Langen, and
# Spindler (2022), with cross-fitted mlr3 nuisances.

#' Double machine learning for causal mediation
#'
#' Estimates `E[Y(d, M(d'))]` for the four combinations of a binary treatment
#' `d` and the mediator's counterfactual `M(d')`, hence the total effect, the
#' natural direct and indirect effects (both decompositions), and optionally
#' the controlled direct effect at a mediator level, using the efficient
#' influence function under sequential ignorability:
#' \deqn{\psi_{d,d'} = \frac{1\{D = d\}\, p_{d'}(M, X)}{p_d(M, X)\, p_{d'}(X)}\,(Y - \mu_d(M, X))
#'   + \frac{1\{D = d'\}}{p_{d'}(X)}\,(\mu_d(M, X) - \nu_{d,d'}(X)) + \nu_{d,d'}(X),}
#' with `mu_d(M, X) = E[Y | D = d, M, X]`, `p_d(X) = P(D = d | X)`,
#' `p_d(M, X) = P(D = d | M, X)`, and the nested regression
#' `nu_{d,d'}(X) = E[mu_d(M, X) | D = d', X]`. For `d = d'` the score is the
#' AIPW score of the average potential outcome. The mediator-density ratio
#' `f(M | d', X) / f(M | d, X)` is written through the two propensities, so
#' the mediator may be binary, continuous, or a vector. All nuisances are
#' cross-fitted; the nested regression is trained on the units with `D = d'`
#' of the training folds using the outcome model fitted on the units with
#' `D = d` of the same folds, so no observation's label enters its own
#' nuisance.
#'
#' @param data A data frame.
#' @param y,d Outcome and binary treatment column names.
#' @param m Character vector of mediator names.
#' @param x Character vector of covariate names.
#' @param learner_y `mlr3` learner for the outcome regressions (default
#'   linear regression, or logistic regression for a binary outcome).
#' @param learner_d `mlr3` classification learner for `P(D = 1 | X)` and
#'   `P(D = 1 | M, X)` (default logistic regression).
#' @param learner_nu `mlr3` regression learner for the nested regression
#'   (default `learner_y`, or linear regression when that is a classifier).
#' @param learner_m `mlr3` classification learner for `P(M = 1 | D, X)`,
#'   used only for the controlled direct effect with a binary mediator.
#' @param m_ref Mediator level for the controlled direct effect (binary
#'   mediator only); `NULL` skips it.
#' @param folds Number of cross-fitting folds.
#' @param n_rep Repeated cross-fitting partitions (median aggregation as in
#'   [est_dml()]).
#' @param trim Observations whose estimated propensities `p_1(X)` or
#'   `p_1(M, X)` fall outside `[trim, 1 - trim]` are dropped from the score
#'   (their share is reported).
#' @param conf_level Confidence level.
#' @param seed Optional seed for the folds.
#'
#' @return A list of class `cm_med_dml` with `effects` (total, nde, nie,
#'   nde_total, nie_pure, and cde when requested), `potential` (the four
#'   `E[Y(d, M(d'))]` with standard errors), `psi` (the influence functions
#'   of the representative repetition), `diagnostics` (nuisance fit,
#'   propensity summaries, share trimmed), `learners`, and settings.
#'
#' @references
#' Tchetgen Tchetgen, E. J. and Shpitser, I. (2012). Semiparametric theory
#' for causal mediation analysis. *Annals of Statistics*, 40(3), 1816-1845.
#'
#' Farbmacher, H., Huber, M., Laffers, L., Langen, H., and Spindler, M.
#' (2022). Causal mediation analysis with double machine learning.
#' *The Econometrics Journal*, 25(2), 277-300.
#'
#' @examples
#' dat <- sim_mediation(1500, dgp = "linear", seed = 1)
#' fit <- mediate_dml(dat, "y", "d", "m", x = c("x1", "x2"), seed = 1)
#' fit
#' @seealso [mediate_reg()], [est_dml()]
#' @export
mediate_dml <- function(data, y, d, m, x, learner_y = NULL, learner_d = NULL, learner_nu = NULL,
                        learner_m = NULL, m_ref = NULL, folds = 5L, n_rep = 1L, trim = 0.01,
                        conf_level = 0.95, seed = NULL) {
  .cm_require_mlr3()
  data <- as.data.frame(data)
  for (v in c(y, d, m, x)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, m, x)])
  data <- data[keep, , drop = FALSE]
  n <- nrow(data)
  yv <- as.numeric(data[[y]])
  dv <- .cm_as_binary(data[[d]], d)
  binary_y <- all(yv %in% c(0, 1))
  if (is.null(learner_y)) learner_y <- if (binary_y) .cm_default_learner("classif") else .cm_default_learner("regr")
  if (is.null(learner_d)) learner_d <- .cm_default_learner("classif")
  if (is.null(learner_nu)) learner_nu <- if (identical(learner_y$task_type, "regr")) learner_y else .cm_default_learner("regr")
  do_cde <- !is.null(m_ref)
  if (do_cde) {
    if (length(m) != 1L || !all(data[[m]] %in% c(0, 1))) stop("The controlled direct effect needs a single binary mediator.", call. = FALSE)
    if (is.null(learner_m)) learner_m <- .cm_default_learner("classif")
  }
  work <- data
  work$.cm_y <- yv
  work$.cm_d <- dv
  fold_sets <- .cm_make_fold_sets(n, folds, n_rep, seed, strata = dv)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  positive_y <- if (binary_y) "1" else NULL

  reps <- lapply(fold_sets, function(fold_id) {
    mu <- list(); nu <- list()
    p_x <- .cm_crossfit_predict(work, ".cm_d", x, learner_d, fold_id, TRUE, positive = "1", task_hint = "p_x")$pred
    p_mx <- .cm_crossfit_predict(work, ".cm_d", c(m, x), learner_d, fold_id, TRUE, positive = "1", task_hint = "p_mx")$pred
    for (dd in c(0, 1)) {
      mu[[dd + 1]] <- .cm_crossfit_predict(work, ".cm_y", c(m, x), learner_y, fold_id, TRUE, subset = dv == dd,
                                           positive = positive_y, task_hint = paste0("mu", dd))$pred
      for (dp in c(0, 1)) {
        nu[[paste(dd, dp)]] <- .cm_crossfit_nested(work, m, x, learner_y, learner_nu, fold_id, dd, dp, positive_y)
      }
    }
    pm <- NULL
    if (do_cde) {
      work$.cm_m <- as.integer(data[[m]])
      pm <- .cm_crossfit_predict(work, ".cm_m", c(d, x), learner_m, fold_id, TRUE, positive = "1", task_hint = "p_m")$pred
    }
    ok <- p_x >= trim & p_x <= 1 - trim & p_mx >= trim & p_mx <= 1 - trim
    pd <- function(dd, p) if (dd == 1) p else 1 - p
    score <- function(dd, dp) {
      mu_d <- mu[[dd + 1]]
      as.numeric(dv == dd) * pd(dp, p_mx) / (pd(dd, p_mx) * pd(dp, p_x)) * (yv - mu_d) +
        as.numeric(dv == dp) / pd(dp, p_x) * (mu_d - nu[[paste(dd, dp)]]) + nu[[paste(dd, dp)]]
    }
    # d = d': AIPW with mu_d(M, X) and its regression on X among D = d
    score_same <- function(dd) {
      mu_d <- mu[[dd + 1]]
      as.numeric(dv == dd) / pd(dd, p_x) * (yv - mu_d) + as.numeric(dv == dd) / pd(dd, p_x) * (mu_d - nu[[paste(dd, dd)]]) + nu[[paste(dd, dd)]]
    }
    psi <- cbind(`E[Y(1,M(1))]` = score_same(1), `E[Y(0,M(0))]` = score_same(0),
                 `E[Y(1,M(0))]` = score(1, 0), `E[Y(0,M(1))]` = score(0, 1))
    cde_psi <- NULL
    if (do_cde) {
      mv <- as.integer(data[[m]])
      cde_psi <- sapply(c(1, 0), function(dd) {
        mu_ref <- .cm_crossfit_at(work, ".cm_y", c(m, x), learner_y, fold_id, subset = dv == dd,
                                  modify = function(df) { df[[m]] <- m_ref; df }, positive = positive_y)
        pmr <- if (m_ref == 1) pm else 1 - pm
        as.numeric(dv == dd & mv == m_ref) / (pd(dd, p_x) * pmr) * (yv - mu_ref) + mu_ref
      })
    }
    psi <- psi[ok, , drop = FALSE]
    if (do_cde) cde_psi <- cde_psi[ok, , drop = FALSE]
    list(psi = psi, cde_psi = cde_psi, share_trimmed = mean(!ok), p_x = p_x, p_mx = p_mx, mu = mu, nu = nu)
  })

  summarise_rep <- function(r) {
    est <- unname(colMeans(r$psi))
    se <- apply(r$psi, 2, stats::sd) / sqrt(nrow(r$psi))
    eff <- c(total = est[1] - est[2], nde = est[3] - est[2], nie = est[1] - est[3],
             nde_total = est[1] - est[4], nie_pure = est[4] - est[2])
    contrasts <- rbind(total = c(1, -1, 0, 0), nde = c(0, -1, 1, 0), nie = c(1, 0, -1, 0),
                       nde_total = c(1, 0, 0, -1), nie_pure = c(0, -1, 0, 1))
    eff_se <- apply(contrasts, 1, function(w) {
      v <- r$psi %*% w
      stats::sd(v) / sqrt(length(v))
    })
    if (!is.null(r$cde_psi)) {
      v <- r$cde_psi[, 1] - r$cde_psi[, 2]
      eff <- c(eff, cde = mean(v)); eff_se <- c(eff_se, cde = stats::sd(v) / sqrt(length(v)))
    }
    list(potential = data.frame(term = colnames(r$psi), estimate = as.numeric(est), std.error = as.numeric(se)),
         effects = data.frame(term = names(eff), estimate = as.numeric(eff), std.error = as.numeric(eff_se)))
  }
  sums <- lapply(reps, summarise_rep)
  agg <- function(field) {
    terms <- sums[[1]][[field]]$term
    est_mat <- sapply(sums, function(s) s[[field]]$estimate)
    se_mat <- sapply(sums, function(s) s[[field]]$std.error)
    est_mat <- matrix(est_mat, nrow = length(terms)); se_mat <- matrix(se_mat, nrow = length(terms))
    out <- t(sapply(seq_along(terms), function(i) unlist(.cm_aggregate_reps(est_mat[i, ], se_mat[i, ]))))
    data.frame(term = terms, estimate = out[, 1], std.error = out[, 2],
               conf.low = out[, 1] - crit * out[, 2], conf.high = out[, 1] + crit * out[, 2], stringsAsFactors = FALSE)
  }
  effects <- agg("effects")
  potential <- agg("potential")
  rep_sel <- which.min(abs(sapply(sums, function(s) s$effects$estimate[3]) - effects$estimate[3]))
  r <- reps[[rep_sel]]
  diagnostics <- list(
    share_trimmed = mean(sapply(reps, `[[`, "share_trimmed")),
    propensity_x = .cm_summary(r$p_x), propensity_mx = .cm_summary(r$p_mx),
    nuisance = rbind(p_x = .cm_fit_quality(dv, r$p_x), p_mx = .cm_fit_quality(dv, r$p_mx),
                     mu0 = .cm_fit_quality(yv[dv == 0], r$mu[[1]][dv == 0]), mu1 = .cm_fit_quality(yv[dv == 1], r$mu[[2]][dv == 1]))
  )
  structure(list(effects = effects, potential = potential, psi = r$psi, cde_psi = r$cde_psi,
                 diagnostics = diagnostics, n = n, n_used = nrow(r$psi), folds = folds, n_rep = n_rep,
                 learners = list(y = .cm_learner_label(learner_y), d = .cm_learner_label(learner_d),
                                 nu = .cm_learner_label(learner_nu)),
                 spec = list(y = y, d = d, m = m, x = x, m_ref = m_ref, trim = trim), conf_level = conf_level,
                 call = match.call()), class = "cm_med_dml")
}

# Cross-fitted predictions on held-out rows after modifying them (for
# counterfactual mediator levels).
.cm_crossfit_at <- function(work, target, features, learner, fold_id, subset, modify, positive = NULL) {
  n <- nrow(work)
  pred <- rep(NA_real_, n)
  for (k in sort(unique(fold_id))) {
    te <- which(fold_id == k)
    tr <- which(fold_id != k & subset)
    fit <- .cm_fit_mlr3(learner, work[tr, , drop = FALSE], target, features, positive = positive, task_hint = "mu_at")
    pred[te] <- .cm_predict_fit(fit, modify(work[te, , drop = FALSE]))
  }
  pred
}

# Nested regression nu_{d,d'}(X) = E[mu_d(M, X) | D = d', X], cross-fitted: for
# each held-out fold, fit mu_d on the D = d units of the training folds,
# evaluate it on the D = d' units of the training folds, regress those values
# on X, and predict the held-out fold.
.cm_crossfit_nested <- function(work, m, x, learner_y, learner_nu, fold_id, dd, dp, positive_y) {
  n <- nrow(work)
  pred <- rep(NA_real_, n)
  for (k in sort(unique(fold_id))) {
    te <- which(fold_id == k)
    tr_d <- which(fold_id != k & work$.cm_d == dd)
    tr_dp <- which(fold_id != k & work$.cm_d == dp)
    fit_mu <- .cm_fit_mlr3(learner_y, work[tr_d, , drop = FALSE], ".cm_y", c(m, x), positive = positive_y, task_hint = "mu_nested")
    lab <- .cm_predict_fit(fit_mu, work[tr_dp, , drop = FALSE])
    nu_data <- work[tr_dp, , drop = FALSE]
    nu_data$.cm_lab <- lab
    fit_nu <- .cm_fit_mlr3(learner_nu, nu_data, ".cm_lab", x, task_hint = "nu")
    pred[te] <- .cm_predict_fit(fit_nu, work[te, , drop = FALSE])
  }
  pred
}

#' @export
print.cm_med_dml <- function(x, ...) {
  cat("Double machine learning mediation (", x$folds, " folds, ", x$n_rep, " repetition(s), n = ", x$n,
      if (x$diagnostics$share_trimmed > 0) paste0(", ", round(100 * x$diagnostics$share_trimmed, 1), "% trimmed"), ")\n", sep = "")
  cat("  learners: outcome = ", x$learners$y, ", propensities = ", x$learners$d, ", nested = ", x$learners$nu, "\n", sep = "")
  print(x$effects, digits = 4, row.names = FALSE)
  invisible(x)
}
