# R/mediate-iv.R
#
# Mediation with an endogenous mediator instrumented by Z (Cattan, Salvanes,
# and Tominey 2025; Dippel, Ferrara, Heblich, and Pinto 2020; Huber 2019).

#' Mediation analysis with an instrumented mediator
#'
#' Decomposes the total effect of a treatment `D` on `Y` into a direct effect
#' and an indirect effect through a mediator `M` that is endogenous, using an
#' instrument `Z` for the mediator. Three regressions, all with the same
#' controls and fixed effects (through `fixest`):
#'
#' * first stage: `M ~ D + Z + X`, giving the effect of `D` on `M`
#'   (`delta_1`) and the instrument's strength;
#' * second stage: `Y ~ D + M + X` with `M` instrumented by `Z`, giving the
#'   direct effect `gamma_1` and the mediator effect `gamma_2`;
#' * reduced form: `Y ~ D + X` (the total effect) and `Y ~ D + Z + X`.
#'
#' The indirect effect is `delta_1 * gamma_2`, and the direct plus indirect
#' effect equals the coefficient on `D` in the reduced form that includes
#' `Z` exactly. Beyond the usual instrument validity conditional on `D` and
#' `X`, the decomposition assumes a homogeneous mediator effect (the
#' compliers' `gamma_2` equals the population's; Huber 2019). The
#' `homogeneity_by` argument re-estimates `gamma_2` within subgroups as a
#' check. Standard errors come from a bootstrap of the whole system,
#' clustered when `cluster` is given.
#'
#' @param data A data frame.
#' @param y,d,m,z Outcome, treatment, mediator, and instrument column names
#'   (`z` may be a vector of instruments).
#' @param x Optional control column names.
#' @param fe Optional fixed-effect column names.
#' @param cluster Optional cluster column for the bootstrap.
#' @param weights Optional column of regression weights (passed to every
#'   `fixest` regression and carried through the bootstrap).
#' @param homogeneity_by Optional column defining subgroups for the
#'   homogeneity check.
#' @param n_boot Bootstrap replications.
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_med_iv` with `effects` (total from the reduced
#'   form without `Z`, `total_z` with `Z`, direct, indirect, first_stage,
#'   mediator_effect, and the naive OLS direct effect), `first_stage_F`,
#'   `homogeneity` (subgroup mediator effects), the fitted `fixest` models,
#'   and settings.
#'
#' @references
#' Cattan, S., Salvanes, K. G., and Tominey, E. (2025). First-generation
#' elite: the role of school social networks. *American Economic Review*,
#' 115(12), 4369-4403.
#'
#' Dippel, C., Ferrara, A., and Heblich, S. (2020). Causal mediation analysis
#' in instrumental-variables regressions. *Stata Journal*, 20(3), 613-626.
#'
#' Huber, M. (2019). A review of causal mediation analysis for assessing
#' direct and indirect treatment effects. IZA Discussion Paper 12660.
#'
#' @examples
#' dat <- sim_mediation(2000, dgp = "iv_mediator", seed = 1)
#' fit <- mediate_iv(dat, "y", "d", "m", z = "z", x = c("x1", "x2"), n_boot = 99, seed = 1)
#' fit
#' attr(dat, "truth")[c("direct", "indirect")]
#' @export
mediate_iv <- function(data, y, d, m, z, x = NULL, fe = NULL, cluster = NULL, weights = NULL,
                       homogeneity_by = NULL, n_boot = 499L, conf_level = 0.95, seed = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, m, z, x, fe, cluster, weights, homogeneity_by)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, m, z, x, fe, cluster, weights)])
  data <- data[keep, , drop = FALSE]
  wf <- if (is.null(weights)) NULL else stats::as.formula(paste0("~", weights))
  fe_part <- if (!is.null(fe)) paste(" |", paste(fe, collapse = " + ")) else ""
  f_first <- stats::as.formula(paste(m, "~", paste(c(d, z, x), collapse = " + "), fe_part))
  f_rf <- stats::as.formula(paste(y, "~", paste(c(d, x), collapse = " + "), fe_part))
  f_rfz <- stats::as.formula(paste(y, "~", paste(c(d, z, x), collapse = " + "), fe_part))
  f_ols <- stats::as.formula(paste(y, "~", paste(c(d, m, x), collapse = " + "), fe_part))
  f_iv <- stats::as.formula(paste(y, "~", paste(c(d, x), collapse = " + "),
                                  if (is.null(fe)) " | " else fe_part, if (!is.null(fe)) " | " else "",
                                  m, " ~ ", paste(z, collapse = " + ")))
  fit_all <- function(dd) {
    first <- fixest::feols(f_first, data = dd, weights = wf, warn = FALSE, notes = FALSE)
    iv <- fixest::feols(f_iv, data = dd, weights = wf, warn = FALSE, notes = FALSE)
    rf <- fixest::feols(f_rf, data = dd, weights = wf, warn = FALSE, notes = FALSE)
    rfz <- fixest::feols(f_rfz, data = dd, weights = wf, warn = FALSE, notes = FALSE)
    ols <- fixest::feols(f_ols, data = dd, weights = wf, warn = FALSE, notes = FALSE)
    delta1 <- stats::coef(first)[[d]]
    gamma1 <- stats::coef(iv)[[d]]
    gamma2 <- stats::coef(iv)[[paste0("fit_", m)]]
    list(stat = c(total = stats::coef(rf)[[d]], total_z = stats::coef(rfz)[[d]], direct = gamma1,
                  indirect = delta1 * gamma2, first_stage = delta1, mediator_effect = gamma2,
                  naive_direct = stats::coef(ols)[[d]], naive_mediator_effect = stats::coef(ols)[[m]]),
         models = list(first_stage = first, iv = iv, reduced_form = rf, reduced_form_z = rfz, ols = ols))
  }
  point <- fit_all(data)
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  effects <- .cm_boot(data, function(dd) fit_all(dd)$stat, n_boot, cluster = cl, seed = seed, conf_level = conf_level)
  fs_F <- tryCatch({
    fs <- fixest::fitstat(point$models$iv, "ivf")
    fs[[1]]$stat
  }, error = function(e) NA_real_)
  homogeneity <- NULL
  if (!is.null(homogeneity_by)) {
    g <- data[[homogeneity_by]]
    homogeneity <- do.call(rbind, lapply(sort(unique(g)), function(lev) {
      sub <- data[g == lev, , drop = FALSE]
      r <- tryCatch(fit_all(sub)$stat, error = function(e) NULL)
      if (is.null(r)) return(NULL)
      data.frame(group = lev, n = nrow(sub), mediator_effect = r[["mediator_effect"]], direct = r[["direct"]],
                 first_stage = r[["first_stage"]])
    }))
  }
  structure(list(effects = effects, first_stage_F = fs_F, homogeneity = homogeneity, models = point$models,
                 n = nrow(data), spec = list(y = y, d = d, m = m, z = z, x = x, fe = fe, weights = weights), conf_level = conf_level,
                 call = match.call()), class = "cm_med_iv")
}

#' @export
print.cm_med_iv <- function(x, ...) {
  cat("Mediation with an instrumented mediator (", x$spec$m, " instrumented by ", paste(x$spec$z, collapse = ", "),
      "; n = ", x$n, ")\n", sep = "")
  print(x$effects, digits = 4, row.names = FALSE)
  cat("  first-stage F for the instrument(s): ", format(round(x$first_stage_F, 2)), "\n", sep = "")
  if (!is.null(x$homogeneity)) {
    cat("  mediator effect by subgroup (homogeneity check):\n")
    print(x$homogeneity, digits = 4, row.names = FALSE)
  }
  invisible(x)
}
