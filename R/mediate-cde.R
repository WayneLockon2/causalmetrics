# R/mediate-cde.R
#
# Controlled direct effects by sequential g-estimation (Acharya, Blackwell,
# and Sen 2016), which allows mediator-outcome confounders that are
# themselves affected by the treatment.

#' Controlled direct effect by sequential g-estimation
#'
#' Estimates the controlled direct effect `E[Y(1, m) - Y(0, m)]` of a binary
#' treatment holding the mediator at `m_ref`, when some confounders of the
#' mediator-outcome relationship (`x_post`) are affected by the treatment, so
#' that conditioning on them in one regression would bias the treatment
#' coefficient (they are post-treatment) and omitting them would bias the
#' mediator coefficient. The two-step procedure of Acharya, Blackwell, and
#' Sen (2016):
#'
#' 1. regress `Y` on `D`, `M`, `x_pre`, and `x_post` (with `D:M` when
#'    `interaction = TRUE`), and estimate the mediator's effect `delta`;
#' 2. form the demediated outcome `Y - delta * (M - m_ref)`;
#' 3. regress the demediated outcome on `D` and `x_pre` only; the coefficient
#'    on `D` is the controlled direct effect.
#'
#' Standard errors come from a bootstrap of both steps (clustered when
#' `cluster` is given), which accounts for the estimated `delta`.
#'
#' @param data A data frame.
#' @param y,d,m Outcome, binary treatment, and mediator column names.
#' @param x_pre Pre-treatment covariates.
#' @param x_post Post-treatment covariates that confound the mediator and the
#'   outcome (intermediate confounders).
#' @param m_ref Mediator level for the controlled direct effect.
#' @param interaction Allow a `D:M` interaction in the first step (the
#'   demediation then uses `delta + delta_int * D`).
#' @param cluster Optional cluster column for the bootstrap.
#' @param n_boot Bootstrap replications.
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_med_cde` with `effects` (`cde` and the naive
#'   single-regression coefficient on `D` for comparison), `delta` (the
#'   mediator effect from step 1), the two fitted models, and settings.
#'
#' @references
#' Acharya, A., Blackwell, M., and Sen, M. (2016). Explaining causal findings
#' without bias: detecting and assessing direct effects. *American Political
#' Science Review*, 110(3), 512-529.
#'
#' @examples
#' dat <- sim_mediation(1500, dgp = "post_treatment_confounder", seed = 1)
#' fit <- mediate_cde(dat, "y", "d", "m", x_pre = c("x1", "x2"), x_post = "z", n_boot = 99, seed = 1)
#' fit
#' attr(dat, "truth")$cde
#' @export
mediate_cde <- function(data, y, d, m, x_pre = NULL, x_post = NULL, m_ref = 0, interaction = FALSE,
                        cluster = NULL, n_boot = 499L, conf_level = 0.95, seed = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, m, x_pre, x_post, cluster)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, m, x_pre, x_post, cluster)])
  data <- data[keep, , drop = FALSE]
  data[[d]] <- .cm_as_binary(data[[d]], d)
  stat <- function(dd) {
    rhs1 <- c(d, m, x_pre, x_post)
    if (interaction) rhs1 <- c(rhs1, paste0(d, ":", m))
    step1 <- stats::lm(.cm_med_formula(y, rhs1), data = dd)
    cf <- stats::coef(step1)
    delta <- cf[[m]] + if (interaction) cf[[paste0(d, ":", m)]] * dd[[d]] else 0
    dd$.cm_ytilde <- dd[[y]] - delta * (dd[[m]] - m_ref)
    step2 <- stats::lm(.cm_med_formula(".cm_ytilde", c(d, x_pre)), data = dd)
    naive <- stats::lm(.cm_med_formula(y, c(d, m, x_pre, x_post)), data = dd)
    c(cde = stats::coef(step2)[[d]], naive_direct = stats::coef(naive)[[d]], delta = cf[[m]],
      total = stats::coef(stats::lm(.cm_med_formula(y, c(d, x_pre)), data = dd))[[d]])
  }
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  tab <- .cm_boot(data, stat, n_boot, cluster = cl, seed = seed, conf_level = conf_level)
  point <- stat(data)
  rhs1 <- c(d, m, x_pre, x_post)
  if (interaction) rhs1 <- c(rhs1, paste0(d, ":", m))
  step1 <- stats::lm(.cm_med_formula(y, rhs1), data = data)
  data$.cm_ytilde <- data[[y]] - (stats::coef(step1)[[m]] + if (interaction) stats::coef(step1)[[paste0(d, ":", m)]] * data[[d]] else 0) * (data[[m]] - m_ref)
  step2 <- stats::lm(.cm_med_formula(".cm_ytilde", c(d, x_pre)), data = data)
  structure(list(effects = tab, delta = point[["delta"]], models = list(step1 = step1, step2 = step2),
                 n = nrow(data), spec = list(y = y, d = d, m = m, x_pre = x_pre, x_post = x_post, m_ref = m_ref, interaction = interaction),
                 conf_level = conf_level, call = match.call()), class = "cm_med_cde")
}

#' @export
print.cm_med_cde <- function(x, ...) {
  cat("Controlled direct effect by sequential g-estimation (mediator ", x$spec$m, " fixed at ", x$spec$m_ref,
      "; n = ", x$n, ")\n", sep = "")
  print(x$effects, digits = 4, row.names = FALSE)
  cat("  `naive_direct` is the treatment coefficient in one regression that conditions on the post-treatment covariates.\n")
  invisible(x)
}
