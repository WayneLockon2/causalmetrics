# R/mediate-reg.R
#
# Regression-based causal mediation: product of coefficients, the regression
# with treatment-mediator interaction, and g-computation with quasi-Bayesian
# or bootstrap inference, for one or several parallel mediators.

#' Regression-based causal mediation analysis
#'
#' Estimates the total effect, the natural direct effect (NDE), the natural
#' indirect effect (NIE), and the controlled direct effect (CDE) of a binary
#' treatment through one or several parallel mediators from a mediator model
#' `M ~ D + X` and an outcome model `Y ~ D + M + X` (optionally with `D:M`),
#' under sequential ignorability (Imai, Keele, and Yamamoto 2010): no
#' unmeasured treatment-outcome, treatment-mediator, or mediator-outcome
#' confounding given `X`, and no mediator-outcome confounder affected by the
#' treatment.
#'
#' Effects are computed by g-computation over the sample: for each unit the
#' mediator model gives the distribution of `M(d')` and the outcome model the
#' mean of `Y(d, M(d'))`; averaging over units gives `E[Y(d, M(d'))]` for the
#' four combinations of `d, d'`. With linear models and no interaction this
#' reduces to the product of coefficients (Baron and Kenny 1986; Sobel 1982),
#' and with the interaction to the closed forms of VanderWeele (2015). Binary
#' mediators or outcomes use logistic models and the effects are on the
#' probability scale.
#'
#' @param data A data frame.
#' @param y,d Outcome and binary treatment column names.
#' @param m Character vector of mediator names (one mediator, or several
#'   parallel mediators that do not cause each other).
#' @param x Optional covariate names.
#' @param interaction Include treatment-mediator interactions in the outcome
#'   model.
#' @param outcome,mediator `"linear"` or `"logit"` model for the outcome and
#'   for the mediators (recycled across mediators).
#' @param method Inference: `"simulation"` (quasi-Bayesian draws of the
#'   coefficients from their robust covariance, Imai et al. 2010),
#'   `"bootstrap"` (nonparametric, clustered when `cluster` is given), or
#'   `"delta"` (Sobel's formula for a single mediator with linear models and
#'   no interaction, treating the two regressions as independent).
#' @param m_ref Mediator value(s) at which the controlled direct effect is
#'   evaluated (default 0 for each mediator).
#' @param cluster Optional cluster column for the covariance or the bootstrap.
#' @param weights Optional column of regression weights (used in every
#'   model and in the averaging over units).
#' @param n_boot Bootstrap replications.
#' @param n_sim Coefficient draws for `method = "simulation"`, and Monte Carlo
#'   draws of a continuous mediator inside a nonlinear outcome model.
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_mediation` with `effects` (term, estimate,
#'   std.error, conf.low, conf.high for total, nde, nie, nde_total, nie_pure,
#'   cde, and per-mediator `nie_<m>` when several mediators are given),
#'   `shares` (proportion mediated with a standard error from the draws or
#'   the bootstrap; with linear models also the two shares
#'   of Wheeler et al. 2022, `S1 = alpha * gamma / beta` with `alpha` the
#'   mediator coefficient in the outcome model, `gamma` the treatment effect
#'   on the mediator, `beta` the total effect, and `S2 = delta * gamma / beta`
#'   with `delta` the mediator-outcome association in the control group,
#'   whose standard errors come from the delta method on the stacked
#'   influence functions of the regressions, so cross-equation covariances
#'   are carried, as in Stata's `suest`),
#'   `potential` (the four `E[Y(d, M(d'))]`), `models`, and the settings.
#'
#' @references
#' Imai, K., Keele, L., and Yamamoto, T. (2010). Identification, inference and
#' sensitivity analysis for causal mediation effects. *Statistical Science*,
#' 25(1), 51-71.
#'
#' VanderWeele, T. J. (2015). *Explanation in Causal Inference*. Oxford.
#'
#' Wheeler, L., Garlick, R., Johnson, E., Shaw, P., and Gargano, M. (2022).
#' LinkedIn(to) job opportunities. *AEJ: Applied Economics*, 14(2), 101-125.
#'
#' @examples
#' dat <- sim_mediation(1000, dgp = "linear", seed = 1)
#' fit <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), seed = 1)
#' fit
#' attr(dat, "truth")[c("total", "nde", "nie")]
#' @seealso [mediate_dml()], [mediate_sensitivity()], [plot_mediation()]
#' @export
mediate_reg <- function(data, y, d, m, x = NULL, interaction = FALSE,
                        outcome = c("linear", "logit"), mediator = c("linear", "logit"),
                        method = c("simulation", "bootstrap", "delta"), m_ref = NULL,
                        cluster = NULL, weights = NULL, n_boot = 499L, n_sim = 1000L, conf_level = 0.95, seed = NULL) {
  outcome <- match.arg(outcome)
  mediator <- match.arg(mediator)
  method <- match.arg(method)
  data <- as.data.frame(data)
  for (v in c(y, d, m, x, cluster, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, m, x, cluster, weights)])
  data <- data[keep, , drop = FALSE]
  data[[d]] <- .cm_as_binary(data[[d]], d)
  data$.cm_w <- if (is.null(weights)) rep(1, nrow(data)) else as.numeric(data[[weights]])
  if (any(data$.cm_w < 0)) stop("`weights` must be non-negative.", call. = FALSE)
  J <- length(m)
  if (is.null(m_ref)) m_ref <- rep(0, J)
  if (length(m_ref) != J) stop("`m_ref` must have one value per mediator.", call. = FALSE)
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  if (method == "delta" && (J > 1L || interaction || outcome != "linear" || mediator != "linear")) {
    stop("`method = \"delta\"` is for one mediator with linear models and no interaction; use \"simulation\" or \"bootstrap\".", call. = FALSE)
  }
  spec <- list(y = y, d = d, m = m, x = x, interaction = interaction, outcome = outcome,
               mediator = mediator, m_ref = m_ref, n_sim = n_sim, weights = weights)

  fit <- .cm_med_fit(data, spec, cl)
  point <- .cm_med_effects(fit, data, spec, draws = NULL)
  shares <- .cm_med_shares(fit, point, spec, data = data, cluster = cl)

  if (method == "simulation") {
    draws <- .cm_with_seed(seed, .cm_med_simulate(fit, data, spec, n_sim))
    effects <- .cm_sim_table(point, draws, conf_level)
    if (is.na(shares$std.error[1])) shares$std.error[1] <- stats::sd(draws[, "nie"] / draws[, "total"], na.rm = TRUE)
  } else if (method == "bootstrap") {
    stat <- function(dd) {
      e <- .cm_med_effects(.cm_med_fit(dd, spec, NULL), dd, spec, draws = NULL)
      c(e, `proportion mediated` = unname(e["nie"] / e["total"]))
    }
    tab <- .cm_boot(data, stat, n_boot, cluster = cl, seed = seed, conf_level = conf_level)
    effects <- tab[tab$term != "proportion mediated", ]
    rownames(effects) <- NULL
    shares$std.error[1] <- tab$std.error[tab$term == "proportion mediated"]
  } else {
    effects <- .cm_med_delta(fit, point, conf_level)
  }
  is_po <- grepl("^E\\[", effects$term)
  potential <- effects[is_po, ]
  effects <- effects[!is_po, ]
  rownames(potential) <- rownames(effects) <- NULL
  structure(list(effects = effects, shares = shares, potential = potential,
                 models = fit$models, spec = spec, method = method, n = nrow(data), conf_level = conf_level,
                 call = match.call()), class = "cm_mediation")
}

# Fit the mediator, outcome, total, and control-group models.
.cm_med_fit <- function(data, spec, cluster) {
  J <- length(spec$m)
  med_models <- lapply(spec$m, function(mj) {
    f <- .cm_med_formula(mj, c(spec$d, spec$x))
    if (spec$mediator == "logit") stats::glm(f, data = data, family = stats::binomial(), weights = .cm_w) else stats::lm(f, data = data, weights = .cm_w)
  })
  rhs <- c(spec$d, spec$m, spec$x)
  if (spec$interaction) rhs <- c(rhs, paste0(spec$d, ":", spec$m))
  f_y <- .cm_med_formula(spec$y, rhs)
  out_model <- if (spec$outcome == "logit") stats::glm(f_y, data = data, family = stats::binomial(), weights = .cm_w) else stats::lm(f_y, data = data, weights = .cm_w)
  f_t <- .cm_med_formula(spec$y, c(spec$d, spec$x))
  tot_model <- if (spec$outcome == "logit") stats::glm(f_t, data = data, family = stats::binomial(), weights = .cm_w) else stats::lm(f_t, data = data, weights = .cm_w)
  ctrl <- data[data[[spec$d]] == 0, , drop = FALSE]
  f_c <- .cm_med_formula(spec$y, c(spec$m, spec$x))
  ctrl_model <- if (spec$outcome == "logit") stats::glm(f_c, data = ctrl, family = stats::binomial(), weights = .cm_w) else stats::lm(f_c, data = ctrl, weights = .cm_w)
  vc <- function(fit, cl) .cm_vcov_robust(fit, cl)
  cl_ctrl <- if (is.null(cluster)) NULL else cluster[data[[spec$d]] == 0]
  list(models = list(mediator = med_models, outcome = out_model, total = tot_model, control = ctrl_model),
       coef = list(mediator = lapply(med_models, stats::coef), outcome = stats::coef(out_model),
                   total = stats::coef(tot_model), control = stats::coef(ctrl_model)),
       vcov = list(mediator = lapply(med_models, vc, cluster), outcome = vc(out_model, cluster),
                   total = vc(tot_model, cluster), control = vc(ctrl_model, cl_ctrl)),
       sigma = vapply(med_models, function(f) if (inherits(f, "glm")) NA_real_ else stats::sigma(f), numeric(1)))
}

# Design matrix of the covariates as the regressions expand them (factors
# become dummies in formula order), without the intercept.
.cm_med_X <- function(data, x) {
  if (!length(x)) return(matrix(0, nrow(data), 0))
  f <- stats::as.formula(paste("~", paste(x, collapse = " + ")))
  stats::model.matrix(f, data = data)[, -1, drop = FALSE]
}

# Potential-outcome means E[Y(d, M(d'))] by g-computation, from coefficient
# vectors (`draws` replaces the point estimates when given).
.cm_med_effects <- function(fit, data, spec, draws = NULL) {
  J <- length(spec$m)
  n <- nrow(data)
  needs_mc <- spec$outcome == "logit" && !(spec$mediator == "logit" && J == 1L)
  if (needs_mc) {
    R <- min(spec$n_sim, 200L)
    mc_norm <- .cm_with_seed(20240905, matrix(stats::rnorm(n * R), n, R))
    mc_unif <- .cm_with_seed(20240906, matrix(stats::runif(n * R), n, R))
  }
  w <- if (".cm_w" %in% names(data)) data$.cm_w else rep(1, n)
  wmean <- function(v) sum(w * v) / sum(w)
  cf_m <- if (is.null(draws)) fit$coef$mediator else draws$mediator
  cf_y <- if (is.null(draws)) fit$coef$outcome else draws$outcome
  X <- .cm_med_X(data, spec$x)
  # mediator means (or probabilities) under d'
  med_mean <- function(j, dd) {
    b <- cf_m[[j]]
    eta <- b[1] + b[2] * dd + if (ncol(X)) as.numeric(X %*% b[-(1:2)]) else 0
    if (spec$mediator == "logit") stats::plogis(eta) else eta
  }
  # outcome mean given d and a mediator matrix (n x J)
  out_mean <- function(dd, M) {
    b <- cf_y
    eta <- b[1] + b[2] * dd
    for (j in seq_len(J)) eta <- eta + b[2 + j] * M[, j]
    if (ncol(X)) eta <- eta + as.numeric(X %*% b[2 + J + seq_len(ncol(X))])
    if (spec$interaction) for (j in seq_len(J)) eta <- eta + b[2 + J + ncol(X) + j] * dd * M[, j]
    if (spec$outcome == "logit") stats::plogis(eta) else eta
  }
  # E[Y(d, M(d'))] with mediators j in `set1` drawn under d1 and the others under d0
  po <- function(dd, d_set) {
    means <- sapply(seq_len(J), function(j) med_mean(j, d_set[j]))
    means <- matrix(means, n, J)
    if (spec$outcome == "linear" || (spec$mediator == "logit" && J == 1L && spec$outcome == "logit")) {
      if (spec$outcome == "logit" && spec$mediator == "logit") {
        p1 <- means[, 1]
        return(wmean(p1 * out_mean(dd, matrix(1, n, 1)) + (1 - p1) * out_mean(dd, matrix(0, n, 1))))
      }
      return(wmean(out_mean(dd, means)))
    }
    # nonlinear outcome: integrate over the mediator distribution by Monte Carlo
    acc <- 0
    for (r in seq_len(R)) {
      Mr <- means
      for (j in seq_len(J)) {
        Mr[, j] <- if (spec$mediator == "logit") as.numeric(mc_unif[, r] < means[, j]) else means[, j] + fit$sigma[j] * mc_norm[, r]
      }
      acc <- acc + wmean(out_mean(dd, Mr))
    }
    acc / R
  }
  e11 <- po(1, rep(1, J)); e00 <- po(0, rep(0, J)); e10 <- po(1, rep(0, J)); e01 <- po(0, rep(1, J))
  M_ref <- matrix(spec$m_ref, n, J, byrow = TRUE)
  cde <- wmean(out_mean(1, M_ref)) - wmean(out_mean(0, M_ref))
  out <- c(total = e11 - e00, nde = e10 - e00, nie = e11 - e10, nde_total = e11 - e01, nie_pure = e01 - e00, cde = cde)
  if (J > 1L) {
    for (j in seq_len(J)) {
      d_set <- rep(0, J); d_set[j] <- 1
      out[paste0("nie_", spec$m[j])] <- po(1, d_set) - e10
    }
  }
  c(out, `E[Y(1,M(1))]` = e11, `E[Y(0,M(0))]` = e00, `E[Y(1,M(0))]` = e10, `E[Y(0,M(1))]` = e01)
}

.cm_med_simulate <- function(fit, data, spec, n_sim) {
  J <- length(spec$m)
  dm <- lapply(seq_len(J), function(j) .cm_draw_coefs(fit$coef$mediator[[j]], fit$vcov$mediator[[j]], n_sim))
  dy <- .cm_draw_coefs(fit$coef$outcome, fit$vcov$outcome, n_sim)
  t(sapply(seq_len(n_sim), function(s) {
    .cm_med_effects(fit, data, spec, draws = list(mediator = lapply(dm, function(D) D[s, ]), outcome = dy[s, ]))
  }))
}

.cm_med_shares <- function(fit, point, spec, data = NULL, cluster = NULL) {
  J <- length(spec$m)
  out <- data.frame(share = "proportion mediated", estimate = unname(point["nie"] / point["total"]),
                    std.error = NA_real_, stringsAsFactors = FALSE)
  if (spec$outcome == "linear" && spec$mediator == "linear" && !spec$interaction) {
    beta <- fit$coef$total[[spec$d]]
    # stacked influence functions of the four regressions for delta-method
    # standard errors that carry the cross-equation covariance
    V <- NULL
    if (!is.null(data)) {
      n <- nrow(data)
      ctrl_rows <- which(data[[spec$d]] == 0)
      infl <- cbind(.cm_infl(fit$models$total, n)[, spec$d, drop = FALSE],
                    do.call(cbind, lapply(seq_len(J), function(j) .cm_infl(fit$models$mediator[[j]], n)[, spec$d, drop = FALSE])),
                    .cm_infl(fit$models$outcome, n)[, spec$m, drop = FALSE],
                    .cm_infl(fit$models$control, n, ctrl_rows)[, spec$m, drop = FALSE])
      colnames(infl) <- c("beta", paste0("gamma", seq_len(J)), paste0("alpha", seq_len(J)), paste0("delta", seq_len(J)))
      V <- .cm_stacked_vcov(infl, cluster)
    }
    for (j in seq_len(J)) {
      gamma <- fit$coef$mediator[[j]][[spec$d]]
      alpha <- fit$coef$outcome[[spec$m[j]]]
      delta <- fit$coef$control[[spec$m[j]]]
      s1 <- alpha * gamma / beta; s2 <- delta * gamma / beta
      se1 <- se2 <- NA_real_
      if (!is.null(V)) {
        g1 <- stats::setNames(numeric(ncol(V)), colnames(V))
        g1[c("beta", paste0("gamma", j), paste0("alpha", j))] <- c(-s1 / beta, s1 / gamma, s1 / alpha)
        g2 <- stats::setNames(numeric(ncol(V)), colnames(V))
        g2[c("beta", paste0("gamma", j), paste0("delta", j))] <- c(-s2 / beta, s2 / gamma, s2 / delta)
        se1 <- sqrt(as.numeric(t(g1) %*% V %*% g1)); se2 <- sqrt(as.numeric(t(g2) %*% V %*% g2))
      }
      out <- rbind(out,
                   data.frame(share = paste0("S1 (alpha*gamma/beta)", if (J > 1) paste0(": ", spec$m[j])), estimate = s1, std.error = se1),
                   data.frame(share = paste0("S2 (delta*gamma/beta)", if (J > 1) paste0(": ", spec$m[j])), estimate = s2, std.error = se2))
    }
    if (J == 1L) out$std.error[1] <- out$std.error[2]
  }
  rownames(out) <- NULL
  out
}

# Sobel-type delta method for one mediator, linear models, no interaction.
.cm_med_delta <- function(fit, point, conf_level) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  a1 <- fit$coef$mediator[[1]][[2]]; va1 <- fit$vcov$mediator[[1]][2, 2]
  b1 <- fit$coef$outcome[[2]]; vb1 <- fit$vcov$outcome[2, 2]
  b2 <- fit$coef$outcome[[3]]; vb2 <- fit$vcov$outcome[3, 3]
  beta <- fit$coef$total[[2]]; vbeta <- fit$vcov$total[2, 2]
  se <- c(total = sqrt(vbeta), nde = sqrt(vb1), nie = sqrt(a1^2 * vb2 + b2^2 * va1),
          nde_total = sqrt(vb1), nie_pure = sqrt(a1^2 * vb2 + b2^2 * va1), cde = sqrt(vb1))
  po <- point[grepl("^E\\[", names(point))]
  se <- c(se, stats::setNames(rep(NA_real_, length(po)), names(po)))
  est <- unname(point[names(se)])
  data.frame(term = names(se), estimate = as.numeric(est), std.error = as.numeric(se),
             conf.low = as.numeric(est - crit * se), conf.high = as.numeric(est + crit * se),
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @export
print.cm_mediation <- function(x, ...) {
  cat("Causal mediation by regression (", x$spec$outcome, " outcome, ", x$spec$mediator, " mediator",
      if (length(x$spec$m) > 1) paste0("s (", length(x$spec$m), ")"), if (x$spec$interaction) ", with D:M interaction",
      "; inference: ", x$method, ")\n", sep = "")
  cat("  n = ", x$n, "; treatment ", x$spec$d, "; mediator", if (length(x$spec$m) > 1) "s", ": ", paste(x$spec$m, collapse = ", "), "\n", sep = "")
  print(x$effects, digits = 4, row.names = FALSE)
  sh <- x$shares
  cat("  proportion mediated = ", format(round(sh$estimate[1], 3)),
      if (!is.na(sh$std.error[1])) paste0(" (", format(round(sh$std.error[1], 3)), ")"), "\n", sep = "")
  if (nrow(sh) > 1) {
    for (i in 2:nrow(sh)) cat("  ", sh$share[i], " = ", format(round(sh$estimate[i], 3)),
                              if (!is.na(sh$std.error[i])) paste0(" (", format(round(sh$std.error[i], 3)), ")"), "\n", sep = "")
  }
  invisible(x)
}

#' Sensitivity of the natural indirect effect to an unmeasured mediator-outcome confounder
#'
#' For a `cm_mediation` object fitted with one mediator, linear models, and no
#' interaction, computes the natural indirect effect as a function of the
#' correlation `rho` between the errors of the mediator model and the outcome
#' model (Imai, Keele, and Yamamoto 2010, Theorem 3):
#' \deqn{NIE(\rho) = \beta_2 \frac{\sigma_1}{\sigma_2}
#'   \Big(\tilde\rho - \rho \sqrt{\frac{1 - \tilde\rho^2}{1 - \rho^2}}\Big),}
#' where `beta_2` is the treatment effect on the mediator, `sigma_1` and
#' `sigma_2` the residual standard deviations of the total-effect regression
#' `Y ~ D + X` and of the mediator regression, and `tilde rho` the
#' correlation of their residuals. `rho = 0` is sequential ignorability.
#'
#' @param fit A `cm_mediation` object from [mediate_reg()].
#' @param rho Grid of error correlations.
#' @return A list of class `cm_med_sensitivity` with `curve` (rho, nie), the
#'   value `rho_zero` at which the indirect effect crosses zero, and the
#'   implied `rho_tilde`.
#' @examples
#' dat <- sim_mediation(1000, dgp = "correlated_errors", rho = 0.4, seed = 1)
#' fit <- mediate_reg(dat, "y", "d", "m", x = c("x1", "x2"), seed = 1)
#' sens <- mediate_sensitivity(fit)
#' sens$rho_zero
#' @export
mediate_sensitivity <- function(fit, rho = seq(-0.9, 0.9, by = 0.05)) {
  if (!inherits(fit, "cm_mediation")) stop("`fit` must come from mediate_reg().", call. = FALSE)
  sp <- fit$spec
  if (length(sp$m) > 1L || sp$interaction || sp$outcome != "linear" || sp$mediator != "linear") {
    stop("The sensitivity formula needs one mediator, linear models, and no interaction.", call. = FALSE)
  }
  e1 <- stats::residuals(fit$models$total)
  e2 <- stats::residuals(fit$models$mediator[[1]])
  s1 <- stats::sd(e1); s2 <- stats::sd(e2)
  rho_tilde <- stats::cor(e1, e2)
  beta2 <- fit$models$mediator[[1]]$coefficients[[sp$d]]
  nie <- beta2 * (s1 / s2) * (rho_tilde - rho * sqrt((1 - rho_tilde^2) / (1 - rho^2)))
  curve <- data.frame(rho = rho, nie = nie)
  f <- function(r) rho_tilde - r * sqrt((1 - rho_tilde^2) / (1 - r^2))
  rho_zero <- tryCatch(stats::uniroot(f, c(-0.999, 0.999))$root, error = function(e) NA_real_)
  structure(list(curve = curve, rho_zero = rho_zero, rho_tilde = rho_tilde, nie_at_zero = nie[which.min(abs(rho))],
                 beta2 = beta2, sigma = c(total = s1, mediator = s2)), class = "cm_med_sensitivity")
}

#' @export
print.cm_med_sensitivity <- function(x, ...) {
  cat("Sensitivity of the natural indirect effect to the mediator-outcome error correlation\n")
  cat("  NIE at rho = 0: ", format(round(x$nie_at_zero, 4)), "; rho at which NIE = 0: ",
      format(round(x$rho_zero, 3)), " (residual correlation rho_tilde = ", format(round(x$rho_tilde, 3)), ")\n", sep = "")
  invisible(x)
}

#' Plot mediation effects
#'
#' @param x A `cm_mediation`, `cm_med_dml`, `cm_med_iv`, or `cm_med_cde`
#'   object, or a data frame with `term`, `estimate`, `conf.low`, `conf.high`.
#' @param terms Which terms to draw (default total, nde, nie, and the
#'   per-mediator indirect effects).
#' @return A ggplot object.
#' @export
plot_mediation <- function(x, terms = NULL) {
  eff <- if (is.data.frame(x)) x else x$effects
  if (is.null(terms)) terms <- eff$term[eff$term %in% c("total", "direct", "indirect", "nde", "nie") | grepl("^nie_", eff$term)]
  eff <- eff[eff$term %in% terms, ]
  eff$term <- factor(eff$term, levels = terms)
  ggplot2::ggplot(eff, ggplot2::aes(x = .data$term, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_col(fill = "#0072B2", alpha = 0.7, width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0.2) +
    ggplot2::labs(x = NULL, y = "Effect") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot a mediation sensitivity curve
#'
#' @param x A `cm_med_sensitivity` object.
#' @return A ggplot object of the indirect effect against the error correlation.
#' @export
plot_mediate_sensitivity <- function(x) {
  ggplot2::ggplot(x$curve, ggplot2::aes(x = .data$rho, y = .data$nie)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.9, colour = "#0072B2") +
    ggplot2::labs(x = "Correlation of mediator and outcome errors (rho)", y = "Natural indirect effect") +
    ggplot2::theme_minimal(base_size = 11)
}
