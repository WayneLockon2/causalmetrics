# R/front-door.R
#
# Front-door adjustment: the effect of D on Y identified through a mediator M
# when the D-Y confounding is unobserved (Pearl 1995; Bellemare, Bloem, and
# Wexler 2020), with a regression estimator and, for a binary mediator, the
# efficient influence-function estimator of Fulcher et al. (2020).

#' Front-door adjustment
#'
#' Under the front-door criterion (the mediator `M` is unconfounded given
#' `D` and `X`; `D` affects `Y` only through `M`; the `M`-`Y` relation is
#' unconfounded given `D` and `X`), the mean potential outcome is
#' \deqn{E[Y(d)] = E_X\Big[\sum_m P(m \mid d, X) \sum_{d'} E[Y \mid d', m, X]\, P(d' \mid X)\Big],}
#' even when unobserved variables confound `D` and `Y`. The regression
#' estimator fits `M ~ D + X` (logistic for a binary mediator) and
#' `Y ~ D + M + X` and plugs in; for a binary mediator `method = "aipw"`
#' uses the influence function
#' \deqn{\varphi_d = \frac{f(M \mid d, X)}{f(M \mid D, X)}\,(Y - \mu(D, M, X))
#'  + \frac{1\{D = d\}}{p(d \mid X)}\,(\xi(M, X) - \eta_d(X)) + \sum_m f(m \mid d, X)\,\mu(D, m, X),}
#' with `xi(m, X) = sum_d' mu(d', m, X) p(d' | X)` and
#' `eta_d(X) = sum_m f(m | d, X) xi(m, X)`, whose nuisances are cross-fitted
#' `mlr3` learners. Inference by bootstrap (regression) or influence function
#' (aipw).
#'
#' @param data A data frame.
#' @param y,d,m Outcome, binary treatment, and mediator column names.
#' @param x Optional covariate names.
#' @param method `"regression"` (any mediator; continuous mediators require a
#'   linear outcome model in `M`) or `"aipw"` (binary mediator).
#' @param learner_y,learner_d,learner_m `mlr3` learners for `method = "aipw"`
#'   (defaults: linear/logistic regression).
#' @param folds Cross-fitting folds for `method = "aipw"`.
#' @param cluster Optional cluster column for the bootstrap.
#' @param n_boot Bootstrap replications for `method = "regression"`.
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_front_door` with `effects` (`total` through
#'   the front door, `backdoor` for the naive regression adjustment on `X`,
#'   and for the regression method the two components `d_on_m` and
#'   `m_on_y`), `potential` (`E[Y(1)]`, `E[Y(0)]`), and settings.
#'
#' @references
#' Pearl, J. (1995). Causal diagrams for empirical research. *Biometrika*,
#' 82(4), 669-688.
#'
#' Bellemare, M. F., Bloem, J. R., and Wexler, N. (2020). The paper of how:
#' estimating treatment effects using the front-door criterion. Working paper.
#'
#' Fulcher, I. R., Shpitser, I., Marealle, S., and Tchetgen Tchetgen, E. J.
#' (2020). Robust inference on population indirect causal effects: the
#' generalized front door criterion. *JRSS-B*, 82(1), 199-214.
#'
#' @examples
#' dat <- sim_mediation(2000, dgp = "front_door", seed = 1)
#' fit <- front_door(dat, "y", "d", "m", x = c("x1", "x2"), n_boot = 99, seed = 1)
#' fit
#' attr(dat, "truth")$total
#' @export
front_door <- function(data, y, d, m, x = NULL, method = c("regression", "aipw"),
                       learner_y = NULL, learner_d = NULL, learner_m = NULL, folds = 5L,
                       cluster = NULL, n_boot = 499L, conf_level = 0.95, seed = NULL) {
  method <- match.arg(method)
  data <- as.data.frame(data)
  for (v in c(y, d, m, x, cluster)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, m, x, cluster)])
  data <- data[keep, , drop = FALSE]
  data[[d]] <- .cm_as_binary(data[[d]], d)
  binary_m <- all(data[[m]] %in% c(0, 1))
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  n <- nrow(data)
  X <- if (length(x)) as.matrix(data[, x, drop = FALSE]) else matrix(0, n, 0)

  if (method == "regression") {
    stat <- function(dd) {
      Xd <- if (length(x)) as.matrix(dd[, x, drop = FALSE]) else matrix(0, nrow(dd), 0)
      f_m <- .cm_med_formula(m, c(d, x))
      f_y <- .cm_med_formula(y, c(d, m, x))
      f_d <- .cm_med_formula(d, x)
      fit_m <- if (binary_m) stats::glm(f_m, data = dd, family = stats::binomial()) else stats::lm(f_m, data = dd)
      fit_y <- stats::lm(f_y, data = dd)
      fit_d <- stats::glm(f_d, data = dd, family = stats::binomial())
      p1 <- stats::fitted(fit_d)
      pred_y <- function(dd_val, m_val) {
        nd <- dd; nd[[d]] <- dd_val; nd[[m]] <- m_val
        as.numeric(stats::predict(fit_y, newdata = nd))
      }
      pred_m <- function(dd_val) {
        nd <- dd; nd[[d]] <- dd_val
        as.numeric(stats::predict(fit_m, newdata = nd, type = "response"))
      }
      po <- function(dd_val) {
        if (binary_m) {
          pm1 <- pred_m(dd_val)
          xi1 <- pred_y(1, 1) * p1 + pred_y(0, 1) * (1 - p1)
          xi0 <- pred_y(1, 0) * p1 + pred_y(0, 0) * (1 - p1)
          mean(pm1 * xi1 + (1 - pm1) * xi0)
        } else {
          mm <- pred_m(dd_val)
          mean(pred_y(1, mm) * p1 + pred_y(0, mm) * (1 - p1))
        }
      }
      e1 <- po(1); e0 <- po(0)
      backdoor <- stats::coef(stats::lm(.cm_med_formula(y, c(d, x)), data = dd))[[d]]
      c(total = e1 - e0, backdoor = backdoor, d_on_m = stats::coef(fit_m)[[d]], m_on_y = stats::coef(fit_y)[[m]],
        `E[Y(1)]` = e1, `E[Y(0)]` = e0)
    }
    cl <- if (is.null(cluster)) NULL else data[[cluster]]
    tab <- .cm_boot(data, stat, n_boot, cluster = cl, seed = seed, conf_level = conf_level)
    effects <- tab[tab$term %in% c("total", "backdoor", "d_on_m", "m_on_y"), ]
    potential <- tab[tab$term %in% c("E[Y(1)]", "E[Y(0)]"), ]
    psi <- NULL
  } else {
    if (!binary_m) stop("`method = \"aipw\"` needs a binary mediator.", call. = FALSE)
    .cm_require_mlr3()
    yv <- as.numeric(data[[y]]); dv <- data[[d]]; mv <- as.integer(data[[m]])
    binary_y <- all(yv %in% c(0, 1))
    if (is.null(learner_y)) learner_y <- if (binary_y) .cm_default_learner("classif") else .cm_default_learner("regr")
    if (is.null(learner_d)) learner_d <- .cm_default_learner("classif")
    if (is.null(learner_m)) learner_m <- .cm_default_learner("classif")
    work <- data; work$.cm_y <- yv; work$.cm_d <- dv; work$.cm_m <- mv
    fold_id <- .cm_make_folds(dv, folds, seed)
    positive_y <- if (binary_y) "1" else NULL
    p1 <- .cm_crossfit_predict(work, ".cm_d", x, learner_d, fold_id, TRUE, positive = "1", task_hint = "p_d")$pred
    pm <- list()
    for (dd in c(0, 1)) {
      pm[[dd + 1]] <- .cm_crossfit_at(work, ".cm_m", c(d, x), learner_m, fold_id, subset = rep(TRUE, n),
                                      modify = function(df) { df[[d]] <- dd; df }, positive = "1")
    }
    mu <- list()
    for (dd in c(0, 1)) for (mm in c(0, 1)) {
      mu[[paste(dd, mm)]] <- .cm_crossfit_at(work, ".cm_y", c(d, m, x), learner_y, fold_id, subset = rep(TRUE, n),
                                             modify = function(df) { df[[d]] <- dd; df[[m]] <- mm; df }, positive = positive_y)
    }
    mu_obs <- ifelse(dv == 1, ifelse(mv == 1, mu[["1 1"]], mu[["1 0"]]), ifelse(mv == 1, mu[["0 1"]], mu[["0 0"]]))
    xi <- function(mm) mu[[paste(1, mm)]] * p1 + mu[[paste(0, mm)]] * (1 - p1)
    f_m_given <- function(dd, mm) if (mm == 1) pm[[dd + 1]] else 1 - pm[[dd + 1]]
    f_obs <- ifelse(mv == 1, ifelse(dv == 1, pm[[2]], pm[[1]]), ifelse(dv == 1, 1 - pm[[2]], 1 - pm[[1]]))
    phi <- sapply(c(1, 0), function(dd) {
      pd <- if (dd == 1) p1 else 1 - p1
      f_d <- ifelse(mv == 1, pm[[dd + 1]], 1 - pm[[dd + 1]])
      xi_obs <- ifelse(mv == 1, xi(1), xi(0))
      eta <- pm[[dd + 1]] * xi(1) + (1 - pm[[dd + 1]]) * xi(0)
      mu_d_obs <- ifelse(dv == 1, mu[["1 1"]] * pm[[dd + 1]] + mu[["1 0"]] * (1 - pm[[dd + 1]]),
                         mu[["0 1"]] * pm[[dd + 1]] + mu[["0 0"]] * (1 - pm[[dd + 1]]))
      f_d / f_obs * (yv - mu_obs) + as.numeric(dv == dd) / pd * (xi_obs - eta) + mu_d_obs
    })
    e1 <- mean(phi[, 1]); e0 <- mean(phi[, 2]); diffv <- phi[, 1] - phi[, 2]
    backdoor <- stats::coef(stats::lm(.cm_med_formula(y, c(d, x)), data = data))[[d]]
    effects <- data.frame(term = c("total", "backdoor"), estimate = c(e1 - e0, backdoor),
                          std.error = c(stats::sd(diffv) / sqrt(n), NA_real_), stringsAsFactors = FALSE)
    effects$conf.low <- effects$estimate - crit * effects$std.error
    effects$conf.high <- effects$estimate + crit * effects$std.error
    potential <- data.frame(term = c("E[Y(1)]", "E[Y(0)]"), estimate = c(e1, e0),
                            std.error = c(stats::sd(phi[, 1]), stats::sd(phi[, 2])) / sqrt(n))
    potential$conf.low <- potential$estimate - crit * potential$std.error
    potential$conf.high <- potential$estimate + crit * potential$std.error
    psi <- phi
  }
  structure(list(effects = effects, potential = potential, psi = psi, method = method, n = n,
                 spec = list(y = y, d = d, m = m, x = x, binary_m = binary_m), conf_level = conf_level,
                 call = match.call()), class = "cm_front_door")
}

#' @export
print.cm_front_door <- function(x, ...) {
  cat("Front-door adjustment (", x$method, "; mediator ", x$spec$m, if (x$spec$binary_m) ", binary", "; n = ", x$n, ")\n", sep = "")
  print(x$effects, digits = 4, row.names = FALSE)
  cat("  `backdoor` is the regression of the outcome on the treatment and the covariates, biased when the treatment-outcome confounding is unobserved.\n")
  invisible(x)
}
