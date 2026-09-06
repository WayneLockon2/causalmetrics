# R/cate-learner.R
#
# Meta-learners for the conditional average treatment effect. Every learner
# ends with one regression of a pseudo-label on the heterogeneity variables,
# fitted on all rows; the pseudo-labels differ by method and are always
# cross-fitted.

#' Meta-learners for the conditional average treatment effect
#'
#' Estimates the CATE `tau(x) = E[Y(1) - Y(0) | X = x]` by reducing the
#' problem to regressions solved by `mlr3` learners (Künzel et al. 2019;
#' Nie and Wager 2021; Kennedy 2023; Chernozhukov et al. 2026, ch. 15).
#' Each method produces a cross-fitted pseudo-label `L_i` (and, for the
#' R-learner, weights `W_i`) and then fits `learner_final` to `L_i` on the
#' heterogeneity variables `x_het` using all rows. The fitted final stage
#' predicts on new data.
#'
#' * `"dr"` (DR-learner): `L = Y(eta)`, the doubly robust pseudo-outcome of
#'   [dr_scores()]. Converges to the best approximation of the CATE in the
#'   final learner's class; error bounded by the product of the two nuisance
#'   errors.
#' * `"r"` (R-learner): weighted regression of `(Y - l(Z)) / (D - p(Z))` on
#'   `x_het` with weights `(D - p(Z))^2`, which minimizes the Robinson loss
#'   `E[(Y - l - tau(X)(D - p))^2]`. Never divides by the propensity; targets
#'   the overlap-weighted projection of the CATE.
#' * `"x"` (X-learner): imputed effects `Y - mu0(Z)` on the treated and
#'   `mu1(Z) - Y` on the controls, each regressed on `Z` within its arm, blended
#'   with weights `1 - p(Z)` and `p(Z)`. With `adapt = TRUE` the arm-specific
#'   stages are trained with the covariate-shift weights `(1 - p)^2 / p` and
#'   `p^2 / (1 - p)` (the domain-adapted X-learner).
#' * `"t"` (T-learner): `L = mu1(Z) - mu0(Z)` from two outcome regressions.
#' * `"s"` (S-learner): `L = g(1, Z) - g(0, Z)` from one regression of `Y`
#'   on `(D, Z)`; `interactions = TRUE` adds `D * Z` columns.
#'
#' @param data,y,d,x Data frame, outcome, treatment, and covariate names, as
#'   in [dr_scores()]. Not needed when `scores` is supplied.
#' @param x_het Character vector of heterogeneity variables (default `x`).
#' @param method One of `"dr"`, `"r"`, `"x"`, `"t"`, `"s"`.
#' @param learner `mlr3` regression learner for every first-stage outcome
#'   regression (default linear regression).
#' @param learner_p `mlr3` classification learner for the propensity score
#'   (default logistic regression).
#' @param learner_final `mlr3` regression learner for the final stage
#'   (default `learner`). The R-learner needs a learner that accepts
#'   observation weights (linear models, ranger, glmnet, rpart do).
#' @param scores Optional `cm_scores` object; its nuisances, folds, and data
#'   are reused and `data`, `y`, `d`, `x` are taken from it.
#' @param folds,seed,p_clip Cross-fitting settings passed to [dr_scores()].
#' @param interactions S-learner only: add treatment-covariate interactions.
#' @param adapt X-learner only: use the covariate-shift weights.
#'
#' @return A list of class `cm_cate` with `tau_hat` (fitted CATE on the
#'   estimation rows), `label` and `weights` (the final-stage pseudo-label and
#'   weights), `final` (the fitted final stage), `x_het`, `method`, `scores`,
#'   and `learners`. Methods: [predict()] on new data, [print()], [tidy()].
#'
#' @references
#' Künzel, S. R., Sekhon, J. S., Bickel, P. J., and Yu, B. (2019).
#' Metalearners for estimating heterogeneous treatment effects using machine
#' learning. *PNAS*, 116(10), 4156-4165.
#'
#' Nie, X. and Wager, S. (2021). Quasi-oracle estimation of heterogeneous
#' treatment effects. *Biometrika*, 108(2), 299-319.
#'
#' Kennedy, E. H. (2023). Towards optimal doubly robust estimation of
#' heterogeneous causal effects. *Electronic Journal of Statistics*, 17(2).
#'
#' Foster, D. J. and Syrgkanis, V. (2023). Orthogonal statistical learning.
#' *Annals of Statistics*, 51(3), 879-908.
#'
#' @examples
#' dat <- sim_hte(800, dgp = "smooth", seed = 1)
#' x <- paste0("x", 1:5)
#' fit <- cate_learner(dat, "y", "d", x, x_het = c("x1", "x2"), method = "dr", seed = 1)
#' fit
#' head(predict(fit, dat))
#' cor(fit$tau_hat, dat$tau_true)
#' @seealso [dr_scores()], [cate_score()], [cate_validate()]
#' @export
cate_learner <- function(data = NULL, y = NULL, d = NULL, x = NULL, x_het = x,
                         method = c("dr", "r", "x", "t", "s"),
                         learner = NULL, learner_p = NULL, learner_final = NULL,
                         scores = NULL, folds = 5L, seed = NULL, p_clip = c(0.01, 0.99),
                         interactions = FALSE, adapt = FALSE) {
  method <- match.arg(method)
  .cm_require_mlr3()
  if (is.null(learner)) learner <- .cm_default_learner("regr")
  if (is.null(learner_final)) learner_final <- learner
  if (is.null(scores)) {
    if (is.null(data) || is.null(y) || is.null(d) || is.null(x)) {
      stop("Supply `data`, `y`, `d`, and `x`, or a `scores` object.", call. = FALSE)
    }
    scores <- dr_scores(data, y, d, x, learner_p = learner_p, learner_mu = learner,
                        folds = folds, seed = seed, p_clip = p_clip)
  } else {
    .cm_check_scores(scores)
    if (scores$type != "dr") stop("`scores` must be of type \"dr\".", call. = FALSE)
    if (is.null(x)) x <- scores$x
  }
  if (is.null(x_het)) x_het <- x
  if (is.null(x_het)) stop("`x_het` is required.", call. = FALSE)
  work <- scores$data
  for (v in x_het) .cm_check_column(v, work)
  n <- scores$n
  yv <- scores$y
  dv <- scores$d
  p <- scores$nuisance$p
  mu0 <- scores$nuisance$mu0
  mu1 <- scores$nuisance$mu1
  fold_id <- scores$fold_id
  work$.cm_y <- yv
  work$.cm_d <- dv
  weights <- NULL
  stages <- list()

  if (method == "dr") {
    label <- scores$score
  } else if (method == "t") {
    label <- mu1 - mu0
  } else if (method == "r") {
    y_tilde <- scores$residuals$y_tilde
    d_tilde <- scores$residuals$d_tilde
    weights <- d_tilde^2
    label <- y_tilde / d_tilde
    label[!is.finite(label)] <- 0
  } else if (method == "s") {
    if (is.null(x)) stop("The S-learner needs `x`.", call. = FALSE)
    if (length(unique(fold_id)) < 2L) stop("The S-learner needs cross-fitting folds; supply `scores` built with folds.", call. = FALSE)
    cf <- .cm_crossfit_counterfactual(work, ".cm_y", c(".cm_d", x), ".cm_d", learner, fold_id,
                                      interactions = interactions)
    label <- cf$g1 - cf$g0
    stages <- list(g1 = cf$g1, g0 = cf$g0)
  } else if (method == "x") {
    if (is.null(x)) stop("The X-learner needs `x`.", call. = FALSE)
    if (length(unique(fold_id)) < 2L) stop("The X-learner needs cross-fitting folds; supply `scores` built with folds.", call. = FALSE)
    work$.cm_imp_t <- yv - mu0
    work$.cm_imp_c <- mu1 - yv
    w_t <- if (adapt) (1 - p)^2 / p else NULL
    w_c <- if (adapt) p^2 / (1 - p) else NULL
    delta_t <- .cm_crossfit_weighted(work, ".cm_imp_t", x, learner, fold_id,
                                     subset = dv == 1L, weights = w_t, task_hint = "x_treated")
    delta_c <- .cm_crossfit_weighted(work, ".cm_imp_c", x, learner, fold_id,
                                     subset = dv == 0L, weights = w_c, task_hint = "x_control")
    label <- delta_t * (1 - p) + delta_c * p
    stages <- list(delta_treated = delta_t, delta_control = delta_c)
  }

  work$.cm_label <- label
  final <- .cm_fit_mlr3(learner_final, work, ".cm_label", x_het, weights = weights,
                        task_hint = paste0(method, "_final"))
  tau_hat <- .cm_predict_fit(final, work)

  structure(list(
    tau_hat = tau_hat, label = label, weights = weights, final = final,
    x_het = x_het, x = x, method = method, scores = scores, stages = stages,
    learners = list(first_stage = .cm_learner_label(learner),
                    propensity = scores$learners$p,
                    final = .cm_learner_label(learner_final)),
    n = n, interactions = interactions, adapt = adapt, call = match.call()
  ), class = "cm_cate")
}

#' @export
predict.cm_cate <- function(object, newdata, ...) {
  if (identical(object$method, "ensemble")) return(.cm_predict_ensemble(object, newdata))
  if (identical(object$method, "budget")) return(.cm_predict_fit(object$final, newdata))
  .cm_predict_fit(object$final, newdata)
}

#' @export
print.cm_cate <- function(x, ...) {
  nm <- c(dr = "DR-learner", r = "R-learner", x = "X-learner", t = "T-learner",
          s = "S-learner", ensemble = "Ensemble")[x$method]
  cat("CATE model: ", nm, "\n", sep = "")
  if (x$method == "ensemble") {
    cat("  rule = ", x$ensemble$rule, ", base models: ",
        paste(sprintf("%s (%.3f)", names(x$ensemble$weights), x$ensemble$weights), collapse = ", "), "\n", sep = "")
  } else {
    cat("  heterogeneity variables: ", paste(x$x_het, collapse = ", "), "\n", sep = "")
    cat("  first stage = ", x$learners$first_stage, ", propensity = ", x$learners$propensity,
        ", final = ", x$learners$final, "\n", sep = "")
  }
  s <- .cm_summary(x$tau_hat)
  cat("  fitted CATE: mean = ", format(round(s["mean"], 3)), ", sd = ",
      format(round(stats::sd(x$tau_hat), 3)), ", range = [", format(round(s["min"], 3)), ", ",
      format(round(s["max"], 3)), "]\n", sep = "")
  invisible(x)
}
