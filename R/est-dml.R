#' Double/debiased machine learning for a treatment effect
#'
#' `est_dml()` implements the double/debiased machine learning (DML) recipe of
#' Chernozhukov et al. (2018): a Neyman-orthogonal score, cross-fitted
#' nuisance predictions, and influence-function inference. Two models are
#' available.
#'
#' * `model = "plr"`, the partially linear regression
#'   `Y = theta * D + g(X) + e`, estimated by the partialling-out score
#'   (Robinson 1988): residualize `Y` and `D` on `X`, then regress the
#'   residuals on each other. `D` may be continuous or binary.
#' * `model = "irm"`, the interactive regression model with a binary `D` and
#'   fully heterogeneous effects, estimated by the doubly robust (AIPW) score
#'   for the ATE or the ATT. This is [est_aipw()] with a choice of fold
#'   solution and repeated cross-fitting added.
#'
#' The function owns the causal-estimation layer only. Nuisance predictions
#' come either from `mlr3` learners fitted internally with cross-fitting, or
#' from outside (any R package, Python, a CSV) through the `*_hat` arguments,
#' which must then be out-of-fold predictions.
#'
#' @param data A data frame or `data.table`.
#' @param y Character scalar. Outcome column name.
#' @param d Character scalar. Treatment column name. Must be binary (0/1,
#'   logical, or coercible) for `model = "irm"`; numeric for `model = "plr"`.
#' @param x Character vector of pre-treatment covariate column names used when
#'   nuisance predictions are estimated internally.
#' @param model `"plr"` (partially linear regression, default) or `"irm"`
#'   (interactive regression model).
#' @param estimand For `model = "irm"`: `"ATE"` (default) or `"ATT"`. Ignored
#'   for `model = "plr"`, which estimates the coefficient `theta`.
#' @param l_hat,m_hat Optional out-of-fold predictions of `E[Y | X]` and
#'   `E[D | X]` for `model = "plr"`. Each is a numeric vector of length
#'   `nrow(data)` or a character scalar naming a column in `data`.
#' @param p_hat,mu0_hat,mu1_hat Optional out-of-fold predictions of
#'   `P(D = 1 | X)`, `E[Y | D = 0, X]`, and `E[Y | D = 1, X]` for
#'   `model = "irm"`. `mu1_hat` is not used for the ATT.
#' @param fold_id Optional fold identifier (vector of length `nrow(data)` or a
#'   column name). Use it to record the folds of externally produced
#'   nuisances or to impose a fixed partition; implies `n_rep = 1`.
#' @param learner_l,learner_m Optional `mlr3` learners for `E[Y | X]` and
#'   `E[D | X]` in the partially linear model. Defaults are `regr.lm`, or
#'   `classif.log_reg` with probabilities when the target is binary.
#' @param learner_p,learner_mu0,learner_mu1 Optional `mlr3` learners for the
#'   interactive model, as in [est_aipw()].
#' @param folds Number of cross-fitting folds (default 5).
#' @param n_rep Number of repeated cross-fitting splits. With `n_rep > 1` the
#'   reported estimate and standard error are the median aggregates of
#'   Chernozhukov et al. (2018, Section 3.4). Requires internally estimated
#'   nuisances.
#' @param cross_fit Logical. `FALSE` fits the nuisances on the full sample and
#'   evaluates the score on the same observations; useful only to demonstrate
#'   the overfitting bias that cross-fitting removes.
#' @param seed Optional integer seed for the fold partitions. Learners that use
#'   randomness (forests, boosting) additionally respond to `set.seed()` in the
#'   calling session.
#' @param solve How the orthogonal score is solved for `theta` across the
#'   folds. `"pooled"` (default; DML2 in Chernozhukov et al. 2018) stacks the
#'   scores of all folds and solves once. `"fold_average"` (DML1) solves within
#'   each fold and averages the fold estimates.
#' @param p_clip Clipping bounds for the propensity score (`model = "irm"`).
#' @param trim Optional propensity trimming bounds (`model = "irm"`); trimming
#'   changes the target population.
#' @param outcome_type `"auto"`, `"continuous"`, or `"binary"`; matters only
#'   when a classification learner predicts the outcome.
#' @param conf_level Confidence level for the Wald interval.
#' @param na_action `"fail"` or `"omit"` for rows with missing required values.
#'
#' @return A list of class `cm_dml` with the aggregated `estimate`,
#'   `std.error`, `conf.low`, `conf.high`, the per-repetition results
#'   (`repetitions`), per-fold estimates (`fold_estimates`), the score of the
#'   representative repetition (`score`), residuals (`residuals`, partially
#'   linear model), nuisance predictions (`nuisance`), `fold_id`, learner
#'   labels, diagnostics, and the call. [tidy()] and [glance()] methods make
#'   the object usable with `modelsummary`.
#'
#' @details
#' **Scores.** Both scores are linear in `theta`, `psi = psi_b - psi_a *
#' theta`, so `theta_hat = sum(psi_b) / sum(psi_a)` and the standard error is
#' `sqrt(sum(psi_hat^2) / (n - 1)) / (abs(mean(psi_a)) * sqrt(n))` with
#' `psi_hat` evaluated at `theta_hat`. For the partially linear model
#' `psi_a = (D - m(X))^2` and `psi_b = (Y - l(X)) (D - m(X))`; the standard
#' error equals the HC1 sandwich of the residual-on-residual regression. For
#' the interactive model the ATE score is the AIPW score (`psi_a = 1`) and the
#' ATT score is `psi_b = (D - (1 - D) p / (1 - p)) (Y - mu0) / mean(D)`,
#' `psi_a = D / mean(D)`.
#'
#' **Cross-fitting.** With `folds = K`, each nuisance is fitted on `K - 1`
#' folds and predicted on the held-out fold, so every observation's nuisance
#' prediction comes from a model that never saw it. The pooled solution
#' solves the stacked score once; the fold-average solution averages the
#' fold-specific solutions. Both use the same variance estimator.
#'
#' **Repeated cross-fitting.** With `n_rep = S`, the whole procedure is
#' repeated on `S` independent partitions. The reported estimate is the median
#' of the `S` estimates and the reported variance is the median of
#' `se_s^2 + (theta_s - theta_median)^2`, which accounts for the variability
#' across partitions. The `score`, `residuals`, and `nuisance` components refer
#' to the repetition whose estimate is closest to the median
#' (`rep_selected`).
#'
#' **Diagnostics.** `diagnostics$nuisance` reports the cross-fitted RMSE and
#' an R-squared-type fit for each nuisance; large RMSE relative to the
#' outcome's scale means the score is doing little debiasing.
#' `diagnostics$identification` (partially linear model) reports the mean
#' squared treatment residual `mean((D - m(X))^2)`, the denominator of the
#' estimator: when it is small relative to `var(D)`, the covariates predict
#' the treatment almost perfectly and the coefficient is weakly identified.
#' For the interactive model the propensity, weight, and effective sample size
#' summaries of [est_aipw()] are reported.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#' Newey, W., and Robins, J. (2018). Double/debiased machine learning for
#' treatment and structural parameters. *The Econometrics Journal*, 21(1),
#' C1-C68.
#'
#' Robinson, P. M. (1988). Root-N-consistent semiparametric regression.
#' *Econometrica*, 56(4), 931-954.
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' d <- 0.5 * x1 + sin(x2) + rnorm(n)
#' y <- 1.5 * d + x1^2 + exp(x2 / 2) + rnorm(n)
#' dat <- data.frame(y = y, d = d, x1 = x1, x2 = x2)
#'
#' # Partially linear model with the default linear nuisances (biased here,
#' # because g(X) is nonlinear) and with random forests.
#' if (requireNamespace("mlr3", quietly = TRUE) &&
#'     requireNamespace("mlr3learners", quietly = TRUE)) {
#'   est_dml(dat, y = "y", d = "d", x = c("x1", "x2"), seed = 1)
#'   if (requireNamespace("ranger", quietly = TRUE)) {
#'     est_dml(dat, y = "y", d = "d", x = c("x1", "x2"),
#'             learner_l = mlr3::lrn("regr.ranger", num.trees = 200),
#'             learner_m = mlr3::lrn("regr.ranger", num.trees = 200),
#'             folds = 5, n_rep = 2, seed = 1)
#'   }
#' }
#'
#' # Supplied out-of-fold nuisances (from any engine) need no learners.
#' dat$l_hat <- x1^2 + exp(x2 / 2) + 1.5 * (0.5 * x1 + sin(x2))
#' dat$m_hat <- 0.5 * x1 + sin(x2)
#' est_dml(dat, y = "y", d = "d", l_hat = "l_hat", m_hat = "m_hat")
#' @seealso [est_aipw()], [tidy.cm_dml()]
#' @export
est_dml <- function(data,
                    y,
                    d,
                    x = NULL,
                    model = c("plr", "irm"),
                    estimand = c("ATE", "ATT"),
                    l_hat = NULL,
                    m_hat = NULL,
                    p_hat = NULL,
                    mu0_hat = NULL,
                    mu1_hat = NULL,
                    fold_id = NULL,
                    learner_l = NULL,
                    learner_m = NULL,
                    learner_p = NULL,
                    learner_mu0 = NULL,
                    learner_mu1 = NULL,
                    folds = 5L,
                    n_rep = 1L,
                    cross_fit = TRUE,
                    seed = NULL,
                    solve = c("pooled", "fold_average"),
                    p_clip = c(0.01, 0.99),
                    trim = NULL,
                    outcome_type = c("auto", "continuous", "binary"),
                    conf_level = 0.95,
                    na_action = c("fail", "omit")) {
  call <- match.call()
  estimand_supplied <- !missing(estimand)
  model <- match.arg(model)
  if (!is.character(estimand) || length(estimand) == 0L) {
    stop("`estimand` must be \"ATE\" or \"ATT\".", call. = FALSE)
  }
  estimand <- match.arg(toupper(estimand), c("ATE", "ATT"))
  solve <- match.arg(solve)
  outcome_type <- match.arg(outcome_type)
  na_action <- match.arg(na_action)

  if (model == "plr" && estimand_supplied && estimand == "ATT") {
    stop("`estimand = \"ATT\"` is only available for `model = \"irm\"`. The partially linear model estimates a single coefficient theta.", call. = FALSE)
  }
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("`data` must be a data frame or data.table.", call. = FALSE)
  }
  if (!.cm_is_string(y) || !.cm_is_string(d)) {
    stop("`y` and `d` must be character scalars naming columns in `data`.", call. = FALSE)
  }
  if (!is.logical(cross_fit) || length(cross_fit) != 1L || is.na(cross_fit)) {
    stop("`cross_fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }
  n_rep <- .cm_check_count(n_rep, "n_rep", min = 1L)

  dt0 <- as.data.frame(data)
  n0 <- nrow(dt0)
  if (n0 < 2L) {
    stop("`data` must contain at least two rows.", call. = FALSE)
  }
  .cm_check_column(y, dt0)
  .cm_check_column(d, dt0)

  if (is.null(x)) {
    x <- character(0L)
  }
  if (!is.character(x)) {
    stop("`x` must be a character vector of covariate column names.", call. = FALSE)
  }
  if (length(x) > 0L) {
    if (anyDuplicated(x)) {
      stop("`x` contains duplicated column names.", call. = FALSE)
    }
    if (any(x %in% c(y, d))) {
      stop("`x` must not include the outcome or treatment column.", call. = FALSE)
    }
    reserved <- c(".cm_y", ".cm_d", ".cm_row_id")
    if (any(x %in% reserved)) {
      stop("`x` contains names reserved for internal use.", call. = FALSE)
    }
    invisible(lapply(x, .cm_check_column, data = dt0))
  }

  y_vec <- dt0[[y]]
  if (is.logical(y_vec)) y_vec <- as.integer(y_vec)
  if (!is.numeric(y_vec)) {
    stop("`y` must be numeric or logical.", call. = FALSE)
  }
  d_raw <- dt0[[d]]
  if (model == "irm") {
    d_vec <- as.numeric(.cm_as_binary(d_raw, "d"))
  } else if (is.numeric(d_raw)) {
    d_vec <- as.numeric(d_raw)
  } else {
    d_vec <- as.numeric(.cm_as_binary(d_raw, "d"))
  }
  d_binary <- all(d_vec[!is.na(d_vec)] %in% c(0, 1))

  # Supplied nuisances -------------------------------------------------------
  if (model == "plr") {
    sup <- list(
      l_hat = .cm_get_optional_numeric(l_hat, dt0, n0, "l_hat"),
      m_hat = .cm_get_optional_numeric(m_hat, dt0, n0, "m_hat")
    )
    if (!is.null(p_hat) || !is.null(mu0_hat) || !is.null(mu1_hat)) {
      stop("`p_hat`, `mu0_hat`, and `mu1_hat` belong to `model = \"irm\"`; use `l_hat` and `m_hat` for the partially linear model.", call. = FALSE)
    }
  } else {
    sup <- list(
      p_hat = .cm_get_optional_numeric(p_hat, dt0, n0, "p_hat"),
      mu0_hat = .cm_get_optional_numeric(mu0_hat, dt0, n0, "mu0_hat"),
      mu1_hat = if (estimand == "ATE") .cm_get_optional_numeric(mu1_hat, dt0, n0, "mu1_hat") else NULL
    )
    if (!is.null(l_hat) || !is.null(m_hat)) {
      stop("`l_hat` and `m_hat` belong to `model = \"plr\"`; use `p_hat`, `mu0_hat`, and `mu1_hat` for the interactive model.", call. = FALSE)
    }
  }
  fid_supplied <- .cm_get_optional_vector(fold_id, dt0, n0, "fold_id")

  need <- if (model == "plr") {
    c(l_hat = is.null(sup$l_hat), m_hat = is.null(sup$m_hat))
  } else {
    c(p_hat = is.null(sup$p_hat), mu0_hat = is.null(sup$mu0_hat),
      mu1_hat = estimand == "ATE" && is.null(sup$mu1_hat))
  }
  need_learners <- any(need)
  any_supplied <- any(!need[names(need) != "mu1_hat" | estimand == "ATE"])

  if (need_learners && length(x) == 0L) {
    stop("`x` must be supplied when any nuisance prediction is estimated internally.", call. = FALSE)
  }
  if (n_rep > 1L) {
    if (!need_learners || any_supplied) {
      stop("`n_rep > 1` repeats the cross-fitting of internally estimated nuisances and cannot be combined with supplied nuisance predictions.", call. = FALSE)
    }
    if (!cross_fit || !is.null(fid_supplied)) {
      stop("`n_rep > 1` requires `cross_fit = TRUE` and no user-supplied `fold_id`.", call. = FALSE)
    }
  }

  # Missing values -----------------------------------------------------------
  keep <- is.finite(y_vec) & is.finite(d_vec)
  for (xj in x) keep <- keep & !is.na(dt0[[xj]])
  for (nm in names(sup)) if (!is.null(sup[[nm]])) keep <- keep & is.finite(sup[[nm]])
  if (!is.null(fid_supplied)) keep <- keep & !is.na(fid_supplied)

  n_missing <- sum(!keep)
  if (n_missing > 0L && identical(na_action, "fail")) {
    stop("Missing or non-finite values found in required variables or supplied nuisance predictions. Use na_action = 'omit' to drop them.", call. = FALSE)
  }
  if (n_missing > 0L && identical(na_action, "omit")) {
    warning(n_missing, " row(s) omitted because of missing or non-finite required values.", call. = FALSE)
  }

  dt <- dt0[keep, , drop = FALSE]
  y_vec <- as.numeric(y_vec[keep])
  d_vec <- d_vec[keep]
  for (nm in names(sup)) if (!is.null(sup[[nm]])) sup[[nm]] <- sup[[nm]][keep]
  if (!is.null(fid_supplied)) fid_supplied <- fid_supplied[keep]
  n <- length(y_vec)
  if (n < 2L) {
    stop("Fewer than two complete observations remain.", call. = FALSE)
  }

  if (outcome_type == "auto") {
    y_unique <- unique(y_vec)
    outcome_type <- if (length(y_unique) <= 2L && all(y_unique %in% c(0, 1))) "binary" else "continuous"
  }
  if (model == "irm" && (sum(d_vec == 1) == 0L || sum(d_vec == 0) == 0L)) {
    stop("Both treated and control observations are required.", call. = FALSE)
  }
  if (model == "plr" && stats::var(d_vec) == 0) {
    stop("`d` has no variation.", call. = FALSE)
  }

  work <- data.frame(
    .cm_y = y_vec,
    .cm_d = if (d_binary) as.integer(d_vec) else d_vec,
    .cm_row_id = seq_len(n),
    check.names = FALSE
  )
  for (xj in x) work[[xj]] <- dt[[xj]]

  # Folds ---------------------------------------------------------------------
  if (!is.null(fid_supplied)) {
    fold_sets <- list(as.integer(as.factor(fid_supplied)))
    if (need_learners && cross_fit && length(unique(fold_sets[[1L]])) < 2L) {
      stop("`fold_id` must contain at least two folds when cross_fit = TRUE and learners are used.", call. = FALSE)
    }
  } else if (need_learners && cross_fit) {
    strata <- if (d_binary) as.integer(d_vec) else NULL
    fold_sets <- .cm_make_fold_sets(n, folds, n_rep, seed, strata = strata)
  } else {
    fold_sets <- list(rep.int(1L, n))
  }

  # Learners ------------------------------------------------------------------
  learners <- list()
  if (need_learners) {
    .cm_require_mlr3()
    if (model == "plr") {
      if (need[["l_hat"]] && is.null(learner_l)) {
        learner_l <- .cm_default_learner(if (outcome_type == "binary") "classif" else "regr")
      }
      if (need[["m_hat"]] && is.null(learner_m)) {
        learner_m <- .cm_default_learner(if (d_binary) "classif" else "regr")
      }
      learners <- list(
        l_hat = if (need[["l_hat"]]) .cm_learner_label(learner_l) else "supplied",
        m_hat = if (need[["m_hat"]]) .cm_learner_label(learner_m) else "supplied"
      )
    } else {
      if (need[["p_hat"]] && is.null(learner_p)) learner_p <- .cm_default_learner("classif")
      if (need[["mu0_hat"]] && is.null(learner_mu0)) learner_mu0 <- .cm_default_learner("regr")
      if (need[["mu1_hat"]] && is.null(learner_mu1)) learner_mu1 <- .cm_default_learner("regr")
      learners <- list(
        p_hat = if (need[["p_hat"]]) .cm_learner_label(learner_p) else "supplied",
        mu0_hat = if (need[["mu0_hat"]]) .cm_learner_label(learner_mu0) else "supplied",
        mu1_hat = if (estimand == "ATT") "not needed for ATT" else if (need[["mu1_hat"]]) .cm_learner_label(learner_mu1) else "supplied"
      )
    }
  } else {
    learners <- if (model == "plr") {
      list(l_hat = "supplied", m_hat = "supplied")
    } else {
      list(p_hat = "supplied", mu0_hat = "supplied",
           mu1_hat = if (estimand == "ATT") "not needed for ATT" else "supplied")
    }
  }

  # Repetitions ----------------------------------------------------------------
  reps <- vector("list", length(fold_sets))
  for (r in seq_along(fold_sets)) {
    reps[[r]] <- .cm_dml_fit_rep(
      work = work, x = x, fid = fold_sets[[r]], model = model, estimand = estimand,
      sup = sup, need = need, cross_fit = cross_fit, solve = solve,
      learner_l = learner_l, learner_m = learner_m,
      learner_p = learner_p, learner_mu0 = learner_mu0, learner_mu1 = learner_mu1,
      outcome_type = outcome_type, d_binary = d_binary, p_clip = p_clip, trim = trim
    )
  }

  repetitions <- data.frame(
    rep = seq_along(reps),
    estimate = vapply(reps, function(r) r$estimate, numeric(1)),
    std.error = vapply(reps, function(r) r$std.error, numeric(1)),
    n = vapply(reps, function(r) r$n, integer(1))
  )
  if (length(reps) == 1L) {
    estimate <- repetitions$estimate[1L]
    std_error <- repetitions$std.error[1L]
    rep_selected <- 1L
  } else {
    agg <- .cm_aggregate_reps(repetitions$estimate, repetitions$std.error)
    estimate <- agg$estimate
    std_error <- agg$std.error
    rep_selected <- which.min(abs(repetitions$estimate - estimate))
  }
  sel <- reps[[rep_selected]]

  fold_estimates <- do.call(rbind, lapply(seq_along(reps), function(r) {
    cbind(rep = r, reps[[r]]$fold_estimates)
  }))

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  diagnostics <- list(
    call = list(model = model, estimand = if (model == "plr") "theta" else estimand,
                solve = solve, cross_fit = cross_fit, folds = length(unique(sel$fold_id)),
                n_rep = length(reps), aggregation = if (length(reps) > 1L) "median" else "none"),
    sample = list(
      n_before_missing = n0,
      n_after_missing = n,
      n_analysis = sel$n,
      omitted_missing = n_missing,
      n_treated = if (d_binary) sum(sel$d == 1) else NA_integer_,
      n_control = if (d_binary) sum(sel$d == 0) else NA_integer_
    ),
    nuisance = sel$nuisance_quality,
    prediction_mode = if (need_learners && cross_fit) "out_of_fold" else if (need_learners) "full_sample" else "supplied",
    learners = learners,
    fold_summary = sel$fold_summary,
    repetitions = repetitions,
    fold_estimates = fold_estimates
  )
  if (model == "plr") {
    diagnostics$identification <- sel$identification
  } else {
    diagnostics$propensity <- sel$propensity
    diagnostics$trimming <- sel$trimming
    diagnostics$weights <- sel$weights
  }

  out <- list(
    estimate = unname(estimate),
    std.error = unname(std_error),
    conf.low = unname(estimate - z * std_error),
    conf.high = unname(estimate + z * std_error),
    conf.level = conf_level,
    model = model,
    estimand = if (model == "plr") "theta" else estimand,
    score_type = if (model == "plr") "partialling_out" else paste0("AIPW_", estimand),
    solve = solve,
    outcome = y,
    treatment = d,
    n = sel$n,
    n_treated = if (d_binary) sum(sel$d == 1) else NA_integer_,
    n_control = if (d_binary) sum(sel$d == 0) else NA_integer_,
    folds = length(unique(sel$fold_id)),
    n_rep = length(reps),
    rep_selected = rep_selected,
    repetitions = repetitions,
    fold_estimates = fold_estimates,
    score = sel$score,
    residuals = sel$residuals,
    nuisance = sel$nuisance,
    fold_id = sel$fold_id,
    learners = learners,
    diagnostics = diagnostics,
    call = call
  )
  class(out) <- "cm_dml"
  out
}

# One complete DML fit on a given fold partition: nuisances, score, solution,
# and diagnostics. Returns a list consumed by est_dml().
.cm_dml_fit_rep <- function(work, x, fid, model, estimand, sup, need, cross_fit, solve,
                            learner_l, learner_m, learner_p, learner_mu0, learner_mu1,
                            outcome_type, d_binary, p_clip, trim) {
  y <- work$.cm_y
  d <- as.numeric(work$.cm_d)
  n <- length(y)
  positive_y <- if (outcome_type == "binary") "1" else NULL
  fold_summary <- NULL

  if (model == "plr") {
    if (need[["l_hat"]]) {
      fit <- .cm_crossfit_predict(work, ".cm_y", x, learner_l, fid, cross_fit,
                                  positive = positive_y, task_hint = "outcome_l",
                                  what = "E[Y | X]")
      l_hat <- fit$pred
      fold_summary <- fit$fold_summary
    } else {
      l_hat <- sup$l_hat
    }
    if (need[["m_hat"]]) {
      fit <- .cm_crossfit_predict(work, ".cm_d", x, learner_m, fid, cross_fit,
                                  positive = if (d_binary) "1" else NULL,
                                  task_hint = "treatment_m", what = "E[D | X]")
      m_hat <- fit$pred
      if (is.null(fold_summary)) fold_summary <- fit$fold_summary
    } else {
      m_hat <- sup$m_hat
    }
    .cm_check_finite(l_hat, "l_hat")
    .cm_check_finite(m_hat, "m_hat")

    y_tilde <- y - l_hat
    d_tilde <- d - m_hat
    psi <- .cm_score_plr(y_tilde, d_tilde)
    solution <- .cm_solve_linear_score(psi$a, psi$b, fid, solve = solve)

    var_d <- stats::var(d)
    return(list(
      estimate = solution$estimate,
      std.error = solution$std.error,
      fold_estimates = solution$fold_estimates,
      score = solution$psi,
      residuals = data.frame(y_tilde = y_tilde, d_tilde = d_tilde),
      nuisance = data.frame(l_hat = l_hat, m_hat = m_hat),
      fold_id = fid,
      d = d,
      n = n,
      fold_summary = fold_summary,
      nuisance_quality = rbind(
        data.frame(nuisance = "l_hat", target = "E[Y | X]", t(.cm_fit_quality(y, l_hat))),
        data.frame(nuisance = "m_hat", target = "E[D | X]", t(.cm_fit_quality(d, m_hat)))
      ),
      identification = list(
        mean_sq_treatment_residual = mean(d_tilde^2),
        share_treatment_variance_unexplained = if (is.finite(var_d) && var_d > 0) mean(d_tilde^2) / var_d else NA_real_,
        jacobian = solution$jacobian
      )
    ))
  }

  # Interactive regression model ------------------------------------------------
  d_int <- as.integer(d)
  if (need[["p_hat"]]) {
    fit <- .cm_crossfit_predict(work, ".cm_d", x, learner_p, fid, cross_fit,
                                positive = "1", task_hint = "propensity",
                                what = "the propensity score")
    p_raw <- fit$pred
    fold_summary <- fit$fold_summary
  } else {
    p_raw <- sup$p_hat
  }
  if (need[["mu0_hat"]]) {
    fit <- .cm_crossfit_predict(work, ".cm_y", x, learner_mu0, fid, cross_fit,
                                subset = d_int == 0L, positive = positive_y,
                                task_hint = "outcome",
                                what = "mu0 estimation (no controls in a training fold)")
    mu0 <- fit$pred
    if (is.null(fold_summary)) fold_summary <- fit$fold_summary
  } else {
    mu0 <- sup$mu0_hat
  }
  if (estimand == "ATE") {
    if (need[["mu1_hat"]]) {
      fit <- .cm_crossfit_predict(work, ".cm_y", x, learner_mu1, fid, cross_fit,
                                  subset = d_int == 1L, positive = positive_y,
                                  task_hint = "outcome",
                                  what = "mu1 estimation (no treated observations in a training fold)")
      mu1 <- fit$pred
      if (is.null(fold_summary)) fold_summary <- fit$fold_summary
    } else {
      mu1 <- sup$mu1_hat
    }
  } else {
    mu1 <- NULL
  }

  .cm_check_finite(p_raw, "p_hat")
  .cm_check_finite(mu0, "mu0_hat")
  if (!is.null(mu1)) .cm_check_finite(mu1, "mu1_hat")
  if (any(p_raw < 0 | p_raw > 1)) {
    stop("`p_hat` must be between 0 and 1 before clipping.", call. = FALSE)
  }

  if (!is.null(trim)) {
    trim <- .cm_check_bounds(trim, "trim", strict = FALSE)
    trim_keep <- p_raw >= trim[1L] & p_raw <= trim[2L]
  } else {
    trim_keep <- rep.int(TRUE, n)
  }
  n_trimmed <- sum(!trim_keep)
  if (n_trimmed > 0L) {
    y <- y[trim_keep]
    d <- d[trim_keep]
    d_int <- d_int[trim_keep]
    p_raw <- p_raw[trim_keep]
    mu0 <- mu0[trim_keep]
    if (!is.null(mu1)) mu1 <- mu1[trim_keep]
    fid <- fid[trim_keep]
  }
  if (length(y) < 2L || sum(d_int == 1L) == 0L || sum(d_int == 0L) == 0L) {
    stop("Trimming removed too many observations; both treatment groups must remain.", call. = FALSE)
  }

  p_bounds <- if (is.null(p_clip)) c(0, 1) else .cm_check_bounds(p_clip, "p_clip", strict = TRUE)
  p <- pmin(pmax(p_raw, p_bounds[1L]), p_bounds[2L])
  if (any(p <= 0 | p >= 1)) {
    stop("Propensity scores must be strictly between 0 and 1 after clipping. Use nonzero clipping bounds.", call. = FALSE)
  }

  psi <- .cm_score_irm(y, d, p, mu0, mu1, estimand)
  solution <- .cm_solve_linear_score(psi$a, psi$b, fid, solve = solve)

  if (estimand == "ATE") {
    w_treated <- d / p
    w_control <- (1 - d) / (1 - p)
    g_obs <- d * mu1 + (1 - d) * mu0
    quality <- rbind(
      data.frame(nuisance = "p_hat", target = "P(D = 1 | X)", t(.cm_fit_quality(d, p_raw))),
      data.frame(nuisance = "mu_hat", target = "E[Y | D, X] (own arm)", t(.cm_fit_quality(y, g_obs)))
    )
  } else {
    w_treated <- d
    w_control <- (1 - d) * p / (1 - p)
    quality <- rbind(
      data.frame(nuisance = "p_hat", target = "P(D = 1 | X)", t(.cm_fit_quality(d, p_raw))),
      data.frame(nuisance = "mu0_hat", target = "E[Y | D = 0, X] (controls)",
                 t(.cm_fit_quality(y[d_int == 0L], mu0[d_int == 0L])))
    )
  }

  nuisance <- data.frame(p_hat = p, mu0_hat = mu0)
  if (!is.null(mu1)) nuisance$mu1_hat <- mu1

  list(
    estimate = solution$estimate,
    std.error = solution$std.error,
    fold_estimates = solution$fold_estimates,
    score = solution$psi,
    residuals = NULL,
    nuisance = nuisance,
    fold_id = fid,
    d = d,
    n = length(y),
    fold_summary = fold_summary,
    nuisance_quality = quality,
    propensity = list(
      raw_summary = .cm_summary(p_raw),
      clipped_summary = .cm_summary(p),
      p_clip = p_bounds,
      n_clipped_low = sum(p_raw < p_bounds[1L]),
      n_clipped_high = sum(p_raw > p_bounds[2L]),
      common_support = .cm_common_support(p, d_int)
    ),
    trimming = list(trim = trim, n_before_trim = n, n_trimmed = n_trimmed, n_after_trim = length(y)),
    weights = list(
      ess_treated = .cm_ess(w_treated[d_int == 1L]),
      ess_control = .cm_ess(w_control[d_int == 0L]),
      treated_weight_summary = .cm_summary(w_treated[d_int == 1L]),
      control_weight_summary = .cm_summary(w_control[d_int == 0L])
    )
  )
}

#' Print method for `cm_dml` objects
#'
#' @param x A `cm_dml` object returned by [est_dml()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.cm_dml <- function(x, ...) {
  label <- if (x$model == "plr") "partially linear regression, theta" else paste0("interactive regression model, ", x$estimand)
  cat("Double ML estimate (", label, ")\n", sep = "")
  solution <- if (x$solve == "pooled") {
    paste0("pooled across ", x$folds, " fold(s)")
  } else {
    paste0("average of ", x$folds, " fold estimate(s)")
  }
  cat("  Score: ", x$score_type, "; solution: ", solution, "; ",
      x$n_rep, " repetition(s)",
      if (x$n_rep > 1L) " (median aggregation)" else "", "\n", sep = "")
  cat("  Estimate:   ", formatC(x$estimate, digits = 4, format = "f"), "\n", sep = "")
  cat("  Std. Error: ", formatC(x$std.error, digits = 4, format = "f"), "\n", sep = "")
  ci_label <- round(100 * x$conf.level, 1)
  cat("  ", ci_label, "% CI:     [",
      formatC(x$conf.low, digits = 4, format = "f"), ", ",
      formatC(x$conf.high, digits = 4, format = "f"), "]\n", sep = "")
  if (!is.na(x$n_treated)) {
    cat("  N:          ", x$n, " (treated: ", x$n_treated, ", control: ", x$n_control, ")\n", sep = "")
  } else {
    cat("  N:          ", x$n, "\n", sep = "")
  }
  q <- x$diagnostics$nuisance
  if (!is.null(q)) {
    cat("  Nuisance RMSE: ",
        paste(sprintf("%s = %s", q$nuisance, formatC(q$rmse, digits = 3, format = "g")), collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}
