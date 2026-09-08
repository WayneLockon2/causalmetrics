# R/est-dml-iv.R
#
# The instrumental-variable models of est_dml(): the partially linear IV
# model (orthogonal partialling-out score with an instrument) and the
# interactive IV model (the LATE as a ratio of two AIPW-type scores). Called
# by est_dml() when model is "pliv" or "iivm"; returns a `cm_dml` object.

.cm_est_dml_iv <- function(data, y, d, z, x, model, z_hat, m_hat, l_hat, p_hat,
                           mu0_hat, mu1_hat, m0_hat, m1_hat, fold_id,
                           learner_l, learner_m, learner_z, learner_p,
                           learner_mu0, learner_mu1, learner_m0, learner_m1,
                           folds, n_rep, cross_fit, seed, solve, p_clip,
                           outcome_type, conf_level, na_action, weak_iv, theta_grid,
                           call) {
  if (is.null(z)) stop("`z` (the instrument column(s)) is required for model = \"", model, "\".", call. = FALSE)
  if (!is.character(z)) stop("`z` must be a character vector of instrument column names.", call. = FALSE)
  dt0 <- as.data.frame(data)
  n0 <- nrow(dt0)
  .cm_check_column(y, dt0); .cm_check_column(d, dt0)
  for (v in z) .cm_check_column(v, dt0)
  if (is.null(x)) x <- character(0)
  for (v in x) .cm_check_column(v, dt0)
  if (any(x %in% c(y, d, z))) stop("`x` must not include the outcome, treatment, or instrument columns.", call. = FALSE)
  if (model == "iivm" && length(z) != 1L) stop("The interactive IV model uses one binary instrument.", call. = FALSE)
  n_rep <- .cm_check_count(n_rep, "n_rep", min = 1L)

  y_vec <- as.numeric(dt0[[y]])
  d_vec <- if (model == "iivm") as.numeric(.cm_as_binary(dt0[[d]], "d")) else as.numeric(dt0[[d]])
  Zm <- as.matrix(dt0[, z, drop = FALSE])
  storage.mode(Zm) <- "numeric"
  if (model == "iivm") Zm[, 1] <- as.numeric(.cm_as_binary(dt0[[z]], "z"))
  d_binary <- all(d_vec[!is.na(d_vec)] %in% c(0, 1))
  z_binary <- all(Zm[!is.na(Zm)] %in% c(0, 1)) && ncol(Zm) == 1L

  # supplied nuisances
  if (model == "pliv") {
    sup <- list(l_hat = .cm_get_optional_numeric(l_hat, dt0, n0, "l_hat"),
                m_hat = .cm_get_optional_numeric(m_hat, dt0, n0, "m_hat"),
                z_hat = if (ncol(Zm) == 1L) .cm_get_optional_numeric(z_hat, dt0, n0, "z_hat") else NULL)
    if (ncol(Zm) > 1L && !is.null(z_hat)) stop("`z_hat` can be supplied only with one instrument; with several, learners estimate E[Z | X].", call. = FALSE)
    need <- c(l_hat = is.null(sup$l_hat), m_hat = is.null(sup$m_hat), z_hat = is.null(sup$z_hat))
  } else {
    sup <- list(p_hat = .cm_get_optional_numeric(p_hat, dt0, n0, "p_hat"),
                mu0_hat = .cm_get_optional_numeric(mu0_hat, dt0, n0, "mu0_hat"),
                mu1_hat = .cm_get_optional_numeric(mu1_hat, dt0, n0, "mu1_hat"),
                m0_hat = .cm_get_optional_numeric(m0_hat, dt0, n0, "m0_hat"),
                m1_hat = .cm_get_optional_numeric(m1_hat, dt0, n0, "m1_hat"))
    need <- vapply(sup, is.null, logical(1))
  }
  fid_supplied <- .cm_get_optional_vector(fold_id, dt0, n0, "fold_id")
  need_learners <- any(need)
  if (need_learners && length(x) == 0L) stop("`x` must be supplied when any nuisance prediction is estimated internally.", call. = FALSE)
  if (n_rep > 1L && (!need_learners || any(!need) || !cross_fit || !is.null(fid_supplied))) {
    stop("`n_rep > 1` requires internally estimated nuisances with cross_fit = TRUE and no supplied `fold_id`.", call. = FALSE)
  }

  keep <- is.finite(y_vec) & is.finite(d_vec) & stats::complete.cases(Zm)
  for (xj in x) keep <- keep & !is.na(dt0[[xj]])
  for (nm in names(sup)) if (!is.null(sup[[nm]])) keep <- keep & is.finite(sup[[nm]])
  if (!is.null(fid_supplied)) keep <- keep & !is.na(fid_supplied)
  n_missing <- sum(!keep)
  if (n_missing > 0L && na_action == "fail") stop("Missing or non-finite values found in required variables. Use na_action = 'omit' to drop them.", call. = FALSE)
  if (n_missing > 0L) warning(n_missing, " row(s) omitted because of missing or non-finite required values.", call. = FALSE)
  dt <- dt0[keep, , drop = FALSE]
  y_vec <- y_vec[keep]; d_vec <- d_vec[keep]; Zm <- Zm[keep, , drop = FALSE]
  for (nm in names(sup)) if (!is.null(sup[[nm]])) sup[[nm]] <- sup[[nm]][keep]
  if (!is.null(fid_supplied)) fid_supplied <- fid_supplied[keep]
  n <- length(y_vec)
  if (outcome_type == "auto") outcome_type <- if (all(y_vec %in% c(0, 1))) "binary" else "continuous"

  work <- data.frame(.cm_y = y_vec, .cm_d = if (d_binary) as.integer(d_vec) else d_vec, .cm_row_id = seq_len(n))
  for (j in seq_len(ncol(Zm))) work[[paste0(".cm_z", j)]] <- if (z_binary) as.integer(Zm[, j]) else Zm[, j]
  for (xj in x) work[[xj]] <- dt[[xj]]

  # folds
  if (!is.null(fid_supplied)) {
    fold_sets <- list(as.integer(as.factor(fid_supplied)))
  } else if (need_learners && cross_fit) {
    strata <- if (model == "iivm") as.integer(Zm[, 1]) * 2L + as.integer(d_vec) else if (d_binary) as.integer(d_vec) else NULL
    fold_sets <- .cm_make_fold_sets(n, folds, n_rep, seed, strata = strata)
  } else {
    fold_sets <- list(rep.int(1L, n))
  }

  # learners
  learners <- list()
  if (need_learners) {
    .cm_require_mlr3()
    if (model == "pliv") {
      if (need[["l_hat"]] && is.null(learner_l)) learner_l <- .cm_default_learner(if (outcome_type == "binary") "classif" else "regr")
      if (need[["m_hat"]] && is.null(learner_m)) learner_m <- .cm_default_learner(if (d_binary) "classif" else "regr")
      if (need[["z_hat"]] && is.null(learner_z)) learner_z <- .cm_default_learner(if (z_binary) "classif" else "regr")
      learners <- list(l_hat = if (need[["l_hat"]]) .cm_learner_label(learner_l) else "supplied",
                       m_hat = if (need[["m_hat"]]) .cm_learner_label(learner_m) else "supplied",
                       z_hat = if (need[["z_hat"]]) .cm_learner_label(learner_z) else "supplied")
    } else {
      if (need[["p_hat"]] && is.null(learner_p)) learner_p <- .cm_default_learner("classif")
      if (need[["mu0_hat"]] && is.null(learner_mu0)) learner_mu0 <- .cm_default_learner(if (outcome_type == "binary") "classif" else "regr")
      if (need[["mu1_hat"]] && is.null(learner_mu1)) learner_mu1 <- .cm_default_learner(if (outcome_type == "binary") "classif" else "regr")
      if (need[["m0_hat"]] && is.null(learner_m0)) learner_m0 <- .cm_default_learner("classif")
      if (need[["m1_hat"]] && is.null(learner_m1)) learner_m1 <- .cm_default_learner("classif")
      learners <- list(p_hat = if (need[["p_hat"]]) .cm_learner_label(learner_p) else "supplied",
                       mu0_hat = if (need[["mu0_hat"]]) .cm_learner_label(learner_mu0) else "supplied",
                       mu1_hat = if (need[["mu1_hat"]]) .cm_learner_label(learner_mu1) else "supplied",
                       m0_hat = if (need[["m0_hat"]]) .cm_learner_label(learner_m0) else "supplied",
                       m1_hat = if (need[["m1_hat"]]) .cm_learner_label(learner_m1) else "supplied")
    }
  } else {
    learners <- as.list(stats::setNames(rep("supplied", length(need)), names(need)))
  }

  reps <- lapply(fold_sets, function(fid) {
    .cm_dml_iv_fit_rep(work, x, fid, model, sup, need, cross_fit, solve,
                       learner_l, learner_m, learner_z, learner_p, learner_mu0, learner_mu1,
                       learner_m0, learner_m1, outcome_type, d_binary, z_binary, p_clip, ncol(Zm))
  })
  repetitions <- data.frame(rep = seq_along(reps),
                            estimate = vapply(reps, function(r) r$estimate, numeric(1)),
                            std.error = vapply(reps, function(r) r$std.error, numeric(1)),
                            n = vapply(reps, function(r) r$n, integer(1)))
  if (length(reps) == 1L) {
    estimate <- repetitions$estimate[1]; std_error <- repetitions$std.error[1]; rep_selected <- 1L
  } else {
    agg <- .cm_aggregate_reps(repetitions$estimate, repetitions$std.error)
    estimate <- agg$estimate; std_error <- agg$std.error
    rep_selected <- which.min(abs(repetitions$estimate - estimate))
  }
  sel <- reps[[rep_selected]]
  fold_estimates <- do.call(rbind, lapply(seq_along(reps), function(r) cbind(rep = r, reps[[r]]$fold_estimates)))
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)

  weak <- NULL
  if (isTRUE(weak_iv)) {
    weak <- .cm_dml_iv_weak(sel, model, theta_grid, conf_level, estimate, std_error)
  }

  diagnostics <- list(
    call = list(model = model, estimand = if (model == "pliv") "theta" else "LATE", solve = solve,
                cross_fit = cross_fit, folds = length(unique(sel$fold_id)), n_rep = length(reps),
                aggregation = if (length(reps) > 1L) "median" else "none"),
    sample = list(n_before_missing = n0, n_after_missing = n, n_analysis = sel$n, omitted_missing = n_missing,
                  n_treated = if (d_binary) sum(sel$d == 1) else NA_integer_,
                  n_control = if (d_binary) sum(sel$d == 0) else NA_integer_),
    nuisance = sel$nuisance_quality,
    prediction_mode = if (need_learners && cross_fit) "out_of_fold" else if (need_learners) "full_sample" else "supplied",
    learners = learners, fold_summary = sel$fold_summary, repetitions = repetitions,
    fold_estimates = fold_estimates, first_stage = sel$first_stage, weak_iv = weak
  )
  if (model == "iivm") diagnostics$compliance <- sel$compliance

  out <- list(
    estimate = unname(estimate), std.error = unname(std_error),
    conf.low = unname(estimate - zq * std_error), conf.high = unname(estimate + zq * std_error),
    conf.level = conf_level, model = model, estimand = if (model == "pliv") "theta" else "LATE",
    score_type = if (model == "pliv") "partialling_out_IV" else "AIPW_ratio_LATE",
    solve = solve, outcome = y, treatment = d, instrument = z, n = sel$n,
    n_treated = if (d_binary) sum(sel$d == 1) else NA_integer_,
    n_control = if (d_binary) sum(sel$d == 0) else NA_integer_,
    folds = length(unique(sel$fold_id)), n_rep = length(reps), rep_selected = rep_selected,
    repetitions = repetitions, fold_estimates = fold_estimates,
    score = sel$score, score_components = sel$score_components, residuals = sel$residuals,
    nuisance = sel$nuisance, fold_id = sel$fold_id, learners = learners,
    first_stage = sel$first_stage, weak_iv = weak, compliance = sel$compliance,
    data = dt, covariates = x, diagnostics = diagnostics, call = call
  )
  class(out) <- "cm_dml"
  out
}

.cm_dml_iv_fit_rep <- function(work, x, fid, model, sup, need, cross_fit, solve,
                               learner_l, learner_m, learner_z, learner_p, learner_mu0, learner_mu1,
                               learner_m0, learner_m1, outcome_type, d_binary, z_binary, p_clip, k) {
  y <- work$.cm_y
  d <- as.numeric(work$.cm_d)
  n <- length(y)
  positive_y <- if (outcome_type == "binary") "1" else NULL
  fold_summary <- NULL
  zcols <- paste0(".cm_z", seq_len(k))

  if (model == "pliv") {
    if (need[["l_hat"]]) {
      fit <- .cm_crossfit_predict(work, ".cm_y", x, learner_l, fid, cross_fit, positive = positive_y,
                                  task_hint = "outcome_l", what = "E[Y | X]")
      l_hat <- fit$pred; fold_summary <- fit$fold_summary
    } else l_hat <- sup$l_hat
    if (need[["m_hat"]]) {
      fit <- .cm_crossfit_predict(work, ".cm_d", x, learner_m, fid, cross_fit,
                                  positive = if (d_binary) "1" else NULL, task_hint = "treatment_m", what = "E[D | X]")
      m_hat <- fit$pred; if (is.null(fold_summary)) fold_summary <- fit$fold_summary
    } else m_hat <- sup$m_hat
    Zt <- matrix(NA_real_, n, k)
    if (need[["z_hat"]]) {
      for (j in seq_len(k)) {
        fit <- .cm_crossfit_predict(work, zcols[j], x, learner_z, fid, cross_fit,
                                    positive = if (z_binary) "1" else NULL, task_hint = "instrument_z", what = "E[Z | X]")
        Zt[, j] <- work[[zcols[j]]] - fit$pred
        if (is.null(fold_summary)) fold_summary <- fit$fold_summary
      }
      z_hat <- work[[zcols[1]]] - Zt[, 1]
    } else {
      z_hat <- sup$z_hat
      Zt[, 1] <- work[[zcols[1]]] - z_hat
    }
    .cm_check_finite(l_hat, "l_hat"); .cm_check_finite(m_hat, "m_hat"); .cm_check_finite(Zt, "z_hat")
    y_tilde <- y - l_hat
    d_tilde <- d - m_hat
    # instrument for the linear score: the residualized instrument (k = 1) or
    # the first-stage projection of d_tilde on the residualized instruments
    if (k == 1L) {
      w_inst <- Zt[, 1]
    } else {
      pi_hat <- stats::lm.fit(Zt, d_tilde)$coefficients
      pi_hat[is.na(pi_hat)] <- 0
      w_inst <- as.numeric(Zt %*% pi_hat)
    }
    psi_a <- d_tilde * w_inst
    psi_b <- y_tilde * w_inst
    solution <- .cm_solve_linear_score(psi_a, psi_b, fid, solve = solve)
    fs <- .cm_first_stage_strength(d_tilde, Zt)
    return(list(
      estimate = solution$estimate, std.error = solution$std.error, fold_estimates = solution$fold_estimates,
      score = solution$psi, score_components = data.frame(psi_a = psi_a, psi_b = psi_b),
      residuals = data.frame(y_tilde = y_tilde, d_tilde = d_tilde, z_tilde = Zt),
      nuisance = data.frame(l_hat = l_hat, m_hat = m_hat, z_hat = z_hat),
      fold_id = fid, d = d, n = n, fold_summary = fold_summary,
      nuisance_quality = rbind(
        data.frame(nuisance = "l_hat", target = "E[Y | X]", t(.cm_fit_quality(y, l_hat))),
        data.frame(nuisance = "m_hat", target = "E[D | X]", t(.cm_fit_quality(d, m_hat))),
        data.frame(nuisance = "z_hat", target = "E[Z | X]", t(.cm_fit_quality(work[[zcols[1]]], z_hat)))
      ),
      first_stage = list(t = fs$t, F = fs$F, F_robust = fs$F_robust, F_effective = fs$F_effective,
                         coefficients = fs$coefficients, jacobian = solution$jacobian),
      compliance = NULL
    ))
  }

  # interactive IV model -------------------------------------------------------
  zi <- as.integer(work$.cm_z1)
  di <- as.integer(d)
  get_nuis <- function(nm, target, learner, subset, positive, what) {
    if (need[[nm]]) {
      fit <- .cm_crossfit_predict(work, target, x, learner, fid, cross_fit, subset = subset,
                                  positive = positive, task_hint = nm, what = what)
      if (is.null(fold_summary)) fold_summary <<- fit$fold_summary
      fit$pred
    } else sup[[nm]]
  }
  p_raw <- get_nuis("p_hat", ".cm_z1", learner_p, NULL, "1", "the instrument propensity P(Z = 1 | X)")
  mu0 <- get_nuis("mu0_hat", ".cm_y", learner_mu0, zi == 0L, positive_y, "E[Y | Z = 0, X]")
  mu1 <- get_nuis("mu1_hat", ".cm_y", learner_mu1, zi == 1L, positive_y, "E[Y | Z = 1, X]")
  m0 <- get_nuis("m0_hat", ".cm_d", learner_m0, zi == 0L, "1", "E[D | Z = 0, X]")
  m1 <- get_nuis("m1_hat", ".cm_d", learner_m1, zi == 1L, "1", "E[D | Z = 1, X]")
  for (v in list(p_raw, mu0, mu1, m0, m1)) .cm_check_finite(v, "nuisance")
  p_bounds <- if (is.null(p_clip)) c(0, 1) else .cm_check_bounds(p_clip, "p_clip", strict = TRUE)
  p <- pmin(pmax(p_raw, p_bounds[1]), p_bounds[2])
  H <- zi / p - (1 - zi) / (1 - p)
  mu_z <- ifelse(zi == 1L, mu1, mu0)
  m_z <- ifelse(zi == 1L, m1, m0)
  psi_b <- mu1 - mu0 + H * (y - mu_z)
  psi_a <- m1 - m0 + H * (d - m_z)
  solution <- .cm_solve_linear_score(psi_a, psi_b, fid, solve = solve)
  comp_share <- mean(psi_a)
  comp_se <- stats::sd(psi_a) / sqrt(n)
  itt <- mean(psi_b)
  itt_se <- stats::sd(psi_b) / sqrt(n)
  fs <- .cm_first_stage_strength(d - m_z + (m1 - m0) * 0, matrix(zi - p, ncol = 1))
  list(
    estimate = solution$estimate, std.error = solution$std.error, fold_estimates = solution$fold_estimates,
    score = solution$psi, score_components = data.frame(psi_a = psi_a, psi_b = psi_b),
    residuals = NULL,
    nuisance = data.frame(p_hat = p, mu0_hat = mu0, mu1_hat = mu1, m0_hat = m0, m1_hat = m1),
    fold_id = fid, d = d, n = n, fold_summary = fold_summary,
    nuisance_quality = rbind(
      data.frame(nuisance = "p_hat", target = "P(Z = 1 | X)", t(.cm_fit_quality(zi, p_raw))),
      data.frame(nuisance = "mu_hat", target = "E[Y | Z, X] (own arm)", t(.cm_fit_quality(y, mu_z))),
      data.frame(nuisance = "m_hat", target = "E[D | Z, X] (own arm)", t(.cm_fit_quality(d, m_z)))
    ),
    first_stage = list(t = comp_share / comp_se, F = (comp_share / comp_se)^2, F_robust = (comp_share / comp_se)^2,
                       F_effective = (comp_share / comp_se)^2, jacobian = solution$jacobian),
    compliance = list(complier_share = comp_share, std.error = comp_se, itt = itt, itt_se = itt_se,
                      propensity = .cm_summary(p_raw), clipped_share = mean(p_raw < p_bounds[1] | p_raw > p_bounds[2]))
  )
}

# Weak-IV-robust confidence set from the cross-fitted score components:
# C(theta) = n M(theta)^2 / Var_n(psi_b - theta psi_a), inverted over a grid
# (analytic for one instrument).
.cm_dml_iv_weak <- function(sel, model, theta_grid, conf_level, estimate, std_error) {
  a <- sel$score_components$psi_a
  b <- sel$score_components$psi_b
  n <- length(a)
  if (is.null(theta_grid)) theta_grid <- seq(estimate - 10 * std_error, estimate + 10 * std_error, length.out = 2001L)
  ar <- .cm_ar_statistic(b, a, matrix(1, n, 1), theta_grid)
  crit <- stats::qchisq(conf_level, df = 1)
  set <- .cm_ar_set_analytic(ar, crit)
  list(statistic = data.frame(theta = theta_grid, statistic = ar$statistic),
       crit_val = crit, conf_level = conf_level, set = set$intervals, type = set$type,
       note = "Anderson-Rubin (C(alpha)) set on the orthogonal score; valid under weak identification.")
}
