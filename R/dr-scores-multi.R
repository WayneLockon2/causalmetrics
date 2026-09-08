# R/dr-scores-multi.R
#
# Doubly robust scores with more than two treatment arms: an n x K matrix
# Gamma[i, k] = mu_k(x_i) + 1{W_i = k} (Y_i - mu_k(x_i)) / e_k(x_i), one
# column per arm, as used by multi-action policy learning (policytree) and
# off-policy evaluation. Built by dr_scores() when `d` has more than two
# levels.

.cm_dr_scores_multi <- function(data, y, d, x, arms, p_hat, mu_hat, learner_p, learner_mu,
                                type, folds, fold_id, seed, p_clip, trim, outcome_type, na_action) {
  w_raw <- data[[d]]
  if (is.null(arms)) arms <- if (is.factor(w_raw)) levels(droplevels(w_raw)) else sort(unique(w_raw))
  arms <- as.character(arms)
  K <- length(arms)
  if (K < 2L) stop("`d` must have at least two arms.", call. = FALSE)
  n_all <- nrow(data)
  # supplied nuisances: matrices or K column names
  get_mat <- function(arg, nm) {
    if (is.null(arg)) return(NULL)
    if (is.character(arg)) {
      if (length(arg) != K) stop("`", nm, "` must name ", K, " columns, one per arm.", call. = FALSE)
      for (v in arg) .cm_check_column(v, data)
      M <- as.matrix(data[, arg])
    } else {
      M <- as.matrix(arg)
      if (nrow(M) != n_all || ncol(M) != K) stop("`", nm, "` must be an n x K matrix.", call. = FALSE)
    }
    storage.mode(M) <- "double"
    colnames(M) <- arms
    M
  }
  e_in <- get_mat(p_hat, "p_hat")
  mu_in <- get_mat(mu_hat, "mu_hat")
  fold_in <- .cm_get_optional_vector(fold_id, data, n_all, "fold_id")

  keep <- stats::complete.cases(data[, c(y, d, x), drop = FALSE]) & !is.na(match(as.character(w_raw), arms))
  if (!is.null(e_in)) keep <- keep & stats::complete.cases(e_in)
  if (!is.null(mu_in)) keep <- keep & stats::complete.cases(mu_in)
  if (!all(keep)) {
    if (na_action == "fail") stop("Missing values in the required columns (or arm labels outside `arms`); use na_action = \"omit\".", call. = FALSE)
    data <- data[keep, , drop = FALSE]; w_raw <- w_raw[keep]
    e_in <- e_in[keep, , drop = FALSE]; mu_in <- mu_in[keep, , drop = FALSE]; fold_in <- fold_in[keep]
  }
  n <- nrow(data)
  yv <- as.numeric(data[[y]])
  w_index <- match(as.character(w_raw), arms)
  if (any(tabulate(w_index, K) < 2L)) stop("Every arm needs at least two observations.", call. = FALSE)

  need_e <- is.null(e_in) && type != "reg"
  need_mu <- is.null(mu_in) && type != "ipw"
  if ((need_e || need_mu) && is.null(x)) stop("`x` is required when nuisances are estimated internally.", call. = FALSE)
  if (!is.null(fold_in)) {
    fold_id <- as.integer(factor(fold_in))
  } else if (need_e || need_mu) {
    fold_id <- .cm_make_fold_sets(n, folds, 1L, seed, strata = w_index)[[1]]
  } else {
    fold_id <- rep(1L, n)
  }
  if (outcome_type == "auto") outcome_type <- if (all(yv %in% c(0, 1))) "binary" else "continuous"
  labels <- list(p = "supplied", mu = "supplied")
  work <- data
  work$.cm_y <- yv
  if (need_e || need_mu) .cm_require_mlr3()
  if (need_e) {
    if (is.null(learner_p)) learner_p <- .cm_default_learner("classif")
    labels$p <- paste0(.cm_learner_label(learner_p), " (one-vs-rest, normalized)")
    e_in <- sapply(seq_len(K), function(k) {
      work$.cm_wk <- as.integer(w_index == k)
      .cm_crossfit_predict(work, ".cm_wk", x, learner_p, fold_id, cross_fit = TRUE,
                           positive = "1", task_hint = "propensity", what = paste("the propensity of arm", arms[k]))$pred
    })
    e_in <- matrix(e_in, nrow = n)
    e_in <- e_in / rowSums(e_in)
    colnames(e_in) <- arms
  }
  if (need_mu) {
    if (is.null(learner_mu)) learner_mu <- .cm_default_learner("regr")
    labels$mu <- .cm_learner_label(learner_mu)
    positive_y <- if (outcome_type == "binary") "1" else NULL
    mu_in <- sapply(seq_len(K), function(k) {
      .cm_crossfit_predict(work, ".cm_y", x, learner_mu, fold_id, cross_fit = TRUE,
                           subset = w_index == k, positive = positive_y, task_hint = "outcome",
                           what = paste("the outcome model of arm", arms[k]))$pred
    })
    mu_in <- matrix(mu_in, nrow = n)
    colnames(mu_in) <- arms
  }
  if (type == "reg") e_in <- matrix(rep(tabulate(w_index, K) / n, each = n), n, dimnames = list(NULL, arms))
  if (type == "ipw") mu_in <- matrix(0, n, K, dimnames = list(NULL, arms))
  trimmed_share <- 0
  if (!is.null(trim)) {
    if (length(trim) != 2L || trim[1] >= trim[2]) stop("`trim` must be two increasing bounds.", call. = FALSE)
    e_obs <- e_in[cbind(seq_len(n), w_index)]
    keep_t <- e_obs > trim[1] & e_obs < trim[2]
    trimmed_share <- mean(!keep_t)
    if (sum(keep_t) < 4L) stop("Trimming removed almost every observation.", call. = FALSE)
    data <- data[keep_t, , drop = FALSE]; yv <- yv[keep_t]; w_index <- w_index[keep_t]; fold_id <- fold_id[keep_t]
    e_in <- e_in[keep_t, , drop = FALSE]; mu_in <- mu_in[keep_t, , drop = FALSE]
    n <- length(yv)
  }
  e_raw <- e_in
  e <- .cm_clip(e_in, p_clip[1], p_clip[2])
  I <- matrix(0, n, K)
  I[cbind(seq_len(n), w_index)] <- 1
  gamma <- switch(type,
    dr = mu_in + I * (yv - mu_in) / e,
    ipw = I * yv / e,
    reg = mu_in
  )
  colnames(gamma) <- arms
  crit <- stats::qnorm(0.975)
  arm_tab <- data.frame(arm = arms, n = tabulate(w_index, K), estimate = colMeans(gamma),
                        std.error = apply(gamma, 2, stats::sd) / sqrt(n), stringsAsFactors = FALSE)
  contrasts <- do.call(rbind, lapply(seq_len(K)[-1], function(k) {
    dlt <- gamma[, k] - gamma[, 1]
    data.frame(contrast = paste0(arms[k], " - ", arms[1]), estimate = mean(dlt),
               std.error = stats::sd(dlt) / sqrt(n), stringsAsFactors = FALSE)
  }))
  nuis_fit <- t(sapply(seq_len(K), function(k) {
    i <- w_index == k
    c(e_rmse = if (type != "reg") .cm_fit_quality(as.numeric(I[, k]), e_raw[, k])[["rmse"]] else NA_real_,
      mu_rmse = if (type != "ipw") .cm_fit_quality(yv[i], mu_in[i, k])[["rmse"]] else NA_real_)
  }))
  rownames(nuis_fit) <- arms
  structure(list(
    gamma = gamma, nuisance = list(e = e, mu = mu_in), y = yv, w = arms[w_index], w_index = w_index,
    arms = arms, x = x, y_name = y, d_name = d, data = data, fold_id = fold_id, n = n, type = type,
    ate = list(arms = arm_tab, contrasts = contrasts), learners = labels, outcome_type = outcome_type,
    p_clip = p_clip, trim = trim,
    diagnostics = list(propensity = apply(e_raw, 2, .cm_summary), trimmed_share = trimmed_share,
                       nuisance = nuis_fit, score = apply(gamma, 2, .cm_summary)),
    call = match.call()
  ), class = "cm_scores_multi")
}

#' @export
print.cm_scores_multi <- function(x, ...) {
  cat("Doubly robust scores (", x$type, ") for ", length(x$arms), " arms: ", paste(x$arms, collapse = ", "), "\n", sep = "")
  cat("  n = ", x$n, ", folds = ", length(unique(x$fold_id)),
      if (!is.null(x$trim)) paste0(", trimmed share = ", round(x$diagnostics$trimmed_share, 3)), "\n", sep = "")
  cat("  nuisances: e = ", x$learners$p, ", mu = ", x$learners$mu, "\n", sep = "")
  cat("  expected outcome by arm (mean of the score column):\n")
  print(x$ate$arms, digits = 4, row.names = FALSE)
  cat("  contrasts against the first arm:\n")
  print(x$ate$contrasts, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Pairwise scores from a multi-arm score object
#'
#' Builds the binary `cm_scores` object of one arm against another from the
#' nuisances stored in a multi-arm object: the rows of the two arms, the
#' pairwise propensity `e_t / (e_t + e_c)`, the two outcome models, and the
#' doubly robust pseudo-outcome of [dr_scores()]. The result works with
#' every binary function ([cate_learner()], [cate_blp()],
#' [cate_validate()], [policy_learn()]).
#'
#' @param x A `cm_scores_multi` object.
#' @param treat,control Arm labels.
#' @return A `cm_scores` object on the rows of the two arms, with `d = 1`
#'   for `treat`.
#' @examples
#' dat <- sim_hte(900, dgp = "smooth", seed = 1)
#' dat$arm <- sample(c("a", "b", "c"), 900, replace = TRUE)
#' dat$y3 <- dat$y + (dat$arm == "b") * dat$tau_true
#' sc <- dr_scores(dat, "y3", "arm", paste0("x", 1:5), seed = 1)
#' contrast_scores(sc, "b", "a")
#' @export
contrast_scores <- function(x, treat, control) {
  if (!inherits(x, "cm_scores_multi")) stop("`x` must be a `cm_scores_multi` object.", call. = FALSE)
  t_i <- match(as.character(treat), x$arms); c_i <- match(as.character(control), x$arms)
  if (is.na(t_i) || is.na(c_i) || t_i == c_i) stop("`treat` and `control` must be two different arms of `x`.", call. = FALSE)
  rows <- x$w_index %in% c(t_i, c_i)
  yv <- x$y[rows]
  dv <- as.integer(x$w_index[rows] == t_i)
  e_t <- x$nuisance$e[rows, t_i]; e_c <- x$nuisance$e[rows, c_i]
  p <- .cm_clip(e_t / (e_t + e_c), x$p_clip[1], x$p_clip[2])
  mu1 <- x$nuisance$mu[rows, t_i]; mu0 <- x$nuisance$mu[rows, c_i]
  H <- dv / p - (1 - dv) / (1 - p)
  mu_d <- ifelse(dv == 1L, mu1, mu0)
  score <- switch(x$type, dr = (mu1 - mu0) + H * (yv - mu_d), ipw = H * yv, reg = mu1 - mu0)
  l <- p * mu1 + (1 - p) * mu0
  n <- sum(rows)
  data <- x$data[rows, , drop = FALSE]
  structure(list(
    score = as.numeric(score),
    nuisance = data.frame(p = p, mu0 = mu0, mu1 = mu1, l = l),
    residuals = data.frame(y_tilde = yv - l, d_tilde = dv - p),
    y = yv, d = dv, x = x$x, y_name = x$y_name, d_name = x$d_name,
    data = data, fold_id = x$fold_id[rows], n = n, type = x$type,
    ate = list(estimate = mean(score), std.error = stats::sd(score) / sqrt(n)),
    learners = list(p = x$learners$p, mu0 = x$learners$mu, mu1 = x$learners$mu),
    outcome_type = x$outcome_type, p_clip = x$p_clip, trim = x$trim,
    diagnostics = list(propensity = .cm_summary(p), clipped_share = 0, trimmed_share = 0,
                       common_support = .cm_common_support(p, dv),
                       nuisance = rbind(p = .cm_fit_quality(dv, p), mu0 = .cm_fit_quality(yv[dv == 0L], mu0[dv == 0L]),
                                        mu1 = .cm_fit_quality(yv[dv == 1L], mu1[dv == 1L])),
                       score = .cm_summary(score)),
    contrast = c(treat = as.character(treat), control = as.character(control)), call = match.call()
  ), class = "cm_scores")
}
