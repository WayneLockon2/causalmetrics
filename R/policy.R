# R/policy.R
#
# Policy evaluation and learning from doubly robust scores: the value of a
# fixed rule with a confidence interval, empirical welfare maximization over
# depth-limited trees (exact search), linear rules, or weighted classifiers,
# budgeted rules from a CATE model, and the welfare-curvature targeting
# frontier of Haushofer et al. (2025).

#' Value of a treatment policy
#'
#' The gain of the rule `pi(X)` over treating no one is
#' `V(pi) = E[pi(X) (Y(1) - Y(0) - cost)] = E[pi(X) (Y(eta) - cost)]`, an
#' average of the doubly robust pseudo-outcome, so it is estimated with a
#' standard error from the influence function. With `baseline = "all"` the
#' gain over treating everyone is `E[(pi(X) - 1)(Y(eta) - cost)]`. Several
#' policies can be evaluated at once and are compared with the first one.
#'
#' With a multi-arm score object (see [dr_scores()] with more than two
#' arms) a policy maps units to arms, its value is the mean of
#' `Gamma[i, pi(X_i)]`, the expected outcome under the policy, and
#' `baseline` names the arm of a uniform comparison policy: the gain is then
#' `E[Gamma[i, pi(X_i)] - Gamma[i, baseline]]` and the percentage gain over
#' the uniform policy is reported. With `type = "ipw"` scores this is the
#' inverse-propensity-score reward of the off-policy evaluation literature.
#'
#' @param scores A `cm_scores` or `cm_scores_multi` object (held out from
#'   the policy's estimation sample for an honest value).
#' @param policy Binary case: a 0/1 (or logical) vector, a `cm_policy`
#'   object, a function of the data returning 0/1, or a named list of such
#'   objects. Multi-arm case: a vector of arm labels, a `cm_policy` object,
#'   a function returning arm labels, or a named list of them.
#' @param baseline Binary case: `"none"` (value over treating no one) or
#'   `"all"`. Multi-arm case: `"none"` (the expected outcome under the
#'   policy) or an arm label (gain over the uniform policy that assigns that
#'   arm to everyone).
#' @param cost Cost of treating one unit, in outcome units (binary case), or
#'   a vector of per-arm costs (multi-arm case).
#' @param conf_level Confidence level.
#' @param n_boot Number of nonparametric bootstrap draws (units resampled
#'   with replacement); 0 (default) reports influence-function standard
#'   errors only.
#' @param seed Optional seed for the bootstrap.
#' @param weights Optional positive unit weights (a numeric vector of length
#'   `nrow(scores$data)` or a column name): survey weights, or the
#'   time-decaying weights `h_t` of policy learning with adaptively collected
#'   data (Zhan et al. 2024). The value is `sum(w_i c_i) / sum(w_i)` with
#'   `c_i` the unit contribution and its standard error
#'   `sqrt(sum(w_i^2 (c_i - value)^2)) / sum(w_i)`, the adaptively weighted
#'   form of Hadad et al. (2021), which reduces to the unweighted formula for
#'   constant weights. Shares treated are unweighted.
#'
#' @return A data frame with one row per policy: `policy`, the share
#'   treated (binary) or the share assigned to each arm (multi-arm),
#'   `estimate`, `std.error`, `conf.low`, `conf.high`, and, when several
#'   policies are given, `diff_vs_first` and `diff_se`. In the multi-arm case
#'   with an arm baseline, `baseline_value` and `gain_pct` are added. With
#'   `n_boot > 0`, `boot.se`, `boot.low`, `boot.high` (percentile interval)
#'   are added and the draws are attached as attributes `"boot"` (values) and
#'   `"boot_diff"` (differences against the first policy).
#' @examples
#' dat <- sim_hte(1000, dgp = "policy", seed = 1)
#' sc <- dr_scores(dat, "y", "d", paste0("x", 1:5), seed = 1)
#' policy_value(sc, list(oracle = dat$tau_true > 0, all = rep(1, 1000)))
#' @export
policy_value <- function(scores, policy, baseline = "none", cost = 0, conf_level = 0.95,
                         n_boot = 0L, seed = NULL, weights = NULL) {
  if (inherits(scores, "cm_scores_multi")) {
    return(.cm_policy_value_multi(scores, policy, baseline, cost, conf_level, n_boot, seed, weights))
  }
  baseline <- match.arg(baseline, c("none", "all"))
  .cm_check_scores(scores)
  pols <- if (is.list(policy) && !inherits(policy, "cm_policy")) policy else list(policy = policy)
  if (is.null(names(pols)) || any(names(pols) == "")) names(pols) <- paste0("policy", seq_along(pols))
  g <- scores$score - cost
  n <- scores$n
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  A <- vapply(names(pols), function(nm) .cm_policy_assignment(pols[[nm]], scores$data, nm), numeric(n))
  A <- matrix(A, nrow = n)
  contrib <- if (baseline == "none") A * g else (A - 1) * g
  w <- .cm_policy_weights(weights, scores)
  wv <- .cm_weighted_value(contrib, w)
  est <- wv$est
  se <- wv$se
  out <- data.frame(policy = names(pols), share_treated = colMeans(A), estimate = est, std.error = se,
                    conf.low = est - crit * se, conf.high = est + crit * se, stringsAsFactors = FALSE)
  if (length(pols) > 1L) {
    d <- .cm_weighted_value(contrib - contrib[, 1], w)
    out$diff_vs_first <- d$est
    out$diff_se <- d$se
  }
  rownames(out) <- NULL
  attr(out, "baseline") <- baseline
  if (n_boot > 0L) out <- .cm_policy_bootstrap(out, contrib, n_boot, seed, conf_level, w)
  out
}

# Unit weights for policy values: NULL, a numeric vector, or a column name.
.cm_policy_weights <- function(weights, scores) {
  if (is.null(weights)) return(NULL)
  w <- .cm_get_optional_numeric(weights, scores$data, scores$n, "weights")
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) stop("`weights` must be finite, non-negative, and not all zero.", call. = FALSE)
  as.numeric(w)
}

# Weighted mean of each column of `contrib` with the standard error
# sqrt(sum w^2 (c - est)^2) / sum w (unweighted: sd / sqrt(n)).
.cm_weighted_value <- function(contrib, w = NULL) {
  contrib <- as.matrix(contrib)
  n <- nrow(contrib)
  if (is.null(w)) {
    return(list(est = colMeans(contrib), se = apply(contrib, 2, stats::sd) / sqrt(n)))
  }
  sw <- sum(w)
  est <- colSums(contrib * w) / sw
  se <- sqrt(colSums(sweep(contrib, 2, est)^2 * w^2)) / sw
  list(est = est, se = se)
}

# Nonparametric bootstrap of the policy values: resample units, recompute
# the means of the per-unit contributions, and attach `boot.se`,
# `boot.low`, `boot.high` (percentile interval) plus the draws as attributes
# "boot" (values) and "boot_diff" (differences against the first policy).
.cm_policy_bootstrap <- function(out, contrib, n_boot, seed, conf_level, w = NULL) {
  n <- nrow(contrib)
  B <- .cm_with_seed(seed, {
    t(vapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      .cm_weighted_value(contrib[idx, , drop = FALSE], if (is.null(w)) NULL else w[idx])$est
    }, numeric(ncol(contrib))))
  })
  B <- matrix(B, ncol = ncol(contrib))
  colnames(B) <- out$policy
  alpha <- (1 - conf_level) / 2
  out$boot.se <- apply(B, 2, stats::sd)
  out$boot.low <- apply(B, 2, stats::quantile, probs = alpha, names = FALSE)
  out$boot.high <- apply(B, 2, stats::quantile, probs = 1 - alpha, names = FALSE)
  attr(out, "boot") <- B
  attr(out, "boot_diff") <- B - B[, 1]
  out
}

.cm_policy_assignment <- function(policy, data, nm = "policy") {
  n <- nrow(data)
  if (inherits(policy, "cm_policy")) return(as.numeric(stats::predict(policy, data)))
  if (is.function(policy)) return(as.numeric(policy(data)))
  if (is.logical(policy) || is.numeric(policy)) {
    if (length(policy) != n) stop("`", nm, "` must have length ", n, ".", call. = FALSE)
    v <- as.numeric(policy)
    if (any(!v %in% c(0, 1))) stop("`", nm, "` must be 0/1.", call. = FALSE)
    return(v)
  }
  stop("`", nm, "` must be a 0/1 vector, a `cm_policy` object, or a function of the data.", call. = FALSE)
}

# Multi-arm assignment as arm indices (1..K).
.cm_policy_arm_index <- function(policy, scores, nm = "policy") {
  n <- scores$n
  lab <- if (inherits(policy, "cm_policy")) stats::predict(policy, scores$data) else if (is.function(policy)) policy(scores$data) else policy
  if (length(lab) != n) stop("`", nm, "` must have length ", n, ".", call. = FALSE)
  idx <- match(as.character(lab), as.character(scores$arms))
  if (anyNA(idx)) stop("`", nm, "` contains values that are not arms of `scores`.", call. = FALSE)
  idx
}

.cm_policy_value_multi <- function(scores, policy, baseline = "none", cost = 0, conf_level = 0.95,
                                   n_boot = 0L, seed = NULL, weights = NULL) {
  pols <- if (is.list(policy) && !inherits(policy, "cm_policy")) policy else list(policy = policy)
  if (is.null(names(pols)) || any(names(pols) == "")) names(pols) <- paste0("policy", seq_along(pols))
  K <- length(scores$arms)
  cost <- rep_len(cost, K)
  G <- sweep(scores$gamma, 2, cost)
  n <- scores$n
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  base_idx <- NULL
  if (!identical(baseline, "none")) {
    base_idx <- match(as.character(baseline), as.character(scores$arms))
    if (is.na(base_idx)) stop("`baseline` must be \"none\" or one of the arms.", call. = FALSE)
  }
  contrib <- sapply(names(pols), function(nm) {
    idx <- .cm_policy_arm_index(pols[[nm]], scores, nm)
    v <- G[cbind(seq_len(n), idx)]
    if (!is.null(base_idx)) v <- v - G[, base_idx]
    v
  })
  contrib <- matrix(contrib, nrow = n)
  shares <- t(sapply(names(pols), function(nm) {
    idx <- .cm_policy_arm_index(pols[[nm]], scores, nm)
    tabulate(idx, K) / n
  }))
  colnames(shares) <- paste0("share_", scores$arms)
  w <- .cm_policy_weights(weights, scores)
  wv <- .cm_weighted_value(contrib, w)
  est <- wv$est
  se <- wv$se
  out <- data.frame(policy = names(pols), shares, estimate = est, std.error = se,
                    conf.low = est - crit * se, conf.high = est + crit * se, stringsAsFactors = FALSE)
  if (!is.null(base_idx)) {
    bv <- if (is.null(w)) mean(G[, base_idx]) else sum(w * G[, base_idx]) / sum(w)
    out$baseline_value <- bv
    out$gain_pct <- 100 * est / bv
  }
  if (length(pols) > 1L) {
    d <- .cm_weighted_value(contrib - contrib[, 1], w)
    out$diff_vs_first <- d$est
    out$diff_se <- d$se
  }
  rownames(out) <- NULL
  attr(out, "baseline") <- baseline
  if (n_boot > 0L) out <- .cm_policy_bootstrap(out, contrib, n_boot, seed, conf_level, w)
  out
}

#' Learn a treatment policy from doubly robust scores
#'
#' Empirical welfare maximization (Kitagawa and Tetenov 2018; Athey and
#' Wager 2021): choose the rule in a restricted class that maximizes the
#' sample analogue of `E[(2 pi(X) - 1)(Y(eta) - cost)]`. Because
#' `E[(2 pi - 1) Gamma] = E[1{2 pi - 1 = sign(Gamma)} |Gamma|]` up to a
#' constant, this is a weighted classification problem with labels
#' `sign(Gamma_i)` and weights `|Gamma_i|`, `Gamma_i = Y_i(eta) - cost`.
#' With a multi-arm score object the objective is the mean of
#' `Gamma[i, pi(X_i)]` over the arms, as in `policytree`.
#'
#' * `"tree"`: exact search over depth-1 or depth-2 trees for any number of
#'   arms. Depth 1 scans every split of every variable in one pass over
#'   sorted cumulative sums. Depth 2 scans root splits over a grid of at most
#'   `max_root_splits` thresholds per variable (all thresholds when there are
#'   fewer) and solves both children exactly. With `engine = "policytree"`
#'   the search is delegated to that package (any depth).
#' * `"linear"`: rules `1{x'b > 0}`; fitted by the weighted logistic
#'   surrogate of the classification problem (binary only).
#' * `"classifier"`: any `mlr3` classification learner accepting weights
#'   (default `classif.rpart`; binary only).
#' * `"budget"`: treat the top `budget` share ranked by a CATE model
#'   `tau_hat` (a `cm_cate` object; binary only).
#'
#' A share `holdout` of the rows is set aside before fitting so the reported
#' held-out value is not contaminated by the search.
#'
#' @param scores A `cm_scores` or `cm_scores_multi` object.
#' @param x Character vector of policy variables (default the `x` of
#'   `scores`).
#' @param method `"tree"`, `"linear"`, `"classifier"`, or `"budget"`.
#' @param depth Tree depth, 1 or 2 (any depth with `engine = "policytree"`).
#' @param cost Cost of treating one unit (binary) or per-arm costs
#'   (multi-arm), in outcome units.
#' @param budget Share of the population to treat (`method = "budget"`).
#' @param tau_hat CATE model for `method = "budget"`.
#' @param learner `mlr3` classification learner for `method = "classifier"`.
#' @param holdout Share of rows held out for an honest value (0 disables).
#' @param engine `"exhaustive"` (default) or `"policytree"`.
#' @param max_root_splits Maximum root thresholds per variable in the
#'   depth-2 exhaustive search.
#' @param baseline Multi-arm case: the arm of the uniform policy the
#'   held-out value is compared with (default the first arm).
#' @param weights Optional positive unit weights (vector or column name).
#'   The objective becomes `sum(w_i Gamma[i, pi(X_i)]) / sum(w_i)`: the
#'   scores are scaled by the normalized weights in the tree searches, and
#'   the classification weights are multiplied by them. Use the
#'   time-decaying weights `h_t` of Zhan et al. (2024) with adaptively
#'   collected data, or survey weights. The reported values use the same
#'   weights.
#' @param split_step policytree engine only: consider every `split_step`-th
#'   split point (`policytree`'s `split.step`); larger values are faster
#'   and approximate.
#' @param seed Optional seed for the holdout split.
#' @param conf_level Confidence level for the values.
#'
#' @return A list of class `cm_policy`: `rule` (tree nodes, coefficients,
#'   fitted classifier with the underlying model in `rule$model`, or
#'   threshold; also available as `model`), `assign` (0/1 or arm labels on all
#'   rows), `value` (data frame with the in-sample and held-out values;
#'   binary: over no treatment and over treating everyone; multi-arm: the
#'   expected outcome and the gain over the uniform baseline arm),
#'   `holdout_id`, and settings. [predict()] applies the rule to new data;
#'   [plot_policy_tree()] draws a tree.
#' @examples
#' dat <- sim_hte(2000, dgp = "policy", seed = 1)
#' sc <- dr_scores(dat, "y", "d", paste0("x", 1:5), seed = 1)
#' pol <- policy_learn(sc, method = "tree", depth = 2, seed = 1)
#' pol
#' attr(dat, "optimal_tree")
#' @references
#' Athey, S. and Wager, S. (2021). Policy learning with observational data.
#' *Econometrica*, 89(1), 133-161.
#'
#' Kitagawa, T. and Tetenov, A. (2018). Who should be treated? Empirical
#' welfare maximization methods for treatment choice. *Econometrica*, 86(2),
#' 591-616.
#' @seealso [policy_value()], [policy_frontier()]
#' @export
policy_learn <- function(scores, x = NULL, method = c("tree", "linear", "classifier", "budget"),
                         depth = 2L, cost = 0, budget = NULL, tau_hat = NULL, learner = NULL,
                         holdout = 0.5, engine = c("exhaustive", "policytree"),
                         max_root_splits = 50L, baseline = NULL, weights = NULL, split_step = 1L,
                         seed = NULL, conf_level = 0.95) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  multi <- inherits(scores, "cm_scores_multi")
  if (!multi) .cm_check_scores(scores)
  if (multi && method != "tree") stop("Multi-arm policies support method = \"tree\" only.", call. = FALSE)
  x <- x %||% scores$x
  if (is.null(x) && method != "budget") stop("`x` is required.", call. = FALSE)
  data <- scores$data
  for (v in x) .cm_check_column(v, data)
  n <- scores$n
  if (multi) {
    K <- length(scores$arms)
    G <- sweep(scores$gamma, 2, rep_len(cost, K))
    labels <- as.character(scores$arms)
  } else {
    gamma <- scores$score - cost
    G <- cbind(0, gamma)
    labels <- c("do not treat", "treat")
  }
  w <- .cm_policy_weights(weights, scores)
  w_norm <- if (is.null(w)) rep(1, n) else w / mean(w)
  Gw <- G * w_norm
  if (holdout < 0 || holdout >= 1) stop("`holdout` must lie in [0, 1).", call. = FALSE)
  test <- rep(FALSE, n)
  if (holdout > 0) {
    m <- floor(holdout * n)
    test[.cm_with_seed(seed, sample.int(n, m))] <- TRUE
  }
  train <- !test
  rule <- NULL
  if (method == "tree") {
    depth <- .cm_check_count(depth, "depth", min = 1L)
    X <- as.matrix(data[, x, drop = FALSE])
    if (!is.numeric(X)) stop("Policy variables must be numeric for trees.", call. = FALSE)
    if (engine == "policytree") {
      .cm_check_package("policytree")
      fit <- policytree::policy_tree(X[train, , drop = FALSE], Gw[train, , drop = FALSE], depth = depth,
                                     split.step = split_step)
      rule <- list(type = "policytree", fit = fit, x = x)
    } else {
      if (depth > 2L) stop("The exhaustive search supports depth 1 or 2; use engine = \"policytree\" for deeper trees.", call. = FALSE)
      rule <- .cm_tree_search(X[train, , drop = FALSE], Gw[train, , drop = FALSE], depth, max_root_splits)
      rule$type <- "tree"
      rule$x <- x
    }
  } else if (method == "linear") {
    df <- data[train, x, drop = FALSE]
    df$.cm_lab <- as.integer(gamma[train] > 0)
    fit <- suppressWarnings(stats::glm(stats::reformulate(x, ".cm_lab"), data = df,
                                       family = stats::binomial(), weights = abs(gamma[train]) * w_norm[train]))
    rule <- list(type = "linear", coefficients = stats::coef(fit), x = x, fit = fit)
  } else if (method == "classifier") {
    .cm_require_mlr3()
    if (is.null(learner)) learner <- mlr3::lrn("classif.rpart", predict_type = "prob")
    df <- data[train, x, drop = FALSE]
    df$.cm_lab <- as.integer(gamma[train] > 0)
    fit <- .cm_fit_mlr3(learner, df, ".cm_lab", x, weights = abs(gamma[train]) * w_norm[train], task_hint = "policy")
    rule <- list(type = "classifier", fit = fit, x = x, model = fit$learner$model)
  } else {
    if (is.null(budget) || budget <= 0 || budget > 1) stop("`budget` must lie in (0, 1].", call. = FALSE)
    if (is.null(tau_hat)) stop("`tau_hat` (a `cm_cate` object) is required for method = \"budget\".", call. = FALSE)
    tau <- .cm_cate_predictions(tau_hat, scores, "tau_hat")
    thr <- .cm_quantile_threshold(tau[train], budget)
    rule <- list(type = "budget", tau_model = tau_hat, threshold = thr, budget = budget)
  }
  obj <- structure(list(rule = rule, method = method, depth = if (method == "tree") depth else NA_integer_,
                        cost = cost, x = x, n = n, holdout_id = test, engine = engine,
                        multi = multi, arms = if (multi) scores$arms else NULL, labels = labels,
                        weights = w, conf_level = conf_level, call = match.call()), class = "cm_policy")
  obj$assign <- predict(obj, data)
  obj$model <- switch(rule$type, linear = rule$fit, classifier = rule$fit$learner, policytree = rule$fit,
                      budget = rule$tau_model, rule)
  if (multi) {
    base_arm <- baseline %||% scores$arms[1]
    val <- function(rows, label) {
      if (!any(rows)) return(NULL)
      sub <- .cm_subset_scores(scores, rows)
      wr <- if (is.null(w)) NULL else w[rows]
      v0 <- policy_value(sub, obj$assign[rows], baseline = "none", cost = cost, conf_level = conf_level, weights = wr)
      v1 <- policy_value(sub, obj$assign[rows], baseline = base_arm, cost = cost, conf_level = conf_level, weights = wr)
      out <- data.frame(sample = label, n = sum(rows), value = v0$estimate, se = v0$std.error,
                        gain_vs_baseline = v1$estimate, se_gain = v1$std.error, gain_pct = v1$gain_pct)
      cbind(out, v0[, grep("^share_", names(v0)), drop = FALSE])
    }
  } else {
    val <- function(rows, label) {
      if (!any(rows)) return(NULL)
      sub <- .cm_subset_scores(scores, rows)
      wr <- if (is.null(w)) NULL else w[rows]
      v0 <- policy_value(sub, obj$assign[rows], baseline = "none", cost = cost, conf_level = conf_level, weights = wr)
      v1 <- policy_value(sub, obj$assign[rows], baseline = "all", cost = cost, conf_level = conf_level, weights = wr)
      data.frame(sample = label, n = sum(rows), share_treated = v0$share_treated,
                 value_vs_none = v0$estimate, se_vs_none = v0$std.error,
                 value_vs_all = v1$estimate, se_vs_all = v1$std.error)
    }
  }
  obj$value <- rbind(val(train, "estimation"), val(test, "holdout"))
  rownames(obj$value) <- NULL
  obj
}

# A cm_scores (or cm_scores_multi) object restricted to `rows` (used for
# held-out values).
.cm_subset_scores <- function(scores, rows) {
  out <- scores
  if (inherits(scores, "cm_scores_multi")) {
    out$gamma <- scores$gamma[rows, , drop = FALSE]
    out$nuisance <- lapply(scores$nuisance, function(M) M[rows, , drop = FALSE])
    out$y <- scores$y[rows]; out$w <- scores$w[rows]; out$w_index <- scores$w_index[rows]
    out$data <- scores$data[rows, , drop = FALSE]
    out$fold_id <- scores$fold_id[rows]
    out$n <- sum(rows)
    return(out)
  }
  out$score <- scores$score[rows]
  out$nuisance <- scores$nuisance[rows, , drop = FALSE]
  out$residuals <- scores$residuals[rows, , drop = FALSE]
  out$y <- scores$y[rows]; out$d <- scores$d[rows]
  out$data <- scores$data[rows, , drop = FALSE]
  out$fold_id <- scores$fold_id[rows]
  out$n <- sum(rows)
  out
}

# Exact tree search -----------------------------------------------------------

# Best depth-1 rule for rows `idx` given pre-sorted orders. Returns a leaf
# (action 0/1) or a split with actions on each side. Value = sum of gamma
# over treated rows.
# Best depth-1 rule for rows `idx` given a reward matrix G (n x K, one
# column per action) and pre-sorted orders. Returns a leaf (the action with
# the largest total reward) or a split with the best action on each side.
# The binary case passes G = cbind(0, gamma): action 1 is "do not treat".
.cm_row_max <- function(M) {
  if (ncol(M) == 1L) return(as.numeric(M[, 1]))
  do.call(pmax, as.data.frame(M))
}

.cm_split1 <- function(X, G, idx, orders) {
  totals <- colSums(G[idx, , drop = FALSE])
  best <- list(value = max(totals), split = FALSE, action = which.max(totals))
  n_idx <- length(idx)
  if (n_idx < 2L) return(best)
  member <- logical(nrow(X))
  member[idx] <- TRUE
  for (j in seq_len(ncol(X))) {
    o <- orders[[j]][member[orders[[j]]]]
    xs <- X[o, j]
    cs <- apply(G[o, , drop = FALSE], 2, cumsum)
    if (is.null(dim(cs))) cs <- matrix(cs, ncol = ncol(G))
    ok <- which(xs[-n_idx] < xs[-1L])
    if (length(ok) == 0L) next
    left <- cs[ok, , drop = FALSE]
    right <- sweep(-left, 2, totals, "+")
    v <- .cm_row_max(left) + .cm_row_max(right)
    m <- which.max(v)
    if (v[m] > best$value + 1e-12) {
      k <- ok[m]
      best <- list(value = v[m], split = TRUE, var = j,
                   threshold = (xs[k] + xs[k + 1L]) / 2,
                   left_action = which.max(left[m, ]), right_action = which.max(right[m, ]))
    }
  }
  best
}

.cm_tree_search <- function(X, G, depth, max_root_splits = 50L) {
  n <- nrow(X)
  p <- ncol(X)
  orders <- lapply(seq_len(p), function(j) order(X[, j]))
  idx <- seq_len(n)
  vars <- colnames(X)
  if (depth == 1L) {
    b <- .cm_split1(X, G, idx, orders)
    return(.cm_tree_node(b, vars))
  }
  best <- .cm_split1(X, G, idx, orders)
  best_tree <- .cm_tree_node(best, vars)
  for (j in seq_len(p)) {
    xs_sorted <- X[orders[[j]], j]
    u <- unique(xs_sorted)
    if (length(u) < 2L) next
    cands <- if (length(u) - 1L <= max_root_splits) {
      (u[-length(u)] + u[-1L]) / 2
    } else {
      unique(stats::quantile(xs_sorted, probs = seq_len(max_root_splits) / (max_root_splits + 1L), names = FALSE, type = 1))
    }
    for (t in cands) {
      L <- idx[X[, j] <= t]
      R <- idx[X[, j] > t]
      if (length(L) == 0L || length(R) == 0L) next
      bl <- .cm_split1(X, G, L, orders)
      br <- .cm_split1(X, G, R, orders)
      v <- bl$value + br$value
      if (v > best_tree$value + 1e-12) {
        best_tree <- list(value = v, split = TRUE, var = vars[j], threshold = t,
                          left = .cm_tree_node(bl, vars), right = .cm_tree_node(br, vars))
      }
    }
  }
  best_tree
}

.cm_tree_node <- function(b, vars) {
  if (!b$split) return(list(value = b$value, split = FALSE, action = b$action))
  list(value = b$value, split = TRUE, var = vars[b$var], threshold = b$threshold,
       left = list(value = NA_real_, split = FALSE, action = b$left_action),
       right = list(value = NA_real_, split = FALSE, action = b$right_action))
}

# Actions are column indices of G (1-based).
.cm_tree_predict <- function(node, data) {
  n <- nrow(data)
  if (!node$split) return(rep(node$action, n))
  left <- data[[node$var]] <= node$threshold
  out <- numeric(n)
  out[left] <- .cm_tree_predict(node$left, data[left, , drop = FALSE])
  out[!left] <- .cm_tree_predict(node$right, data[!left, , drop = FALSE])
  out
}

.cm_tree_text <- function(node, indent = 0, labels = c("do not treat", "treat")) {
  pad <- strrep("  ", indent)
  if (!node$split) return(paste0(pad, "action = ", labels[node$action], "\n"))
  paste0(pad, "if ", node$var, " <= ", format(signif(node$threshold, 4)), ":\n",
         .cm_tree_text(node$left, indent + 1, labels),
         pad, "else:\n",
         .cm_tree_text(node$right, indent + 1, labels))
}

#' @export
predict.cm_policy <- function(object, newdata, ...) {
  newdata <- as.data.frame(newdata)
  r <- object$rule
  multi <- isTRUE(object$multi)
  idx <- switch(r$type,
    tree = .cm_tree_predict(r, newdata),
    policytree = as.numeric(predict(r$fit, as.matrix(newdata[, r$x, drop = FALSE]))),
    linear = {
      X <- stats::model.matrix(stats::reformulate(r$x), newdata)
      as.numeric(X %*% r$coefficients > 0) + 1
    },
    classifier = as.numeric(.cm_predict_fit(r$fit, newdata) > 0.5) + 1,
    budget = as.numeric(predict(r$tau_model, newdata) >= r$threshold) + 1
  )
  if (multi) return(object$arms[idx])
  as.numeric(idx) - 1
}

#' @export
print.cm_policy <- function(x, ...) {
  cat("Treatment policy (", x$method, if (x$method == "tree") paste0(", depth ", x$depth, ", ", x$engine),
      if (isTRUE(x$multi)) paste0(", ", length(x$arms), " arms"), ")\n", sep = "")
  if (isTRUE(x$multi)) {
    sh <- round(100 * prop.table(table(factor(x$assign, levels = x$arms))), 1)
    cat("  share assigned to each arm: ", paste0(names(sh), ": ", sh, "%", collapse = ", "), "\n", sep = "")
  } else {
    cat("  cost = ", x$cost, ", share treated = ", format(round(mean(x$assign), 3)), "\n", sep = "")
  }
  if (x$rule$type == "tree") cat(.cm_tree_text(x$rule, labels = x$labels))
  if (x$rule$type == "policytree") print(x$rule$fit)
  if (x$rule$type == "linear") { cat("  treat if x'b > 0 with b:\n"); print(round(x$rule$coefficients, 4)) }
  if (x$rule$type == "budget") cat("  treat if predicted CATE >= ", format(signif(x$rule$threshold, 4)),
                                   " (top ", 100 * x$rule$budget, "%)\n", sep = "")
  cat(if (isTRUE(x$multi)) "  expected outcome and gain over the uniform baseline arm:\n" else "  value over no treatment / over treating everyone:\n")
  print(x$value, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Draw a depth-limited policy tree
#'
#' @param x A `cm_policy` object with `method = "tree"` and the exhaustive
#'   engine.
#' @return A ggplot object.
#' @export
plot_policy_tree <- function(x) {
  if (!inherits(x, "cm_policy") || x$rule$type != "tree") stop("`x` must be a tree policy from the exhaustive engine.", call. = FALSE)
  labels <- x$labels
  nodes <- list(); edges <- list()
  walk <- function(node, xpos, ypos, width, id) {
    lab <- if (node$split) paste0(node$var, " <= ", format(signif(node$threshold, 3))) else labels[node$action]
    nodes[[length(nodes) + 1L]] <<- data.frame(x = xpos, y = ypos, label = lab, type = if (node$split) "split" else "leaf")
    if (node$split) {
      for (side in c("left", "right")) {
        dx <- if (side == "left") -width else width
        edges[[length(edges) + 1L]] <<- data.frame(x = xpos, y = ypos, xend = xpos + dx, yend = ypos - 1,
                                                   label = if (side == "left") "yes" else "no")
        walk(node[[side]], xpos + dx, ypos - 1, width / 2, paste0(id, side))
      }
    }
  }
  walk(x$rule, 0, 0, 2, "r")
  nd <- do.call(rbind, nodes); ed <- do.call(rbind, edges)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), colour = "grey50") +
    ggplot2::geom_label(data = ed, ggplot2::aes(x = (.data$x + .data$xend) / 2, y = (.data$y + .data$yend) / 2, label = .data$label), size = 3, label.size = 0) +
    ggplot2::geom_label(data = nd, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, fill = .data$type), size = 3.5) +
    ggplot2::scale_fill_manual(values = c(split = "grey90", leaf = "#BFE3D0"), guide = "none") +
    ggplot2::theme_void()
  p
}

#' Targeting frontier between impact and deprivation
#'
#' Implements the targeting exercise of Haushofer, Niehaus, Paramo, Miguel,
#' and Walker (2025): a planner with utility `u` allocates a transfer to the
#' `budget` share of units with the largest predicted utility gain
#' `u(y0 + tau) - u(y0)`, where `y0` is the predicted untreated outcome
#' (deprivation) and `tau` the predicted effect. With linear utility the
#' planner targets impact only; with strong curvature it targets deprivation.
#' For each curvature value the function reports the overlap of the chosen set
#' with impact-only and deprivation-only targeting, and the doubly robust
#' value of the rule.
#'
#' @param scores A `cm_scores` object (test sample).
#' @param tau_hat Predicted effects: a `cm_cate` object, a numeric vector,
#'   or a column name.
#' @param y0_hat Predicted untreated outcomes (default the cross-fitted
#'   `mu0` of `scores`).
#' @param budget Share of units to treat.
#' @param utility `"cara"` (`u(y) = (1 - exp(-a y)) / a`, linear at `a = 0`)
#'   or `"crra"` (`u(y) = (y^(1 - r) - 1) / (1 - r)`, `log(y)` at `r = 1`,
#'   requires positive outcomes).
#' @param curvature Numeric vector of curvature parameters.
#' @param cost Cost per treated unit for the value calculation.
#' @param conf_level Confidence level.
#'
#' @return A list of class `cm_frontier` with `table` (curvature, overlap
#'   with impact-only and deprivation-only targeting, mean predicted effect
#'   and untreated outcome of the selected units, doubly robust value and
#'   standard error), `selected` (logical matrix, one column per curvature),
#'   `impact_only`, `deprivation_only`, and the inputs.
#'   [plot_policy_frontier()] draws the selected region.
#' @references
#' Haushofer, J., Niehaus, P., Paramo, C., Miguel, E., and Walker, M. (2025).
#' Targeting impact versus deprivation. *American Economic Review*.
#' @export
policy_frontier <- function(scores, tau_hat, y0_hat = NULL, budget = 0.3,
                            utility = c("cara", "crra"), curvature = c(0, 0.5, 1, 2, 5),
                            cost = 0, conf_level = 0.95) {
  utility <- match.arg(utility)
  .cm_check_scores(scores)
  tau <- .cm_cate_predictions(tau_hat, scores, "tau_hat")
  y0 <- if (is.null(y0_hat)) scores$nuisance$mu0 else .cm_cate_predictions(y0_hat, scores, "y0_hat")
  n <- scores$n
  k <- max(1L, round(budget * n))
  u <- function(y, a) {
    if (utility == "cara") {
      if (a == 0) y else (1 - exp(-a * y)) / a
    } else {
      if (any(y <= 0)) stop("CRRA utility needs positive outcomes; rescale `y0_hat`.", call. = FALSE)
      if (a == 1) log(y) else (y^(1 - a) - 1) / (1 - a)
    }
  }
  top_k <- function(v) {
    sel <- logical(n)
    sel[order(v, decreasing = TRUE)[seq_len(k)]] <- TRUE
    sel
  }
  impact_only <- top_k(tau)
  deprivation_only <- top_k(-y0)
  selected <- vapply(curvature, function(a) top_k(u(y0 + tau, a) - u(y0, a)), logical(n))
  selected <- matrix(selected, nrow = n, dimnames = list(NULL, paste0("a=", curvature)))
  rows <- lapply(seq_along(curvature), function(i) {
    s <- selected[, i]
    v <- policy_value(scores, as.numeric(s), cost = cost, conf_level = conf_level)
    data.frame(curvature = curvature[i],
               overlap_impact = mean(s[impact_only]),
               overlap_deprivation = mean(s[deprivation_only]),
               mean_tau_selected = mean(tau[s]), mean_y0_selected = mean(y0[s]),
               value = v$estimate, std.error = v$std.error)
  })
  tab <- do.call(rbind, rows)
  structure(list(table = tab, selected = selected, impact_only = impact_only,
                 deprivation_only = deprivation_only, tau = tau, y0 = y0, budget = budget,
                 utility = utility, curvature = curvature, n = n), class = "cm_frontier")
}

#' @export
print.cm_frontier <- function(x, ...) {
  cat("Targeting frontier (", x$utility, " utility, budget = ", x$budget, ", n = ", x$n, ")\n", sep = "")
  print(x$table, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Plot the region selected by a welfare-maximizing planner
#'
#' @param x A `cm_frontier` object.
#' @param curvature One of the curvature values in `x`; default the first.
#' @return A ggplot object of predicted effect against predicted untreated
#'   outcome, selected units highlighted.
#' @export
plot_policy_frontier <- function(x, curvature = x$curvature[1]) {
  j <- match(curvature, x$curvature)
  if (is.na(j)) stop("`curvature` must be one of the values in `x`.", call. = FALSE)
  df <- data.frame(y0_hat = x$y0, tau_hat = x$tau, selected = ifelse(x$selected[, j], "selected", "not selected"))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$y0_hat, y = .data$tau_hat, colour = .data$selected)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_colour_manual(values = c(selected = "#B22222", `not selected` = "grey70")) +
    ggplot2::labs(x = "Predicted untreated outcome (deprivation)", y = "Predicted treatment effect (impact)",
                  colour = NULL, subtitle = paste0(x$utility, " curvature = ", curvature, ", budget = ", x$budget)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
