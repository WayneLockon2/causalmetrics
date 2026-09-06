#' Efficient estimation for staggered rollouts with random timing
#'
#' When treatment timing is (as good as) randomly assigned across units,
#' Roth and Sant'Anna (2023) show that difference-in-differences estimators
#' are unbiased but inefficient, and derive the efficient estimator in the
#' class `theta_hat(beta) = theta_hat_0 - beta' X_hat`, where `theta_hat_0`
#' is the plug-in estimator that compares cohort means after treatment and
#' `X_hat` is the same comparison in the last pre-treatment period, which
#' has mean zero under no anticipation. `beta = 1` gives the Callaway-Sant'Anna
#' estimator with not-yet-treated comparisons, `beta = 0` the simple
#' difference in means, and the efficient `beta*` is estimated from the
#' within-cohort covariance matrices of the outcome path. The function
#' reproduces the `staggered` package.
#'
#' @inheritParams att_gt
#' @param estimand `"simple"` (average over all post-treatment cohort-period
#'   cells, weighted by cohort size), `"cohort"`, `"calendar"`, or
#'   `"eventstudy"` (effect `event_time` periods after adoption).
#' @param event_time Event time(s) for `estimand = "eventstudy"`.
#' @param beta `NULL` for the efficient plug-in `beta*` (default), `1` for the
#'   DiD (Callaway-Sant'Anna) weights, `0` for the difference in means, or
#'   another scalar.
#' @param control `"notyet"` (all cohorts not yet treated at the comparison
#'   period, default) or `"last"` (the last cohort only, as in Sun and
#'   Abraham).
#' @param n_perm Number of permutations for the Fisher randomization test with
#'   the studentized statistic (0 skips it).
#'
#' @return An object of class `cm_staggered` with `estimates` (one row per
#'   event time: `estimate`, `std.error` (the adjusted, less conservative
#'   estimator), `se_neyman` (the conservative Neyman estimator), `beta`,
#'   `theta0`, `xhat`, and the permutation p-value when requested), the
#'   cohort sizes, and the settings.
#'
#' @details Under random timing the finite population of potential outcome
#'   paths is fixed and only the assignment of cohorts is random. With
#'   cohort means `Ybar_g` (a vector over periods) and within-cohort
#'   covariance matrices `S_g`, any estimator in the class is linear in the
#'   cohort means, `theta_hat_0 = sum_g A_g Ybar_g` and
#'   `X_hat = sum_g A0_g Ybar_g`, so `Var(theta_hat_0, X_hat)` is estimable
#'   from `S_g / N_g`, and `beta* = Var(X_hat)^{-1} Cov(X_hat, theta_hat_0)`.
#'   The Neyman variance ignores a non-estimable negative term (the variance
#'   of individual effects) and is conservative; the adjusted variance
#'   subtracts the part of that term that the pre-treatment outcomes
#'   explain, following Section 4 of the paper. The permutation test
#'   reassigns cohorts across units and recomputes the studentized
#'   statistic; it is exact for the sharp null of no effect and
#'   asymptotically valid for the weak null.
#'
#' @references Roth, J. and Sant'Anna, P. H. C. (2023). Efficient estimation
#'   for staggered rollout designs. *Journal of Political Economy:
#'   Microeconomics*, 1(4), 669-709.
#' @examples
#' dat <- sim_did_panel(n_units = 300, n_periods = 8, groups = c(3, 5, 7),
#'                      never_share = 0.3, x_select = 0, seed = 1)
#' staggered_efficient(dat, id = "id", time = "time", group = "g", y = "y", estimand = "simple")
#' staggered_efficient(dat, id = "id", time = "time", group = "g", y = "y",
#'                     estimand = "eventstudy", event_time = 0:2)
#' @export
staggered_efficient <- function(data, id, time, group, y,
                                estimand = c("simple", "cohort", "calendar", "eventstudy"),
                                event_time = 0, beta = NULL, control = c("notyet", "last"),
                                n_perm = 0L, seed = NULL) {
  estimand <- match.arg(estimand)
  control <- match.arg(control)
  for (v in c(id, time, group, y)) .cm_check_column(v, data)
  n_perm <- .cm_check_count(n_perm, "n_perm", 0L)
  if (!is.null(beta) && (!is.numeric(beta) || length(beta) != 1L)) stop("`beta` must be NULL or a single number.", call. = FALSE)
  dt <- data.table::as.data.table(data)[, c(id, time, group, y), with = FALSE]
  data.table::setnames(dt, c(".id", ".t", ".g", ".y"))
  gg <- dt$.g
  gg[is.na(gg) | gg == 0] <- Inf
  dt[, .g := gg]
  if (anyNA(dt$.y)) stop("`y` contains missing values.", call. = FALSE)
  tlist <- sort(unique(dt$.t))
  nT <- length(tlist)
  cnt <- dt[, .N, by = .id]
  if (any(cnt$N != nT) || anyDuplicated(dt[, list(.id, .t)])) stop("A balanced panel is required.", call. = FALSE)
  early <- dt[.g <= tlist[1L], unique(.id)]
  if (length(early)) {
    warning(length(early), " units treated in the first period or earlier were dropped (no pre-treatment period).", call. = FALSE)
    dt <- dt[!.id %in% early]
  }
  units <- dt[.t == tlist[1L], list(.id, .g)]
  sizes <- units[, .N, by = .g][order(.g)]
  single <- sizes[N == 1L, .g]
  if (length(single)) {
    warning("Cohorts with a single unit were dropped: ", paste(single, collapse = ", "), ".", call. = FALSE)
    dt <- dt[!.g %in% single]
    units <- units[!.g %in% single]
    sizes <- sizes[!.g %in% single]
  }
  glist <- sizes$.g
  N_g <- sizes$N
  G <- length(glist)
  if (G < 2L) stop("At least two cohorts (including a never-treated or last-treated cohort) are required.", call. = FALSE)
  Ymat <- matrix(NA_real_, nrow(units), nT)
  Ymat[cbind(match(dt$.id, units$.id), match(dt$.t, tlist))] <- dt$.y
  g_unit <- units$.g

  summaries <- function(g_unit) {
    Ybar <- t(sapply(glist, function(g) colMeans(Ymat[g_unit == g, , drop = FALSE])))
    S <- lapply(glist, function(g) stats::var(Ymat[g_unit == g, , drop = FALSE]))
    list(Ybar = Ybar, S = S)
  }
  mats <- .cm_stag_matrices(estimand, event_time, glist, tlist, N_g, control)
  compute <- function(g_unit, beta_in) {
    sm <- summaries(g_unit)
    lapply(mats, function(m) .cm_stag_estimate(m$A, m$A0, sm$Ybar, sm$S, N_g, glist, tlist, beta_in))
  }
  obs <- compute(g_unit, beta)
  tab <- do.call(rbind, lapply(seq_along(obs), function(k) {
    data.frame(event_time = if (estimand == "eventstudy") event_time[k] else NA_real_,
               estimate = obs[[k]]$estimate, std.error = obs[[k]]$se, se_neyman = obs[[k]]$se_neyman,
               beta = obs[[k]]$beta, theta0 = obs[[k]]$theta0, xhat = obs[[k]]$xhat)
  }))
  if (n_perm > 0L) {
    perm <- .cm_with_seed(seed, lapply(seq_len(n_perm), function(b) {
      gp <- sample(g_unit)
      r <- tryCatch(compute(gp, beta), error = function(e) NULL)
      if (is.null(r)) return(rep(NA_real_, 2L * length(obs)))
      unlist(lapply(r, function(o) c(o$estimate / o$se, o$estimate / o$se_neyman)))
    }))
    perm <- do.call(rbind, perm)
    for (k in seq_along(obs)) {
      t_obs <- tab$estimate[k] / tab$std.error[k]
      t_ney <- tab$estimate[k] / tab$se_neyman[k]
      tab$fisher_p.value[k] <- mean(abs(t_obs) < abs(perm[, 2L * k - 1L]), na.rm = TRUE)
      tab$fisher_p.value_neyman[k] <- mean(abs(t_ney) < abs(perm[, 2L * k]), na.rm = TRUE)
    }
  }
  if (estimand != "eventstudy") tab$event_time <- NULL
  out <- list(estimates = tab, estimand = estimand, control = control, beta_input = beta,
              cohorts = data.frame(group = glist, n = N_g), n_units = nrow(units), periods = tlist,
              n_perm = n_perm, call = match.call())
  class(out) <- "cm_staggered"
  out
}

# Contrast vector over cohorts for ATE(t, g): cohort g against the cohorts
# not yet treated at max(t, g) (or the last cohort only).
.cm_stag_contrast <- function(t, g, glist, N_g, control) {
  if (t >= max(glist)) stop("Period ", t, " has no comparison cohort.", call. = FALSE)
  ctrl <- if (control == "notyet") which(glist > max(g, t)) else which(glist > t & glist == max(glist))
  a <- numeric(length(glist))
  a[ctrl] <- -N_g[ctrl] / sum(N_g[ctrl])
  a[glist == g] <- a[glist == g] + 1
  a
}

# Build the G x T matrices A_theta (weights on cohort means defining the
# estimand) and A0 (the same contrasts in the period before adoption).
.cm_stag_matrices <- function(estimand, event_time, glist, tlist, N_g, control) {
  G <- length(glist); nT <- length(tlist)
  col_of <- function(t) match(t, tlist)
  tg_pair <- function(t, g, w) {
    A <- matrix(0, G, nT); A0 <- matrix(0, G, nT)
    A[, col_of(t)] <- w * .cm_stag_contrast(t, g, glist, N_g, control)
    if (any(tlist == g - 1)) A0[, col_of(g - 1)] <- w * .cm_stag_contrast(t, g, glist, N_g, control)
    list(A = A, A0 = A0)
  }
  add <- function(x, y) list(A = x$A + y$A, A0 = x$A0 + y$A0)
  zero <- list(A = matrix(0, G, nT), A0 = matrix(0, G, nT))
  one <- function(e) {
    if (estimand == "simple") {
      pairs <- expand.grid(g = glist, t = tlist)
      pairs <- pairs[pairs$t >= pairs$g & pairs$t < max(glist), ]
      pairs$N <- N_g[match(pairs$g, glist)]
      w <- pairs$N / sum(pairs$N)
      Reduce(add, Map(function(t, g, w) tg_pair(t, g, w), pairs$t, pairs$g, w), zero)
    } else if (estimand == "cohort") {
      elig <- which(glist < max(glist) & glist <= max(tlist))
      Ne <- sum(N_g[elig])
      Reduce(add, lapply(elig, function(k) {
        g <- glist[k]
        ts <- tlist[tlist >= g & tlist < max(glist)]
        Reduce(add, lapply(ts, function(t) tg_pair(t, g, (N_g[k] / Ne) / length(ts))), zero)
      }), zero)
    } else if (estimand == "calendar") {
      ts <- tlist[tlist >= min(glist) & tlist < max(glist)]
      Reduce(add, lapply(ts, function(t) {
        ks <- which(glist <= t)
        Nt <- sum(N_g[ks])
        Reduce(add, lapply(ks, function(k) tg_pair(t, glist[k], (N_g[k] / Nt) / length(ts))), zero)
      }), zero)
    } else {
      maxG <- max(glist)
      elig <- which(pmax(glist + e, glist) < maxG & pmax(glist + e, glist) <= max(tlist))
      if (length(elig) == 0L) stop("No cohort is observed at event time ", e, " with a comparison cohort.", call. = FALSE)
      Ne <- sum(N_g[elig])
      Reduce(add, lapply(elig, function(k) {
        g <- glist[k]
        if (!any(tlist == g + e)) return(zero)
        tg_pair(g + e, g, N_g[k] / Ne)
      }), zero)
    }
  }
  if (estimand == "eventstudy") lapply(event_time, one) else list(one(NA))
}

.cm_ginv <- function(M) {
  s <- svd(M)
  d <- s$d
  pos <- d > max(dim(M)) * max(d) * .Machine$double.eps
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / d[pos])
}

.cm_stag_estimate <- function(A_theta, A0, Ybar, S, N_g, glist, tlist, beta) {
  G <- length(glist)
  theta0 <- sum(sapply(seq_len(G), function(g) sum(A_theta[g, ] * Ybar[g, ])))
  xhat <- sum(sapply(seq_len(G), function(g) sum(A0[g, ] * Ybar[g, ])))
  quad <- function(a, b) sum(sapply(seq_len(G), function(g) as.numeric(t(a[g, ]) %*% S[[g]] %*% b[g, ]) / N_g[g]))
  Xvar <- quad(A0, A0)
  Xtheta <- quad(A0, A_theta)
  thetaVar <- quad(A_theta, A_theta)
  if (is.null(beta)) beta <- if (Xvar > 0) Xtheta / Xvar else 0
  estimate <- theta0 - xhat * beta
  var_neyman <- thetaVar + beta^2 * Xvar - 2 * Xtheta * beta
  if (var_neyman < 0) { warning("The conservative variance is negative; set to zero.", call. = FALSE); var_neyman <- 0 }
  # adjustment using pre-treatment periods common to all cohorts with nonzero weight
  nonzero <- which(apply(abs(A_theta), 1, max) > 0)
  gMin <- if (length(nonzero)) glist[min(nonzero)] else glist[1L]
  n_pre <- sum(tlist < gMin)
  adj <- 0
  if (n_pre > 0L && gMin > min(tlist)) {
    M <- cbind(diag(n_pre), matrix(0, n_pre, length(tlist) - n_pre))
    ks <- which(glist >= gMin)
    bsum <- Reduce(`+`, lapply(ks, function(g) .cm_ginv(M %*% S[[g]] %*% t(M)) %*% M %*% S[[g]] %*% A_theta[g, ]))
    avgMSM <- Reduce(`+`, lapply(ks, function(g) M %*% S[[g]] %*% t(M))) / length(ks)
    adj <- as.numeric(t(bsum) %*% avgMSM %*% bsum) / sum(N_g)
  }
  se_adj <- if (var_neyman - adj < 0) 0 else sqrt(var_neyman - adj)
  list(estimate = estimate, se = se_adj, se_neyman = sqrt(var_neyman), beta = beta, theta0 = theta0, xhat = xhat)
}

#' @export
print.cm_staggered <- function(x, digits = 4, ...) {
  cat("Efficient staggered-rollout estimator (Roth and Sant'Anna 2023), estimand: ", x$estimand,
      "; comparison: ", if (x$control == "notyet") "not-yet-treated cohorts" else "last cohort",
      "; beta: ", if (is.null(x$beta_input)) "efficient plug-in" else x$beta_input, "\n", sep = "")
  cat("  ", x$n_units, " units in ", nrow(x$cohorts), " cohorts (", paste(x$cohorts$n, collapse = ", "), ")\n", sep = "")
  print(format(x$estimates, digits = digits), row.names = FALSE)
  if (x$n_perm > 0) cat("Fisher p-values from ", x$n_perm, " permutations of cohort assignment (studentized statistic).\n", sep = "")
  invisible(x)
}

#' @rdname cm_tidiers
#' @export
tidy.cm_staggered <- function(x, ...) {
  tab <- x$estimates
  data.frame(term = if ("event_time" %in% names(tab)) paste0("event_time::", tab$event_time) else x$estimand,
             event_time = if ("event_time" %in% names(tab)) tab$event_time else NA_real_,
             estimate = tab$estimate, std.error = tab$std.error,
             statistic = tab$estimate / tab$std.error,
             p.value = 2 * stats::pnorm(-abs(tab$estimate / tab$std.error)),
             conf.low = tab$estimate - stats::qnorm(0.975) * tab$std.error,
             conf.high = tab$estimate + stats::qnorm(0.975) * tab$std.error,
             stringsAsFactors = FALSE)
}
