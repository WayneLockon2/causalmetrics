# R/synth-inference.R
#
# Inference and robustness for synthetic controls: in-space and in-time
# placebos, conformal inference (Chernozhukov, Wuthrich, Zhu 2021), the
# specification test of Ferman and Pinto (2021), and leave-one-donor-out.

# Refit a cm_synth object on a modified block; keeps every setting.
.cm_synth_refit <- function(x, blk, s = x$settings) {
  f <- .cm_synth_fit(blk, s)
  f$mode <- "average"; f$columns <- x$columns
  class(f) <- "cm_synth"
  f
}

# Placebo block: unit `j` (a donor position) becomes the treated unit.
.cm_synth_swap_block <- function(blk, j, include_treated = FALSE) {
  donors_idx <- seq_len(blk$N0); treated_idx <- blk$N0 + seq_len(blk$N1)
  new_donors <- setdiff(donors_idx, j)
  if (include_treated) new_donors <- c(new_donors, treated_idx)
  ord <- c(new_donors, j)
  b <- blk
  b$Y <- blk$Y[ord, , drop = FALSE]
  b$X <- lapply(blk$X, function(M) M[ord, , drop = FALSE])
  b$N0 <- length(new_donors); b$N1 <- 1L
  b$unit_names <- blk$unit_names[ord]
  b$treated_names <- blk$unit_names[j]
  b
}

#' Placebo tests for a synthetic control
#'
#' In-space placebos reassign the treatment to each donor and re-estimate
#' the synthetic control with the same predictors, giving the permutation
#' distribution of the gap and of the post/pre mean squared prediction error
#' (MSPE) ratio (Abadie, Diamond, and Hainmueller 2010). In-time placebos
#' move the treatment date to an earlier period, using only the data before
#' the real treatment (backdating; Abadie, Diamond, and Hainmueller 2015).
#'
#' @param x A `cm_synth` object in the average mode.
#' @param type `"space"` (default) or `"time"`.
#' @param placebo_time For `type = "time"`, the pseudo adoption period; must
#'   be a pre-treatment period with at least two periods before it.
#' @param mspe_limit For `type = "space"`, placebos whose pre-treatment MSPE
#'   exceeds `mspe_limit` times the treated unit's are dropped from the
#'   p-values (Andersson 2019 uses 20). The full set is kept in the object.
#' @param include_treated For `type = "space"`, keep the real treated unit in
#'   the donor pool of each placebo (default `FALSE`, as in Abadie et al.).
#' @param lags,pre_window For `type = "time"`, the predictor lags and the
#'   covariate window relative to the placebo date; default all pre-placebo
#'   periods.
#'
#' @return For `type = "space"`, an object of class `cm_synth_placebo` with
#'   `gaps` (one column per unit, the treated unit first), `mspe` (pre and
#'   post MSPE, their ratio, and the rank of each unit), `p_value` (share of
#'   kept units, the treated one included, whose ratio is at least the
#'   treated unit's), `p_value_by_period` (share of kept units whose absolute
#'   gap in each post period is at least the treated unit's), `kept`, and
#'   the treated fit. For `type = "time"`, a `cm_synth` object fitted at the
#'   placebo date.
#' @references Abadie, A., Diamond, A., and Hainmueller, J. (2015).
#'   Comparative politics and the synthetic control method. *American
#'   Journal of Political Science*, 59(2), 495-510.
#' @examples
#' dat <- sim_synth_panel(n_donors = 12, t_pre = 12, t_post = 4, effect = 3, seed = 2)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
#' pl <- synth_placebo(fit)
#' pl$p_value
#' @export
synth_placebo <- function(x, type = c("space", "time"), placebo_time = NULL, mspe_limit = Inf,
                          include_treated = FALSE, lags = "all", pre_window = NULL) {
  type <- match.arg(type)
  if (!inherits(x, "cm_synth") || x$mode != "average") stop("`x` must be a `cm_synth` object in the average mode.", call. = FALSE)
  blk <- x$block
  if (type == "time") {
    if (is.null(placebo_time)) stop("`placebo_time` is required for in-time placebos.", call. = FALSE)
    pre_t <- blk$time_values[seq_len(blk$T0)]
    if (!(placebo_time %in% pre_t) || sum(pre_t < placebo_time) < 2L) stop("`placebo_time` must be a pre-treatment period with at least two periods before it.", call. = FALSE)
    keep <- seq_len(blk$T0)
    b <- blk
    b$Y <- blk$Y[, keep, drop = FALSE]
    b$X <- lapply(blk$X, function(M) M[, keep, drop = FALSE])
    b$time_values <- pre_t
    b$T0 <- sum(pre_t < placebo_time); b$T1 <- length(pre_t) - b$T0
    s <- x$settings; s$lags <- lags; s$pre_window <- pre_window
    f <- .cm_synth_refit(x, b, s)
    f$placebo_time <- placebo_time
    f$call <- match.call()
    return(f)
  }
  donors_idx <- seq_len(blk$N0)
  fits <- lapply(donors_idx, function(j) .cm_synth_refit(x, .cm_synth_swap_block(blk, j, include_treated)))
  gaps <- cbind(x$effects$gap, sapply(fits, function(f) f$effects$gap))
  colnames(gaps) <- c(as.character(x$treated_names[1L]), as.character(blk$unit_names[donors_idx]))
  if (x$N1 > 1L) colnames(gaps)[1L] <- "treated (average)"
  post <- x$effects$post
  pre_mspe <- colMeans(gaps[!post, , drop = FALSE]^2)
  post_mspe <- colMeans(gaps[post, , drop = FALSE]^2)
  ratio <- post_mspe / pre_mspe
  kept <- pre_mspe <= mspe_limit * pre_mspe[1L]
  kept[1L] <- TRUE
  mspe <- data.frame(id = colnames(gaps), pre_mspe = pre_mspe, post_mspe = post_mspe, ratio = ratio,
                     kept = kept, stringsAsFactors = FALSE)
  mspe$rank <- NA_integer_
  mspe$rank[kept] <- rank(-ratio[kept], ties.method = "min")
  rownames(mspe) <- NULL
  names(mspe)[1L] <- x$columns[["id"]]
  p_value <- mean(ratio[kept] >= ratio[1L])
  p_by_period <- apply(abs(gaps[post, kept, drop = FALSE]), 1, function(r) mean(r >= r[1L]))
  p_value_by_period <- data.frame(time = x$effects$time[post], gap = x$effects$gap[post], p_value = p_by_period)
  out <- list(gaps = gaps, time = x$effects$time, post = post, mspe = mspe, p_value = p_value,
              p_value_by_period = p_value_by_period, kept = kept, mspe_limit = mspe_limit,
              include_treated = include_treated, fit = x, fits = fits, call = match.call())
  class(out) <- "cm_synth_placebo"
  out
}

#' @export
print.cm_synth_placebo <- function(x, digits = 3, ...) {
  cat("In-space placebo test: ", ncol(x$gaps) - 1L, " placebo units", if (is.finite(x$mspe_limit)) paste0(", ", sum(x$kept) - 1L, " kept (pre-MSPE within ", x$mspe_limit, "x the treated unit's)") else "", "\n", sep = "")
  cat("  Treated post/pre MSPE ratio: ", formatC(x$mspe$ratio[1L], digits = digits, format = "g"),
      "; rank ", x$mspe$rank[1L], " of ", sum(x$kept), "; p-value ", formatC(x$p_value, digits = digits, format = "g"), "\n", sep = "")
  invisible(x)
}

.cm_plot_synth_placebo <- function(x, type, mspe_limit = NULL) {
  type <- match.arg(type, c("placebos", "ratio"))
  base <- ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  kept <- if (is.null(mspe_limit)) x$kept else { k <- x$mspe$pre_mspe <= mspe_limit * x$mspe$pre_mspe[1L]; k[1L] <- TRUE; k }
  cut <- x$time[sum(!x$post)] + 0.5 * (x$time[sum(!x$post) + 1L] - x$time[sum(!x$post)])
  if (type == "placebos") {
    g <- x$gaps[, kept, drop = FALSE]
    d <- data.frame(time = rep(x$time, ncol(g)), gap = as.vector(g),
                    id = rep(colnames(g), each = length(x$time)),
                    series = rep(c("Treated", rep("Placebo", ncol(g) - 1L)), each = length(x$time)))
    return(ggplot2::ggplot() +
             ggplot2::geom_line(data = d[d$series == "Placebo", ], ggplot2::aes(x = .data$time, y = .data$gap, group = .data$id), colour = "grey75") +
             ggplot2::geom_line(data = d[d$series == "Treated", ], ggplot2::aes(x = .data$time, y = .data$gap), colour = "black", linewidth = 0.9) +
             ggplot2::geom_hline(yintercept = 0, linetype = "dotted") + ggplot2::geom_vline(xintercept = cut, linetype = "dotted") +
             ggplot2::labs(x = NULL, y = "Gap: unit minus its synthetic control",
                           subtitle = paste0("Treated unit in black; ", sum(kept) - 1L, " placebo units in grey")) + base)
  }
  m <- x$mspe[kept, ]
  m$id <- factor(m[[1L]], levels = m[[1L]][order(m$ratio)])
  m$series <- ifelse(seq_len(nrow(m)) == 1L, "Treated", "Placebo")
  ggplot2::ggplot(m, ggplot2::aes(x = .data$id, y = .data$ratio, fill = .data$series)) + ggplot2::geom_col() +
    ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = c(Treated = "black", Placebo = "grey65")) +
    ggplot2::labs(x = NULL, y = "Post-period MSPE / pre-period MSPE", fill = NULL) + base
}

# Moving-block permutation p-value of a statistic of the residual vector u.
.cm_block_pvalue <- function(u, T1, stat, permutations = "moving_block", n_perm = 999L) {
  Tt <- length(u)
  s_obs <- stat(u[(Tt - T1 + 1L):Tt])
  if (permutations == "moving_block") {
    s_perm <- vapply(0:(Tt - 1L), function(j) {
      idx <- ((seq_len(Tt) + j - 1L) %% Tt) + 1L
      stat(u[idx][(Tt - T1 + 1L):Tt])
    }, numeric(1))
  } else {
    s_perm <- c(s_obs, vapply(seq_len(n_perm), function(b) stat(sample(u)[(Tt - T1 + 1L):Tt]), numeric(1)))
  }
  list(p_value = mean(s_perm >= s_obs - 1e-12), statistic = s_obs, permuted = s_perm)
}

# Refit with every period as a predictor lag (outcome only) on adjusted data.
.cm_synth_null_fit <- function(x, y_adj, keep_periods = NULL) {
  blk <- x$block
  if (!is.null(keep_periods)) {
    blk$Y <- blk$Y[, keep_periods, drop = FALSE]
    blk$time_values <- blk$time_values[keep_periods]
    y_adj <- y_adj[keep_periods]
  }
  treated_idx <- blk$N0 + seq_len(blk$N1)
  # replace the treated rows by the adjusted treated average
  blk$Y <- rbind(blk$Y[seq_len(blk$N0), , drop = FALSE], matrix(y_adj, 1L, ncol(blk$Y)))
  blk$N1 <- 1L
  blk$X <- list()
  blk$unit_names <- c(blk$unit_names[seq_len(blk$N0)], "treated")
  blk$treated_names <- "treated"
  blk$T0 <- ncol(blk$Y); blk$T1 <- 0L
  s <- x$settings
  s$lags <- "all"; s$pre_window <- NULL; s$x <- NULL; s$v <- "equal"; s$augment <- "none"
  f <- .cm_synth_fit(blk, s)
  y_adj - f$effects$synthetic
}

#' Conformal inference for a synthetic control
#'
#' The test of Chernozhukov, Wuthrich, and Zhu (2021): under a sharp null
#' `H0: tau_t = theta_t` for the post-treatment periods, subtract the null
#' effect from the treated outcome, refit the synthetic control on all
#' periods with every period's outcome as a predictor, and compare the
#' post-treatment residuals with their moving-block permutations. Inverting
#' the test over `grid` gives a confidence set for a constant effect;
#' `per_period = TRUE` repeats the test period by period, dropping the
#' other post-treatment periods.
#'
#' @param x A `cm_synth` object in the average mode.
#' @param null The null effect: a scalar (constant over post periods) or a
#'   vector with one entry per post period. Default 0.
#' @param q Norm of the test statistic, `(mean |u_t|^q)^(1/q)` over the
#'   post periods (default 1).
#' @param permutations `"moving_block"` (all cyclic shifts of the residual
#'   vector; the default, valid under weak dependence) or `"iid"` (random
#'   permutations).
#' @param n_perm Number of random permutations for `permutations = "iid"`.
#' @param grid Values of a constant effect over which to invert the test;
#'   by default 61 points around the estimate, widened until the confidence
#'   set no longer touches an end of the grid.
#' @param level Confidence level of the reported set (default 0.95).
#' @param per_period Also compute a p-value and a confidence set for each
#'   post period separately.
#' @param seed Optional seed for `"iid"` permutations.
#'
#' @return A list of class `cm_synth_conformal` with `p_value` (for `null`),
#'   `statistic`, `conf_set` (`grid` values not rejected at `level`, with
#'   `conf.low` and `conf.high` as its range), `grid_p` (p-value at each
#'   grid value), and, if requested, `per_period` (time, gap, p-value at
#'   zero, and interval).
#' @references Chernozhukov, V., Wuthrich, K., and Zhu, Y. (2021). An exact
#'   and robust conformal inference method for counterfactual and synthetic
#'   controls. *JASA*, 116(536), 1849-1864.
#' @examples
#' dat <- sim_synth_panel(n_donors = 12, t_pre = 15, t_post = 4, effect = 3, seed = 3)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
#' ci <- synth_conformal(fit)
#' c(p_value = ci$p_value, ci$conf.low, ci$conf.high)
#' @export
synth_conformal <- function(x, null = 0, q = 1, permutations = c("moving_block", "iid"), n_perm = 999L,
                            grid = NULL, level = 0.95, per_period = FALSE, seed = NULL) {
  permutations <- match.arg(permutations)
  if (!inherits(x, "cm_synth") || x$mode != "average") stop("`x` must be a `cm_synth` object in the average mode.", call. = FALSE)
  if (x$T1 < 1L) stop("No post-treatment periods.", call. = FALSE)
  post <- x$effects$post
  y_tr <- x$effects$treated
  T1 <- x$T1
  stat <- function(u) mean(abs(u)^q)^(1 / q)
  test_at <- function(theta) {
    th <- if (length(theta) == 1L) rep(theta, T1) else theta
    y_adj <- y_tr; y_adj[post] <- y_adj[post] - th
    u <- .cm_synth_null_fit(x, y_adj)
    .cm_block_pvalue(u, T1, stat, permutations, n_perm)
  }
  res <- .cm_with_seed(seed, {
    main <- test_at(null)
    spread <- max(4 * x$pre_rmspe, 2 * stats::sd(x$effects$gap[post]), abs(x$estimate), 1e-8)
    user_grid <- !is.null(grid)
    if (!user_grid) grid <- seq(x$estimate - spread, x$estimate + spread, length.out = 61L)
    else spread <- max(abs(range(grid) - x$estimate))
    grid_p <- vapply(grid, function(th) test_at(th)$p_value, numeric(1))
    # widen the default grid until the confidence set does not touch its ends
    if (!user_grid) {
      step <- grid[2L] - grid[1L]
      for (k in seq_len(6L)) {
        lo <- grid_p[1L] > 1 - level; hi <- grid_p[length(grid_p)] > 1 - level
        if (!lo && !hi) break
        if (lo) { ext <- seq(grid[1L] - 30 * step, grid[1L] - step, by = step); grid <- c(ext, grid)
                  grid_p <- c(vapply(ext, function(th) test_at(th)$p_value, numeric(1)), grid_p) }
        if (hi) { ext <- seq(grid[length(grid)] + step, grid[length(grid)] + 30 * step, by = step); grid <- c(grid, ext)
                  grid_p <- c(grid_p, vapply(ext, function(th) test_at(th)$p_value, numeric(1))) }
      }
    }
    pp <- NULL
    if (per_period) {
      pre_idx <- which(!post)
      pp <- do.call(rbind, lapply(which(post), function(t_idx) {
        keep <- c(pre_idx, t_idx)
        one <- function(th) {
          y_adj <- y_tr; y_adj[t_idx] <- y_adj[t_idx] - th
          u <- .cm_synth_null_fit(x, y_adj, keep)
          .cm_block_pvalue(u, 1L, stat, permutations, n_perm)$p_value
        }
        g_t <- x$effects$gap[t_idx]
        gr <- seq(g_t - spread, g_t + spread, length.out = 61L)
        pg <- vapply(gr, one, numeric(1))
        acc <- gr[pg > 1 - level]
        data.frame(time = x$effects$time[t_idx], gap = g_t, p_value = one(0),
                   conf.low = if (length(acc)) min(acc) else NA_real_, conf.high = if (length(acc)) max(acc) else NA_real_)
      }))
    }
    list(main = main, grid = grid, grid_p = grid_p, pp = pp)
  })
  acc <- res$grid[res$grid_p > 1 - level]
  out <- list(p_value = res$main$p_value, statistic = res$main$statistic, null = null,
              conf_set = acc, conf.low = if (length(acc)) min(acc) else NA_real_,
              conf.high = if (length(acc)) max(acc) else NA_real_, level = level,
              grid_p = data.frame(theta = res$grid, p_value = res$grid_p), per_period = res$pp,
              q = q, permutations = permutations, estimate = x$estimate, call = match.call())
  class(out) <- "cm_synth_conformal"
  out
}

#' @export
print.cm_synth_conformal <- function(x, digits = 3, ...) {
  cat("Conformal inference (", x$permutations, " permutations, q = ", x$q, ")\n", sep = "")
  cat("  Average post-treatment effect ", formatC(x$estimate, digits = digits, format = "g"),
      "; p-value for effect = ", paste(formatC(x$null, digits = digits, format = "g"), collapse = ", "),
      ": ", formatC(x$p_value, digits = digits, format = "g"), "\n", sep = "")
  cat("  ", round(100 * x$level), "% confidence set for a constant effect: [",
      formatC(x$conf.low, digits = digits, format = "g"), ", ", formatC(x$conf.high, digits = digits, format = "g"), "]\n", sep = "")
  invisible(x)
}

#' Specification test of Ferman and Pinto (2021)
#'
#' Contrasts the demeaned synthetic control with difference-in-differences.
#' Under a linear factor model both are unbiased when the common factors
#' have mean zero after treatment; a large difference between them signals
#' that both are biased. The test imposes the null, fits the demeaned
#' synthetic control with every period's outcome as a predictor, forms
#' `u_t = w'(y_t - ybar) - (1/J) 1'(y_t - ybar)` (Ferman and Pinto, eq. 11),
#' and compares the post-treatment average of `u_t` with its moving-block
#' permutations (their Proposition 4).
#'
#' @param x A `cm_synth` object in the average mode.
#' @param q Norm for the secondary statistic `(mean |u_t|^q)^(1/q)`.
#' @return A list of class `cm_synth_spec` with `p_value` (statistic
#'   `|mean of u_t over post periods|`), `p_value_norm` (the `q`-norm
#'   statistic), `u` (the contrast series), and `paths` (the demeaned
#'   synthetic control and DiD counterfactuals fitted on the pre-treatment
#'   periods, for the plot).
#' @examples
#' dat <- sim_synth_panel(n_donors = 12, t_pre = 15, t_post = 4, effect = 3, seed = 4)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
#' synth_spec_test(fit)
#' @export
synth_spec_test <- function(x, q = 2) {
  if (!inherits(x, "cm_synth") || x$mode != "average") stop("`x` must be a `cm_synth` object in the average mode.", call. = FALSE)
  blk <- x$block
  Tt <- ncol(blk$Y); T1 <- x$T1; T0 <- x$T0
  donors_idx <- seq_len(blk$N0); treated_idx <- blk$N0 + seq_len(blk$N1)
  y0 <- colMeans(blk$Y[treated_idx, , drop = FALSE])
  Y0 <- blk$Y[donors_idx, , drop = FALSE]
  # demeaned SC on all periods, all lags, simplex weights
  b <- blk; b$T0 <- Tt; b$T1 <- 0L; b$X <- list()
  s <- x$settings; s$lags <- "all"; s$pre_window <- NULL; s$x <- NULL; s$v <- "equal"
  s$demean <- TRUE; s$augment <- "none"; s$constraints <- "simplex"
  f_all <- .cm_synth_fit(b, s)
  w <- f_all$weights$weight
  sc_eff <- f_all$effects$gap
  did_eff <- (y0 - colMeans(Y0)) - (mean(y0) - mean(Y0))
  u <- sc_eff - did_eff
  stat_mean <- function(v) abs(mean(v))
  stat_norm <- function(v) mean(abs(v)^q)^(1 / q)
  r1 <- .cm_block_pvalue(u, T1, stat_mean)
  r2 <- .cm_block_pvalue(u, T1, stat_norm)
  # counterfactual paths fitted on the pre-treatment periods, for the plot
  f_pre <- if (x$settings$demean && x$settings$augment == "none" && identical(x$settings$lags, "all") && is.null(x$settings$x)) x else {
    s2 <- s; b2 <- blk; b2$X <- list()
    .cm_synth_fit(b2, s2)
  }
  did_cf <- mean(y0[seq_len(T0)]) + (colMeans(Y0) - mean(Y0[, seq_len(T0)]))
  paths <- data.frame(time = blk$time_values, treated = y0, demeaned_sc = f_pre$effects$synthetic, did = did_cf,
                      post = seq_len(Tt) > T0)
  out <- list(p_value = r1$p_value, statistic = r1$statistic, p_value_norm = r2$p_value, statistic_norm = r2$statistic,
              u = data.frame(time = blk$time_values, u = u, post = seq_len(Tt) > T0), weights = f_all$weights,
              paths = paths, estimates = c(demeaned_sc = mean(f_pre$effects$gap[(T0 + 1L):Tt]), did = mean((y0 - did_cf)[(T0 + 1L):Tt])),
              q = q, call = match.call())
  class(out) <- "cm_synth_spec"
  out
}

#' @export
print.cm_synth_spec <- function(x, digits = 3, ...) {
  cat("Ferman-Pinto specification test: demeaned synthetic control versus DiD\n")
  cat("  Estimates: demeaned SC ", formatC(x$estimates[["demeaned_sc"]], digits = digits, format = "g"),
      ", DiD ", formatC(x$estimates[["did"]], digits = digits, format = "g"), "\n", sep = "")
  cat("  p-value (mean contrast): ", formatC(x$p_value, digits = digits, format = "g"),
      "; p-value (", x$q, "-norm): ", formatC(x$p_value_norm, digits = digits, format = "g"), "\n", sep = "")
  invisible(x)
}

.cm_plot_synth_spec <- function(x) {
  p <- x$paths
  d <- data.frame(time = rep(p$time, 3), y = c(p$treated, p$demeaned_sc, p$did),
                  series = rep(c("Treated", "Demeaned synthetic control", "DiD counterfactual"), each = nrow(p)))
  cut <- p$time[sum(!p$post)] + 0.5 * (p$time[sum(!p$post) + 1L] - p$time[sum(!p$post)])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$y, colour = .data$series, linetype = .data$series)) +
    ggplot2::geom_line() + ggplot2::geom_vline(xintercept = cut, linetype = "dotted") +
    ggplot2::scale_colour_manual(values = c(Treated = "black", `Demeaned synthetic control` = "grey45", `DiD counterfactual` = "#B22222")) +
    ggplot2::scale_linetype_manual(values = c(Treated = "solid", `Demeaned synthetic control` = "longdash", `DiD counterfactual` = "dotdash")) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, linetype = NULL,
                  subtitle = paste0("Specification test p-value ", formatC(x$p_value, digits = 3, format = "g"))) +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}

#' Leave-one-donor-out synthetic controls
#'
#' Refits the synthetic control after dropping, one at a time, each donor
#' with weight above `min_weight`.
#'
#' @param x A `cm_synth` object in the average mode.
#' @param min_weight Donors with weight at or below this value are not
#'   dropped (default 0.001, as in Andersson 2019).
#' @return An object of class `cm_synth_loo` with `table` (dropped donor,
#'   its weight, the new estimate and pre-treatment RMSPE) and `paths` (the
#'   synthetic paths, one column per refit).
#' @examples
#' dat <- sim_synth_panel(n_donors = 10, t_pre = 12, t_post = 4, effect = 2, seed = 5)
#' fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
#' synth_loo(fit)$table
#' @export
synth_loo <- function(x, min_weight = 0.001) {
  if (!inherits(x, "cm_synth") || x$mode != "average") stop("`x` must be a `cm_synth` object in the average mode.", call. = FALSE)
  blk <- x$block
  drop <- which(x$weights$weight > min_weight)
  if (length(drop) == 0L) stop("No donor has weight above `min_weight`.", call. = FALSE)
  if (blk$N0 - 1L < 2L) stop("At least three donors are needed.", call. = FALSE)
  fits <- lapply(drop, function(j) {
    keep <- c(setdiff(seq_len(blk$N0), j), blk$N0 + seq_len(blk$N1))
    b <- blk
    b$Y <- blk$Y[keep, , drop = FALSE]; b$X <- lapply(blk$X, function(M) M[keep, , drop = FALSE])
    b$N0 <- blk$N0 - 1L; b$unit_names <- blk$unit_names[keep]
    .cm_synth_refit(x, b)
  })
  tab <- data.frame(dropped = x$weights[[1L]][drop], weight = x$weights$weight[drop],
                    estimate = vapply(fits, function(f) f$estimate, numeric(1)),
                    pre_rmspe = vapply(fits, function(f) f$pre_rmspe, numeric(1)), stringsAsFactors = FALSE)
  paths <- sapply(fits, function(f) f$effects$synthetic)
  colnames(paths) <- paste0("without ", tab$dropped)
  out <- list(table = tab, paths = paths, time = x$effects$time, treated = x$effects$treated,
              synthetic = x$effects$synthetic, post = x$effects$post, estimate = x$estimate, fits = fits, call = match.call())
  class(out) <- "cm_synth_loo"
  out
}

#' @export
print.cm_synth_loo <- function(x, digits = 3, ...) {
  cat("Leave-one-donor-out: full estimate ", formatC(x$estimate, digits = digits, format = "g"),
      "; range over refits [", formatC(min(x$table$estimate), digits = digits, format = "g"), ", ",
      formatC(max(x$table$estimate), digits = digits, format = "g"), "]\n", sep = "")
  print(x$table, digits = digits, row.names = FALSE)
  invisible(x)
}

.cm_plot_synth_loo <- function(x) {
  d <- data.frame(time = rep(x$time, ncol(x$paths)), y = as.vector(x$paths),
                  id = rep(colnames(x$paths), each = length(x$time)))
  main <- data.frame(time = rep(x$time, 2), y = c(x$treated, x$synthetic),
                     series = rep(c("Treated", "Synthetic (all donors)"), each = length(x$time)))
  cut <- x$time[sum(!x$post)] + 0.5 * (x$time[sum(!x$post) + 1L] - x$time[sum(!x$post)])
  ggplot2::ggplot() +
    ggplot2::geom_line(data = d, ggplot2::aes(x = .data$time, y = .data$y, group = .data$id), colour = "grey75") +
    ggplot2::geom_line(data = main, ggplot2::aes(x = .data$time, y = .data$y, colour = .data$series, linetype = .data$series)) +
    ggplot2::geom_vline(xintercept = cut, linetype = "dotted") +
    ggplot2::scale_colour_manual(values = c(Treated = "black", `Synthetic (all donors)` = "grey30")) +
    ggplot2::scale_linetype_manual(values = c(Treated = "solid", `Synthetic (all donors)` = "longdash")) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, linetype = NULL, subtitle = "Grey: synthetic control with one donor left out") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}
