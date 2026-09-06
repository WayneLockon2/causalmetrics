#' Difference-in-differences with a continuous or multi-valued treatment
#'
#' `att_dose()` implements the two-period framework of Callaway,
#' Goodman-Bacon, and Sant'Anna (2024) for a dose `D >= 0` with untreated
#' units at `D = 0`. It reports level effects `ATT(d | d)`, the effect of
#' dose `d` among units that received `d`, which parallel trends against the
#' untreated identifies, and causal responses, the slope of the dose-response
#' curve, which need the stronger assumption that every dose group would have
#' followed the same trend at every dose ("strong parallel trends"). It also
#' returns the two weight decompositions of the two-way fixed effects
#' coefficient, on level effects and on causal-response increments.
#'
#' @inheritParams att_gt
#' @param dose Column with the dose received in the post-treatment period
#'   (constant within unit; 0 for untreated units).
#' @param dose_type `"binned"` (default; doses grouped by `breaks` or into
#'   `n_bins` quantile bins) or `"spline"` (a B-spline regression of the
#'   outcome change on the dose among treated units).
#' @param n_bins,breaks Bins for `dose_type = "binned"`.
#' @param df Degrees of freedom of the B-spline for `dose_type = "spline"`.
#' @param grid Dose values at which the spline curve is evaluated (default:
#'   50 points over the treated dose range).
#' @param n_boot Bootstrap replications for the spline curve (units resampled
#'   with replacement).
#'
#' @return An object of class `cm_att_dose` with `by_dose` (level effects by
#'   bin with standard errors), `acr` (causal responses between consecutive
#'   bins, valid under strong parallel trends), `curve` (spline level effects
#'   and derivatives on the grid with bootstrap bands, when
#'   `dose_type = "spline"`), `twfe` (the regression coefficient of the
#'   outcome change on the dose), and `twfe_weights` (the two decompositions).
#'
#' @details The panel must have exactly two periods. The two-way fixed
#'   effects coefficient equals `Cov(dY, D) / Var(D)`, which is a weighted
#'   sum of level effects with weights `(d_k - E[D]) p_k / Var(D)` that can
#'   be negative, and a weighted sum of causal-response increments with
#'   non-negative weights `P(D >= d_k) (E[D | D >= d_k] - E[D]) / Var(D)`
#'   that integrate to one. The two readings of the same coefficient rely on
#'   different assumptions; the function reports both weight sets so the
#'   reader can see which interpretation is being invoked.
#'
#' @references Callaway, B., Goodman-Bacon, A., and Sant'Anna, P. H. C.
#'   (2024). Difference-in-differences with a continuous treatment. Working
#'   paper, arXiv:2107.02637.
#' @examples
#' set.seed(1)
#' n <- 600
#' dose <- ifelse(runif(n) < 0.3, 0, runif(n, 0.5, 3))
#' dat <- data.frame(id = rep(1:n, each = 2), time = rep(1:2, n),
#'                   dose = rep(dose, each = 2))
#' dat$y <- rnorm(2 * n) + rep(rnorm(n), each = 2) + 0.5 * dat$time +
#'   (dat$time == 2) * (2 * dat$dose - 0.3 * dat$dose^2)
#' att_dose(dat, id = "id", time = "time", y = "y", dose = "dose", n_bins = 4)
#' @export
att_dose <- function(data, id, time, y, dose, dose_type = c("binned", "spline"),
                     n_bins = 4L, breaks = NULL, df = 4L, grid = NULL, n_boot = 199L,
                     conf_level = 0.95, seed = NULL) {
  dose_type <- match.arg(dose_type)
  for (v in c(id, time, y, dose)) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, c(id, time, y, dose), with = FALSE]
  data.table::setnames(dt, c(".id", ".t", ".y", ".dose"))
  tlist <- sort(unique(dt$.t))
  if (length(tlist) != 2L) stop("`att_dose()` uses a two-period panel (before and after treatment).", call. = FALSE)
  cnt <- dt[, .N, by = .id]
  if (any(cnt$N != 2L)) stop("Every unit must be observed in both periods.", call. = FALSE)
  chk <- dt[, data.table::uniqueN(.dose), by = .id]
  if (any(chk$V1 > 1L)) stop("`dose` must be constant within unit.", call. = FALSE)
  u <- dt[, list(dy = .y[.t == tlist[2L]] - .y[.t == tlist[1L]], dose = .dose[1L]), by = .id]
  if (any(u$dose < 0)) stop("`dose` must be non-negative.", call. = FALSE)
  if (!any(u$dose == 0)) stop("Untreated units (dose 0) are required.", call. = FALSE)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  n <- nrow(u)
  dy0 <- u$dy[u$dose == 0]
  m0 <- mean(dy0)
  v0 <- stats::var(dy0) / length(dy0)

  # ---- binned level effects and causal responses -------------------------
  treated <- u[u$dose > 0, ]
  if (is.null(breaks)) {
    n_bins <- .cm_check_count(n_bins, "n_bins", 1L)
    breaks <- unique(stats::quantile(treated$dose, probs = seq(0, 1, length.out = n_bins + 1L)))
  }
  treated$bin <- cut(treated$dose, breaks = breaks, include.lowest = TRUE)
  by_dose <- do.call(rbind, lapply(levels(treated$bin), function(b) {
    s <- treated[treated$bin == b, ]
    data.frame(bin = b, dose = mean(s$dose), n = nrow(s), att = mean(s$dy) - m0,
               std.error = sqrt(stats::var(s$dy) / nrow(s) + v0), stringsAsFactors = FALSE)
  }))
  by_dose$conf.low <- by_dose$att - z * by_dose$std.error
  by_dose$conf.high <- by_dose$att + z * by_dose$std.error
  d_prev <- c(0, by_dose$dose[-nrow(by_dose)])
  att_prev <- c(0, by_dose$att[-nrow(by_dose)])
  # variance of the increment: bins are independent samples; the untreated mean cancels
  bin_var <- by_dose$std.error^2 - v0
  var_prev <- c(0, bin_var[-nrow(by_dose)])
  acr <- data.frame(from = d_prev, to = by_dose$dose,
                    increment = by_dose$att - att_prev,
                    acr = (by_dose$att - att_prev) / (by_dose$dose - d_prev))
  acr$std.error <- sqrt(bin_var + var_prev + ifelse(seq_len(nrow(acr)) == 1L, v0, 0)) / (by_dose$dose - d_prev)
  acr$conf.low <- acr$acr - z * acr$std.error
  acr$conf.high <- acr$acr + z * acr$std.error

  # ---- TWFE coefficient and its two weight decompositions ----------------
  beta <- stats::cov(u$dy, u$dose) / stats::var(u$dose)
  ED <- mean(u$dose); VD <- stats::var(u$dose) * (n - 1) / n
  dk <- c(0, by_dose$dose)
  pk <- c(mean(u$dose == 0), by_dose$n / n)
  # use bin means as the discrete dose values; recompute beta on the binned dose for exact identities
  dose_binned <- ifelse(u$dose == 0, 0, by_dose$dose[as.integer(cut(u$dose, breaks = breaks, include.lowest = TRUE))])
  ED_b <- mean(dose_binned); VD_b <- mean((dose_binned - ED_b)^2)
  beta_b <- stats::cov(u$dy, dose_binned) * (n - 1) / n / VD_b
  w_level <- (dk - ED_b) * pk / VD_b
  w_slope <- vapply(seq_along(dk), function(k) {
    if (k == 1L) return(NA_real_)
    ge <- dose_binned >= dk[k]
    mean(ge) * (mean(dose_binned[ge]) - ED_b) / VD_b
  }, numeric(1))
  twfe_weights <- data.frame(dose = dk, share = pk, w_level = w_level, w_slope = w_slope)
  level_effects <- c(0, by_dose$att)
  checks <- c(twfe_binned = beta_b,
              sum_w_level_x_att = sum(w_level * level_effects),
              sum_w_slope_x_increment = sum(w_slope[-1L] * acr$increment),
              sum_w_slope_x_dose_gap = sum(w_slope[-1L] * (dk[-1L] - dk[-length(dk)])))

  # ---- spline curve --------------------------------------------------------
  curve <- NULL
  if (dose_type == "spline") {
    if (is.null(grid)) grid <- seq(min(treated$dose), max(treated$dose), length.out = 50L)
    fit_curve <- function(dd, dyy, dy0m) {
      B <- splines::bs(dd, df = df)
      fit <- stats::lm.fit(cbind(1, B), dyy)
      predict_at <- function(g) {
        g <- pmin(pmax(g, min(dd)), max(dd))
        as.numeric(cbind(1, suppressWarnings(stats::predict(B, g))) %*% fit$coefficients)
      }
      m <- predict_at(grid)
      h <- 1e-4 * diff(range(dd))
      deriv <- (predict_at(pmin(grid + h, max(dd))) - predict_at(pmax(grid - h, min(dd)))) /
        (pmin(grid + h, max(dd)) - pmax(grid - h, min(dd)))
      list(att = m - dy0m, acrt = deriv)
    }
    point <- fit_curve(treated$dose, treated$dy, m0)
    boot <- .cm_with_seed(seed, replicate(n_boot, {
      idx <- sample(n, replace = TRUE)
      ub <- u[idx, ]
      tb <- ub[ub$dose > 0, ]
      if (nrow(tb) < df + 2L || sum(ub$dose == 0) < 2L) return(c(rep(NA_real_, 2L * length(grid))))
      r <- tryCatch(fit_curve(tb$dose, tb$dy, mean(ub$dy[ub$dose == 0])), error = function(e) NULL)
      if (is.null(r)) return(c(rep(NA_real_, 2L * length(grid))))
      c(r$att, r$acrt)
    }))
    boot <- boot[, colSums(is.na(boot)) == 0, drop = FALSE]
    L <- length(grid)
    se_att <- apply(boot[seq_len(L), , drop = FALSE], 1, stats::sd)
    se_acrt <- apply(boot[L + seq_len(L), , drop = FALSE], 1, stats::sd)
    sup_t <- function(est, se, draws) {
      tt <- abs(sweep(sweep(draws, 1, est, "-"), 1, se, "/"))
      as.numeric(stats::quantile(apply(tt, 2, max, na.rm = TRUE), conf_level, na.rm = TRUE))
    }
    cv_att <- sup_t(point$att, se_att, boot[seq_len(L), , drop = FALSE])
    cv_acrt <- sup_t(point$acrt, se_acrt, boot[L + seq_len(L), , drop = FALSE])
    curve <- data.frame(dose = grid, att = point$att, att_se = se_att,
                        att_band.low = point$att - cv_att * se_att, att_band.high = point$att + cv_att * se_att,
                        acrt = point$acrt, acrt_se = se_acrt,
                        acrt_band.low = point$acrt - cv_acrt * se_acrt, acrt_band.high = point$acrt + cv_acrt * se_acrt)
  }
  out <- list(by_dose = by_dose, acr = acr, curve = curve, twfe = beta, twfe_weights = twfe_weights,
              checks = checks, dose_type = dose_type, n = n, n_untreated = length(dy0),
              breaks = breaks, conf_level = conf_level, call = match.call())
  class(out) <- "cm_att_dose"
  out
}

#' @export
print.cm_att_dose <- function(x, digits = 4, ...) {
  cat("Difference-in-differences with a continuous treatment (Callaway, Goodman-Bacon, Sant'Anna)\n")
  cat("  Units: ", x$n, " (", x$n_untreated, " untreated); TWFE coefficient on the dose: ",
      formatC(x$twfe, digits = digits, format = "g"), "\n", sep = "")
  cat("Level effects ATT(d | d) by dose bin (parallel trends vs untreated):\n")
  print(format(x$by_dose[, c("bin", "dose", "n", "att", "std.error", "conf.low", "conf.high")], digits = digits), row.names = FALSE)
  cat("Causal responses between bins (strong parallel trends):\n")
  print(format(x$acr, digits = digits), row.names = FALSE)
  cat("TWFE weights: on level effects (w_level, can be negative) and on increments (w_slope, sum to one over dose gaps)\n")
  print(format(x$twfe_weights, digits = digits), row.names = FALSE)
  invisible(x)
}

#' Plot dose-response estimates
#'
#' @param x An object from [att_dose()].
#' @param what `"att"` (level effects) or `"acr"` (causal responses).
#' @return A `ggplot` object.
#' @export
plot_att_dose <- function(x, what = c("att", "acr")) {
  what <- match.arg(what)
  if (what == "att") {
    p <- ggplot2::ggplot(x$by_dose, ggplot2::aes(x = .data$dose, y = .data$att)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Dose", y = "ATT(d | d)")
    if (!is.null(x$curve)) {
      p <- p + ggplot2::geom_ribbon(data = x$curve, ggplot2::aes(x = .data$dose, ymin = .data$att_band.low, ymax = .data$att_band.high), inherit.aes = FALSE, alpha = 0.15) +
        ggplot2::geom_line(data = x$curve, ggplot2::aes(x = .data$dose, y = .data$att), colour = "steelblue")
    }
  } else {
    p <- ggplot2::ggplot(x$acr, ggplot2::aes(x = .data$to, y = .data$acr)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Dose", y = "Causal response (per unit of dose)")
    if (!is.null(x$curve)) {
      p <- p + ggplot2::geom_ribbon(data = x$curve, ggplot2::aes(x = .data$dose, ymin = .data$acrt_band.low, ymax = .data$acrt_band.high), inherit.aes = FALSE, alpha = 0.15) +
        ggplot2::geom_line(data = x$curve, ggplot2::aes(x = .data$dose, y = .data$acrt), colour = "steelblue")
    }
  }
  p + ggplot2::theme_minimal(base_size = 11)
}
