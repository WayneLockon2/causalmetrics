# R/ssiv.R
#
# Shift-share instrument diagnostics: Rotemberg weights (exogenous-shares
# view) and the shock-level regression (exogenous-shocks view).

#' Rotemberg weights of a shift-share instrument
#'
#' A shift-share (Bartik) instrument `B_i = sum_k s_ik g_k` combines unit
#' exposure shares `s_ik` with shocks `g_k`. Goldsmith-Pinkham, Sorkin, and
#' Swift (2020) show that the two-stage least squares coefficient equals a
#' weighted average of the just-identified estimates that use each share as
#' the instrument, `beta = sum_k alpha_k beta_k`, with Rotemberg weights
#' `alpha_k = g_k (Z_k' D^perp) / sum_j g_j (Z_j' D^perp)`, where `Z_k` is
#' the residualized share and `D^perp` the residualized treatment. The
#' weights say which industries (shocks) drive the estimate and are the
#' object to inspect when the identifying assumption is exogeneity of the
#' shares.
#'
#' @param data A data frame with one row per unit.
#' @param y,d Outcome and treatment column names.
#' @param shares Character vector of share column names (one per shock).
#' @param shocks Numeric vector of shocks, in the order of `shares`.
#' @param x Optional control column names.
#' @param weights Optional weights column name.
#' @return A list of class `cm_ssiv_rotemberg` with `table` (shock,
#'   share mean, shock value, `alpha_k`, `beta_k`, first-stage t), the 2SLS
#'   `estimate` and `std.error`, `check` (`sum(alpha_k beta_k)`), and
#'   `summary` (share of negative weights, sum of negative weights, top
#'   five weights, correlation of weights with shocks).
#' @references
#' Goldsmith-Pinkham, P., Sorkin, I., and Swift, H. (2020). Bartik
#' instruments: what, when, why, and how. *American Economic Review*,
#' 110(8), 2586-2624.
#' @examples
#' dat <- sim_iv(400, dgp = "shift_share", n_industries = 10, seed = 1)
#' ssiv_rotemberg(dat, "y", "d", shares = paste0("s", 1:10), shocks = attr(dat, "shocks"), x = "x1")
#' @export
ssiv_rotemberg <- function(data, y, d, shares, shocks, x = NULL, weights = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, shares, x, weights)) .cm_check_column(v, data)
  if (length(shocks) != length(shares)) stop("`shocks` must have one value per share column.", call. = FALSE)
  keep <- stats::complete.cases(data[, c(y, d, shares, x, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  n <- nrow(data)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(data[[weights]])
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  S <- as.matrix(data[, shares, drop = FALSE])
  B <- as.numeric(S %*% shocks)
  ts <- .cm_tsls(as.numeric(data[[y]]), as.numeric(data[[d]]), matrix(B, ncol = 1), X, w)
  yt <- ts$residuals$y; dt <- ts$residuals$d
  St <- .cm_residualize(S, X, w)
  zd <- colSums(w * St * dt)
  zy <- colSums(w * St * yt)
  alpha <- shocks * zd / sum(shocks * zd)
  beta_k <- zy / zd
  fs <- vapply(seq_along(shares), function(k) {
    f <- .cm_first_stage_strength(dt, St[, k, drop = FALSE], w)
    f$t
  }, numeric(1))
  tab <- data.frame(shock = shares, share_mean = colMeans(S), shock_value = shocks,
                    alpha_k = unname(alpha), beta_k = unname(beta_k), first_stage_t = fs)
  rownames(tab) <- NULL
  ord <- order(-abs(tab$alpha_k))
  summ <- list(share_negative = mean(alpha < 0), sum_negative = sum(alpha[alpha < 0]),
               top_five = tab[ord[seq_len(min(5L, nrow(tab)))], c("shock", "alpha_k", "beta_k")],
               cor_alpha_shock = suppressWarnings(stats::cor(alpha, shocks)),
               herfindahl = sum(alpha^2))
  structure(list(table = tab, estimate = ts$estimate, std.error = ts$std.error, check = sum(alpha * beta_k),
                 summary = summ, n = n, instrument = B), class = "cm_ssiv_rotemberg")
}

#' @export
print.cm_ssiv_rotemberg <- function(x, ...) {
  cat("Shift-share IV: 2SLS estimate = ", format(round(x$estimate, 4)), " (SE ", format(round(x$std.error, 4)),
      "), n = ", x$n, "; sum(alpha_k beta_k) = ", format(round(x$check, 4)), "\n", sep = "")
  cat("  share of negative Rotemberg weights = ", format(round(x$summary$share_negative, 3)),
      "; sum of negative weights = ", format(round(x$summary$sum_negative, 3)),
      "; Herfindahl of weights = ", format(round(x$summary$herfindahl, 3)), "\n", sep = "")
  cat("  top weights:\n"); print(x$summary$top_five, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Shock-level regression for a shift-share instrument
#'
#' Borusyak, Hull, and Jaravel (2022) show that the shift-share IV estimate
#' equals a shock-level instrumental-variables regression: aggregate the
#' residualized outcome and treatment to the shock level with exposure
#' weights `s_k = sum_i w_i s_ik`, `ybar_k = sum_i w_i s_ik y_i / s_k`, and
#' regress `ybar_k` on `dbar_k` instrumented by `g_k` with weights `s_k`
#' (no intercept when shares sum to one; otherwise the sum of shares is
#' added to the controls). Inference at the shock level is valid when the
#' shocks are as good as randomly assigned, which is the exogenous-shocks
#' identification story.
#'
#' @inheritParams ssiv_rotemberg
#' @param cluster Optional vector (length of `shares`) assigning shocks to
#'   clusters for the shock-level standard errors.
#' @return A list of class `cm_ssiv_shock` with `estimate`, `std.error`
#'   (heteroskedasticity-robust at the shock level), `unit_level` (the
#'   unit-level 2SLS estimate, which should coincide), the shock-level
#'   `table`, and the number of shocks.
#' @references
#' Borusyak, K., Hull, P., and Jaravel, X. (2022). Quasi-experimental
#' shift-share research designs. *Review of Economic Studies*, 89(1),
#' 181-213.
#' @examples
#' dat <- sim_iv(400, dgp = "shift_share", n_industries = 10, seed = 1)
#' ssiv_shock_level(dat, "y", "d", shares = paste0("s", 1:10), shocks = attr(dat, "shocks"), x = "x1")
#' @export
ssiv_shock_level <- function(data, y, d, shares, shocks, x = NULL, weights = NULL, cluster = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, shares, x, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, shares, x, weights), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  n <- nrow(data)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(data[[weights]])
  S <- as.matrix(data[, shares, drop = FALSE])
  share_sum <- rowSums(S)
  incomplete <- any(abs(share_sum - 1) > 1e-8)
  X <- if (is.null(x)) NULL else as.matrix(data[, x, drop = FALSE])
  if (incomplete) X <- cbind(X, share_sum = share_sum)
  B <- as.numeric(S %*% shocks)
  unit <- .cm_tsls(as.numeric(data[[y]]), as.numeric(data[[d]]), matrix(B, ncol = 1), X, w)
  yt <- unit$residuals$y; dt <- unit$residuals$d
  s_k <- colSums(w * S)
  ybar <- colSums(w * S * yt) / s_k
  dbar <- colSums(w * S * dt) / s_k
  # shock-level IV without intercept, weights s_k
  est <- sum(s_k * shocks * ybar) / sum(s_k * shocks * dbar)
  e_k <- ybar - est * dbar
  psi <- s_k * shocks * e_k
  jac <- sum(s_k * shocks * dbar)
  K <- length(shocks)
  if (is.null(cluster)) {
    v <- sum(psi^2) * K / (K - 1) / jac^2
    G <- K
  } else {
    cs <- rowsum(psi, as.integer(factor(cluster)))
    G <- nrow(cs)
    v <- sum(cs^2) * G / (G - 1) / jac^2
  }
  tab <- data.frame(shock = shares, exposure = s_k, shock_value = shocks, ybar = ybar, dbar = dbar)
  rownames(tab) <- NULL
  structure(list(estimate = est, std.error = sqrt(v), unit_level = c(estimate = unit$estimate, std.error = unit$std.error),
                 table = tab, n_shocks = K, n_clusters = G, incomplete_shares = incomplete, n = n), class = "cm_ssiv_shock")
}

#' @export
print.cm_ssiv_shock <- function(x, ...) {
  cat("Shift-share IV, shock-level regression (", x$n_shocks, " shocks", if (x$n_clusters < x$n_shocks) paste0(", ", x$n_clusters, " clusters"), ")\n", sep = "")
  cat("  estimate = ", format(round(x$estimate, 4)), " (shock-level SE ", format(round(x$std.error, 4)), ")\n", sep = "")
  cat("  unit-level 2SLS: ", format(round(x$unit_level[1], 4)), " (SE ", format(round(x$unit_level[2], 4)), ")",
      if (x$incomplete_shares) "; incomplete shares: sum of shares added as a control", "\n", sep = "")
  invisible(x)
}
