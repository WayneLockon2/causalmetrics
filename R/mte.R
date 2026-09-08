# R/mte.R
#
# Marginal treatment effects from a continuous or multi-valued instrument:
# local instrumental variables on the propensity score, with the weights that
# turn the MTE curve into policy-relevant averages.

#' Marginal treatment effect curve and policy-relevant averages
#'
#' In the selection model `D = 1{p(Z, X) > V}` with `V ~ U(0, 1)`, the
#' marginal treatment effect `MTE(u, x) = E[Y(1) - Y(0) | X = x, V = u]` is
#' the derivative of `E[Y | X, p(Z, X) = p]` in `p` at `p = u` (local
#' instrumental variables; Heckman and Vytlacil 2005). The function estimates
#' the propensity by logistic regression of `D` on `(Z, X)`, then fits the
#' separable model `E[Y | X, p] = X'b0 + p X'(b1 - b0) + K(p)` with `K` a
#' polynomial in `p` (or a local quadratic fit of the residual on `p`), and
#' reports `MTE(u) = mean(X)'(b1 - b0) + K'(u)` on a grid. The curve is then
#' integrated against the Heckman-Vytlacil weights: `1` for the ATE,
#' `P(p > u) / E[p]` for the ATT, `P(p <= u) / E[1 - p]` for the ATU, the
#' uniform weight on `[p_low, p_high]` for the LATE of an instrument shift,
#' and `(F_p(u) - F_p'(u)) / (E[p'] - E[p])` for a policy that moves the
#' propensity from `p` to `p'`. Bootstrap replications re-estimate both
#' stages.
#'
#' @param data A data frame.
#' @param y,d Outcome and binary treatment column names.
#' @param z Character vector of instrument column names.
#' @param x Optional covariate column names.
#' @param method `"polynomial"` (default) or `"local"` (local quadratic
#'   derivative of the residualized outcome in `p`).
#' @param degree Polynomial degree of `K(p)`.
#' @param bandwidth Bandwidth for the local method (default
#'   `1.06 sd(p) n^-1/5`).
#' @param grid Values of `u` at which to report the curve.
#' @param late_range Two propensity values defining the LATE weight.
#' @param policy Optional function mapping the estimated propensity to the
#'   counterfactual propensity under a policy (for the PRTE).
#' @param n_boot Bootstrap replications (0 to skip).
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_mte` with `curve` (u, mte, std.error,
#'   conf.low, conf.high, and the support indicator), `effects` (ATE, ATT,
#'   ATU, LATE, and PRTE with standard errors and the share of each weight
#'   outside the support of the estimated propensity), `support` (range and
#'   quantiles of the propensity), the fitted models, and the weight
#'   functions on the grid. [plot_mte()] draws curve and weights.
#' @references
#' Heckman, J. J. and Vytlacil, E. (2005). Structural equations, treatment
#' effects, and econometric policy evaluation. *Econometrica*, 73(3),
#' 669-738.
#'
#' Cornelissen, T., Dustmann, C., Raute, A., and Schönberg, U. (2016). From
#' LATE to MTE: alternative methods for the evaluation of policy
#' interventions. *Labour Economics*, 41, 47-60.
#' @examples
#' dat <- sim_iv(3000, dgp = "mte", seed = 1)
#' m <- mte_curve(dat, "y", "d", "z", x = "x1", degree = 3, n_boot = 20, seed = 1)
#' m
#' @export
mte_curve <- function(data, y, d, z, x = NULL, method = c("polynomial", "local"), degree = 3L,
                      bandwidth = NULL, grid = seq(0.05, 0.95, by = 0.05), late_range = NULL,
                      policy = NULL, n_boot = 99L, conf_level = 0.95, seed = NULL) {
  method <- match.arg(method)
  data <- as.data.frame(data)
  for (v in c(y, d, z, x)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, z, x), drop = FALSE])
  data <- data[keep, , drop = FALSE]
  n <- nrow(data)
  degree <- .cm_check_count(degree, "degree", min = 1L)
  dv <- as.numeric(.cm_as_binary(data[[d]], d))
  yv <- as.numeric(data[[y]])
  fine <- seq(0.005, 0.995, by = 0.005)

  fit_once <- function(idx) {
    df <- data[idx, , drop = FALSE]
    df$.d <- dv[idx]; df$.y <- yv[idx]
    ps <- stats::glm(stats::reformulate(c(z, x), ".d"), data = df, family = stats::binomial())
    p <- pmin(pmax(stats::fitted(ps), 1e-4), 1 - 1e-4)
    Xm <- if (is.null(x)) NULL else as.matrix(df[, x, drop = FALSE])
    Xc <- if (is.null(x)) NULL else sweep(Xm, 2, colMeans(Xm))
    P <- outer(p, seq_len(degree), "^")
    colnames(P) <- paste0("p", seq_len(degree))
    design <- cbind(1, Xc, if (!is.null(x)) Xc * p, P)
    fit <- stats::lm.fit(design, df$.y)
    cf <- fit$coefficients; cf[is.na(cf)] <- 0
    kx <- if (is.null(x)) 0L else ncol(Xc)
    b_int <- if (kx > 0) cf[(1 + kx + 1):(1 + 2 * kx)] else numeric(0)
    alpha <- cf[(length(cf) - degree + 1):length(cf)]
    k_prime <- function(u) as.numeric(outer(u, seq_len(degree), function(uu, j) j * uu^(j - 1)) %*% alpha)
    if (method == "polynomial") {
      mte_fun <- function(u) k_prime(u)   # covariates centered, so mean(X)'(b1 - b0) = 0 at the mean
    } else {
      # local quadratic derivative of the outcome net of the covariate part
      resid <- df$.y - as.numeric(cbind(1, Xc, if (!is.null(x)) Xc * p) %*% cf[seq_len(1 + 2 * kx)])
      h <- if (is.null(bandwidth)) 1.06 * stats::sd(p) * n^(-1 / 5) else bandwidth
      mte_fun <- function(u) vapply(u, function(u0) {
        w <- stats::dnorm((p - u0) / h)
        if (sum(w > 1e-8) < 10) return(NA_real_)
        Xl <- cbind(1, p - u0, (p - u0)^2)
        stats::lm.wfit(Xl, resid, w)$coefficients[2]
      }, numeric(1))
    }
    list(p = p, mte = mte_fun, coef = cf, ps = ps)
  }
  weights_for <- function(p, u, pol) {
    Fp <- stats::ecdf(p)
    w <- list(ATE = rep(1, length(u)),
              ATT = (1 - Fp(u)) / mean(p),
              ATU = Fp(u) / mean(1 - p))
    if (!is.null(late_range)) {
      lo <- min(late_range); hi <- max(late_range)
      w$LATE <- as.numeric(u > lo & u <= hi) / (hi - lo)
    }
    if (!is.null(pol)) {
      p2 <- pmin(pmax(pol(p), 1e-4), 1 - 1e-4)
      Fp2 <- stats::ecdf(p2)
      w$PRTE <- (Fp(u) - Fp2(u)) / (mean(p2) - mean(p))
    }
    w
  }
  effects_from <- function(fit) {
    m <- fit$mte(fine)
    w <- weights_for(fit$p, fine, policy)
    du <- fine[2] - fine[1]
    inside <- fine >= min(fit$p) & fine <= max(fit$p)
    vapply(w, function(wk) sum(m * wk * du, na.rm = TRUE), numeric(1))
  }
  point <- fit_once(seq_len(n))
  curve <- point$mte(grid)
  eff <- effects_from(point)
  w_fine <- weights_for(point$p, fine, policy)
  du <- fine[2] - fine[1]
  inside <- fine >= min(point$p) & fine <= max(point$p)
  outside_share <- vapply(w_fine, function(wk) sum(abs(wk[!inside]) * du) / max(sum(abs(wk) * du), 1e-12), numeric(1))

  se_curve <- rep(NA_real_, length(grid)); se_eff <- rep(NA_real_, length(eff))
  if (n_boot > 0) {
    boots <- .cm_with_seed(seed, lapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      f <- tryCatch(fit_once(idx), error = function(e) NULL)
      if (is.null(f)) return(NULL)
      list(curve = f$mte(grid), eff = effects_from(f))
    }))
    boots <- boots[!vapply(boots, is.null, logical(1))]
    if (length(boots) > 1L) {
      se_curve <- apply(sapply(boots, `[[`, "curve"), 1, stats::sd, na.rm = TRUE)
      se_eff <- apply(sapply(boots, `[[`, "eff"), 1, stats::sd, na.rm = TRUE)
    }
  }
  zq <- stats::qnorm(1 - (1 - conf_level) / 2)
  curve_df <- data.frame(u = grid, mte = curve, std.error = se_curve, conf.low = curve - zq * se_curve,
                         conf.high = curve + zq * se_curve, in_support = grid >= min(point$p) & grid <= max(point$p))
  effects <- data.frame(estimand = names(eff), estimate = unname(eff), std.error = unname(se_eff),
                        weight_outside_support = unname(outside_share))
  effects$conf.low <- effects$estimate - zq * effects$std.error
  effects$conf.high <- effects$estimate + zq * effects$std.error
  structure(list(curve = curve_df, effects = effects,
                 support = list(range = range(point$p), quantiles = stats::quantile(point$p, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))),
                 weights = data.frame(u = fine, as.data.frame(w_fine)), propensity = point$p,
                 propensity_model = point$ps, coefficients = point$coef, method = method, degree = degree,
                 n = n, n_boot = n_boot, conf_level = conf_level), class = "cm_mte")
}

#' @export
print.cm_mte <- function(x, ...) {
  cat("Marginal treatment effects (", x$method, if (x$method == "polynomial") paste0(", degree ", x$degree), "), n = ", x$n, "\n", sep = "")
  cat("  propensity support: [", format(round(x$support$range[1], 3)), ", ", format(round(x$support$range[2], 3)), "]\n", sep = "")
  print(x$effects, digits = 3, row.names = FALSE)
  invisible(x)
}

#' Plot a marginal treatment effect curve and its policy weights
#'
#' @param x A `cm_mte` object.
#' @param what `"curve"` or `"weights"`.
#' @return A ggplot object.
#' @export
plot_mte <- function(x, what = c("curve", "weights")) {
  what <- match.arg(what)
  if (what == "curve") {
    df <- x$curve
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$u, y = .data$mte))
    if (any(is.finite(df$std.error))) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), fill = "grey80", alpha = 0.7)
    p + ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::annotate("rect", xmin = x$support$range[1], xmax = x$support$range[2], ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "steelblue") +
      ggplot2::labs(x = "Unobserved resistance to treatment (u)", y = "Marginal treatment effect",
                    subtitle = "Shaded: support of the estimated propensity score") +
      ggplot2::theme_minimal(base_size = 11)
  } else {
    w <- x$weights
    long <- do.call(rbind, lapply(setdiff(names(w), "u"), function(nm) data.frame(u = w$u, weight = w[[nm]], estimand = nm)))
    ggplot2::ggplot(long, ggplot2::aes(x = .data$u, y = .data$weight, colour = .data$estimand)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::labs(x = "u", y = "Weight on MTE(u)", colour = NULL) +
      ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  }
}
