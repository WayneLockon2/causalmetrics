# R/cate-blp.R
#
# Inference on low-dimensional summaries of the CATE: the best linear
# predictor in a dictionary of covariate functions, and group average
# treatment effects, both with pointwise and simultaneous confidence bands.

#' Best linear predictor of the CATE
#'
#' Regresses the doubly robust pseudo-outcome of [dr_scores()] on a
#' dictionary `p(X)` given by `formula`, yielding the best linear
#' approximation `p(x)' beta_0` to the CATE (Semenova and Chernozhukov 2021;
#' Chernozhukov et al. 2026, ch. 14). Because the pseudo-outcome is
#' conditionally Neyman orthogonal, ordinary least squares with a
#' heteroskedasticity-robust (HC1) covariance is valid as if the nuisances
#' were known, provided the product of their errors vanishes faster than
#' `n^-1/2`. Simultaneous bands for the whole curve use the multiplier
#' bootstrap on the influence functions.
#'
#' @param scores A `cm_scores` object.
#' @param formula One-sided formula in the covariates of `scores$data` (or
#'   `data`), e.g. `~ x1 + I(x1^2)` or `~ 0 + group` for group effects.
#'   `NULL` fits an intercept only, which returns the ATE.
#' @param data Optional data frame with `nrow(scores$data)` rows providing
#'   extra regressors (for instance CATE predictions from another sample).
#' @param conf_level Confidence level.
#' @param uniform Compute the simultaneous critical value.
#' @param n_boot Multiplier bootstrap draws.
#' @param seed Optional seed for the bootstrap.
#'
#' @return A list of class `cm_blp`: `coefficients` (term, estimate,
#'   std.error, statistic, p.value, conf.low, conf.high, band.low,
#'   band.high), `vcov`, `crit_val` (simultaneous), `inffunc`, the design
#'   information, `r.squared`, `adj.r.squared`, and a Wald test that all
#'   non-intercept coefficients are zero (`test_heterogeneity`). [predict()]
#'   evaluates the curve on new data with pointwise and simultaneous bands;
#'   [plot_cate_blp()] draws it. [tidy()] and [glance()] make the object
#'   usable with `modelsummary`.
#'
#' @references
#' Semenova, V. and Chernozhukov, V. (2021). Debiased machine learning of
#' conditional average treatment effects and other causal functions.
#' *The Econometrics Journal*, 24(2), 264-289.
#'
#' @examples
#' dat <- sim_hte(1000, dgp = "smooth", seed = 1)
#' sc <- dr_scores(dat, "y", "d", paste0("x", 1:5), seed = 1)
#' blp <- cate_blp(sc, ~ x1 + I(x2^2), seed = 1)
#' blp
#' grid <- data.frame(x1 = seq(-2, 2, by = 0.5), x2 = 0)
#' predict(blp, grid)
#' @seealso [cate_gate()], [plot_cate_blp()]
#' @export
cate_blp <- function(scores, formula = NULL, data = NULL, conf_level = 0.95,
                     uniform = TRUE, n_boot = 999L, seed = NULL) {
  .cm_check_scores(scores)
  df <- .cm_blp_frame(scores, data)
  if (is.null(formula)) formula <- ~ 1
  if (!inherits(formula, "formula")) stop("`formula` must be a one-sided formula.", call. = FALSE)
  formula <- stats::as.formula(formula)
  if (length(formula) == 3L) stop("`formula` must be one-sided, e.g. `~ x1 + x2`.", call. = FALSE)
  X <- stats::model.matrix(formula, df)
  if (nrow(X) != nrow(df)) stop("Missing values in the regressors; complete the data first.", call. = FALSE)
  fit <- .cm_ols_if(X, scores$score, conf_level = conf_level, uniform = uniform,
                    n_boot = n_boot, seed = seed)
  est <- fit$beta
  se <- fit$se
  z <- est / se
  crit <- fit$crit
  cu <- if (is.finite(fit$crit_unif)) fit$crit_unif else crit
  coefs <- data.frame(
    term = colnames(X), estimate = est, std.error = se, statistic = z,
    p.value = 2 * stats::pnorm(-abs(z)),
    conf.low = est - crit * se, conf.high = est + crit * se,
    band.low = est - cu * se, band.high = est + cu * se,
    stringsAsFactors = FALSE
  )
  rownames(coefs) <- NULL
  slope <- which(colnames(X) != "(Intercept)")
  test <- if (length(slope) > 0L && length(slope) < ncol(X)) {
    .cm_if_wald(est[slope], fit$inffunc[, slope, drop = FALSE], fit$n)
  } else if (length(slope) == ncol(X) && ncol(X) > 1L) {
    # no intercept (group dummies): test equality of the coefficients
    C <- cbind(-1, diag(ncol(X) - 1L))
    .cm_if_wald(as.numeric(C %*% est), fit$inffunc %*% t(C), fit$n)
  } else {
    list(statistic = NA_real_, df = 0L, p.value = NA_real_)
  }
  structure(list(
    coefficients = coefs, beta = est, vcov = fit$vcov, inffunc = fit$inffunc,
    Qinv = fit$Qinv, crit_val = cu, crit_pointwise = crit, conf_level = conf_level,
    formula = formula, terms = stats::terms(formula), xlevels = .cm_xlevels(formula, df),
    test_heterogeneity = test, n = fit$n, n_boot = n_boot, seed = seed,
    residual_sd = stats::sd(fit$residuals),
    r.squared = 1 - sum(fit$residuals^2) / sum((scores$score - mean(scores$score))^2),
    adj.r.squared = 1 - (sum(fit$residuals^2) / max(fit$n - ncol(X), 1)) /
      (sum((scores$score - mean(scores$score))^2) / (fit$n - 1)),
    df = ncol(X), call = match.call()
  ), class = "cm_blp")
}

.cm_blp_frame <- function(scores, data = NULL) {
  df <- scores$data
  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (nrow(data) != nrow(df)) stop("`data` must have the same rows as `scores$data`.", call. = FALSE)
    for (v in setdiff(names(data), names(df))) df[[v]] <- data[[v]]
    for (v in intersect(names(data), names(df))) df[[v]] <- data[[v]]
  }
  df
}

.cm_xlevels <- function(formula, df) {
  vars <- all.vars(formula)
  lv <- lapply(vars, function(v) if (is.factor(df[[v]]) || is.character(df[[v]])) levels(factor(df[[v]])) else NULL)
  names(lv) <- vars
  lv[!vapply(lv, is.null, logical(1))]
}

#' @export
predict.cm_blp <- function(object, newdata, uniform = TRUE, n_boot = NULL, seed = NULL, ...) {
  newdata <- as.data.frame(newdata)
  for (v in names(object$xlevels)) newdata[[v]] <- factor(newdata[[v]], levels = object$xlevels[[v]])
  Xn <- stats::model.matrix(stats::delete.response(object$terms), newdata)
  if (ncol(Xn) != length(object$beta)) stop("`newdata` does not match the regressors of the model.", call. = FALSE)
  est <- as.numeric(Xn %*% object$beta)
  se <- sqrt(rowSums((Xn %*% object$vcov) * Xn))
  crit <- object$crit_pointwise
  cu <- crit
  if (uniform && nrow(Xn) > 1L) {
    if_pred <- object$inffunc %*% t(Xn)
    sup <- .cm_sup_crit(if_pred, object$n, n_boot = n_boot %||% object$n_boot,
                        conf_level = object$conf_level, seed = seed %||% object$seed)
    if (is.finite(sup$crit)) cu <- sup$crit
  }
  out <- data.frame(estimate = est, std.error = se,
                    conf.low = est - crit * se, conf.high = est + crit * se,
                    band.low = est - cu * se, band.high = est + cu * se)
  attr(out, "crit_val") <- cu
  cbind(newdata, out)
}

#' @export
print.cm_blp <- function(x, ...) {
  cat("Best linear predictor of the CATE: score ~ ", deparse(x$formula[[2]]), "\n", sep = "")
  cat("  n = ", x$n, ", ", 100 * x$conf_level, "% pointwise |z| = ",
      format(round(x$crit_pointwise, 3)), ", simultaneous = ", format(round(x$crit_val, 3)), "\n", sep = "")
  print(x$coefficients[, c("term", "estimate", "std.error", "conf.low", "conf.high", "band.low", "band.high")],
        digits = 4, row.names = FALSE)
  if (is.finite(x$test_heterogeneity$p.value)) {
    cat("  Wald test of no heterogeneity: chi2(", x$test_heterogeneity$df, ") = ",
        format(round(x$test_heterogeneity$statistic, 2)), ", p = ",
        format(signif(x$test_heterogeneity$p.value, 3)), "\n", sep = "")
  }
  invisible(x)
}

#' Group average treatment effects
#'
#' Doubly robust group average treatment effects: the best linear predictor
#' of [cate_blp()] with group indicators as the dictionary, so each
#' coefficient is the mean of the pseudo-outcome in its group. Returns
#' pointwise and simultaneous intervals and a Wald test that all groups share
#' the same effect.
#'
#' @param scores A `cm_scores` object.
#' @param groups A column name in `scores$data` or a vector of length
#'   `nrow(scores$data)` defining the groups (factor, character, or numeric).
#' @inheritParams cate_blp
#'
#' @return A list of class `cm_gate` with `table` (group, n, share, estimate,
#'   std.error, conf.low, conf.high, band.low, band.high), `test_equal`,
#'   `crit_val`, and the underlying `cm_blp` object.
#' @examples
#' dat <- sim_hte(1000, dgp = "smooth", seed = 1)
#' sc <- dr_scores(dat, "y", "d", paste0("x", 1:5), seed = 1)
#' g <- cut(dat$x1, quantile(dat$x1, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
#' cate_gate(sc, g, seed = 1)
#' @seealso [cate_blp()], [plot_cate_gate()]
#' @export
cate_gate <- function(scores, groups, conf_level = 0.95, n_boot = 999L, seed = NULL) {
  .cm_check_scores(scores)
  g <- if (.cm_is_string(groups)) {
    .cm_check_column(groups, scores$data)
    scores$data[[groups]]
  } else groups
  if (length(g) != scores$n) stop("`groups` must have length nrow(scores$data).", call. = FALSE)
  g <- factor(g)
  if (nlevels(g) < 2L) stop("`groups` must define at least two groups.", call. = FALSE)
  blp <- cate_blp(scores, ~ 0 + .cm_group, data = data.frame(.cm_group = g),
                  conf_level = conf_level, uniform = TRUE, n_boot = n_boot, seed = seed)
  tab <- blp$coefficients
  tab$term <- levels(g)
  names(tab)[1] <- "group"
  tab$n <- as.integer(table(g))
  tab$share <- tab$n / scores$n
  tab <- tab[, c("group", "n", "share", "estimate", "std.error", "conf.low", "conf.high", "band.low", "band.high")]
  structure(list(table = tab, test_equal = blp$test_heterogeneity, crit_val = blp$crit_val,
                 conf_level = conf_level, blp = blp, groups = g), class = "cm_gate")
}

#' @export
print.cm_gate <- function(x, ...) {
  cat("Group average treatment effects (doubly robust)\n")
  print(x$table, digits = 4, row.names = FALSE)
  cat("  Wald test of equal effects: chi2(", x$test_equal$df, ") = ",
      format(round(x$test_equal$statistic, 2)), ", p = ", format(signif(x$test_equal$p.value, 3)), "\n", sep = "")
  invisible(x)
}

#' Plot a best linear predictor of the CATE along one covariate
#'
#' @param x A `cm_blp` object.
#' @param newdata Data frame on which to evaluate the curve; the column
#'   `x_var` is placed on the horizontal axis.
#' @param x_var Name of the covariate on the horizontal axis (default: the
#'   first column of `newdata`).
#' @param ... Passed to the `predict()` method of `cm_blp` objects (for instance `seed`).
#' @return A ggplot object.
#' @export
plot_cate_blp <- function(x, newdata, x_var = names(newdata)[1], ...) {
  pr <- predict(x, newdata, ...)
  pr$x_var <- pr[[x_var]]
  ggplot2::ggplot(pr, ggplot2::aes(x = .data$x_var, y = .data$estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$band.low, ymax = .data$band.high), fill = "grey80", alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), fill = "grey60", alpha = 0.6) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(x = x_var, y = "Best linear predictor of the CATE",
                  subtitle = "Dark band: pointwise; light band: simultaneous") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot group average treatment effects
#'
#' @param x A `cm_gate` object.
#' @return A ggplot object with pointwise (thick) and simultaneous (thin)
#'   intervals.
#' @export
plot_cate_gate <- function(x) {
  tab <- x$table
  tab$group <- factor(tab$group, levels = tab$group)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$group, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$band.low, ymax = .data$band.high), width = 0.25, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0.12, linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = NULL, y = "Group average treatment effect",
                  subtitle = "Thick: pointwise; thin: simultaneous") +
    ggplot2::theme_minimal(base_size = 11)
}
