# R/gelbach-decomp.R
#
# Gelbach's (2016) decomposition of the change in a treatment coefficient
# between a base regression and a regression with added covariates.

#' Gelbach decomposition of a coefficient change
#'
#' When covariates `x_add` are added to a base regression `Y ~ D + X_base`,
#' the coefficient on `D` moves from `beta_base` to `beta_full`. Gelbach
#' (2016) shows that the change decomposes exactly, and independently of
#' the order in which covariates are added, as
#' \deqn{\hat\beta_{base} - \hat\beta_{full} = \sum_k \hat\Gamma_k \hat\beta_k^{full},}
#' where `Gamma_k` is the coefficient on `D` in the auxiliary regression of
#' the added covariate `k` on `D` and `X_base`, and `beta_k^full` is the
#' coefficient of covariate `k` in the full regression. Each term is the
#' part of the coefficient change attributable to covariate `k`. With
#' post-treatment covariates (candidate mediators) the decomposition is a
#' descriptive accounting of how much of the treatment coefficient the
#' mediators absorb, not an estimate of natural effects.
#'
#' @param data A data frame.
#' @param y,d Outcome and treatment column names (the treatment may be
#'   continuous).
#' @param x_base Covariates in the base regression.
#' @param x_add Covariates added in the full regression.
#' @param groups Optional named list assigning the added covariates to
#'   groups whose contributions are summed.
#' @param fe Optional fixed-effect column names (regressions run with
#'   `fixest::feols`).
#' @param cluster Optional cluster column for the bootstrap.
#' @param weights Optional column of regression weights (every regression,
#'   including the auxiliary ones, is weighted; the identity still holds
#'   exactly).
#' @param n_boot Bootstrap replications for the standard errors (0 skips
#'   inference).
#' @param conf_level Confidence level.
#' @param seed Optional seed.
#'
#' @return A list of class `cm_gelbach` with `coefficients` (base, full, and
#'   their difference), `contributions` (covariate or group, `Gamma`,
#'   `beta_full`, `contribution`, share of the change, and bootstrap
#'   intervals), and the fitted models. [plot_gelbach()] draws the
#'   contributions.
#'
#' @references
#' Gelbach, J. B. (2016). When do covariates matter? And which ones, and how
#' much? *Journal of Labor Economics*, 34(2), 509-543.
#'
#' @examples
#' dat <- sim_mediation(1000, dgp = "parallel", seed = 1)
#' g <- gelbach_decomp(dat, "y", "d", x_base = c("x1", "x2"), x_add = c("m1", "m2"), n_boot = 99)
#' g
#' @export
gelbach_decomp <- function(data, y, d, x_base = NULL, x_add, groups = NULL, fe = NULL,
                           cluster = NULL, weights = NULL, n_boot = 499L, conf_level = 0.95, seed = NULL) {
  data <- as.data.frame(data)
  for (v in c(y, d, x_base, x_add, fe, cluster, weights)) .cm_check_column(v, data)
  keep <- stats::complete.cases(data[, c(y, d, x_base, x_add, fe, cluster, weights)])
  data <- data[keep, , drop = FALSE]
  wf <- if (is.null(weights)) NULL else stats::as.formula(paste0("~", weights))
  if (!is.null(groups)) {
    miss <- setdiff(unlist(groups), x_add)
    if (length(miss)) stop("`groups` names covariates not in `x_add`: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  fml <- function(lhs, rhs) {
    f <- paste(lhs, "~", paste(rhs, collapse = " + "))
    if (!is.null(fe)) f <- paste(f, "|", paste(fe, collapse = " + "))
    stats::as.formula(f)
  }
  compute <- function(dd) {
    base <- fixest::feols(fml(y, c(d, x_base)), data = dd, weights = wf, warn = FALSE, notes = FALSE)
    full <- fixest::feols(fml(y, c(d, x_base, x_add)), data = dd, weights = wf, warn = FALSE, notes = FALSE)
    b_base <- stats::coef(base)[[d]]
    b_full <- stats::coef(full)[[d]]
    Gamma <- vapply(x_add, function(k) stats::coef(fixest::feols(fml(k, c(d, x_base)), data = dd, weights = wf, warn = FALSE, notes = FALSE))[[d]], numeric(1))
    beta_k <- stats::coef(full)[x_add]
    contrib <- Gamma * beta_k
    out <- c(base = b_base, full = b_full, change = b_base - b_full, contrib)
    if (!is.null(groups)) {
      out <- c(out, vapply(groups, function(g) sum(contrib[g]), numeric(1)))
    }
    list(stat = out, Gamma = Gamma, beta_k = beta_k, models = list(base = base, full = full))
  }
  point <- compute(data)
  est <- point$stat
  if (n_boot > 0) {
    cl <- if (is.null(cluster)) NULL else data[[cluster]]
    boot <- .cm_boot(data, function(dd) compute(dd)$stat, n_boot, cluster = cl, seed = seed, conf_level = conf_level)
  } else {
    boot <- data.frame(term = names(est), estimate = as.numeric(est), std.error = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  }
  coefs <- boot[boot$term %in% c("base", "full", "change"), ]
  contributions <- boot[boot$term %in% x_add, ]
  contributions$Gamma <- point$Gamma[contributions$term]
  contributions$beta_full <- point$beta_k[contributions$term]
  contributions$share <- contributions$estimate / est[["change"]]
  names(contributions)[names(contributions) == "term"] <- "covariate"
  names(contributions)[names(contributions) == "estimate"] <- "contribution"
  contributions <- contributions[, c("covariate", "Gamma", "beta_full", "contribution", "std.error", "conf.low", "conf.high", "share")]
  group_table <- NULL
  if (!is.null(groups)) {
    group_table <- boot[boot$term %in% names(groups), ]
    names(group_table)[1:2] <- c("group", "contribution")
    group_table$share <- group_table$contribution / est[["change"]]
  }
  identity_gap <- est[["change"]] - sum(point$Gamma * point$beta_k)
  structure(list(coefficients = coefs, contributions = contributions, groups = group_table,
                 identity_gap = identity_gap, models = point$models, n = nrow(data),
                 spec = list(y = y, d = d, x_base = x_base, x_add = x_add, fe = fe, weights = weights), conf_level = conf_level,
                 call = match.call()), class = "cm_gelbach")
}

#' @export
print.cm_gelbach <- function(x, ...) {
  cat("Gelbach decomposition of the coefficient on ", x$spec$d, " (n = ", x$n, ")\n", sep = "")
  print(x$coefficients, digits = 4, row.names = FALSE)
  cat("Contributions of the added covariates (they sum to the change; identity gap ",
      format(signif(x$identity_gap, 2)), "):\n", sep = "")
  print(x$contributions, digits = 4, row.names = FALSE)
  if (!is.null(x$groups)) {
    cat("Group contributions:\n")
    print(x$groups, digits = 4, row.names = FALSE)
  }
  invisible(x)
}

#' Plot a Gelbach decomposition
#'
#' @param x A `cm_gelbach` object.
#' @param groups Draw the group contributions instead of the covariate ones.
#' @return A ggplot object: contributions with bootstrap intervals.
#' @export
plot_gelbach <- function(x, groups = FALSE) {
  df <- if (groups && !is.null(x$groups)) {
    data.frame(covariate = x$groups$group, contribution = x$groups$contribution,
               conf.low = x$groups$conf.low, conf.high = x$groups$conf.high)
  } else {
    x$contributions[, c("covariate", "contribution", "conf.low", "conf.high")]
  }
  df$covariate <- factor(df$covariate, levels = rev(df$covariate))
  ggplot2::ggplot(df, ggplot2::aes(y = .data$covariate, x = .data$contribution)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_col(fill = "#0072B2", alpha = 0.7, width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high), width = 0.2) +
    ggplot2::labs(x = paste0("Contribution to the change in the ", x$spec$d, " coefficient (total ",
                             format(round(x$coefficients$estimate[x$coefficients$term == "change"], 3)), ")"), y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}
