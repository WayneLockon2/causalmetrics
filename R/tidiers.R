# R/tidiers.R
#
# broom-style tidy() and glance() methods so causalmetrics estimators drop
# into modelsummary tables. The generics come from the generics package,
# which broom (in Depends) re-exports; the methods are registered against
# that generic without re-exporting it.

#' Tidy and glance methods for causalmetrics estimators
#'
#' `tidy()` returns the one-row coefficient table (the treatment effect) and
#' `glance()` returns one-row model information for [est_aipw()] and
#' [est_dml()] objects. Both follow the `broom` conventions used by
#' `modelsummary`.
#'
#' @param x A `cm_aipw`, `cm_dml`, `cm_scores`, `cm_blp`, `cm_gate`, `cm_cate`,
#'   `cm_cate_score`, or `cm_policy` object.
#' @param conf.int Logical. Include confidence limits (default `TRUE`).
#' @param conf.level Confidence level; defaults to the level stored in the
#'   object and is recomputed from the standard error otherwise.
#' @param ... Unused.
#'
#' @return `tidy()`: a data frame with columns `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, and (if requested) `conf.low`,
#'   `conf.high`. `glance()`: a one-row data frame with `nobs` and
#'   method information.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x1 <- rnorm(n)
#' p <- plogis(0.5 * x1)
#' d <- rbinom(n, 1, p)
#' y <- 1 + x1 + 2 * d + rnorm(n)
#' dat <- data.frame(y = y, d = d, x1 = x1, p = p, mu0 = 1 + x1, mu1 = 3 + x1)
#' fit <- est_aipw(dat, "y", "d", p_hat = "p", mu0_hat = "mu0", mu1_hat = "mu1")
#' tidy(fit)
#' glance(fit)
#' @name cm_tidiers
NULL

#' @importFrom generics tidy glance
NULL

.cm_tidy_effect <- function(x, conf.int, conf.level) {
  statistic <- x$estimate / x$std.error
  out <- data.frame(
    term = x$treatment %||% "treatment",
    estimate = x$estimate,
    std.error = x$std.error,
    statistic = statistic,
    p.value = 2 * stats::pnorm(-abs(statistic)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    if (is.null(conf.level) || isTRUE(all.equal(conf.level, x$conf.level))) {
      out$conf.low <- x$conf.low
      out$conf.high <- x$conf.high
    } else {
      z <- stats::qnorm(1 - (1 - conf.level) / 2)
      out$conf.low <- x$estimate - z * x$std.error
      out$conf.high <- x$estimate + z * x$std.error
    }
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @rdname cm_tidiers
#' @export
tidy.cm_aipw <- function(x, conf.int = TRUE, conf.level = x$conf.level, ...) {
  .cm_tidy_effect(x, conf.int, conf.level)
}

#' @rdname cm_tidiers
#' @export
glance.cm_aipw <- function(x, ...) {
  data.frame(
    nobs = x$n,
    n_treated = x$n_treated,
    n_control = x$n_control,
    estimand = x$estimand,
    method = "AIPW",
    folds = x$diagnostics$call$folds,
    prediction_mode = x$diagnostics$nuisance$prediction_mode,
    stringsAsFactors = FALSE
  )
}

#' @rdname cm_tidiers
#' @export
tidy.cm_dml <- function(x, conf.int = TRUE, conf.level = x$conf.level, ...) {
  .cm_tidy_effect(x, conf.int, conf.level)
}

#' @rdname cm_tidiers
#' @export
glance.cm_dml <- function(x, ...) {
  q <- x$diagnostics$nuisance
  rmse <- if (!is.null(q)) stats::setNames(as.list(q$rmse), paste0("rmse_", q$nuisance)) else list()
  out <- data.frame(
    nobs = x$n,
    model = x$model,
    estimand = x$estimand,
    method = paste0("DML (", x$score_type, ")"),
    solve = x$solve,
    folds = x$folds,
    n_rep = x$n_rep,
    learners = paste(sprintf("%s: %s", names(x$learners), unlist(x$learners)), collapse = "; "),
    stringsAsFactors = FALSE
  )
  for (nm in names(rmse)) out[[nm]] <- rmse[[nm]]
  out
}

# Heterogeneous effects and policy learning ------------------------------------

#' @rdname cm_tidiers
#' @export
tidy.cm_scores <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  out <- data.frame(term = "ATE", estimate = x$ate$estimate, std.error = x$ate$std.error,
                    statistic = x$ate$estimate / x$ate$std.error,
                    p.value = 2 * stats::pnorm(-abs(x$ate$estimate / x$ate$std.error)),
                    stringsAsFactors = FALSE)
  if (conf.int) {
    out$conf.low <- out$estimate - crit * out$std.error
    out$conf.high <- out$estimate + crit * out$std.error
  }
  out
}

#' @rdname cm_tidiers
#' @export
tidy.cm_blp <- function(x, conf.int = TRUE, conf.level = x$conf_level, ...) {
  out <- x$coefficients[, c("term", "estimate", "std.error", "statistic", "p.value")]
  if (conf.int) {
    crit <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - crit * out$std.error
    out$conf.high <- out$estimate + crit * out$std.error
  }
  out
}

#' @rdname cm_tidiers
#' @export
tidy.cm_gate <- function(x, ...) {
  out <- x$table
  names(out)[names(out) == "group"] <- "term"
  out
}

#' @rdname cm_tidiers
#' @export
tidy.cm_cate <- function(x, ...) {
  data.frame(row = seq_along(x$tau_hat), tau_hat = x$tau_hat)
}

#' @rdname cm_tidiers
#' @export
tidy.cm_cate_score <- function(x, ...) {
  out <- x$table
  names(out)[names(out) == "model"] <- "term"
  out
}

#' @rdname cm_tidiers
#' @export
tidy.cm_policy <- function(x, ...) {
  x$value
}

#' @rdname cm_tidiers
#' @export
glance.cm_scores <- function(x, ...) {
  data.frame(nobs = x$n, treated_share = mean(x$d), folds = length(unique(x$fold_id)),
             learner_p = x$learners$p, learner_mu0 = x$learners$mu0, learner_mu1 = x$learners$mu1,
             stringsAsFactors = FALSE)
}
