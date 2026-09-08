# R/rd-frame.R
#
# One tidy table across the regression discontinuity packages, and tidy()
# and glance() methods for their objects so that modelsummary works.

#' Tidy rows from regression discontinuity objects
#'
#' Collects estimates from `rdrobust::rdrobust()`, `RDHonest::RDHonest()`,
#' `rdlocrand::rdrandinf()`, `rdmulti::rdmc()`, `rdhte::rdhte()`, and the
#' package's own `cm_rd_adjust` and `cm_rd_weak_iv` objects into one data
#' frame with common columns, so different inference procedures can be
#' compared in one table or figure.
#'
#' @param ... Named objects.
#' @param methods For `rdrobust` objects, which rows to keep (default all
#'   three).
#' @return A data frame with `model`, `term`, `method`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`, `p.value`, `h`, `n_eff`.
#' @examples
#' dat <- sim_rd(1500, "lee", seed = 1)
#' if (requireNamespace("rdrobust", quietly = TRUE)) {
#'   fit <- rdrobust::rdrobust(dat$y, dat$x, c = 0)
#'   rd_frame(rdrobust = fit)
#' }
#' @export
rd_frame <- function(..., methods = c("conventional", "bias_corrected", "robust")) {
  objs <- list(...)
  if (length(objs) == 0L) stop("Supply at least one object.", call. = FALSE)
  if (is.null(names(objs)) || any(names(objs) == "")) names(objs) <- paste0("model", seq_along(objs))
  rows <- lapply(names(objs), function(nm) {
    t <- .cm_rd_tidy_any(objs[[nm]])
    if (is.null(t)) return(NULL)
    t <- t[t$method %in% c(methods, setdiff(unique(t$method), c("conventional", "bias_corrected", "robust"))), , drop = FALSE]
    cbind(model = nm, t)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.cm_rd_tidy_any <- function(obj) {
  std <- c("term", "method", "estimate", "std.error", "conf.low", "conf.high", "p.value", "h", "n_eff")
  fill <- function(df) {
    for (v in std) if (!v %in% names(df)) df[[v]] <- NA
    df[, std, drop = FALSE]
  }
  if (inherits(obj, "rdrobust")) {
    t <- .cm_rd_tidy_rdrobust(obj)
    t$h <- t$h_left; t$n_eff <- t$n_left + t$n_right
    return(fill(t))
  }
  if (inherits(obj, "RDResults")) {
    co <- obj$coefficients
    return(fill(data.frame(term = as.character(co$term), method = "honest", estimate = co$estimate,
                           std.error = co$std.error, conf.low = co$conf.low, conf.high = co$conf.high,
                           p.value = co$p.value, h = co$bandwidth, n_eff = co$eff.obs, stringsAsFactors = FALSE)))
  }
  if (inherits(obj, "rdhte")) {
    est <- as.numeric(obj$Estimate); se <- as.numeric(obj$se.rb); ci <- obj$ci.rb
    return(fill(data.frame(term = names(obj$Estimate) %||% paste0("group", seq_along(est)), method = "robust",
                           estimate = est, std.error = se, conf.low = ci[, 1], conf.high = ci[, 2],
                           p.value = as.numeric(obj$pv.rb), h = obj$h[1], n_eff = sum(obj$Nh), stringsAsFactors = FALSE)))
  }
  if (inherits(obj, "cm_rd_adjust")) {
    if (is.null(obj$comparison)) return(NULL)
    t <- obj$comparison
    t$term <- paste0("RD effect, ", t$outcome); t$h <- t$h_left; t$n_eff <- t$n_left + t$n_right
    return(fill(t))
  }
  if (inherits(obj, "cm_rd_weak_iv")) {
    return(fill(obj$table))
  }
  if (is.list(obj) && all(c("obs.stat", "p.value", "window") %in% names(obj))) {
    ci <- if (!is.null(obj$ci) && is.matrix(obj$ci)) obj$ci[1, ] else c(NA, NA)
    return(fill(data.frame(term = "RD effect", method = "local randomization", estimate = as.numeric(obj$obs.stat),
                           std.error = NA_real_, conf.low = ci[1], conf.high = ci[2], p.value = as.numeric(obj$p.value),
                           h = diff(as.numeric(obj$window)) / 2, n_eff = if (!is.null(obj$sumstats)) sum(obj$sumstats[1, ]) else NA,
                           stringsAsFactors = FALSE)))
  }
  if (is.list(obj) && all(c("tau", "se.rb", "Coefs") %in% names(obj))) {
    co <- obj$Coefs
    labs <- colnames(co)
    is_cut <- !labs %in% c("weighted", "pooled")
    rows <- data.frame(term = ifelse(is_cut, paste0("cutoff ", labs), labs), method = "robust",
                       estimate = as.numeric(co[1, ]), std.error = sqrt(as.numeric(obj$V[1, ])),
                       conf.low = as.numeric(obj$CI[1, ]), conf.high = as.numeric(obj$CI[2, ]),
                       p.value = as.numeric(obj$Pv[1, ]), h = as.numeric(obj$H[1, ]),
                       n_eff = as.numeric(obj$Nh[1, ]), stringsAsFactors = FALSE)
    return(fill(rows))
  }
  NULL
}

#' Tidy methods for regression discontinuity objects
#'
#' `tidy()` methods for `RDHonest` (`RDResults`) and `rdhte` objects, so
#' that these fits work with `modelsummary`. (`rdrobust` ships its own
#' `tidy()` and `glance()` methods; [rd_frame()] uses the same numbers.)
#'
#' @param x An `RDResults` or `rdhte` object.
#' @param conf.int Include confidence limits.
#' @param conf.level Ignored; the level stored in the object is reported.
#' @param ... Unused.
#' @return A data frame in `broom` layout.
#' @name rd_tidiers
NULL

#' @rdname rd_tidiers
#' @export
tidy.RDResults <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  co <- x$coefficients
  out <- data.frame(term = as.character(co$term), estimate = co$estimate, std.error = co$std.error,
                    statistic = co$estimate / co$std.error, p.value = co$p.value, stringsAsFactors = FALSE)
  if (conf.int) { out$conf.low <- co$conf.low; out$conf.high <- co$conf.high }
  out
}

#' @rdname rd_tidiers
#' @export
tidy.rdhte <- function(x, conf.int = TRUE, conf.level = NULL, ...) {
  est <- as.numeric(x$Estimate); se <- as.numeric(x$se.rb)
  out <- data.frame(term = names(x$Estimate) %||% paste0("group", seq_along(est)), estimate = est, std.error = se,
                    statistic = as.numeric(x$t.rb), p.value = as.numeric(x$pv.rb), stringsAsFactors = FALSE)
  if (conf.int) { out$conf.low <- x$ci.rb[, 1]; out$conf.high <- x$ci.rb[, 2] }
  out
}
