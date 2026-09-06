#' Imputation estimator for staggered difference-in-differences
#'
#' The imputation estimator of Borusyak, Jaravel, and Spiess (2024) (see also
#' Gardner 2022) fits unit and period effects on untreated observations only,
#' imputes the untreated potential outcome of every treated observation, and
#' averages the differences `Y - Y_hat(0)` by event time, cohort, or overall.
#' It uses every pre-treatment period as the baseline, which is efficient
#' when parallel trends holds in all periods and errors are homoskedastic
#' and serially uncorrelated, and it never compares treated units to
#' already-treated units.
#'
#' @inheritParams att_gt
#' @param x Optional covariates (character vector or one-sided formula)
#'   entered linearly in the first stage.
#' @param horizons Event times to report (default: all observed).
#' @param pre_window Number of pre-treatment leads for the placebo
#'   regression on untreated observations (0 to skip).
#' @param cluster Optional cluster column for the standard errors (default:
#'   unit).
#'
#' @return An object of class `cm_did_imputation`: `by_event` (event time,
#'   estimate, standard error, number of cells), `overall`, `by_group`,
#'   `pre` (lead coefficients from the untreated sample with clustered
#'   standard errors), `cells` (imputed effects per treated observation),
#'   the first-stage `fixest` model, and the call.
#'
#' @details Standard errors follow Borusyak, Jaravel, and Spiess (2024,
#'   Section 4): the estimator is linear in the outcome, so its variance is
#'   the sum over clusters of squared weighted residuals, with first-stage
#'   residuals on untreated cells and, conservatively, deviations of the
#'   imputed effects from their event-time mean on treated cells. The
#'   imputation weights on untreated cells are computed exactly by solving
#'   the sparse normal equations of the first stage.
#'
#' @references
#' Borusyak, K., Jaravel, X., and Spiess, J. (2024). Revisiting event-study
#' designs: Robust and efficient estimation. *Review of Economic Studies*,
#' 91(6), 3253-3285.
#'
#' Gardner, J. (2022). Two-stage differences in differences. Working paper.
#' @examples
#' dat <- sim_did_panel(n_units = 300, n_periods = 8, groups = c(3, 6), seed = 1)
#' did_imputation(dat, id = "id", time = "time", group = "g", y = "y")
#' @export
did_imputation <- function(data, id, time, group, y, x = NULL, horizons = NULL,
                           pre_window = 3L, cluster = NULL, weights = NULL) {
  for (v in c(id, time, group, y, weights, cluster)) .cm_check_column(v, data)
  x_vars <- .cm_did_x_vars(x)
  for (v in x_vars) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, unique(c(id, time, group, y, x_vars, weights, cluster)), with = FALSE]
  data.table::setnames(dt, c(id, time, group, y), c(".id", ".t", ".g0", ".y"))
  gg <- dt$.g0
  gg[is.na(gg) | gg == 0] <- Inf
  dt[, .g := gg]
  dt[, .treated := is.finite(.g) & .t >= .g]
  dt[, .e := ifelse(is.finite(.g), .t - .g, NA_real_)]
  dt[, .w := if (is.null(weights)) 1 else get(weights)]
  dt[, .cl := if (is.null(cluster)) .id else get(cluster)]
  untreated <- dt[.treated == FALSE]
  treated <- dt[.treated == TRUE]
  if (nrow(treated) == 0L) stop("No treated observations.", call. = FALSE)
  rhs <- if (length(x_vars)) paste(x_vars, collapse = " + ") else "1"
  fml <- stats::as.formula(paste0(".y ~ ", rhs, " | .id + .t"))
  stage1 <- fixest::feols(fml, data = untreated, weights = ~.w)
  pred <- suppressWarnings(stats::predict(stage1, newdata = treated))
  treated[, .y0 := pred]
  dropped <- sum(is.na(treated$.y0))
  if (dropped > 0L) {
    warning(dropped, " treated observations could not be imputed (their unit or period has no untreated observation) and were dropped.", call. = FALSE)
    treated <- treated[!is.na(.y0)]
  }
  treated[, .tau := .y - .y0]
  if (is.null(horizons)) horizons <- sort(unique(treated$.e))
  treated <- treated[.e %in% horizons]

  # ---- exact linear representation for the variance ----------------------
  # tau_hat(w) = sum_T w Y - v' Y_U with v = X_U (X_U'X_U)^{-1} X_T' w.
  u_ids <- sort(unique(dt$.id))
  t_ids <- sort(unique(dt$.t))
  make_X <- function(d) {
    ui <- match(d$.id, u_ids)
    ti <- match(d$.t, t_ids)
    parts <- list(
      Matrix::sparseMatrix(i = seq_len(nrow(d)), j = ui, x = 1, dims = c(nrow(d), length(u_ids))),
      Matrix::sparseMatrix(i = seq_len(nrow(d)), j = ti, x = 1, dims = c(nrow(d), length(t_ids)))[, -1L, drop = FALSE]
    )
    if (length(x_vars)) parts[[3L]] <- Matrix::Matrix(as.matrix(d[, x_vars, with = FALSE]), sparse = TRUE)
    do.call(cbind, parts)
  }
  XU <- make_X(untreated)
  XT <- make_X(treated)
  wU <- untreated$.w
  XtX <- Matrix::crossprod(XU * sqrt(wU))
  # identify: drop columns with no untreated support (all-zero) to keep XtX invertible
  keep_col <- Matrix::diag(XtX) > 0
  XtX <- XtX[keep_col, keep_col, drop = FALSE]
  resid_U <- stats::resid(stage1)
  ev <- as.numeric(treated$.e)
  tau_mean_e <- stats::ave(treated$.tau, ev, FUN = mean)
  resid_T <- treated$.tau - tau_mean_e
  cl_U <- untreated$.cl
  cl_T <- treated$.cl
  agg_se <- function(wT) {
    r <- as.numeric(Matrix::crossprod(XT[, keep_col, drop = FALSE], wT))
    b <- Matrix::solve(XtX, r)
    v <- as.numeric((XU[, keep_col, drop = FALSE] %*% b) * wU)
    contrib <- c(wT * resid_T, -v * resid_U)
    cls <- c(cl_T, cl_U)
    s <- rowsum(contrib, cls)
    sqrt(sum(s^2))
  }
  est_and_se <- function(sel) {
    wT <- numeric(nrow(treated))
    wT[sel] <- treated$.w[sel] / sum(treated$.w[sel])
    c(estimate = sum(wT * treated$.tau), std.error = agg_se(wT), n_cells = sum(sel))
  }
  by_event <- do.call(rbind, lapply(horizons, function(e) {
    r <- est_and_se(ev == e)
    data.frame(event_time = e, estimate = r[["estimate"]], std.error = r[["std.error"]], n_cells = r[["n_cells"]])
  }))
  ov <- est_and_se(rep(TRUE, nrow(treated)))
  overall <- data.frame(estimate = ov[["estimate"]], std.error = ov[["std.error"]], n_cells = ov[["n_cells"]])
  glist <- sort(unique(treated$.g))
  by_group <- do.call(rbind, lapply(glist, function(g) {
    r <- est_and_se(treated$.g == g)
    data.frame(group = g, estimate = r[["estimate"]], std.error = r[["std.error"]], n_cells = r[["n_cells"]])
  }))
  z <- stats::qnorm(0.975)
  for (nm in c("by_event", "by_group", "overall")) {
    obj <- get(nm)
    obj$conf.low <- obj$estimate - z * obj$std.error
    obj$conf.high <- obj$estimate + z * obj$std.error
    assign(nm, obj)
  }

  # ---- pre-trend leads on the untreated sample ----------------------------
  pre <- NULL
  if (pre_window > 0L) {
    untreated[, .lead := ifelse(is.finite(.g) & (.g - .t) <= pre_window, .g - .t, 0)]
    if (any(untreated$.lead > 0)) {
      fml_pre <- stats::as.formula(paste0(".y ~ i(.lead, ref = 0)", if (length(x_vars)) paste0(" + ", rhs) else "", " | .id + .t"))
      fit_pre <- fixest::feols(fml_pre, data = untreated, weights = ~.w, cluster = ~.cl)
      ct <- fixest::coeftable(fit_pre)
      lead_rows <- grepl("^\\.lead::", rownames(ct))
      pre <- data.frame(
        event_time = -as.numeric(sub("^\\.lead::", "", rownames(ct)[lead_rows])),
        estimate = ct[lead_rows, 1L], std.error = ct[lead_rows, 2L]
      )
      pre$conf.low <- pre$estimate - z * pre$std.error
      pre$conf.high <- pre$estimate + z * pre$std.error
      pre <- pre[order(pre$event_time), ]
      rownames(pre) <- NULL
    }
  }
  out <- list(by_event = by_event, overall = overall, by_group = by_group, pre = pre,
              cells = as.data.frame(treated[, list(.id, .t, .g, .e, .y, .y0, .tau)]),
              stage1 = stage1, n_treated = nrow(treated), n_untreated = nrow(untreated),
              cluster = if (is.null(cluster)) id else cluster, call = match.call())
  names(out$cells) <- c(id, time, "group", "event_time", y, "y0_hat", "tau_hat")
  class(out) <- "cm_did_imputation"
  out
}

#' @export
print.cm_did_imputation <- function(x, digits = 4, ...) {
  cat("Imputation difference-in-differences (Borusyak, Jaravel, and Spiess)\n")
  cat("  Treated observations: ", x$n_treated, "; untreated observations used in the first stage: ", x$n_untreated,
      "; standard errors clustered by ", x$cluster, "\n", sep = "")
  cat("Overall ATT: ", formatC(x$overall$estimate, digits = digits, format = "g"), "  (SE ",
      formatC(x$overall$std.error, digits = digits, format = "g"), ")\n", sep = "")
  print(format(x$by_event, digits = digits), row.names = FALSE)
  if (!is.null(x$pre)) {
    cat("Pre-treatment leads (untreated sample):\n")
    print(format(x$pre, digits = digits), row.names = FALSE)
  }
  invisible(x)
}

#' @rdname cm_tidiers
#' @export
tidy.cm_did_imputation <- function(x, ...) {
  tab <- rbind(
    if (!is.null(x$pre)) data.frame(event_time = x$pre$event_time, estimate = x$pre$estimate, std.error = x$pre$std.error) else NULL,
    data.frame(event_time = x$by_event$event_time, estimate = x$by_event$estimate, std.error = x$by_event$std.error)
  )
  data.frame(term = paste0("event_time::", tab$event_time), event_time = tab$event_time,
             estimate = tab$estimate, std.error = tab$std.error,
             statistic = tab$estimate / tab$std.error,
             p.value = 2 * stats::pnorm(-abs(tab$estimate / tab$std.error)),
             conf.low = tab$estimate - stats::qnorm(0.975) * tab$std.error,
             conf.high = tab$estimate + stats::qnorm(0.975) * tab$std.error,
             stringsAsFactors = FALSE)
}
