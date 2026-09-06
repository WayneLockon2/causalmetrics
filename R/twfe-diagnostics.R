#' Goodman-Bacon decomposition of the static two-way fixed effects estimator
#'
#' With staggered adoption, the coefficient on the treatment dummy in
#' `y ~ d | id + time` is a weighted average of all two-by-two
#' difference-in-differences comparisons in the data: each timing group
#' against the never-treated units, earlier-treated against later-treated
#' units (before the later group is treated), and later-treated against
#' earlier-treated units (after the earlier group is treated). The last type
#' uses already-treated units as controls, which is where heterogeneous
#' dynamic effects bias the regression. `bacon_decomp()` returns every
#' comparison, its estimate, and its weight (Goodman-Bacon 2021).
#'
#' @param data A balanced panel in long format.
#' @param id Column name of the unit identifier.
#' @param time Column name of the period.
#' @param y Column name of the outcome.
#' @param d Column name of the treatment indicator (0/1, absorbing: once 1,
#'   stays 1). Units treated in the first period are always treated and are
#'   dropped with a warning, as in the original decomposition.
#'
#' @return An object of class `cm_bacon`: `decomposition` (one row per
#'   comparison: `treated`, `control`, `type`, `weight`, `estimate`),
#'   `by_type` (weights and weighted averages by comparison type),
#'   `twfe` (the static TWFE coefficient), and `check` (the weighted sum of
#'   the two-by-two estimates, which equals `twfe` in a balanced panel).
#'
#' @references Goodman-Bacon, A. (2021). Difference-in-differences with
#'   variation in treatment timing. *Journal of Econometrics*, 225(2), 254-277.
#' @examples
#' dat <- sim_did_panel(n_units = 200, n_periods = 8, groups = c(3, 6), seed = 1)
#' bacon_decomp(dat, id = "id", time = "time", y = "y", d = "treated")
#' @export
bacon_decomp <- function(data, id, time, y, d) {
  for (v in c(id, time, y, d)) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, c(id, time, y, d), with = FALSE]
  data.table::setnames(dt, c("id", "time", "y", "d"))
  dt[, d := .cm_as_binary(d, "d")]
  tlist <- sort(unique(dt$time))
  nT <- length(tlist)
  cnt <- dt[, .N, by = id]
  if (any(cnt$N != nT) || anyDuplicated(dt[, c("id", "time")])) stop("`bacon_decomp()` requires a balanced panel.", call. = FALSE)
  first <- dt[d == 1L, list(g = as.numeric(min(time))), by = id]
  dt <- merge(dt, first, by = "id", all.x = TRUE)
  dt[is.na(g), g := Inf]
  bad <- dt[is.finite(g), any(d[time < g] == 1L) || any(d[time >= g] == 0L), by = id]
  if (any(bad$V1)) stop("`d` must be an absorbing treatment: 0 before the first treated period and 1 afterwards.", call. = FALSE)
  always <- unique(dt[g == tlist[1L], id])
  if (length(always)) {
    warning(length(always), " always-treated units (treated in the first period) were dropped.", call. = FALSE)
    dt <- dt[!id %in% always]
  }
  # collapse to group-by-time means
  n_units <- data.table::uniqueN(dt$id)
  gsize <- dt[time == tlist[1L], .N, by = g][order(g)]
  gsize[, share := N / n_units]
  cell <- dt[, list(ybar = mean(y)), by = list(g, time)]
  dbar <- dt[, list(dbar = mean(d)), by = g]
  # variance of the double-demeaned treatment
  dt[, dtil := d - mean(d), by = id]
  dt[, dtil := dtil - mean(dtil), by = time]
  VD <- mean(dt$dtil^2)
  twfe <- fixest::feols(y ~ d | id + time, data = dt)
  beta <- unname(stats::coef(twfe)[["d"]])
  groups <- gsize$g[is.finite(gsize$g)]
  has_never <- any(!is.finite(gsize$g))
  two_by_two <- function(tr, co, window) {
    tr_pre <- cell[g == tr & time %in% window & time < tr, mean(ybar)]
    tr_post <- cell[g == tr & time %in% window & time >= tr, mean(ybar)]
    co_pre <- cell[g == co & time %in% window & time < tr, mean(ybar)]
    co_post <- cell[g == co & time %in% window & time >= tr, mean(ybar)]
    (tr_post - tr_pre) - (co_post - co_pre)
  }
  rows <- list()
  share <- function(gg) gsize$share[match(gg, gsize$g)]
  dshare <- function(gg) dbar$dbar[match(gg, dbar$g)]
  if (has_never) {
    nU <- share(Inf)
    for (k in groups) {
      nk <- share(k); nkU <- nk / (nk + nU); Dk <- dshare(k)
      rows[[length(rows) + 1L]] <- data.frame(
        treated = k, control = Inf, type = "treated vs never treated",
        weight = (nk + nU)^2 * nkU * (1 - nkU) * Dk * (1 - Dk) / VD,
        estimate = two_by_two(k, Inf, tlist), stringsAsFactors = FALSE
      )
    }
  }
  if (length(groups) > 1L) {
    for (i in seq_len(length(groups) - 1L)) {
      for (j in (i + 1L):length(groups)) {
        k <- groups[i]; l <- groups[j]
        nk <- share(k); nl <- share(l); nkl <- nk / (nk + nl); Dk <- dshare(k); Dl <- dshare(l)
        rows[[length(rows) + 1L]] <- data.frame(
          treated = k, control = l, type = "earlier vs later treated",
          weight = ((nk + nl) * (1 - Dl))^2 * nkl * (1 - nkl) * ((Dk - Dl) / (1 - Dl)) * ((1 - Dk) / (1 - Dl)) / VD,
          estimate = two_by_two(k, l, tlist[tlist < l]), stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
          treated = l, control = k, type = "later vs earlier treated",
          weight = ((nk + nl) * Dk)^2 * nkl * (1 - nkl) * (Dl / Dk) * ((Dk - Dl) / Dk) / VD,
          estimate = two_by_two(l, k, tlist[tlist >= k]), stringsAsFactors = FALSE
        )
      }
    }
  }
  dec <- do.call(rbind, rows)
  by_type <- do.call(rbind, lapply(split(dec, dec$type), function(s) {
    data.frame(type = s$type[1L], weight = sum(s$weight),
               estimate = sum(s$weight * s$estimate) / sum(s$weight), stringsAsFactors = FALSE)
  }))
  rownames(by_type) <- NULL
  out <- list(decomposition = dec, by_type = by_type, twfe = beta,
              check = sum(dec$weight * dec$estimate), sum_weights = sum(dec$weight),
              n_units = n_units, groups = groups, has_never_treated = has_never, call = match.call())
  class(out) <- "cm_bacon"
  out
}

#' @export
print.cm_bacon <- function(x, digits = 4, ...) {
  cat("Goodman-Bacon decomposition of the static TWFE coefficient\n")
  cat("  TWFE estimate: ", formatC(x$twfe, digits = digits, format = "g"),
      "; weighted sum of 2x2 estimates: ", formatC(x$check, digits = digits, format = "g"),
      "; weights sum to ", formatC(x$sum_weights, digits = 4, format = "f"), "\n", sep = "")
  print(format(x$by_type, digits = digits), row.names = FALSE)
  cat("\nAll comparisons:\n")
  print(format(x$decomposition, digits = digits), row.names = FALSE)
  invisible(x)
}

#' Weights that static two-way fixed effects puts on each treated cell
#'
#' de Chaisemartin and D'Haultfoeuille (2020) show that the static TWFE
#' coefficient equals a weighted sum of the treatment effects of the treated
#' unit-period cells, with weights that sum to one but can be negative. The
#' weight of cell `(i, t)` is the residual of `d` on unit and time effects
#' divided by the sum of those residuals over treated cells. `twfe_weights()`
#' computes them, summarizes them by cohort and period, and reports the share
#' and mass of negative weights.
#'
#' @inheritParams bacon_decomp
#' @param y Optional outcome; when given, the function verifies the identity
#'   `beta = sum(eps * y) / sum_{treated} eps` that underlies the weights.
#'
#' @return An object of class `cm_twfe_weights` with `cells` (unit-period
#'   weights on treated cells), `by_cell` (weights summed by first-treatment
#'   period and calendar period), `share_negative`, `sum_negative`, and, when
#'   `y` is supplied, `twfe` and `check`.
#'
#' @references de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way
#'   fixed effects estimators with heterogeneous treatment effects.
#'   *American Economic Review*, 110(9), 2964-2996.
#' @examples
#' dat <- sim_did_panel(n_units = 200, n_periods = 8, groups = c(3, 6), seed = 1)
#' twfe_weights(dat, id = "id", time = "time", d = "treated", y = "y")
#' @export
twfe_weights <- function(data, id, time, d, y = NULL) {
  for (v in c(id, time, d, y)) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, c(id, time, d, y), with = FALSE]
  data.table::setnames(dt, c("id", "time", "d", if (!is.null(y)) "y"))
  dt[, d := .cm_as_binary(d, "d")]
  fit_d <- fixest::feols(d ~ 1 | id + time, data = dt)
  dt[, eps := stats::resid(fit_d)]
  first <- dt[d == 1L, list(g = min(time)), by = id]
  dt <- merge(dt, first, by = "id", all.x = TRUE)
  treated <- dt[d == 1L]
  treated[, weight := eps / sum(eps)]
  by_cell <- treated[, list(weight = sum(weight), n = .N), by = list(g, time)][order(g, time)]
  out <- list(
    cells = as.data.frame(treated[, list(id, time, g, weight)]),
    by_cell = as.data.frame(by_cell),
    share_negative = mean(treated$weight < 0),
    sum_negative = sum(treated$weight[treated$weight < 0]),
    n_treated_cells = nrow(treated),
    call = match.call()
  )
  if (!is.null(y)) {
    fit_y <- fixest::feols(y ~ d | id + time, data = dt)
    out$twfe <- unname(stats::coef(fit_y)[["d"]])
    out$check <- sum(dt$eps * dt$y) / sum(treated$eps)
  }
  class(out) <- "cm_twfe_weights"
  out
}

#' @export
print.cm_twfe_weights <- function(x, digits = 4, ...) {
  cat("Two-way fixed effects weights on treated cells (de Chaisemartin and D'Haultfoeuille)\n")
  cat("  Treated cells: ", x$n_treated_cells, "; share with negative weight: ",
      formatC(x$share_negative, digits = 3, format = "f"), "; mass of negative weights: ",
      formatC(x$sum_negative, digits = 3, format = "f"), "\n", sep = "")
  if (!is.null(x$twfe)) {
    cat("  TWFE estimate ", formatC(x$twfe, digits = digits, format = "g"),
        " = sum(eps * y) / sum_treated(eps) = ", formatC(x$check, digits = digits, format = "g"), "\n", sep = "")
  }
  print(format(x$by_cell, digits = digits), row.names = FALSE)
  invisible(x)
}
