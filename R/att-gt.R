#' Group-time average treatment effects (Callaway and Sant'Anna 2021)
#'
#' `att_gt()` estimates the building block of modern difference-in-differences
#' with staggered adoption: the average treatment effect at calendar time `t`
#' among units first treated in period `g`, `ATT(g, t)`, for every pair
#' `(g, t)`. Each cell is a two-by-two comparison of the long difference
#' `Y_t - Y_{base}` between group `g` and a comparison group (units never
#' treated, or units not yet treated by `t`). With covariates the cell is
#' estimated by outcome regression, inverse probability weighting, or the
#' doubly robust score of Sant'Anna and Zhao (2020). A two-group design is
#' the special case of one cohort, so the function also covers the canonical
#' two-by-two and two-group event studies.
#'
#' The function owns the causal layer only: cell definitions, the scores,
#' influence functions, the multiplier bootstrap, and simultaneous bands.
#' Regressions and logits inside the cells are parametric (as in `DRDID`)
#' unless `mlr3` learners are supplied, in which case the nuisances are
#' cross-fitted within each cell and the plain orthogonal score is used.
#'
#' @param data A data frame or `data.table` in long format (one row per unit
#'   and period for panel data; one row per observation for repeated
#'   cross-sections).
#' @param id Column name of the unit identifier.
#' @param time Column name of the period (numeric).
#' @param group Column name of the first treatment period (numeric;
#'   never-treated units coded `0`, `NA`, or `Inf`).
#' @param y Column name of the outcome.
#' @param x Covariates for conditional parallel trends: `NULL`, a character
#'   vector of column names, or a one-sided formula. Covariates are taken
#'   from the earlier of the two periods in each cell.
#' @param method `"dr"` (doubly robust, default), `"reg"` (outcome
#'   regression), or `"ipw"` (Hajek-normalized inverse probability weighting).
#' @param control_group `"notyet"` (units not yet treated by the later of the
#'   current and base period, default) or `"never"`.
#' @param base_period `"varying"` (default; the period just before `t` for
#'   pre-treatment cells and `g - 1 - anticipation` for post-treatment cells)
#'   or `"universal"` (always `g - 1 - anticipation`, so the base period
#'   itself has `ATT = 0` by construction).
#' @param anticipation Number of periods before `g` in which units may already
#'   respond; shifts the base period back.
#' @param sampling `"panel"` (balanced panel, default) or `"rcs"` (repeated
#'   cross-sections; `id` may then be omitted).
#' @param weights Optional column name of sampling weights.
#' @param cluster Optional column name for clustered inference (must be
#'   constant within unit).
#' @param learner_p,learner_or Optional `mlr3` learners for the propensity
#'   score and the outcome regression `E[Y_t - Y_base | X, comparison]`.
#'   Supplying either switches the cells to cross-fitted machine-learning
#'   nuisances (panel data only).
#' @param folds Cross-fitting folds for the learner path.
#' @param seed Optional seed for the cross-fitting partition and the
#'   bootstrap draws.
#' @param n_boot Number of multiplier-bootstrap draws (default 999). `0`
#'   returns analytic standard errors and pointwise intervals only.
#' @param boot_weights `"mammen"` (default) or `"rademacher"` multipliers.
#' @param conf_level Confidence level for intervals and the uniform band.
#' @param p_trim Comparison units with a propensity score above this value
#'   are dropped from the cell (as in `DRDID`).
#' @param keep_inffunc Keep the sparse influence-function matrix (needed by
#'   [aggregate_att()]).
#'
#' @return An object of class `cm_att_gt`: `att_gt` (one row per cell with
#'   `group`, `time`, `event_time`, `post`, `att`, `std.error`, `conf.low`,
#'   `conf.high`, `band.low`, `band.high`, `n_treated`, `n_control`),
#'   `crit_val` (uniform band critical value), `pretest` (Wald test that all
#'   pre-treatment cells are zero), `inffunc` (sparse `n x K` influence
#'   functions), `units` (unit-level group, weight, cluster), and the settings.
#'
#' @details
#' **Cells.** For group `g` and period `t`, the treated sample is group `g`
#' and the comparison sample is the never-treated units or, with
#' `control_group = "notyet"`, all units whose first treatment period exceeds
#' `max(t, base) + anticipation`. The outcome is the long difference
#' `Y_t - Y_base`. Cells with `t >= g` are post-treatment (`post = 1`); the
#' others are pre-treatment placebos.
#'
#' **Influence functions.** Every cell returns a unit-level influence
#' function, zero for units outside the cell and scaled by `n / n_cell`.
#' Standard errors are `sqrt(sum_c (sum_{i in c} IF_i)^2) / n` analytically,
#' or the interquartile-range estimate from the multiplier bootstrap. The
#' uniform band uses the bootstrap distribution of the maximum studentized
#' deviation across cells.
#'
#' @references
#' Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-differences
#' with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
#'
#' Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
#' difference-in-differences estimators. *Journal of Econometrics*, 219(1),
#' 101-122.
#'
#' @examples
#' dat <- sim_did_panel(n_units = 300, n_periods = 8, seed = 1)
#' fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y",
#'               x = c("x1", "x2"), control_group = "notyet", n_boot = 99)
#' fit
#' aggregate_att(fit, type = "dynamic")
#' @export
att_gt <- function(data, id = NULL, time, group, y, x = NULL,
                   method = c("dr", "reg", "ipw"),
                   control_group = c("notyet", "never"),
                   base_period = c("varying", "universal"),
                   anticipation = 0L,
                   sampling = c("panel", "rcs"),
                   weights = NULL, cluster = NULL,
                   learner_p = NULL, learner_or = NULL, folds = 5L, seed = NULL,
                   n_boot = 999L, boot_weights = c("mammen", "rademacher"),
                   conf_level = 0.95, p_trim = 0.995, keep_inffunc = TRUE) {
  method <- match.arg(method)
  control_group <- match.arg(control_group)
  base_period <- match.arg(base_period)
  sampling <- match.arg(sampling)
  boot_weights <- match.arg(boot_weights)
  cl <- match.call()
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (v in c(time, group, y)) .cm_check_column(v, data)
  if (!is.null(weights)) .cm_check_column(weights, data)
  if (!is.null(cluster)) .cm_check_column(cluster, data)
  if (sampling == "panel" && is.null(id)) stop("`id` is required for panel data.", call. = FALSE)
  if (!is.null(id)) .cm_check_column(id, data)
  anticipation <- .cm_check_count(anticipation, "anticipation", min = 0L)
  n_boot <- .cm_check_count(n_boot, "n_boot", min = 0L)
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) stop("`conf_level` must be in (0, 1).", call. = FALSE)
  use_learner <- !is.null(learner_p) || !is.null(learner_or)
  if (use_learner && sampling == "rcs") stop("Learner-based nuisances are available for panel data only.", call. = FALSE)
  x_vars <- .cm_did_x_vars(x)
  for (v in x_vars) .cm_check_column(v, data)

  dt <- data.table::as.data.table(data)
  keep_cols <- unique(c(id, time, group, y, x_vars, weights, cluster))
  dt <- dt[, keep_cols, with = FALSE]
  if (is.null(id)) {
    dt[, .cm_id := seq_len(.N)]
    id <- ".cm_id"
  }
  tt <- dt[[time]]
  gg <- dt[[group]]
  if (!is.numeric(tt)) stop("`time` must be numeric.", call. = FALSE)
  if (!is.numeric(gg)) stop("`group` must be numeric (first treatment period; 0, NA, or Inf for never treated).", call. = FALSE)
  gg[is.na(gg) | gg == 0] <- Inf
  if (anyNA(dt[[y]])) stop("`y` contains missing values.", call. = FALSE)
  if (length(x_vars) && anyNA(dt[, x_vars, with = FALSE])) stop("Covariates contain missing values.", call. = FALSE)
  tlist <- sort(unique(tt))
  nT <- length(tlist)
  if (nT < 2L) stop("At least two periods are required.", call. = FALSE)
  late <- is.finite(gg) & gg > max(tlist)
  if (any(late)) {
    message("Units first treated after the last period are used as never treated.")
    gg[late] <- Inf
  }
  early <- is.finite(gg) & gg <= tlist[1L + anticipation]
  if (any(early)) {
    warning("Units first treated in or before period ", tlist[1L + anticipation],
            " have no pre-treatment period and are dropped.", call. = FALSE)
    dt <- dt[!early]
    gg <- gg[!early]
    tt <- tt[!early]
  }
  data.table::set(dt, j = ".cm_g", value = gg)
  data.table::setorderv(dt, c(id, time))
  glist <- sort(unique(gg[is.finite(gg)]))
  if (length(glist) == 0L) stop("No treated groups found.", call. = FALSE)
  w_all <- if (is.null(weights)) rep(1, nrow(dt)) else as.numeric(dt[[weights]])
  if (any(!is.finite(w_all)) || any(w_all < 0)) stop("`weights` must be non-negative and finite.", call. = FALSE)
  w_all <- w_all / mean(w_all)
  design_all <- .cm_did_design(dt, x)
  p_dim <- ncol(design_all)

  # ---- unit-level objects ------------------------------------------------
  ids <- dt[[id]]
  uid <- unique(ids)
  n_units <- length(uid)
  unit_index <- match(ids, uid)
  if (sampling == "panel") {
    counts <- tabulate(unit_index, n_units)
    if (any(counts != nT) || anyDuplicated(dt[, c(id, time), with = FALSE])) {
      stop("`sampling = \"panel\"` requires a balanced panel: every unit observed once in every period. Use `sampling = \"rcs\"` otherwise.", call. = FALSE)
    }
    tt <- dt[[time]]
    t_index <- match(tt, tlist)
    Ymat <- matrix(NA_real_, n_units, nT)
    Ymat[cbind(unit_index, t_index)] <- as.numeric(dt[[y]])
    first <- !duplicated(unit_index)
    Gu <- gg[first][order(unit_index[first])]
    wu <- w_all[first][order(unit_index[first])]
    clu <- if (is.null(cluster)) NULL else dt[[cluster]][first][order(unit_index[first])]
    if (!is.null(cluster)) {
      cl_check <- dt[, list(.cm_nu = data.table::uniqueN(get(cluster))), by = id]
      if (any(cl_check$.cm_nu > 1L)) stop("`cluster` must be constant within unit.", call. = FALSE)
    }
    Xlist <- lapply(seq_len(nT), function(tp) {
      rows <- which(t_index == tp)
      m <- design_all[rows, , drop = FALSE][order(unit_index[rows]), , drop = FALSE]
      dimnames(m) <- list(NULL, colnames(design_all))
      m
    })
    x_invariant <- all(vapply(Xlist[-1L], function(m) isTRUE(all.equal(m, Xlist[[1L]], check.attributes = FALSE)), logical(1)))
    Xraw <- if (use_learner) lapply(seq_len(nT), function(tp) {
      rows <- which(t_index == tp)
      as.data.frame(dt[rows, x_vars, with = FALSE])[order(unit_index[rows]), , drop = FALSE]
    }) else NULL
    n_scale <- n_units
  } else {
    tt <- dt[[time]]
    t_index <- match(tt, tlist)
    g_row <- gg
    first <- !duplicated(unit_index)
    Gu <- gg[first][order(unit_index[first])]
    wu <- as.numeric(tapply(w_all, unit_index, mean))
    clu <- if (is.null(cluster)) NULL else dt[[cluster]][first][order(unit_index[first])]
    n_scale <- nrow(dt)
  }

  # ---- cell loop ---------------------------------------------------------
  tfac <- if (base_period == "varying") 1L else 0L
  t_range <- if (base_period == "varying") seq_len(nT - 1L) else seq_len(nT)
  cells <- list()
  if_i <- list()
  if_x <- list()
  ps_summary <- list()
  counter <- 0L
  for (g in glist) {
    idx_g <- which(tlist + anticipation < g)
    if (length(idx_g) == 0L) {
      warning("Group ", g, " has no pre-treatment period and is dropped.", call. = FALSE)
      next
    }
    pret_g <- max(idx_g)
    for (t in t_range) {
      tcur <- t + tfac
      pret <- if (base_period == "universal" || g <= tlist[tcur]) pret_g else t
      post <- as.integer(g <= tlist[tcur])
      counter <- counter + 1L
      if (base_period == "universal" && pret == tcur) {
        cells[[counter]] <- data.frame(group = g, time = tlist[tcur], att = 0, n_treated = sum(Gu == g), n_control = NA_integer_, post = post, fixed = TRUE)
        if_i[[counter]] <- integer(0)
        if_x[[counter]] <- numeric(0)
        next
      }
      if (control_group == "never") {
        C_u <- !is.finite(Gu)
      } else {
        thr <- tlist[max(t, pret) + tfac] + anticipation
        C_u <- !is.finite(Gu) | (Gu > thr & Gu != g)
      }
      G_u <- Gu == g
      keep <- which(G_u | C_u)
      if (sum(G_u) == 0L || sum(C_u) == 0L) {
        cells[[counter]] <- data.frame(group = g, time = tlist[tcur], att = NA_real_, n_treated = sum(G_u), n_control = sum(C_u), post = post, fixed = FALSE)
        if_i[[counter]] <- integer(0)
        if_x[[counter]] <- numeric(0)
        warning("No treated or no comparison units for group ", g, " in period ", tlist[tcur], ".", call. = FALSE)
        next
      }
      earlier <- min(tcur, pret)
      est <- tryCatch({
        if (sampling == "panel") {
          d_cell <- as.numeric(G_u[keep])
          y1 <- Ymat[keep, tcur]
          y0 <- Ymat[keep, pret]
          w_cell <- wu[keep]
          if (use_learner) {
            .cm_did_panel_cell_learner(
              dy = y1 - y0, d = d_cell, work = Xraw[[earlier]][keep, , drop = FALSE],
              features = x_vars, learner_p = learner_p, learner_or = learner_or,
              folds = folds, seed = seed, p_trim = p_trim
            )
          } else {
            Xc <- Xlist[[if (x_invariant) 1L else earlier]][keep, , drop = FALSE]
            .cm_did_panel_cell(y1, y0, d_cell, Xc, w_cell, method = method, p_trim = p_trim)
          }
        } else {
          rows <- which((t_index == tcur | t_index == pret) & (G_u[unit_index] | C_u[unit_index]))
          d_cell <- as.numeric(G_u[unit_index[rows]])
          post_r <- as.numeric(t_index[rows] == tcur)
          if (sum(d_cell * post_r) == 0 || sum(d_cell * (1 - post_r)) == 0 ||
              sum((1 - d_cell) * post_r) == 0 || sum((1 - d_cell) * (1 - post_r)) == 0) {
            stop("empty cell")
          }
          res <- .cm_did_rc_cell(dt[[y]][rows], post_r, d_cell, design_all[rows, , drop = FALSE], w_all[rows], method = method, p_trim = p_trim)
          res$rows <- rows
          res
        }
      }, error = function(e) {
        warning("Cell (g = ", g, ", t = ", tlist[tcur], ") could not be estimated: ", conditionMessage(e), call. = FALSE)
        NULL
      })
      if (is.null(est)) {
        cells[[counter]] <- data.frame(group = g, time = tlist[tcur], att = NA_real_, n_treated = sum(G_u), n_control = sum(C_u), post = post, fixed = FALSE)
        if_i[[counter]] <- integer(0)
        if_x[[counter]] <- numeric(0)
        next
      }
      if (sampling == "panel") {
        n1 <- length(keep)
        if_i[[counter]] <- keep
        if_x[[counter]] <- est$inffunc * (n_scale / n1)
        n_tr <- sum(G_u)
        n_co <- sum(C_u)
      } else {
        n1 <- length(est$rows)
        agg <- rowsum(est$inffunc * (n_scale / n1), unit_index[est$rows], reorder = TRUE)
        if_i[[counter]] <- as.integer(rownames(agg))
        if_x[[counter]] <- as.numeric(agg[, 1L])
        n_tr <- sum(d_cell)
        n_co <- sum(1 - d_cell)
      }
      cells[[counter]] <- data.frame(group = g, time = tlist[tcur], att = est$att, n_treated = n_tr, n_control = n_co, post = post, fixed = FALSE)
      if (!is.null(est$ps)) {
        ps <- est$ps
        ps_summary[[counter]] <- data.frame(group = g, time = tlist[tcur],
                                            ps_min = min(ps), ps_max = max(ps),
                                            share_trimmed = mean(ps >= p_trim & (if (sampling == "panel") d_cell == 0 else d_cell == 0)))
      }
    }
  }
  if (length(cells) == 0L) stop("No (g, t) cells could be formed.", call. = FALSE)
  tab <- do.call(rbind, cells)
  K <- nrow(tab)
  tab$event_time <- tab$time - tab$group
  inffunc <- Matrix::sparseMatrix(
    i = unlist(if_i, use.names = FALSE), j = rep.int(seq_len(K), lengths(if_i)),
    x = unlist(if_x, use.names = FALSE), dims = c(n_units, K)
  )

  # ---- inference ---------------------------------------------------------
  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha / 2)
  est_ok <- is.finite(tab$att) & !tab$fixed
  se <- rep(NA_real_, K)
  crit_val <- z
  boot <- NULL
  if (n_boot > 0L && any(est_ok)) {
    boot <- .cm_multiplier_bootstrap(inffunc[, est_ok, drop = FALSE], n = n_scale, n_boot = n_boot,
                                     cluster = clu, boot_weights = boot_weights,
                                     conf_level = conf_level, seed = seed)
    se[est_ok] <- boot$se
    crit_val <- boot$crit_val
    if (!is.finite(crit_val) || crit_val < z) {
      warning("The simultaneous critical value could not be computed reliably; using the pointwise critical value.", call. = FALSE)
      crit_val <- z
    }
  } else if (any(est_ok)) {
    se[est_ok] <- .cm_if_analytic_se(inffunc[, est_ok, drop = FALSE], n = n_scale, cluster = clu)
  }
  se[tab$fixed] <- NA_real_
  tab$std.error <- se
  tab$conf.low <- tab$att - z * se
  tab$conf.high <- tab$att + z * se
  tab$band.low <- tab$att - crit_val * se
  tab$band.high <- tab$att + crit_val * se
  pre <- which(tab$post == 0 & est_ok)
  pretest <- if (length(pre)) .cm_if_wald(tab$att[pre], inffunc[, pre, drop = FALSE], n_scale, clu) else NULL

  units <- data.frame(id = uid, group = Gu, weight = wu, stringsAsFactors = FALSE)
  if (!is.null(clu)) units$cluster <- clu
  out <- list(
    att_gt = tab[, c("group", "time", "event_time", "post", "att", "std.error", "conf.low", "conf.high", "band.low", "band.high", "n_treated", "n_control")],
    crit_val = crit_val,
    pretest = pretest,
    inffunc = if (keep_inffunc) inffunc else NULL,
    units = units,
    n = n_scale,
    n_units = n_units,
    groups = glist,
    times = tlist,
    method = if (use_learner) "dr_learner" else method,
    control_group = control_group,
    base_period = base_period,
    anticipation = anticipation,
    sampling = sampling,
    outcome = y,
    covariates = x_vars,
    n_boot = n_boot,
    boot_weights = boot_weights,
    conf_level = conf_level,
    cluster = cluster,
    learners = if (use_learner) list(p = .cm_learner_label(learner_p), or = .cm_learner_label(learner_or)) else NULL,
    diagnostics = list(propensity = if (length(ps_summary)) do.call(rbind, ps_summary) else NULL),
    call = cl
  )
  class(out) <- "cm_att_gt"
  out
}

#' @export
print.cm_att_gt <- function(x, digits = 4, ...) {
  cat("Group-time average treatment effects (Callaway and Sant'Anna)\n")
  cat("  Method: ", x$method, "; comparison: ", x$control_group, "-treated; base period: ", x$base_period,
      if (x$anticipation > 0) paste0("; anticipation: ", x$anticipation) else "", "\n", sep = "")
  cat("  Units: ", x$n_units, "; groups: ", length(x$groups), "; periods: ", length(x$times),
      "; sampling: ", x$sampling, "\n", sep = "")
  if (x$n_boot > 0) {
    cat("  Inference: multiplier bootstrap (", x$n_boot, " draws, ", x$boot_weights, " weights",
        if (!is.null(x$cluster)) paste0(", clustered by ", x$cluster) else "", "); ",
        round(100 * x$conf_level), "% uniform band critical value ", formatC(x$crit_val, digits = 3, format = "f"), "\n", sep = "")
  } else {
    cat("  Inference: analytic influence-function standard errors",
        if (!is.null(x$cluster)) paste0(" clustered by ", x$cluster) else "", "\n", sep = "")
  }
  tab <- x$att_gt[, c("group", "time", "event_time", "att", "std.error", "band.low", "band.high")]
  tab$sig <- ifelse(is.finite(tab$band.low) & (tab$band.low > 0 | tab$band.high < 0), "*", "")
  print(format(tab, digits = digits), row.names = FALSE)
  if (!is.null(x$pretest) && is.finite(x$pretest$p.value)) {
    cat("Pre-treatment cells jointly zero: chi-square(", x$pretest$df, ") = ",
        formatC(x$pretest$statistic, digits = 3, format = "f"), ", p = ",
        formatC(x$pretest$p.value, digits = 3, format = "f"), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname cm_tidiers
#' @export
tidy.cm_att_gt <- function(x, ...) {
  tab <- x$att_gt
  data.frame(
    term = paste0("ATT(", tab$group, ",", tab$time, ")"),
    group = tab$group, time = tab$time, event_time = tab$event_time, post = tab$post,
    estimate = tab$att, std.error = tab$std.error,
    statistic = tab$att / tab$std.error,
    p.value = 2 * stats::pnorm(-abs(tab$att / tab$std.error)),
    conf.low = tab$conf.low, conf.high = tab$conf.high,
    band.low = tab$band.low, band.high = tab$band.high,
    stringsAsFactors = FALSE
  )
}

#' @rdname cm_tidiers
#' @export
glance.cm_att_gt <- function(x, ...) {
  data.frame(
    n_units = x$n_units, n_groups = length(x$groups), n_periods = length(x$times),
    n_cells = nrow(x$att_gt), method = x$method, control_group = x$control_group,
    base_period = x$base_period, n_boot = x$n_boot,
    pretest_statistic = if (is.null(x$pretest)) NA_real_ else x$pretest$statistic,
    pretest_p.value = if (is.null(x$pretest)) NA_real_ else x$pretest$p.value,
    stringsAsFactors = FALSE
  )
}
