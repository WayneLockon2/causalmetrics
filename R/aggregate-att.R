#' Aggregate group-time effects into summary parameters
#'
#' `aggregate_att()` turns the `ATT(g, t)` cells of [att_gt()] into the
#' summary parameters of Callaway and Sant'Anna (2021): an event-study curve
#' by length of exposure (`"dynamic"`), average effects by cohort
#' (`"group"`), by calendar period (`"calendar"`), or one weighted average of
#' all post-treatment cells (`"simple"`). Weights are the group shares among
#' treated units, and the influence functions carry the estimation error of
#' those shares, so standard errors and uniform bands are valid for the
#' aggregate.
#'
#' @param x An object returned by [att_gt()] with `keep_inffunc = TRUE`.
#' @param type `"dynamic"` (default), `"group"`, `"calendar"`, or `"simple"`.
#' @param balance_e For `type = "dynamic"`: keep only groups observed for at
#'   least `balance_e` periods after treatment, so every event time averages
#'   the same set of groups.
#' @param min_e,max_e Event-time window for the dynamic aggregation; `max_e`
#'   also caps the post-treatment cells used by the other types.
#' @param n_boot,boot_weights,conf_level,seed Bootstrap settings; default to
#'   those of `x`.
#' @param na.rm Drop cells whose `ATT(g, t)` is missing.
#'
#' @return An object of class `cm_agg_att` with `overall` (estimate,
#'   standard error, interval), `by` (one row per event time, group, or
#'   period with pointwise interval and uniform band), the aggregation
#'   weights, `crit_val`, and the influence functions.
#'
#' @details The dynamic aggregate at event time `e` averages `ATT(g, g + e)`
#'   over the groups observed at `e`, weighted by group size; its overall
#'   value is the unweighted mean of the post-treatment event times. Group
#'   aggregates average the post-treatment cells of each cohort with equal
#'   weight over time; the overall value weights cohorts by size. Calendar
#'   aggregates average the cohorts treated by each period. See Callaway
#'   and Sant'Anna (2021, Section 3) for the definitions and the discussion
#'   of compositional changes across event times, which `balance_e` removes.
#'
#' @examples
#' dat <- sim_did_panel(n_units = 300, n_periods = 8, seed = 1)
#' fit <- att_gt(dat, id = "id", time = "time", group = "g", y = "y", n_boot = 99)
#' aggregate_att(fit, type = "dynamic", min_e = -3, max_e = 3)
#' aggregate_att(fit, type = "group")
#' @export
aggregate_att <- function(x, type = c("dynamic", "group", "calendar", "simple"),
                          balance_e = NULL, min_e = -Inf, max_e = Inf,
                          n_boot = NULL, boot_weights = NULL, conf_level = NULL, seed = NULL,
                          na.rm = FALSE) {
  type <- match.arg(type)
  if (!inherits(x, "cm_att_gt")) stop("`x` must be an `att_gt()` result.", call. = FALSE)
  if (is.null(x$inffunc)) stop("`x` has no influence functions; rerun att_gt() with keep_inffunc = TRUE.", call. = FALSE)
  if (is.null(n_boot)) n_boot <- x$n_boot
  if (is.null(boot_weights)) boot_weights <- x$boot_weights
  if (is.null(conf_level)) conf_level <- x$conf_level
  tab <- x$att_gt
  IF <- x$inffunc
  n <- x$n
  clu <- x$units$cluster
  wu <- x$units$weight
  Gu <- x$units$group
  ok <- is.finite(tab$att)
  if (!na.rm && !all(ok)) stop("Some ATT(g, t) are missing; set `na.rm = TRUE` to drop them.", call. = FALSE)
  if (na.rm) {
    tab <- tab[ok, ]
    IF <- IF[, ok, drop = FALSE]
  }
  att <- tab$att
  group <- tab$group
  time <- tab$time
  etime <- tab$event_time
  glist <- sort(unique(group))
  pg_all <- vapply(glist, function(g) mean(wu * (Gu == g)), numeric(1))
  pg <- pg_all[match(group, glist)]
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  # influence function of a weighted average of cells, with the estimation
  # error of the group shares (the `wif` term of did::compute.aggte)
  agg_if <- function(which, weights, with_wif = TRUE) {
    base <- as.numeric(IF[, which, drop = FALSE] %*% weights)
    if (!with_wif) return(base)
    spg <- sum(pg[which])
    g_k <- group[which]
    a_k <- att[which]
    # sum_k (w_i 1{G_i = g_k} - pg_k) a_k / spg  -  sum_k (w_i 1{G_i = g_k} - pg_k) * sum_k pg_k a_k / spg^2
    ind <- vapply(seq_len(length(Gu)), function(i) 0, numeric(1))
    att_by_g <- tapply(a_k, g_k, sum)
    cnt_by_g <- tapply(a_k, g_k, length)
    m <- match(Gu, names(att_by_g))
    term1 <- ifelse(is.na(m), 0, wu * att_by_g[m]) - sum(pg[which] * a_k)
    rows <- ifelse(is.na(m), 0, wu * cnt_by_g[m]) - spg
    base + as.numeric(term1) / spg - as.numeric(rows) * sum(pg[which] * a_k) / spg^2
  }
  infer <- function(IFmat, ests) {
    IFmat <- as.matrix(IFmat)
    if (n_boot > 0L) {
      b <- .cm_multiplier_bootstrap(IFmat, n = n, n_boot = n_boot, cluster = clu,
                                    boot_weights = boot_weights, conf_level = conf_level, seed = seed)
      cv <- b$crit_val
      if (!is.finite(cv) || cv < z) cv <- z
      list(se = b$se, crit_val = cv)
    } else {
      list(se = .cm_if_analytic_se(IFmat, n = n, cluster = clu), crit_val = z)
    }
  }
  keepers <- which(tab$post == 1 & etime <= max_e)
  if (length(keepers) == 0L) stop("No post-treatment cells are available.", call. = FALSE)

  if (type == "simple") {
    wts <- pg[keepers] / sum(pg[keepers])
    est <- sum(att[keepers] * wts)
    IFo <- agg_if(keepers, wts)
    inf <- infer(IFo, est)
    by <- NULL
    weights_out <- data.frame(group = group[keepers], time = time[keepers], weight = wts)
    IFby <- NULL
  } else if (type == "group") {
    by_est <- vapply(glist, function(g) mean(att[which(group == g & tab$post == 1 & etime <= max_e)]), numeric(1))
    IFby <- sapply(glist, function(g) {
      wh <- which(group == g & tab$post == 1 & etime <= max_e)
      agg_if(wh, rep(1 / length(wh), length(wh)), with_wif = FALSE)
    })
    inf_by <- infer(IFby, by_est)
    wts <- pg_all / sum(pg_all)
    est <- sum(by_est * wts)
    # overall: weights are group shares among the treated, with their wif
    spg <- sum(pg_all)
    m <- match(Gu, glist)
    term1 <- ifelse(is.na(m), 0, wu * by_est[m]) - sum(pg_all * by_est)
    rows <- ifelse(is.na(m), 0, wu) - spg
    IFo <- as.numeric(IFby %*% wts) + term1 / spg - rows * sum(pg_all * by_est) / spg^2
    inf <- infer(IFo, est)
    by <- data.frame(group = glist, estimate = by_est, std.error = inf_by$se)
    weights_out <- data.frame(group = glist, weight = wts)
    inf_by_cv <- inf_by$crit_val
  } else if (type == "dynamic") {
    include <- rep(TRUE, length(group))
    if (!is.null(balance_e)) {
      include <- (max(x$times) - group) >= balance_e
    }
    eseq <- sort(unique(etime[include]))
    if (!is.null(balance_e)) eseq <- eseq[eseq <= balance_e & eseq >= balance_e - (max(x$times) - min(x$times))]
    eseq <- eseq[eseq >= min_e & eseq <= max_e]
    if (length(eseq) == 0L) stop("No event times fall in the requested window.", call. = FALSE)
    by_est <- vapply(eseq, function(e) {
      wh <- which(etime == e & include)
      sum(att[wh] * pg[wh] / sum(pg[wh]))
    }, numeric(1))
    IFby <- sapply(eseq, function(e) {
      wh <- which(etime == e & include)
      agg_if(wh, pg[wh] / sum(pg[wh]))
    })
    inf_by <- infer(IFby, by_est)
    epos <- eseq >= 0
    est <- mean(by_est[epos])
    IFo <- as.numeric(IFby[, epos, drop = FALSE] %*% rep(1 / sum(epos), sum(epos)))
    inf <- infer(IFo, est)
    by <- data.frame(event_time = eseq, estimate = by_est, std.error = inf_by$se)
    weights_out <- do.call(rbind, lapply(eseq, function(e) {
      wh <- which(etime == e & include)
      data.frame(event_time = e, group = group[wh], weight = pg[wh] / sum(pg[wh]))
    }))
    inf_by_cv <- inf_by$crit_val
  } else {
    tl <- sort(unique(time[tab$post == 1]))
    by_est <- vapply(tl, function(t1) {
      wh <- which(time == t1 & tab$post == 1)
      sum(att[wh] * pg[wh] / sum(pg[wh]))
    }, numeric(1))
    IFby <- sapply(tl, function(t1) {
      wh <- which(time == t1 & tab$post == 1)
      agg_if(wh, pg[wh] / sum(pg[wh]))
    })
    inf_by <- infer(IFby, by_est)
    est <- mean(by_est)
    IFo <- as.numeric(IFby %*% rep(1 / length(tl), length(tl)))
    inf <- infer(IFo, est)
    by <- data.frame(time = tl, estimate = by_est, std.error = inf_by$se)
    weights_out <- do.call(rbind, lapply(tl, function(t1) {
      wh <- which(time == t1 & tab$post == 1)
      data.frame(time = t1, group = group[wh], weight = pg[wh] / sum(pg[wh]))
    }))
    inf_by_cv <- inf_by$crit_val
  }
  if (!is.null(by)) {
    cv <- inf_by_cv
    by$conf.low <- by$estimate - z * by$std.error
    by$conf.high <- by$estimate + z * by$std.error
    by$band.low <- by$estimate - cv * by$std.error
    by$band.high <- by$estimate + cv * by$std.error
  } else {
    cv <- z
  }
  out <- list(
    type = type,
    overall = data.frame(estimate = est, std.error = inf$se,
                         conf.low = est - z * inf$se, conf.high = est + z * inf$se),
    by = by,
    weights = weights_out,
    crit_val = cv,
    inffunc = IFby,
    inffunc_overall = IFo,
    n = n,
    cluster = clu,
    balance_e = balance_e, min_e = min_e, max_e = max_e,
    n_boot = n_boot, conf_level = conf_level,
    source = x[c("method", "control_group", "base_period", "anticipation", "outcome")],
    call = match.call()
  )
  class(out) <- "cm_agg_att"
  out
}

#' @export
print.cm_agg_att <- function(x, digits = 4, ...) {
  label <- switch(x$type, dynamic = "event time", group = "group (cohort)", calendar = "calendar time", simple = "all post-treatment cells")
  cat("Aggregated ATT (", label, "); ", x$source$method, " cells, ", x$source$control_group, "-treated comparison\n", sep = "")
  cat("Overall ATT: ", formatC(x$overall$estimate, digits = digits, format = "g"),
      "  (SE ", formatC(x$overall$std.error, digits = digits, format = "g"), ")  ",
      round(100 * x$conf_level), "% CI [", formatC(x$overall$conf.low, digits = digits, format = "g"), ", ",
      formatC(x$overall$conf.high, digits = digits, format = "g"), "]\n", sep = "")
  if (!is.null(x$by)) {
    tab <- x$by
    tab$sig <- ifelse(is.finite(tab$band.low) & (tab$band.low > 0 | tab$band.high < 0), "*", "")
    print(format(tab, digits = digits), row.names = FALSE)
    cat("Uniform band critical value ", formatC(x$crit_val, digits = 3, format = "f"), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname cm_tidiers
#' @export
tidy.cm_agg_att <- function(x, ...) {
  if (is.null(x$by)) {
    return(data.frame(term = "ATT", estimate = x$overall$estimate, std.error = x$overall$std.error,
                      statistic = x$overall$estimate / x$overall$std.error,
                      p.value = 2 * stats::pnorm(-abs(x$overall$estimate / x$overall$std.error)),
                      conf.low = x$overall$conf.low, conf.high = x$overall$conf.high, stringsAsFactors = FALSE))
  }
  key <- names(x$by)[1L]
  data.frame(
    term = paste0(key, "::", x$by[[key]]),
    level = x$by[[key]],
    estimate = x$by$estimate, std.error = x$by$std.error,
    statistic = x$by$estimate / x$by$std.error,
    p.value = 2 * stats::pnorm(-abs(x$by$estimate / x$by$std.error)),
    conf.low = x$by$conf.low, conf.high = x$by$conf.high,
    band.low = x$by$band.low, band.high = x$by$band.high,
    stringsAsFactors = FALSE
  )
}

#' @rdname cm_tidiers
#' @export
glance.cm_agg_att <- function(x, ...) {
  data.frame(type = x$type, estimate = x$overall$estimate, std.error = x$overall$std.error,
             n = x$n, n_boot = x$n_boot, crit_val = x$crit_val, stringsAsFactors = FALSE)
}
