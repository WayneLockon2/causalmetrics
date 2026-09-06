#' Put event-study estimates from different estimators on one frame
#'
#' `event_study_frame()` collects event-time estimates from
#' [aggregate_att()] (dynamic type), [att_gt()] with a single cohort,
#' [did_imputation()], `fixest` regressions with `i(event_time, ref = ...)` or
#' `sunab()` terms, or plain data frames, and returns one tidy table that
#' [plot_event_study()] can draw. Reference periods omitted from a regression
#' are added back as zeros.
#'
#' @param ... Named objects to combine; the names become the `method`
#'   column. Data frames must have columns `event_time`, `estimate`,
#'   `std.error`.
#' @param ref Event time(s) omitted as the reference in the regressions
#'   (default `-1`), added as zero rows.
#' @param conf_level Confidence level for pointwise intervals computed from
#'   standard errors.
#'
#' @return A data frame with `method`, `event_time`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `band.low`, `band.high` (uniform band when the
#'   source provides one, otherwise `NA`).
#' @examples
#' dat <- sim_did_panel(n_units = 300, n_periods = 8, groups = c(3, 6), seed = 1)
#' cs <- aggregate_att(att_gt(dat, "id", "time", "g", "y", n_boot = 99), type = "dynamic")
#' dat$rel <- ifelse(dat$g > 0, dat$time - dat$g, -1000)
#' tw <- fixest::feols(y ~ i(rel, ref = c(-1, -1000)) | id + time, data = dat, cluster = ~id)
#' frame <- event_study_frame("Callaway-Sant'Anna" = cs, "TWFE event study" = tw)
#' plot_event_study(frame)
#' @export
event_study_frame <- function(..., ref = -1, conf_level = 0.95) {
  objs <- list(...)
  if (length(objs) == 0L) stop("Supply at least one object.", call. = FALSE)
  if (is.null(names(objs)) || any(!nzchar(names(objs)))) names(objs) <- paste0("method_", seq_along(objs))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  one <- function(obj, nm) {
    if (inherits(obj, "cm_agg_att")) {
      if (obj$type != "dynamic") stop("aggregate_att() objects must use type = \"dynamic\".", call. = FALSE)
      d <- data.frame(event_time = obj$by$event_time, estimate = obj$by$estimate, std.error = obj$by$std.error,
                      band.low = obj$by$band.low, band.high = obj$by$band.high)
    } else if (inherits(obj, "cm_att_gt")) {
      if (length(obj$groups) != 1L) stop("att_gt() objects with several cohorts must be aggregated first.", call. = FALSE)
      d <- data.frame(event_time = obj$att_gt$event_time, estimate = obj$att_gt$att, std.error = obj$att_gt$std.error,
                      band.low = obj$att_gt$band.low, band.high = obj$att_gt$band.high)
    } else if (inherits(obj, "cm_did_imputation")) {
      d <- rbind(
        if (!is.null(obj$pre)) data.frame(event_time = obj$pre$event_time, estimate = obj$pre$estimate, std.error = obj$pre$std.error) else NULL,
        data.frame(event_time = obj$by_event$event_time, estimate = obj$by_event$estimate, std.error = obj$by_event$std.error)
      )
      d$band.low <- NA_real_; d$band.high <- NA_real_
      d <- rbind(d, data.frame(event_time = ref, estimate = 0, std.error = 0, band.low = NA, band.high = NA))
    } else if (inherits(obj, "fixest")) {
      ct <- fixest::coeftable(obj)
      nm_coef <- rownames(ct)
      hit <- grepl("::-?[0-9]+(\\.[0-9]+)?$", nm_coef)
      if (!any(hit)) stop("No event-time coefficients (name::value) found in the fixest object `", nm, "`.", call. = FALSE)
      et <- as.numeric(sub("^.*::", "", nm_coef[hit]))
      d <- data.frame(event_time = et, estimate = unname(ct[hit, 1L]), std.error = unname(ct[hit, 2L]),
                      band.low = NA_real_, band.high = NA_real_)
      missing_ref <- setdiff(ref, d$event_time)
      if (length(missing_ref)) d <- rbind(d, data.frame(event_time = missing_ref, estimate = 0, std.error = 0, band.low = NA, band.high = NA))
    } else if (is.data.frame(obj)) {
      for (v in c("event_time", "estimate", "std.error")) if (!v %in% names(obj)) stop("Data frames need columns event_time, estimate, std.error.", call. = FALSE)
      d <- data.frame(event_time = obj$event_time, estimate = obj$estimate, std.error = obj$std.error,
                      band.low = if ("band.low" %in% names(obj)) obj$band.low else NA_real_,
                      band.high = if ("band.high" %in% names(obj)) obj$band.high else NA_real_)
    } else {
      stop("Unsupported object of class ", class(obj)[1L], " for `", nm, "`.", call. = FALSE)
    }
    d$conf.low <- d$estimate - z * d$std.error
    d$conf.high <- d$estimate + z * d$std.error
    d <- d[order(d$event_time), ]
    data.frame(method = nm, d[, c("event_time", "estimate", "std.error", "conf.low", "conf.high", "band.low", "band.high")],
               stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, Map(one, objs, names(objs)))
  rownames(out) <- NULL
  out$method <- factor(out$method, levels = names(objs))
  out
}

#' Plot event-study estimates
#'
#' Draws the estimates of [event_study_frame()] with pointwise intervals as
#' error bars and, where available, uniform bands as ribbons. Methods are
#' distinguished by colour and slightly dodged.
#'
#' @param frame Output of [event_study_frame()], or a single object it accepts.
#' @param ref_line Event time at which to draw the vertical treatment line
#'   (default `-0.5`, between the last pre and first post period).
#' @param dodge Horizontal dodge between methods.
#' @param band Draw uniform bands when present.
#' @param ... Passed to [event_study_frame()] when `frame` is not a frame.
#' @return A `ggplot` object.
#' @export
plot_event_study <- function(frame, ref_line = -0.5, dodge = 0.3, band = TRUE, ...) {
  if (!is.data.frame(frame) || !all(c("method", "event_time", "estimate") %in% names(frame))) {
    frame <- event_study_frame(estimate = frame, ...)
  }
  n_m <- length(unique(frame$method))
  pd <- ggplot2::position_dodge(width = if (n_m > 1L) dodge else 0)
  p <- ggplot2::ggplot(frame, ggplot2::aes(x = .data$event_time, y = .data$estimate, colour = .data$method, group = .data$method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_vline(xintercept = ref_line, linetype = "dotted", colour = "grey50")
  if (band && any(is.finite(frame$band.low))) {
    p <- p + ggplot2::geom_ribbon(
      data = frame[is.finite(frame$band.low), ],
      ggplot2::aes(ymin = .data$band.low, ymax = .data$band.high, fill = .data$method),
      alpha = 0.12, colour = NA, position = pd
    )
  }
  p <- p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), width = 0.2, position = pd) +
    ggplot2::geom_point(position = pd, size = 2) +
    ggplot2::geom_line(position = pd, alpha = 0.6) +
    ggplot2::labs(x = "Event time (periods since treatment)", y = "Estimate", colour = NULL, fill = NULL) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  if (n_m == 1L) p <- p + ggplot2::theme(legend.position = "none")
  p
}
