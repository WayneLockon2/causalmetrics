#' Inputs for a HonestDiD sensitivity analysis
#'
#' Rambachan and Roth (2023) sensitivity analysis needs the event-study
#' coefficients with the reference period removed and their covariance
#' matrix. `honestdid_inputs()` extracts both from a dynamic
#' [aggregate_att()] object, with the covariance computed from the influence
#' functions (clustered when the underlying `att_gt()` was), so the result
#' can be passed directly to `HonestDiD::createSensitivityResults()` and
#' `HonestDiD::createSensitivityResults_relativeMagnitudes()`. Use
#' `base_period = "universal"` in [att_gt()], so the reference coefficient
#' at event time `-1` is exactly zero and every other pre-period is a
#' genuine pre-trend estimate.
#'
#' @param x A dynamic [aggregate_att()] object.
#' @param ref Reference event time to drop (default `-1`).
#' @return A list with `betahat`, `sigma`, `numPrePeriods`, `numPostPeriods`,
#'   and `event_time`.
#' @references Rambachan, A. and Roth, J. (2023). A more credible approach
#'   to parallel trends. *Review of Economic Studies*, 90(5), 2555-2591.
#' @export
honestdid_inputs <- function(x, ref = -1) {
  if (!inherits(x, "cm_agg_att") || x$type != "dynamic") stop("`x` must be a dynamic aggregate_att() object.", call. = FALSE)
  et <- x$by$event_time
  keep <- et != ref
  IFm <- as.matrix(x$inffunc)[, keep, drop = FALSE]
  cs <- .cm_if_cluster_sums(IFm, x$cluster)
  sigma <- crossprod(as.matrix(cs$mat)) / x$n^2
  list(betahat = x$by$estimate[keep], sigma = sigma,
       numPrePeriods = sum(et[keep] < ref), numPostPeriods = sum(et[keep] > ref),
       event_time = et[keep])
}
