#' Power of a pre-trends test and the bias conditional on passing it
#'
#' Roth (2022) shows that an event-study pre-test has low power against
#' violations of parallel trends that would badly bias the post-treatment
#' estimates, and that conditioning on passing the test makes the bias
#' worse. `pretrend_power()` quantifies both for a linear violation: given
#' event-study coefficients, their covariance, and a hypothesized slope, it
#' simulates the estimates under the trend, computes the probability that
#' the standard pre-test rejects (any pre-treatment coefficient significant
#' at level `alpha`), and reports the bias of each post-treatment coefficient
#' unconditionally and conditional on not rejecting. It also solves for the
#' slope that the pre-test would detect with a given power.
#'
#' @param estimate Numeric vector of event-study coefficients, or an object
#'   from [aggregate_att()] with `type = "dynamic"`.
#' @param vcov Covariance matrix of the coefficients (ignored when
#'   `estimate` is an `aggregate_att()` object, whose influence functions
#'   supply it).
#' @param event_time Event times of the coefficients.
#' @param ref Reference event time through which the linear trend passes.
#' @param slope Slope(s) of the hypothesized linear violation, in outcome
#'   units per period. `NULL` reports only the slopes detected with the
#'   target powers.
#' @param target_power Powers for which to solve the detectable slope.
#' @param alpha Significance level of the pre-test.
#' @param n_sim Number of simulation draws.
#' @param seed Optional seed.
#'
#' @return A list with `slopes` (for each requested slope: power of the
#'   pre-test, and for each post-treatment event time the unconditional and
#'   conditional bias), `detectable` (slopes detected with the target powers),
#'   and the inputs.
#' @references Roth, J. (2022). Pretest with caution: Event-study estimates
#'   after testing for parallel trends. *American Economic Review: Insights*,
#'   4(3), 305-322.
#' @export
pretrend_power <- function(estimate, vcov = NULL, event_time = NULL, ref = -1, slope = NULL,
                           target_power = c(0.5, 0.8), alpha = 0.05, n_sim = 5000L, seed = NULL) {
  if (inherits(estimate, "cm_agg_att")) {
    obj <- estimate
    if (obj$type != "dynamic") stop("Use an aggregate_att() object with type = \"dynamic\".", call. = FALSE)
    event_time <- obj$by$event_time
    IFm <- as.matrix(obj$inffunc)
    cs <- .cm_if_cluster_sums(IFm, obj$cluster)
    vcov <- crossprod(as.matrix(cs$mat)) / obj$n^2
    estimate <- obj$by$estimate
  }
  if (is.null(event_time) || length(event_time) != length(estimate)) stop("`event_time` must match `estimate`.", call. = FALSE)
  vcov <- as.matrix(vcov)
  if (any(dim(vcov) != length(estimate))) stop("`vcov` must be square with the dimension of `estimate`.", call. = FALSE)
  pre <- which(event_time < ref)
  post <- which(event_time >= 0)
  if (length(pre) == 0L) stop("No pre-treatment coefficients (event_time < ref).", call. = FALSE)
  se <- sqrt(diag(vcov))
  z <- stats::qnorm(1 - alpha / 2)
  L <- t(chol(vcov + diag(1e-12, nrow(vcov))))
  draws <- .cm_with_seed(seed, matrix(stats::rnorm(n_sim * length(estimate)), n_sim))
  base <- draws %*% t(L)  # n_sim x K, mean zero
  evaluate <- function(s) {
    delta <- s * (event_time - ref)
    sim <- sweep(base, 2, delta, "+")
    pass <- rowSums(abs(sim[, pre, drop = FALSE]) > matrix(z * se[pre], n_sim, length(pre), byrow = TRUE)) == 0
    power <- 1 - mean(pass)
    cond <- if (any(pass)) colMeans(sim[pass, post, drop = FALSE]) else rep(NA_real_, length(post))
    list(power = power, bias = data.frame(event_time = event_time[post], unconditional_bias = delta[post],
                                          conditional_bias = cond, stringsAsFactors = FALSE))
  }
  power_of <- function(s) evaluate(s)$power
  detectable <- vapply(target_power, function(tp) {
    hi <- 1
    while (power_of(hi) < tp && hi < 1e6) hi <- hi * 2
    if (power_of(hi) < tp) return(NA_real_)
    stats::uniroot(function(s) power_of(s) - tp, c(0, hi), tol = 1e-4)$root
  }, numeric(1))
  slopes <- if (is.null(slope)) NULL else lapply(slope, function(s) c(list(slope = s), evaluate(s)))
  out <- list(slopes = slopes,
              detectable = data.frame(target_power = target_power, slope = detectable),
              estimate = estimate, event_time = event_time, se = se, ref = ref, alpha = alpha,
              n_pre = length(pre), n_post = length(post))
  class(out) <- "cm_pretrend_power"
  out
}

#' @export
print.cm_pretrend_power <- function(x, digits = 3, ...) {
  cat("Pre-trends test power (Roth 2022): ", x$n_pre, " pre-treatment coefficients tested at level ", x$alpha, "\n", sep = "")
  cat("Linear violations detected with the target power:\n")
  print(format(x$detectable, digits = digits), row.names = FALSE)
  for (s in x$slopes) {
    cat("\nSlope ", formatC(s$slope, digits = digits, format = "g"), " per period: power ",
        formatC(s$power, digits = 2, format = "f"), "\n", sep = "")
    print(format(s$bias, digits = digits), row.names = FALSE)
  }
  invisible(x)
}
