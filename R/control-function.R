# R/control-function.R
#
# Control-function helpers: first-stage residuals (or generalized residuals
# for a binary endogenous regressor), bootstrap inference for the two-step
# estimator, and average partial effects from the average structural
# function. The second stage stays in lm(), glm(), feols(), or mlogit().

#' First-stage residuals for a control-function second stage
#'
#' Fits the first stage of an endogenous regressor `d` on the instruments
#' `z` and controls `x` and appends the control to the data: the OLS residual
#' for a continuous `d` (`family = "gaussian"`), or the generalized residual
#' of a probit or logit first stage for a binary `d`
#' (`D lambda(xb) - (1 - D) lambda(-xb)` with `lambda` the inverse Mills
#' ratio for the probit; `D - Lambda(xb)` for the logit). Adding the control
#' to a linear second stage reproduces two-stage least squares exactly; adding
#' it to a probit, logit, Poisson, or discrete-choice second stage is the
#' control-function estimator of Rivers and Vuong (1988), Smith and Blundell
#' (1986), Wooldridge (2015), and Petrin and Train (2010).
#'
#' @param data A data frame.
#' @param d Endogenous regressor column name.
#' @param z Character vector of instrument column names.
#' @param x Optional control column names.
#' @param family `"gaussian"` (OLS residual), `"probit"`, or `"logit"`
#'   (generalized residual).
#' @param name Name of the control column.
#' @return `data` with the control added; attribute `"first_stage"` holds the
#'   fitted first stage and, for the Gaussian case, the first-stage F
#'   statistics of [iv_first_stage()].
#' @references
#' Wooldridge, J. M. (2015). Control function methods in applied
#' econometrics. *Journal of Human Resources*, 50(2), 420-445.
#' @examples
#' dat <- sim_iv(1000, dgp = "linear", seed = 1)
#' cf <- cf_residuals(dat, d = "d", z = "z", x = paste0("x", 1:5))
#' coef(lm(y ~ d + x1 + x2 + x3 + x4 + x5 + v_hat, data = cf))["d"]
#' @export
cf_residuals <- function(data, d, z, x = NULL, family = c("gaussian", "probit", "logit"), name = "v_hat") {
  family <- match.arg(family)
  data <- as.data.frame(data)
  for (v in c(d, z, x)) .cm_check_column(v, data)
  if (name %in% names(data)) warning("Column `", name, "` is overwritten.", call. = FALSE)
  rhs <- c(z, x)
  if (family == "gaussian") {
    fit <- stats::lm(stats::reformulate(rhs, d), data = data, na.action = stats::na.exclude)
    v <- stats::residuals(fit)
    fs <- iv_first_stage(data[stats::complete.cases(data[, c(d, z, x), drop = FALSE]), , drop = FALSE], d, z, x)
  } else {
    dv <- .cm_as_binary(data[[d]], d)
    df <- data; df[[d]] <- dv
    link <- if (family == "probit") "probit" else "logit"
    fit <- stats::glm(stats::reformulate(rhs, d), data = df, family = stats::binomial(link = link), na.action = stats::na.exclude)
    xb <- stats::predict(fit, type = "link")
    v <- if (family == "probit") {
      dv * stats::dnorm(xb) / stats::pnorm(xb) - (1 - dv) * stats::dnorm(xb) / stats::pnorm(-xb)
    } else {
      dv - stats::plogis(xb)
    }
    fs <- NULL
  }
  data[[name]] <- as.numeric(v)
  attr(data, "first_stage") <- list(fit = fit, family = family, strength = fs, control = name)
  data
}

#' Bootstrap inference for a two-step control-function estimator
#'
#' Re-runs both stages on bootstrap samples (pairs bootstrap, or cluster
#' bootstrap when `cluster` is given) and reports bootstrap standard errors
#' and percentile intervals for the second-stage coefficients and, if
#' requested, for average partial effects. Second-stage standard errors that
#' treat the control as a fixed regressor are wrong (generated-regressor
#' problem); the bootstrap is the simplest valid alternative.
#'
#' @param data A data frame.
#' @param first A function of a data frame returning the data with the
#'   control added (typically a call to [cf_residuals()]).
#' @param second A function of the augmented data returning a fitted model
#'   (anything with a `coef()` method) or a named numeric vector.
#' @param n_boot Bootstrap replications.
#' @param cluster Optional cluster column name for a cluster bootstrap.
#' @param seed Optional seed.
#' @param ape Optional function of `(fit, data)` returning a named numeric
#'   vector of average partial effects (for instance a call to [cf_ape()]).
#' @param conf_level Confidence level for the percentile intervals.
#' @return A list of class `cm_cf_boot` with `table` (term, estimate,
#'   std.error, conf.low, conf.high, component), the matrix of bootstrap
#'   draws, and the number of failed replications.
#' @examples
#' dat <- sim_iv(600, dgp = "probit_cf", seed = 1)
#' first <- function(df) cf_residuals(df, d = "d", z = "z", x = "x1")
#' second <- function(df) glm(y ~ d + x1 + v_hat, data = df, family = binomial("probit"))
#' cf_bootstrap(dat, first, second, n_boot = 30, seed = 1,
#'              ape = function(fit, df) cf_ape(fit, df, d = "d"))
#' @export
cf_bootstrap <- function(data, first, second, n_boot = 499L, cluster = NULL, seed = NULL, ape = NULL, conf_level = 0.95) {
  data <- as.data.frame(data)
  if (!is.null(cluster)) .cm_check_column(cluster, data)
  extract <- function(df) {
    aug <- first(df)
    fit <- second(aug)
    cf <- if (is.numeric(fit)) fit else stats::coef(fit)
    cf <- cf[!is.na(cf)]
    out <- c(cf)
    comp <- rep("coefficient", length(cf))
    if (!is.null(ape)) {
      a <- ape(fit, aug)
      a <- stats::setNames(as.numeric(a), names(a))
      out <- c(out, a)
      comp <- c(comp, rep("ape", length(a)))
    }
    list(values = out, component = comp)
  }
  pt <- extract(data)
  n <- nrow(data)
  cl <- if (is.null(cluster)) NULL else data[[cluster]]
  draws <- .cm_with_seed(seed, lapply(seq_len(n_boot), function(b) {
    idx <- if (is.null(cl)) sample.int(n, n, replace = TRUE) else {
      u <- unique(cl); pick <- sample(u, length(u), replace = TRUE)
      unlist(lapply(pick, function(g) which(cl == g)))
    }
    tryCatch(extract(data[idx, , drop = FALSE])$values, error = function(e) NULL)
  }))
  ok <- !vapply(draws, is.null, logical(1))
  B <- do.call(rbind, draws[ok])
  B <- B[, names(pt$values), drop = FALSE]
  se <- apply(B, 2, stats::sd)
  a <- (1 - conf_level) / 2
  lo <- apply(B, 2, stats::quantile, probs = a, na.rm = TRUE)
  hi <- apply(B, 2, stats::quantile, probs = 1 - a, na.rm = TRUE)
  tab <- data.frame(term = names(pt$values), estimate = unname(pt$values), std.error = unname(se),
                    conf.low = unname(lo), conf.high = unname(hi), component = pt$component, stringsAsFactors = FALSE)
  structure(list(table = tab, draws = B, n_boot = n_boot, n_failed = sum(!ok), conf_level = conf_level,
                 cluster = cluster), class = "cm_cf_boot")
}

#' @export
print.cm_cf_boot <- function(x, ...) {
  cat("Two-step control-function estimator, ", x$n_boot - x$n_failed, " bootstrap replications",
      if (!is.null(x$cluster)) paste0(" (clustered by ", x$cluster, ")"), "\n", sep = "")
  print(x$table, digits = 4, row.names = FALSE)
  invisible(x)
}

#' Average partial effects from a control-function second stage
#'
#' For a fitted nonlinear second stage (a `glm` with a probit, logit, or
#' Poisson family, or any model with a `predict(newdata, type = "response")`
#' method), computes the average partial effect of `d` from the average
#' structural function: the derivative of the mean prediction with respect
#' to `d`, averaged over the observed covariates and the estimated control,
#' or the mean difference of predictions for a discrete change `delta`.
#' Averaging over the control is what makes the quantity structural
#' (Blundell and Powell 2004; Wooldridge 2015).
#'
#' @param fit The fitted second stage.
#' @param data The augmented data used to fit it.
#' @param d Name of the endogenous regressor.
#' @param delta Optional discrete change in `d` (default: derivative).
#' @param h Step for the numerical derivative (default `1e-4 sd(d)`).
#' @return A named numeric vector with the APE (name `paste0("ape_", d)`) and
#'   attribute `"contributions"` holding the per-observation effects.
#' @export
cf_ape <- function(fit, data, d, delta = NULL, h = NULL) {
  data <- as.data.frame(data)
  .cm_check_column(d, data)
  pred <- function(df) as.numeric(stats::predict(fit, newdata = df, type = "response"))
  if (is.null(delta)) {
    if (is.null(h)) h <- 1e-4 * max(stats::sd(data[[d]]), 1e-8)
    up <- data; up[[d]] <- up[[d]] + h
    dn <- data; dn[[d]] <- dn[[d]] - h
    contrib <- (pred(up) - pred(dn)) / (2 * h)
  } else {
    up <- data; up[[d]] <- up[[d]] + delta
    contrib <- pred(up) - pred(data)
  }
  out <- stats::setNames(mean(contrib, na.rm = TRUE), paste0("ape_", d))
  attr(out, "contributions") <- contrib
  out
}
