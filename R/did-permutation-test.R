#' Permutation (randomization) test for difference-in-differences
#'
#' With few treated clusters, cluster-robust standard errors and the
#' bootstrap are unreliable. A Fisher randomization test instead re-assigns
#' the treatment timing across clusters and recomputes a studentized
#' statistic; under random timing it is exact for the sharp null of no effect
#' for any unit, and with a studentized statistic it is asymptotically valid
#' for the weak null of a zero average effect (Roth and Sant'Anna 2023).
#'
#' @inheritParams att_gt
#' @param group Column with the first treatment period (0, `NA`, or `Inf`
#'   for never treated), constant within cluster.
#' @param cluster Column at whose level treatment is assigned and permuted
#'   (default: the unit).
#' @param statistic Function of a data frame returning a scalar. The default
#'   is the t-statistic of the static two-way fixed effects coefficient with
#'   standard errors clustered by `cluster`.
#' @param n_perm Number of permutations.
#'
#' @return An object of class `cm_did_permutation` with the observed
#'   statistic, the two-sided permutation p-value, and the permutation
#'   distribution.
#' @references Roth, J. and Sant'Anna, P. H. C. (2023). Efficient estimation
#'   for staggered rollout designs. *Journal of Political Economy:
#'   Microeconomics*, 1(4), 669-709.
#' @examples
#' dat <- sim_did_panel(n_units = 60, n_periods = 6, groups = 4, seed = 1)
#' did_permutation_test(dat, id = "id", time = "time", group = "g", y = "y", n_perm = 199)
#' @export
did_permutation_test <- function(data, id, time, group, y, cluster = NULL, statistic = NULL,
                                 n_perm = 999L, seed = NULL) {
  for (v in c(id, time, group, y, cluster)) .cm_check_column(v, data)
  dt <- data.table::as.data.table(data)[, unique(c(id, time, group, y, cluster)), with = FALSE]
  data.table::setnames(dt, c(id, time, group, y), c("id", "time", "g", "y"))
  if (is.null(cluster)) { dt[, cl := id] } else { data.table::setnames(dt, cluster, "cl") }
  dt[is.na(g) | g == 0, g := Inf]
  chk <- dt[, data.table::uniqueN(g), by = cl]
  if (any(chk$V1 > 1L)) stop("`group` must be constant within `cluster`.", call. = FALSE)
  if (is.null(statistic)) {
    statistic <- function(d) {
      d$d <- as.numeric(is.finite(d$g) & d$time >= d$g)
      fit <- fixest::feols(y ~ d | id + time, data = d, cluster = ~cl)
      unname(fixest::coeftable(fit)["d", "t value"])
    }
  }
  obs <- statistic(as.data.frame(dt))
  cl_tab <- unique(dt[, list(cl, g)])
  perm <- .cm_with_seed(seed, vapply(seq_len(n_perm), function(b) {
    g_perm <- sample(cl_tab$g)
    d <- data.table::copy(dt)
    d[, g := g_perm[match(cl, cl_tab$cl)]]
    tryCatch(statistic(as.data.frame(d)), error = function(e) NA_real_)
  }, numeric(1)))
  perm <- perm[is.finite(perm)]
  p <- (1 + sum(abs(perm) >= abs(obs))) / (1 + length(perm))
  out <- list(statistic = obs, p.value = p, permutations = perm, n_perm = length(perm),
              n_clusters = nrow(cl_tab), n_treated_clusters = sum(is.finite(cl_tab$g)), call = match.call())
  class(out) <- "cm_did_permutation"
  out
}

#' @export
print.cm_did_permutation <- function(x, digits = 3, ...) {
  cat("Permutation test for difference-in-differences (treatment timing permuted across ",
      x$n_clusters, " clusters, ", x$n_treated_clusters, " treated)\n", sep = "")
  cat("  Observed statistic: ", formatC(x$statistic, digits = digits, format = "f"),
      "; two-sided permutation p-value: ", formatC(x$p.value, digits = digits, format = "f"),
      " (", x$n_perm, " permutations)\n", sep = "")
  invisible(x)
}
