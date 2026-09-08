# R/est-dml-structured.R
#
# Double machine learning for structured outcome models with combinatorial
# treatments (Farrell, Liang, and Misra 2020; Ye et al. 2025). The nuisance is
# a parameter function theta(x) estimated by any engine (out of fold); the
# package owns the influence function, the cross-fitted estimates, and the
# inference.

#' Double machine learning for a structured outcome model
#'
#' Estimates average treatment effects of treatment combinations from a
#' structured model `E[Y | X = x, T = t] = G(theta(x), t)` with a known link
#' `G` and a nuisance parameter function `theta(x)` estimated out of fold by
#' any learner (a deep network, a forest, a regression). The estimator is the
#' cross-fitted mean of the influence function of Farrell, Liang, and Misra
#' (2020), adapted by Ye et al. (2025) to combinatorial experiments:
#' \deqn{\psi_i(t) = H(x_i, \hat\theta; t, t_0) - H_\theta(x_i, \hat\theta; t, t_0)'
#'   \Lambda(x_i)^{-1} \ell_\theta(y_i, \check t_i, \hat\theta),}
#' with `H = G(theta, t) - G(theta, t0)`, `H_theta` its gradient in `theta`,
#' `Lambda(x) = 2 sum_t nu(t) G_theta(theta, t) G_theta(theta, t)'` over the
#' known assignment distribution `nu` of the observed combinations, and
#' `l_theta = 2 G_theta(theta, t_obs) (G(theta, t_obs) - y)` the gradient of the
#' squared loss. The first term alone is the plug-in (structured learning)
#' estimator; the second term removes the first-order effect of nuisance
#' estimation error (Neyman orthogonality), so the mean of `psi` is
#' root-n normal when `theta_hat` converges faster than `n^-1/4`.
#'
#' @param theta_hat Numeric matrix, `n` rows, one column per nuisance
#'   parameter, of out-of-fold estimates of `theta(x_i)`.
#' @param y Outcome vector.
#' @param t_obs Realized treatment combinations: an `n x m` 0/1 matrix (or a
#'   data frame of indicator columns).
#' @param targets Treatment combinations whose effect against `t0` is wanted:
#'   a `K x m` 0/1 matrix, a data frame, or a character vector such as
#'   `c("110", "101")`.
#' @param t0 Reference combination (default all zeros).
#' @param link `"gen_sigmoid"` (generalized sigmoid form II of Ye et al.:
#'   `theta = (theta_0, theta_1..theta_m, c)`, `G = c / (1 + exp(-(theta_0 +
#'   theta' t)))`), `"linear"` (`theta_0 + theta' t`), `"logit"`
#'   (`1 / (1 + exp(-(theta_0 + theta' t)))`, Farrell et al. 2020), or
#'   `"custom"` with `G` and `G_grad`.
#' @param G,G_grad For `link = "custom"`: functions `(theta, t)` where `theta`
#'   is the `n x d` matrix and `t` a single combination (length `m`) or an
#'   `n x m` matrix, returning a length-`n` vector and an `n x d` gradient
#'   matrix.
#' @param t_dist Assignment distribution of the observed combinations used in
#'   `Lambda(x)`: a data frame with the `m` indicator columns and a `prob`
#'   column, a matrix with a `prob` attribute, or `NULL` for the empirical
#'   shares of `t_obs`. Probabilities that do not sum to one are accepted with
#'   a warning and used as weights (this reproduces some published code).
#' @param fold_id Optional fold identifier of the cross-fitting that produced
#'   `theta_hat`; estimates are reported fold by fold and pooled.
#' @param ridge Added to the diagonal of `Lambda(x)` before inversion (the
#'   device used in the authors' code; the theory uses 0).
#' @param conf_level Confidence level.
#' @param contrast_best Also estimate `tau(t) = mu(t*) - mu(t)`, the gap to
#'   the best target, from the difference of the influence functions.
#'
#' @return A list of class `cm_dml_structured`: `estimates` (target, fold
#'   including `"pooled"`, n, `estimate`, `plugin`, `std.error`, `conf.low`,
#'   `conf.high`, and for the pooled rows `ci_low_avg`, `ci_high_avg`, the
#'   averages of the fold-level bounds as in Ye et al.'s Table 7),
#'   `contrasts` (the gap to the best target), `psi` (`n x K` influence
#'   functions, for group averages), `plugin` (`n x K` plug-in terms),
#'   `best`, `targets`, `t_dist`, `fold_id`, `diagnostics` (share of rows whose
#'   `Lambda(x)` was ill conditioned), and the settings. Methods: [print()],
#'   [tidy()].
#'
#' @details The pooled estimate is the mean of the fold estimates. Its
#'   standard error treats the folds as independent, `sqrt(sum se_s^2) / S`,
#'   which coincides with the influence-function standard error of the
#'   stacked `psi` when the folds are of equal size. `ci_low_avg` and
#'   `ci_high_avg` average the fold bounds instead and are wider by about
#'   `sqrt(S)`; they are reported because the paper's tables use them.
#'
#' @references
#' Farrell, M. H., Liang, T., and Misra, S. (2020). Deep learning for
#' individual heterogeneity: an automatic inference framework.
#' arXiv:2010.14694.
#'
#' Ye, Z., Zhang, Z., Zhang, D. J., Zhang, H., and Zhang, R. (2025). Deep
#' learning-based causal inference for large-scale combinatorial experiments:
#' theory and empirical evidence. *Management Science*.
#'
#' @examples
#' set.seed(1)
#' n <- 2000; m <- 2
#' x <- matrix(runif(n * 3), n)
#' theta <- cbind(x[, 1] - 0.5, 1 + x[, 2], -0.5 + x[, 3], 4)
#' t_obs <- cbind(rbinom(n, 1, 0.5), rbinom(n, 1, 0.5))
#' u <- theta[, 1] + theta[, 2] * t_obs[, 1] + theta[, 3] * t_obs[, 2]
#' y <- theta[, 4] * plogis(u) + rnorm(n, sd = 0.1)
#' # with the true theta the correction has mean zero and the estimator is the
#' # plug-in contrast
#' fit <- est_dml_structured(theta, y, t_obs, targets = c("10", "01", "11"))
#' fit
#' @export
est_dml_structured <- function(theta_hat, y, t_obs, targets, t0 = NULL,
                               link = c("gen_sigmoid", "linear", "logit", "custom"),
                               G = NULL, G_grad = NULL, t_dist = NULL, fold_id = NULL,
                               ridge = 0, conf_level = 0.95, contrast_best = TRUE) {
  link <- match.arg(link)
  theta_hat <- as.matrix(theta_hat)
  storage.mode(theta_hat) <- "double"
  y <- as.numeric(y)
  n <- length(y)
  if (nrow(theta_hat) != n) stop("`theta_hat` must have one row per observation.", call. = FALSE)
  t_obs <- .cm_as_t_matrix(t_obs, n = n)
  m <- ncol(t_obs)
  targets <- .cm_as_t_matrix(targets, m = m)
  if (is.null(t0)) t0 <- rep(0, m)
  t0 <- as.numeric(t0)
  if (length(t0) != m) stop("`t0` must have length ", m, ".", call. = FALSE)
  lk <- .cm_structured_link(link, m, ncol(theta_hat), G, G_grad)
  d <- lk$d
  if (ncol(theta_hat) != d) stop("`theta_hat` must have ", d, " columns for link \"", link, "\".", call. = FALSE)
  if (is.null(fold_id)) fold_id <- rep(1L, n)
  fold_id <- as.integer(factor(fold_id))
  if (length(fold_id) != n) stop("`fold_id` must have length n.", call. = FALSE)
  dist <- .cm_t_dist(t_dist, t_obs, m)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Lambda(x) = 2 sum_t nu(t) G_theta G_theta' , one d x d block per row.
  Lmat <- matrix(0, n, d * d)
  for (j in seq_len(nrow(dist$t))) {
    g <- lk$grad(theta_hat, dist$t[j, ])
    Lmat <- Lmat + 2 * dist$prob[j] * g[, rep(seq_len(d), times = d), drop = FALSE] * g[, rep(seq_len(d), each = d), drop = FALSE]
  }
  if (ridge != 0) {
    diag_idx <- (seq_len(d) - 1L) * d + seq_len(d)
    Lmat[, diag_idx] <- Lmat[, diag_idx] + ridge
  }
  # loss gradient at the realized treatment
  ell <- 2 * lk$grad(theta_hat, t_obs) * (lk$G(theta_hat, t_obs) - y)
  v <- .cm_blockwise_solve(Lmat, ell, d)
  rc <- .cm_block_rcond(Lmat, d)

  K <- nrow(targets)
  labels <- apply(targets, 1, function(r) paste0("(", paste(r, collapse = ", "), ")"))
  G0 <- lk$G(theta_hat, t0)
  grad0 <- lk$grad(theta_hat, t0)
  psi <- matrix(NA_real_, n, K, dimnames = list(NULL, labels))
  plug <- matrix(NA_real_, n, K, dimnames = list(NULL, labels))
  for (k in seq_len(K)) {
    H <- lk$G(theta_hat, targets[k, ]) - G0
    Hg <- lk$grad(theta_hat, targets[k, ]) - grad0
    plug[, k] <- H
    psi[, k] <- H - rowSums(Hg * v)
  }

  est <- .cm_structured_fold_table(psi, plug, fold_id, labels, crit)
  best <- NULL
  contrasts <- NULL
  if (contrast_best && K > 1L) {
    pooled <- est[est$fold == "pooled", ]
    best <- pooled$target[which.max(pooled$estimate)]
    kb <- match(best, labels)
    dpsi <- psi[, kb] - psi
    dplug <- plug[, kb] - plug
    contrasts <- .cm_structured_fold_table(dpsi, dplug, fold_id, labels, crit)
    contrasts$p.value <- 2 * stats::pnorm(-abs(contrasts$estimate / contrasts$std.error))
    contrasts$p.value[contrasts$target == best] <- NA_real_
  }

  structure(list(
    estimates = est, contrasts = contrasts, best = best, psi = psi, plugin = plug,
    targets = targets, labels = labels, t0 = t0, t_dist = dist, fold_id = fold_id,
    link = link, d = d, m = m, n = n, ridge = ridge, conf_level = conf_level,
    diagnostics = list(share_ill_conditioned = mean(rc < 1e-8), rcond_summary = .cm_summary(rc),
                       ridge = ridge),
    call = match.call()
  ), class = "cm_dml_structured")
}

# Fold-by-fold and pooled tables for an n x K matrix of scores.
.cm_structured_fold_table <- function(psi, plug, fold_id, labels, crit) {
  folds <- sort(unique(fold_id))
  S <- length(folds)
  rows <- list()
  for (k in seq_along(labels)) {
    fe <- fs <- fp <- fn <- numeric(S)
    for (s in seq_len(S)) {
      i <- fold_id == folds[s]
      fn[s] <- sum(i)
      fe[s] <- mean(psi[i, k])
      fp[s] <- mean(plug[i, k])
      fs[s] <- stats::sd(psi[i, k]) / sqrt(fn[s])
      rows[[length(rows) + 1L]] <- data.frame(
        target = labels[k], fold = as.character(folds[s]), n = fn[s], estimate = fe[s], plugin = fp[s],
        std.error = fs[s], conf.low = fe[s] - crit * fs[s], conf.high = fe[s] + crit * fs[s],
        ci_low_avg = NA_real_, ci_high_avg = NA_real_, stringsAsFactors = FALSE)
    }
    pe <- mean(fe)
    pse <- sqrt(sum(fs^2)) / S
    rows[[length(rows) + 1L]] <- data.frame(
      target = labels[k], fold = "pooled", n = sum(fn), estimate = pe, plugin = mean(fp),
      std.error = pse, conf.low = pe - crit * pse, conf.high = pe + crit * pse,
      ci_low_avg = mean(fe - crit * fs), ci_high_avg = mean(fe + crit * fs), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Links --------------------------------------------------------------------------

.cm_structured_link <- function(link, m, d_given, G = NULL, G_grad = NULL) {
  tmat <- function(t, n) {
    if (is.matrix(t)) return(t)
    matrix(t, n, length(t), byrow = TRUE)
  }
  index <- function(theta, t) {
    tt <- tmat(t, nrow(theta))
    theta[, 1] + rowSums(theta[, 1 + seq_len(ncol(tt)), drop = FALSE] * tt)
  }
  if (link == "gen_sigmoid") {
    return(list(
      d = m + 2L,
      G = function(theta, t) theta[, m + 2L] * stats::plogis(index(theta, t)),
      grad = function(theta, t) {
        tt <- tmat(t, nrow(theta))
        s <- stats::plogis(index(theta, t))
        ds <- theta[, m + 2L] * s * (1 - s)
        cbind(ds, ds * tt, s)
      }
    ))
  }
  if (link == "linear") {
    return(list(
      d = m + 1L,
      G = function(theta, t) index(theta, t),
      grad = function(theta, t) cbind(1, tmat(t, nrow(theta)))
    ))
  }
  if (link == "logit") {
    return(list(
      d = m + 1L,
      G = function(theta, t) stats::plogis(index(theta, t)),
      grad = function(theta, t) {
        s <- stats::plogis(index(theta, t))
        s * (1 - s) * cbind(1, tmat(t, nrow(theta)))
      }
    ))
  }
  if (!is.function(G) || !is.function(G_grad)) stop("`link = \"custom\"` needs `G` and `G_grad` functions.", call. = FALSE)
  list(d = d_given,
       G = function(theta, t) as.numeric(G(theta, tmat(t, nrow(theta)))),
       grad = function(theta, t) as.matrix(G_grad(theta, tmat(t, nrow(theta)))))
}

.cm_as_t_matrix <- function(t, n = NULL, m = NULL) {
  if (is.character(t)) {
    if (is.null(m)) m <- nchar(t[1])
    if (any(nchar(t) != m) || any(grepl("[^01]", t))) stop("Treatment strings must be 0/1 strings of length ", m, ".", call. = FALSE)
    t <- t(vapply(strsplit(t, ""), function(s) as.numeric(s), numeric(m)))
    if (m == 1L) t <- matrix(t, ncol = 1)
  }
  t <- as.matrix(t)
  storage.mode(t) <- "double"
  if (!is.null(n) && nrow(t) != n) stop("Treatment matrix must have one row per observation.", call. = FALSE)
  if (!is.null(m) && ncol(t) != m) stop("Treatment combinations must have ", m, " columns.", call. = FALSE)
  if (any(!(t %in% c(0, 1)))) stop("Treatment combinations must be 0/1.", call. = FALSE)
  t
}

.cm_t_dist <- function(t_dist, t_obs, m) {
  if (is.null(t_dist)) {
    key <- apply(t_obs, 1, paste, collapse = "")
    tab <- table(key)
    tm <- t(vapply(strsplit(names(tab), ""), as.numeric, numeric(m)))
    if (m == 1L) tm <- matrix(tm, ncol = 1)
    return(list(t = tm, prob = as.numeric(tab) / sum(tab)))
  }
  if (is.data.frame(t_dist)) {
    if (!"prob" %in% names(t_dist)) stop("`t_dist` needs a `prob` column.", call. = FALSE)
    prob <- as.numeric(t_dist$prob)
    tm <- as.matrix(t_dist[, setdiff(names(t_dist), "prob"), drop = FALSE])
  } else {
    prob <- attr(t_dist, "prob")
    if (is.null(prob)) stop("A matrix `t_dist` needs a `prob` attribute.", call. = FALSE)
    tm <- as.matrix(t_dist)
  }
  storage.mode(tm) <- "double"
  if (ncol(tm) != m) stop("`t_dist` must have ", m, " treatment columns.", call. = FALSE)
  if (abs(sum(prob) - 1) > 1e-8) warning("`t_dist$prob` does not sum to one; the values are used as weights in Lambda(x).", call. = FALSE)
  list(t = tm, prob = prob)
}

# Solve L_i v_i = b_i for every row, where row i of `Lmat` holds the d x d block
# column-major. Uses one sparse block-diagonal solve per chunk.
.cm_blockwise_solve <- function(Lmat, B, d, chunk = 20000L) {
  n <- nrow(Lmat)
  V <- matrix(NA_real_, n, d)
  starts <- seq(1L, n, by = chunk)
  for (s in starts) {
    idx <- s:min(s + chunk - 1L, n)
    nb <- length(idx)
    block_off <- rep((seq_len(nb) - 1L) * d, each = d * d)
    ii <- block_off + rep(rep(seq_len(d), times = d), times = nb)
    jj <- block_off + rep(rep(seq_len(d), each = d), times = nb)
    xx <- as.numeric(t(Lmat[idx, , drop = FALSE]))
    Ls <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(nb * d, nb * d))
    sol <- tryCatch(Matrix::solve(Ls, as.numeric(t(B[idx, , drop = FALSE]))), error = function(e) {
      stop("Lambda(x) is singular or nearly singular for some observations (", conditionMessage(e),
           "). Add a `ridge`, or check that the observed combinations identify every nuisance parameter.", call. = FALSE)
    })
    V[idx, ] <- matrix(as.numeric(sol), nb, d, byrow = TRUE)
  }
  V
}

.cm_block_rcond <- function(Lmat, d, max_rows = 2000L) {
  n <- nrow(Lmat)
  rows <- if (n > max_rows) sort(sample.int(n, max_rows)) else seq_len(n)
  vapply(rows, function(i) {
    M <- matrix(Lmat[i, ], d, d)
    tryCatch(rcond(M), error = function(e) 0)
  }, numeric(1))
}

#' @export
print.cm_dml_structured <- function(x, ...) {
  cat("Structured double machine learning (link: ", x$link, ", ", x$m, " treatments, ",
      length(unique(x$fold_id)), " fold(s), n = ", x$n, ")\n", sep = "")
  cat("  reference: (", paste(x$t0, collapse = ", "), "); assignment distribution over ",
      nrow(x$t_dist$t), " observed combinations", if (x$ridge != 0) paste0("; ridge = ", x$ridge), "\n", sep = "")
  pooled <- x$estimates[x$estimates$fold == "pooled", c("target", "estimate", "plugin", "std.error", "conf.low", "conf.high")]
  print(pooled, digits = 4, row.names = FALSE)
  if (!is.null(x$best)) {
    cat("  best target: ", x$best, "\n", sep = "")
  }
  if (x$diagnostics$share_ill_conditioned > 0) {
    cat("  note: ", round(100 * x$diagnostics$share_ill_conditioned, 1),
        "% of sampled Lambda(x) blocks are ill conditioned; consider `ridge`.\n", sep = "")
  }
  invisible(x)
}
