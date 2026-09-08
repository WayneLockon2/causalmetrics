# Helpers for the replication of Ye, Zhang, Zhang, Zhang, and Zhang (2025),
# "Deep learning-based causal inference for large-scale combinatorial
# experiments". Sourced by inst/replications/04_Doubly_Robust_and_Double_ML/
# ye2025deep.Rmd. The first stage (PyTorch) runs through
# inst/python/dedl_nuisance.py; everything else is R.

library(dplyr)
library(tidyr)
library(ggplot2)

# Python detection ---------------------------------------------------------------

# A Python interpreter with torch, numpy, and pandas. Candidates: the
# PYTHON_DEDL environment variable, miniconda/anaconda, then the PATH.
.python_has_torch <- function(bin) {
  if (!nzchar(bin) || !file.exists(bin)) return(FALSE)
  status <- suppressWarnings(system(paste(shQuote(bin), "-c", shQuote("import torch, numpy, pandas")),
                                    ignore.stdout = TRUE, ignore.stderr = TRUE))
  identical(status, 0L)
}

find_python_torch <- function() {
  cands <- unique(c(
    Sys.getenv("PYTHON_DEDL"),
    path.expand("~/miniconda3/bin/python"),
    path.expand("~/anaconda3/bin/python"),
    "/opt/anaconda3/bin/python",
    unname(Sys.which("python3")),
    unname(Sys.which("python"))
  ))
  cands <- cands[nzchar(cands)]
  for (b in cands) if (.python_has_torch(b)) return(b)
  ""
}

# Run the first-stage script and read its output. `args` is a character
# vector of command-line flags. Results are cached on disk by `output`.
run_dedl_python <- function(python, script, input, output, args, history = NULL, epoch_theta = NULL,
                            rerun = FALSE, quiet = TRUE) {
  if (!rerun && file.exists(output)) return(invisible(read.csv(output)))
  cmd_args <- c(shQuote(script), "--input", shQuote(input), "--output", shQuote(output), args)
  if (!is.null(history)) cmd_args <- c(cmd_args, "--history", shQuote(history))
  if (!is.null(epoch_theta)) cmd_args <- c(cmd_args, "--epoch-theta", shQuote(epoch_theta))
  status <- system2(python, cmd_args, stdout = if (quiet) FALSE else "", stderr = if (quiet) FALSE else "")
  if (!identical(status, 0L)) stop("dedl_nuisance.py failed (exit status ", status, ").", call. = FALSE)
  invisible(read.csv(output))
}

# Treatment combinations ---------------------------------------------------------

combo_key <- function(t) apply(as.matrix(t), 1, paste, collapse = "")
combo_label <- function(key) paste0("(", gsub("", ", ", key) |> sub("^, ", "", x = _) |> sub(", $", "", x = _), ")")
combo_matrix <- function(keys) {
  m <- t(vapply(strsplit(keys, ""), as.numeric, numeric(nchar(keys[1]))))
  if (nchar(keys[1]) == 1) m <- matrix(m, ncol = 1)
  m
}
paper_order <- c("000", "001", "010", "100", "111", "110", "101", "011")
observed_cells <- c("000", "001", "010", "100", "111")

# Ground truth: cell means against the baseline, with the notebook's test
# (pooled t-test when Levene's test does not reject, Welch otherwise).
ground_truth_table <- function(y, key, combos = paper_order, base = "000") {
  y0 <- y[key == base]
  do.call(rbind, lapply(combos, function(cc) {
    y1 <- y[key == cc]
    if (cc == base) {
      return(data.frame(combo = cc, n = length(y1), diff = 0, relative = 0, sd = sd(y1), se = 0, p.value = 1))
    }
    lev <- car::leveneTest(c(y0, y1), factor(rep(c(0, 1), c(length(y0), length(y1)))), center = mean)$`Pr(>F)`[1]
    tt <- stats::t.test(y1, y0, var.equal = lev >= 0.05)
    data.frame(combo = cc, n = length(y1), diff = mean(y1) - mean(y0),
               relative = (mean(y1) - mean(y0)) / mean(y0), sd = sd(y1),
               se = sqrt(var(y1) / length(y1) + var(y0) / length(y0)), p.value = tt$p.value)
  }))
}

stars <- function(p) ifelse(p < 1e-4, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# Benchmarks ---------------------------------------------------------------------

# Linear addition: single-treatment ATEs (cell means) added up; variances added
# under independence; the notebook's convention (baseline variance counted once
# per contrast).
la_estimator <- function(y, key, combos = paper_order, base = "000", conf_level = 0.95) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  cell_stats <- function(cc) {
    v <- y[key == cc]
    c(mean = mean(v), var = stats::var(v) / length(v))
  }
  b <- cell_stats(base)
  singles <- c("100", "010", "001")
  eff <- sapply(singles, function(s) cell_stats(s)["mean"] - b["mean"])
  vr <- sapply(singles, function(s) cell_stats(s)["var"] + b["var"])
  do.call(rbind, lapply(combos, function(cc) {
    bits <- as.numeric(strsplit(cc, "")[[1]])
    if (cc %in% c(base, singles, "111") && cc != base) {
      s <- cell_stats(cc)
      est <- s["mean"] - b["mean"]; v <- s["var"] + b["var"]
    } else if (cc == base) {
      est <- 0; v <- 0
    } else {
      est <- sum(eff[bits == 1]); v <- sum(vr[bits == 1])
    }
    data.frame(combo = cc, estimate = unname(est), std.error = sqrt(unname(v)),
               conf.low = unname(est) - crit * sqrt(unname(v)), conf.high = unname(est) + crit * sqrt(unname(v)))
  }))
}

# Linear regression benchmark: y on x and t (training rows), contrasts of
# predictions on the prediction rows against the baseline prediction, with the
# notebook's inference (sd of the row-wise contrast / sqrt(n); Welch p-value
# between the two prediction vectors).
lr_estimator <- function(train, predict_rows, y, tvar, xvar, combos = paper_order, conf_level = 0.95) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  f <- stats::reformulate(c(xvar, tvar), y)
  fit <- stats::lm(f, data = train)
  pred_at <- function(cc) {
    nd <- predict_rows
    bits <- as.numeric(strsplit(cc, "")[[1]])
    for (k in seq_along(tvar)) nd[[tvar[k]]] <- bits[k]
    as.numeric(stats::predict(fit, nd))
  }
  base <- strrep("0", nchar(combos[1]))
  p0 <- pred_at(base)
  preds <- sapply(combos, pred_at)
  tab <- do.call(rbind, lapply(seq_along(combos), function(i) {
    est <- preds[, i] - p0
    data.frame(combo = combos[i], estimate = mean(est), std.error = sd(est) / sqrt(length(est)),
               conf.low = mean(est) - crit * sd(est) / sqrt(length(est)),
               conf.high = mean(est) + crit * sd(est) / sqrt(length(est)),
               p.value = if (combos[i] == base) 1 else stats::t.test(preds[, i], p0)$p.value)
  }))
  list(table = tab, predictions = preds)
}

# Predicted-outcome contrasts from the pure network (PDL) output of the script.
pdl_estimator <- function(pdl, combos = paper_order, conf_level = 0.95) {
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  base <- strrep("0", nchar(combos[1]))
  p0 <- pdl[[paste0("pred_", base)]]
  preds <- sapply(combos, function(cc) pdl[[paste0("pred_", cc)]])
  tab <- do.call(rbind, lapply(seq_along(combos), function(i) {
    est <- preds[, i] - p0
    data.frame(combo = combos[i], estimate = mean(est), std.error = sd(est) / sqrt(length(est)),
               conf.low = mean(est) - crit * sd(est) / sqrt(length(est)),
               conf.high = mean(est) + crit * sd(est) / sqrt(length(est)),
               p.value = if (combos[i] == base) 1 else stats::t.test(preds[, i], p0)$p.value)
  }))
  list(table = tab, predictions = preds)
}

# Metrics of Section 4.4: correct direction ratio, MAPE, MSE, MAE. Following
# the authors' code, a combination whose ground truth is insignificant
# contributes 0 to MAPE, MSE, and MAE and counts as correct when the estimate's
# interval covers zero.
combo_metrics <- function(est, truth, alpha = 0.05) {
  stopifnot(all(est$combo == truth$combo))
  sig <- truth$p.value < alpha
  covers0 <- est$conf.low * est$conf.high <= 0
  correct <- ifelse(sig, !covers0 & sign(est$estimate) == sign(truth$diff), covers0)
  ape <- ifelse(sig, 100 * abs(est$estimate - truth$diff) / abs(truth$diff), 0)
  se <- ifelse(sig, (est$estimate - truth$diff)^2, 0)
  ae <- ifelse(sig, abs(est$estimate - truth$diff), 0)
  data.frame(cdr = paste0(sum(correct), "/", length(correct)), cdr_share = mean(correct),
             mape = mean(ape), mse = mean(se), mae = mean(ae))
}

metrics_table <- function(estimates, truth, unobserved = c("110", "101", "011")) {
  do.call(rbind, lapply(names(estimates), function(nm) {
    e <- estimates[[nm]]
    e <- e[match(truth$combo, e$combo), ]
    keep_u <- truth$combo %in% unobserved
    keep_a <- truth$combo != "000"
    mu <- combo_metrics(e[keep_u, ], truth[keep_u, ])
    ma <- combo_metrics(e[keep_a, ], truth[keep_a, ])
    data.frame(estimator = nm, cdr_u = mu$cdr, mape_u = mu$mape, mse_u = mu$mse, mae_u = mu$mae,
               cdr = ma$cdr, mape = ma$mape, mse = ma$mse, mae = ma$mae)
  }))
}

# Figures ------------------------------------------------------------------------

plot_figure1 <- function(df) {
  true_eff <- df$app_DS - df$app___
  la_eff <- df$app__S + df$app_D_ - 2 * df$app___
  slope <- sum(true_eff * la_eff) / sum(true_eff^2)
  pct <- la_eff / true_eff - 1
  qs <- stats::quantile(pct, c(0.1, 0.9))
  d1 <- data.frame(true = true_eff, la = la_eff)
  p1 <- ggplot(d1, aes(true, la)) +
    geom_point(size = 0.8, colour = "#1f77b4") +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    geom_abline(slope = slope, intercept = 0, colour = "blue", linetype = "dotdash") +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
    coord_cartesian(xlim = c(-200, 200), ylim = c(-200, 200)) +
    labs(x = "True effect size", y = "Linear additive effect size",
         subtitle = paste0("Solid: 45-degree line; dot-dash: linear fit (slope ", round(slope, 2), ")")) +
    theme_minimal(base_size = 11)
  d2 <- data.frame(pct = pct[pct > qs[1] & pct < qs[2]])
  p2 <- ggplot(d2, aes(pct)) +
    geom_histogram(aes(y = after_stat(density)), bins = 250, fill = "#1f77b4") +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    labs(x = "Percentage error between linear additive and true effect size", y = "Frequency") +
    theme_minimal(base_size = 11)
  patchwork::wrap_plots(p1, p2, widths = c(8, 20))
}

plot_figure5 <- function(truth, dedl, sdl, combos = paper_order[-1]) {
  df <- bind_rows(
    data.frame(combo = truth$combo, series = "True effect", estimate = truth$diff,
               low = truth$diff - 1.96 * truth$se, high = truth$diff + 1.96 * truth$se),
    data.frame(combo = dedl$combo, series = "DeDL", estimate = dedl$estimate, low = dedl$conf.low, high = dedl$conf.high),
    data.frame(combo = sdl$combo, series = "SDL", estimate = sdl$estimate, low = sdl$conf.low, high = sdl$conf.high)
  ) %>% filter(combo %in% combos) %>%
    mutate(combo = factor(combo_label(combo), levels = combo_label(sort(combos))),
           series = factor(series, levels = c("True effect", "DeDL", "SDL")))
  ggplot(df, aes(combo, estimate, colour = series)) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.25, position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_colour_manual(values = c("True effect" = "red", "DeDL" = "blue", "SDL" = "darkgreen"), name = NULL) +
    labs(x = "Treatment combination", y = "Effect size") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}

plot_figure7 <- function(curve) {
  # curve: epoch, val_mse, mape_dedl, mape_sdl, mape_lr
  scale <- max(curve$val_mse) / max(c(curve$mape_dedl, curve$mape_sdl, curve$mape_lr))
  df <- curve %>%
    select(epoch, DeDL = mape_dedl, SDL = mape_sdl, `Linear regression` = mape_lr) %>%
    pivot_longer(-epoch, names_to = "estimator", values_to = "mape")
  ggplot() +
    geom_line(data = curve, aes(epoch, val_mse, colour = "Validation MSE (left axis)"), linewidth = 0.9) +
    geom_line(data = df, aes(epoch, mape * scale, colour = estimator), linewidth = 0.9) +
    scale_y_continuous(name = "Validation MSE loss", sec.axis = sec_axis(~ . / scale, name = "MAPE (percent)")) +
    scale_colour_manual(values = c("Validation MSE (left axis)" = "red", DeDL = "blue", SDL = "darkgreen",
                                   "Linear regression" = "orange"), name = NULL) +
    labs(x = "Training epoch") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}

# Appendix D data-generating process ----------------------------------------------

# The validation design of Appendix D.1 (notebook cell "Validation of DeDL"):
# m experiments, d_c uniform covariates, generalized sigmoid outcome with a
# cubic index, observed combinations = baseline, the m singles, and (1, 1, 0, ...).
appendix_d_design <- function(m = 4, d_c = 10, n_train = 2000, n_est = 2000, seed = 1,
                              index_power = 3, gamma_linear = 0, noise = 0.05) {
  set.seed(seed)
  coef <- matrix(runif(d_c * (m + 1), -0.5, 0.5), d_c, m + 1)
  coef_lin <- matrix(runif(d_c * (m + 1), -0.5, 0.5), d_c, m + 1)
  c_true <- runif(1, 10, 20)
  all_t <- as.matrix(expand.grid(rep(list(0:1), m)))[, m:1, drop = FALSE]
  colnames(all_t) <- paste0("t", seq_len(m))
  obs_t <- rbind(rep(0, m), diag(m), c(1, 1, rep(0, m - 2)))
  colnames(obs_t) <- colnames(all_t)
  mean_fun <- function(x, t) {
    beta <- (x %*% coef)^index_power
    lin <- if (gamma_linear != 0) gamma_linear * rowSums((x %*% coef_lin) * cbind(1, t)) else 0
    c_true / (1 + exp(-rowSums(beta * cbind(1, t)))) + lin
  }
  gen <- function(n, tset) {
    x <- matrix(runif(n * d_c), n, d_c)
    t <- tset[sample.int(nrow(tset), n, replace = TRUE), , drop = FALSE]
    y <- mean_fun(x, t) + noise * runif(n, -1, 1)
    df <- data.frame(x, t, y = y)
    names(df)[seq_len(d_c)] <- paste0("x", seq_len(d_c))
    df
  }
  train <- gen(n_train, obs_t)
  est <- gen(n_est, obs_t)
  truth <- sapply(seq_len(nrow(all_t)), function(j) {
    x <- as.matrix(est[, paste0("x", seq_len(d_c))])
    mean(mean_fun(x, matrix(all_t[j, ], n_est, m, byrow = TRUE)) - mean_fun(x, matrix(0, n_est, m)))
  })
  list(train = train, est = est, all_t = all_t, obs_t = obs_t, truth = truth,
       truth_keys = combo_key(all_t), c_true = c_true, m = m, d_c = d_c)
}
