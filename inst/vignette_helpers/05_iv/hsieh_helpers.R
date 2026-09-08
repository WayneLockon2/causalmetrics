# Helpers for inst/replications/05_IV_control_fun/hsieh2026leveraging.Rmd
# (Hsieh, Du, and Lu 2026, Marketing Science). Everything here is specific to
# the replication: the authors' simulation scripts rerun verbatim, a
# show-level exposure simulation, the geometric stocks, the cluster bootstrap
# for the two-step estimator, and the table builders.

# ---------------------------------------------------------------------------
# 1. Geometric stocks (ad stock, control-function stock)
# ---------------------------------------------------------------------------

# adstock: S_t = x_t + lambda * S_{t-1}, S_1 = x_1 (Equation 11 of the paper).
adstock <- function(x, lambda) {
  if (lambda == 0) return(as.numeric(x))
  as.numeric(stats::filter(x, lambda, method = "recursive"))
}

# adstock by unit; `x` must be sorted by (id, time) within id.
adstock_by <- function(x, id, lambda) {
  if (lambda == 0) return(as.numeric(x))
  as.numeric(stats::ave(as.numeric(x), id, FUN = function(v) adstock(v, lambda)))
}

# ---------------------------------------------------------------------------
# 2. The authors' scripts, verbatim up to file paths and the parallel plan
#    (future.seed = TRUE gives the same streams for any number of workers)
# ---------------------------------------------------------------------------

authors_step1 <- function(N_hh = 25000, T_day = 150, workers = 1, seed = 1107) {
  requireNamespace("data.table"); requireNamespace("future.apply")
  `:=` <- data.table::`:=`; .N <- NULL
  old <- options(future.globals.maxSize = Inf); on.exit(options(old), add = TRUE)
  set.seed(seed)
  if (workers > 1) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  panel <- data.table::data.table(hh = rep(1:N_hh, each = T_day), t = rep(1:T_day, times = N_hh))
  panel[, u_endo := stats::rnorm(.N, 0, 1)]
  panel[, u_exo := stats::rnorm(.N, 0, 1)]
  panel[, p_target1 := stats::plogis(1.0 * u_endo - 0.5)]
  panel[, p_target2 := stats::plogis(0.7 * u_endo - 0.3)]
  panel[, p_target3 := stats::plogis(0.5 * u_endo - 0.2)]
  panel[, A_b1 := stats::rbinom(.N, 1, p_target1)]
  panel[, A_b2 := stats::rbinom(.N, 1, p_target2)]
  panel[, A_b3 := stats::rbinom(.N, 1, p_target3)]
  make_view_portion <- function(u) {
    zero_mass <- stats::rbinom(length(u), 1, 0.30)
    port <- pmin(pmax(u + stats::runif(length(u), 0, 0.5), 0), 1)
    ifelse(zero_mass == 1, 0, port)
  }
  panel[, port_watch1 := make_view_portion(u_endo)]
  panel[, port_watch2 := make_view_portion(u_endo)]
  panel[, port_watch3 := make_view_portion(u_endo)]
  panel[, adexpo1 := stats::rbinom(.N, 1, port_watch1) * A_b1]
  panel[, adexpo2 := stats::rbinom(.N, 1, port_watch2) * A_b2]
  panel[, adexpo3 := stats::rbinom(.N, 1, port_watch3) * A_b3]
  panel[, A := adexpo1 + adexpo2 + adexpo3]
  panel[, A_exp := port_watch1 * p_target1 + port_watch2 * p_target2 + port_watch3 * p_target3]
  panel[, IV := A - A_exp]
  # the paper's expected treatment conditions on the realized targeting decision
  panel[, A_exp_paper := A_b1 * port_watch1 + A_b2 * port_watch2 + A_b3 * port_watch3]
  panel[, IV_paper := A - A_exp_paper]
  stock_exp <- function(x, lambda) {
    out <- numeric(length(x)); out[1] <- x[1]
    for (i in 2:length(x)) out[i] <- x[i] + lambda * out[i - 1]
    out
  }
  panel_split <- split(panel, panel$hh)
  fs_split <- future.apply::future_lapply(panel_split, function(dt) {
    data.table::setorder(dt, t)
    dt[, A_stock := stock_exp(A, 0.70)]
    dt
  }, future.seed = TRUE)
  first_stage_df <- data.table::rbindlist(fs_split)
  first_stage_df[, list(hh, t, A, A_stock, IV, A_exp, IV_paper, A_exp_paper, u_endo, u_exo)]
}

authors_step2 <- function(fs, workers = 1, seed = 123) {
  requireNamespace("data.table"); requireNamespace("future.apply")
  `:=` <- data.table::`:=`
  sim <- data.table::as.data.table(fs)
  data.table::setorder(sim, hh, t)
  old <- options(future.globals.maxSize = Inf); on.exit(options(old), add = TRUE)
  set.seed(seed)
  N_hh <- length(unique(sim$hh)); T_day <- max(sim$t)
  if (workers > 1) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  a0 <- -2.5; beta0 <- 0.012; delta <- 0.005
  a_post <- 1.500; a_f1 <- 0.450; a_f2 <- -0.060; a_r1 <- -0.100; a_r2 <- -0.030
  b_post <- 0.000; b_f1 <- 0.020; b_f2 <- -0.010; b_r1 <- 0.000; b_r2 <- 0.005
  g_promo <- 0.2; g_comp <- 0.02
  set_comp <- function(T) {
    x <- numeric(T); x[1] <- stats::rnorm(1)
    for (tt in 2:T) x[tt] <- 0.8 * x[tt - 1] + stats::rnorm(1, 0, 0.5)
    rng <- range(x); (x - rng[1]) / (rng[2] - rng[1])
  }
  comp_list <- future.apply::future_lapply(1:N_hh, function(i) set_comp(T_day), future.seed = TRUE)
  sim[, comp := unlist(comp_list)]
  sim[, purchase := 0L]; sim[, freq := 0L]; sim[, rec := 0L]; sim[, post := 0L]; sim[, promo1m := 0L]; sim[, unobs := 0L]
  simulate_household <- function(dt_hh) {
    data.table::setorder(dt_hh, t)
    n <- nrow(dt_hh)
    Freq <- 0; Rec <- 0; Post <- 0; unobs <- 0; first_t <- NA
    pur_vec <- integer(n); freq_vec <- integer(n); rec_vec <- integer(n); post_vec <- integer(n); promo_vec <- integer(n)
    for (k in 1:n) {
      lf <- log1p(Freq); lr <- log1p(Rec)
      unobs <- 0.95 * unobs + dt_hh$u_endo[k]
      beta_eff <- beta0 + (b_post * Post + b_f1 * lf + b_f2 * (lf^2) + b_r1 * lr + b_r2 * (lr^2))
      mu <- a0 + a_post * Post + a_f1 * lf + a_f2 * (lf^2) + a_r1 * lr + a_r2 * (lr^2) +
        beta_eff * dt_hh$A_stock[k] + g_promo * promo_vec[k] + g_comp * dt_hh$comp[k] + delta * unobs + dt_hh$u_exo[k]
      y <- ifelse(mu > 0, 1, 0)
      pur_vec[k] <- y
      if (y == 1) { Freq <- Freq + 1; Rec <- 0; if (is.na(first_t)) first_t <- dt_hh$t[k]; Post <- 1 } else { Rec <- Rec + 1 }
      if (Post == 1 & k < first_t + 30) promo <- 1 else promo <- 0
      if (Post == 0) Rec <- 0
      if (k < n) { freq_vec[k + 1] <- Freq; rec_vec[k + 1] <- Rec; post_vec[k + 1] <- Post; promo_vec[k + 1] <- promo }
    }
    dt_hh[, purchase := pur_vec]; dt_hh[, freq := freq_vec]; dt_hh[, rec := rec_vec]; dt_hh[, post := post_vec]; dt_hh[, promo1m := promo_vec]
    dt_hh
  }
  sim_split <- split(sim, sim$hh)
  simulated_list <- future.apply::future_lapply(sim_split, simulate_household, future.seed = TRUE)
  sim_final <- data.table::rbindlist(simulated_list)
  as.data.frame(sim_final[, list(hh, t, purchase, freq, rec, post, promo1m, A, A_stock, comp, IV, IV_paper, A_exp, A_exp_paper, u_endo, u_exo)])
}

authors_formulas <- list(
  first = A ~ IV + post + lf + I(lf^2) + lr + I(lr^2) + promo1m + comp,
  naive = purchase ~ A_stock + post + lf + I(lf^2) + lr + I(lr^2) + A_stock:post + A_stock:lf + A_stock:I(lf^2) +
    A_stock:lr + A_stock:I(lr^2) + promo1m + comp,
  full = purchase ~ A_stock + CF_stock + post + lf + I(lf^2) + lr + I(lr^2) + A_stock:post + A_stock:lf + A_stock:I(lf^2) +
    A_stock:lr + A_stock:I(lr^2) + promo1m + comp
)

# Step 4 (point estimates) on a data frame with columns hh, t, purchase, freq, rec, post, promo1m, A, A_stock, IV, comp.
authors_step4 <- function(ss, lambda_CF = 0.95) {
  est_df <- ss[order(ss$hh, ss$t), ]
  est_df$lf <- log1p(est_df$freq); est_df$lr <- log1p(est_df$rec)
  est_FS <- stats::lm(authors_formulas$first, data = est_df)
  est_df$CF <- stats::residuals(est_FS)
  est_df$CF_stock <- adstock_by(est_df$CF, est_df$hh, lambda_CF)
  mod_naive <- stats::glm(authors_formulas$naive, data = est_df, family = stats::binomial(link = "probit"))
  mod_full <- stats::glm(authors_formulas$full, data = est_df, family = stats::binomial(link = "probit"))
  list(first = est_FS, naive = mod_naive, full = mod_full, data = est_df)
}

# Step 3, the authors' bootstrap: unique households of each draw (their filter(hh %in% boot_hh)).
authors_step3 <- function(ss, B = 50, seed = 1107, lambda_CF = 0.95, workers = 1) {
  set.seed(seed)
  hh_list <- unique(ss$hh)
  draws <- lapply(seq_len(B), function(b) sample(hh_list, size = length(hh_list), replace = TRUE))
  one <- function(boot_hh) {
    ss_b <- ss[ss$hh %in% boot_hh, ]
    ss_b$lf <- log1p(ss_b$freq); ss_b$lr <- log1p(ss_b$rec)
    est_FS <- stats::lm(authors_formulas$first, data = ss_b)
    ss_b$CF <- stats::residuals(est_FS)
    ss_b <- ss_b[order(ss_b$hh, ss_b$t), ]
    ss_b$CF_stock <- adstock_by(ss_b$CF, ss_b$hh, lambda_CF)
    mod_full <- stats::glm(authors_formulas$full, data = ss_b, family = stats::binomial(link = "probit"))
    stats::coef(mod_full)
  }
  res <- if (workers > 1) parallel::mclapply(draws, one, mc.cores = workers) else lapply(draws, one)
  B_mat <- do.call(rbind, res)
  list(draws = B_mat, se = apply(B_mat, 2, stats::sd))
}

authors_truth <- c(`(Intercept)` = -2.5, A_stock = 0.012, CF_stock = NA, post = 1.5, lf = 0.45, lf2 = -0.06,
                   lr = -0.1, lr2 = -0.03, promo1m = 0.2, comp = 0.02, `A_stock:post` = 0, `A_stock:lf` = 0.02,
                   `A_stock:lf2` = -0.010, `A_stock:lr` = 0, `A_stock:lr2` = 0.005)

# The authors' formulas write the squares as I(lf^2); ours use the columns lf2, lr2.
rename_authors <- function(x) {
  nm <- if (is.matrix(x)) colnames(x) else names(x)
  nm <- gsub("I(lf^2)", "lf2", gsub("I(lr^2)", "lr2", nm, fixed = TRUE), fixed = TRUE)
  if (is.matrix(x)) colnames(x) <- nm else names(x) <- nm
  x
}

# ---------------------------------------------------------------------------
# 3. Show-level exposure simulation (the paper's Equations 1-8)
# ---------------------------------------------------------------------------

# Within-show ad airing-time densities by network: mixtures of betas on [0, 1].
network_densities <- function() {
  mix <- function(w, a, b) list(
    r = function(n) { k <- sample(seq_along(w), n, replace = TRUE, prob = w); stats::rbeta(n, a[k], b[k]) },
    d = function(x) rowSums(sapply(seq_along(w), function(j) w[j] * stats::dbeta(x, a[j], b[j]))),
    p = function(x) rowSums(sapply(seq_along(w), function(j) w[j] * stats::pbeta(x, a[j], b[j])))
  )
  list(
    `Network 1` = mix(c(0.35, 0.35, 0.30), c(6, 12, 20), c(18, 10, 3)),   # trimodal, like MTV in Figure 2
    `Network 2` = mix(c(0.5, 0.5), c(4, 10), c(10, 4)),                    # bimodal (sports)
    `Network 3` = mix(1, 1.2, 1.2)                                         # close to uniform (news)
  )
}

# Purchase-model parameters: the authors' values plus the household traits and
# random effects used in the show-level design (Table 4 columns 5-6, Table 5).
purchase_par <- function(traits = TRUE, random = TRUE, a0 = -2.5, delta = 0.005) {
  list(a0 = a0, beta0 = 0.012, delta = delta,
       a_post = 1.5, a_f1 = 0.45, a_f2 = -0.06, a_r1 = -0.10, a_r2 = -0.03,
       b_post = 0, b_f1 = 0.02, b_f2 = -0.01, b_r1 = 0, b_r2 = 0.005,
       g_promo = 0.2, g_comp = 0.02,
       a_view = if (traits) -0.05 else 0, a_sports = if (traits) 0.03 else 0,
       b_view = if (traits) -0.004 else 0, b_sports = if (traits) 0.004 else 0,
       sigma1 = if (random) 0.2 else 0, sigma2 = if (random) 0.01 else 0, rho = if (random) -0.3 else 0,
       lambda_A = 0.7, lambda_u = 0.95)
}

# violation = NULL, or list(type = "skip", pi = ) : partial viewers with u_endo > 0 skip the
# focal ad with probability pi (strategic viewership); list(type = "order", pi = ): the focal ad
# is put in the first pod (Beta(1, 8)) with probability pi (nonrandom ordering).
simulate_exposure <- function(n_hh = 25000, n_days = 150, seed = 1107, n_nonfocal = 200000L, sample_frac = 0.02,
                              lambda_A = 0.7, par = purchase_par(), violation = NULL, placebo = TRUE) {
  set.seed(seed)
  dens <- network_densities(); K <- length(dens)
  # households: traits (standardized) and random effects
  z_view <- as.numeric(scale(stats::rnorm(n_hh))); z_sports <- as.numeric(scale(stats::rnorm(n_hh)))
  om <- MASS::mvrnorm(n_hh, c(0, 0), matrix(c(par$sigma1^2, par$rho * par$sigma1 * par$sigma2,
                                                par$rho * par$sigma1 * par$sigma2, par$sigma2^2), 2))
  if (par$sigma1 == 0) om[, 1] <- 0
  if (par$sigma2 == 0) om[, 2] <- 0
  households <- data.frame(hh = seq_len(n_hh), z_view = z_view, z_sports = z_sports, om_a = om[, 1], om_b = om[, 2])
  # nonfocal ads: all other brands' airing times, the sample the empirical density comes from
  nonfocal <- lapply(dens, function(d) d$r(n_nonfocal))
  Fhat <- lapply(nonfocal, stats::ecdf)
  a_t <- c(1.0, 0.7, 0.5); b_t <- c(0.5, 0.3, 0.2)   # the authors' targeting rules for the three shows
  N <- n_hh * n_days
  pos <- function(t) (seq_len(n_hh) - 1L) * n_days + t   # hh-major layout
  cols <- c("u_endo", "u_exo", "A", "A_exp", "IV", "IV_unif", "IV_oracle", "n_target", "n_partial", "view_hr",
            "targeted_hr", "A_pl", "IV_pl")
  out <- lapply(cols, function(v) numeric(N)); names(out) <- cols
  # accumulators for Table 2 (household-show rows with A_b = 1 and some viewing) and Proposition 2 checks
  vars_show <- c("A", "A_exp", "IV", "P_hat", "A_b")
  S1 <- setNames(numeric(5), vars_show); S2 <- matrix(0, 5, 5, dimnames = list(vars_show, vars_show)); n_show <- 0
  S1_all <- setNames(numeric(2), c("IV", "A_b")); S2_all <- matrix(0, 2, 2); n_all <- 0
  samp <- vector("list", n_days); focal_times <- vector("list", n_days)
  for (t in seq_len(n_days)) {
    u_endo <- stats::rnorm(n_hh); u_exo <- stats::rnorm(n_hh)
    A <- A_exp <- IV <- IV_unif <- IV_oracle <- n_target <- n_partial <- view <- tview <- A_pl <- IV_pl <- numeric(n_hh)
    day_rows <- vector("list", K); ft <- vector("list", K)
    for (k in seq_len(K)) {
      p_target <- stats::plogis(a_t[k] * u_endo - b_t[k])
      A_b <- stats::rbinom(n_hh, 1, p_target)
      zero <- stats::rbinom(n_hh, 1, 0.30)
      len <- pmin(pmax(u_endo + stats::runif(n_hh, 0, 0.5) + 0.25 * z_view + (k == 2) * 0.25 * z_sports, 0), 1)
      len[zero == 1] <- 0
      start <- stats::runif(n_hh, 0, 1 - len); end <- start + len
      T_f <- dens[[k]]$r(n_hh)
      if (!is.null(violation) && violation$type == "order") {
        sw <- stats::rbinom(n_hh, 1, violation$pi) == 1
        T_f[sw] <- stats::rbeta(sum(sw), 1, 8)
      }
      exposed <- A_b * as.integer(T_f >= start & T_f <= end)
      if (!is.null(violation) && violation$type == "skip") {
        can <- exposed == 1 & len < 1 & u_endo > 0
        skip <- can & stats::rbinom(n_hh, 1, violation$pi) == 1
        exposed[skip] <- 0L
      }
      P_true <- dens[[k]]$p(end) - dens[[k]]$p(start)
      P_hat <- Fhat[[k]](end) - Fhat[[k]](start)
      iv <- exposed - A_b * P_hat
      A <- A + exposed; A_exp <- A_exp + A_b * P_hat; IV <- IV + iv
      IV_unif <- IV_unif + exposed - A_b * len
      IV_oracle <- IV_oracle + exposed - p_target * P_true
      n_target <- n_target + A_b; partial <- A_b == 1 & len > 0 & len < 1; n_partial <- n_partial + partial
      view <- view + len; tview <- tview + A_b * len
      if (placebo) {
        A_b2 <- stats::rbinom(n_hh, 1, stats::plogis(0.5 * u_endo - 0.4))
        T2 <- dens[[k]]$r(n_hh)
        ex2 <- A_b2 * as.integer(T2 >= start & T2 <= end)
        A_pl <- A_pl + ex2; IV_pl <- IV_pl + ex2 - A_b2 * P_hat
      }
      keep <- A_b == 1 & len > 0
      M <- cbind(A = exposed[keep], A_exp = (A_b * P_hat)[keep], IV = iv[keep], P_hat = P_hat[keep], A_b = A_b[keep])
      S1 <- S1 + colSums(M); S2 <- S2 + crossprod(M); n_show <- n_show + nrow(M)
      M2 <- cbind(iv, A_b); S1_all <- S1_all + colSums(M2); S2_all <- S2_all + crossprod(M2); n_all <- n_all + n_hh
      pick <- which(stats::runif(n_hh) < sample_frac)
      day_rows[[k]] <- data.frame(hh = pick, t = rep(t, length(pick)), network = rep(k, length(pick)), A_b = A_b[pick], len = len[pick], start = start[pick],
                                  end = end[pick], T_focal = T_f[pick], A = exposed[pick], P_hat = P_hat[pick],
                                  P_true = P_true[pick], IV = iv[pick], u_endo = u_endo[pick])
      tf <- T_f[A_b == 1][seq_len(min(sum(A_b == 1), 150L))]; ft[[k]] <- data.frame(network = rep(k, length(tf)), time = tf)
    }
    idx <- pos(t)
    out$u_endo[idx] <- u_endo; out$u_exo[idx] <- u_exo; out$A[idx] <- A; out$A_exp[idx] <- A_exp; out$IV[idx] <- IV
    out$IV_unif[idx] <- IV_unif; out$IV_oracle[idx] <- IV_oracle; out$n_target[idx] <- n_target
    out$n_partial[idx] <- n_partial; out$view_hr[idx] <- view; out$targeted_hr[idx] <- tview
    out$A_pl[idx] <- A_pl; out$IV_pl[idx] <- IV_pl
    samp[[t]] <- do.call(rbind, day_rows); focal_times[[t]] <- do.call(rbind, ft)
  }
  panel <- data.frame(hh = rep(seq_len(n_hh), each = n_days), t = rep(seq_len(n_days), times = n_hh), out)
  # competitor ad spend: household AR(1), normalized to [0, 1] as in the authors' Step 2
  comp <- matrix(0, n_hh, n_days); comp[, 1] <- stats::rnorm(n_hh)
  for (t in 2:n_days) comp[, t] <- 0.8 * comp[, t - 1] + stats::rnorm(n_hh, 0, 0.5)
  comp <- (comp - apply(comp, 1, min)) / (apply(comp, 1, max) - apply(comp, 1, min))
  panel$comp <- as.numeric(t(comp))
  panel$A_stock <- adstock_by(panel$A, panel$hh, lambda_A)
  panel$A_pl_stock <- adstock_by(panel$A_pl, panel$hh, lambda_A)
  panel$IV_stock <- adstock_by(panel$IV, panel$hh, lambda_A)
  panel <- merge(panel, households, by = "hh", sort = FALSE)
  panel <- panel[order(panel$hh, panel$t), ]
  rownames(panel) <- NULL
  ad_times <- rbind(
    do.call(rbind, lapply(seq_len(K), function(k) data.frame(network = k, type = "Nonfocal brands", time = nonfocal[[k]][seq_len(20000L)]))),
    transform(do.call(rbind, focal_times), type = "Focal brand")
  )
  list(panel = panel, households = households, show_sample = do.call(rbind, samp), ad_times = ad_times,
       densities = dens, Fhat = Fhat, moments = list(S1 = S1, S2 = S2, n = n_show, S1_all = S1_all, S2_all = S2_all, n_all = n_all),
       par = par, n_hh = n_hh, n_days = n_days, lambda_A = lambda_A)
}

# Purchases: the authors' Step 2 loop, vectorized across households (loop over days).
# `panel` sorted by (hh, t) with A_stock, u_endo, u_exo, comp, and (if used) z_view, z_sports, om_a, om_b.
simulate_purchases <- function(panel, par = purchase_par(FALSE, FALSE)) {
  stopifnot(!is.unsorted(order(panel$hh, panel$t)))
  n_hh <- length(unique(panel$hh)); T_day <- max(panel$t); N <- nrow(panel)
  pos <- function(t) (seq_len(n_hh) - 1L) * T_day + t
  zv <- if ("z_view" %in% names(panel)) panel$z_view[pos(1)] else 0
  zs <- if ("z_sports" %in% names(panel)) panel$z_sports[pos(1)] else 0
  oa <- if ("om_a" %in% names(panel)) panel$om_a[pos(1)] else 0
  ob <- if ("om_b" %in% names(panel)) panel$om_b[pos(1)] else 0
  Freq <- Rec <- Post <- unobs <- promo <- numeric(n_hh); first_t <- rep(NA_real_, n_hh)
  purchase <- freq <- rec <- post <- promo1m <- integer(N); index <- beta_it <- numeric(N)
  for (k in seq_len(T_day)) {
    idx <- pos(k)
    lf <- log1p(Freq); lr <- log1p(Rec)
    unobs <- par$lambda_u * unobs + panel$u_endo[idx]
    beta_eff <- par$beta0 + par$b_post * Post + par$b_f1 * lf + par$b_f2 * lf^2 + par$b_r1 * lr + par$b_r2 * lr^2 +
      par$b_view * zv + par$b_sports * zs + ob
    mu <- par$a0 + par$a_post * Post + par$a_f1 * lf + par$a_f2 * lf^2 + par$a_r1 * lr + par$a_r2 * lr^2 +
      beta_eff * panel$A_stock[idx] + par$g_promo * promo + par$g_comp * panel$comp[idx] + par$delta * unobs +
      par$a_view * zv + par$a_sports * zs + oa
    y <- as.integer(mu + panel$u_exo[idx] > 0)
    purchase[idx] <- y; freq[idx] <- Freq; rec[idx] <- Rec; post[idx] <- Post; promo1m[idx] <- promo
    index[idx] <- mu; beta_it[idx] <- beta_eff
    Freq <- Freq + y
    Rec <- ifelse(y == 1, 0, Rec + 1)
    first_t <- ifelse(y == 1 & is.na(first_t), k, first_t)
    Post <- pmax(Post, y)
    promo <- as.integer(Post == 1 & !is.na(first_t) & k < first_t + 30)
    Rec[Post == 0] <- 0
  }
  panel$purchase <- purchase; panel$freq <- freq; panel$rec <- rec; panel$post <- post; panel$promo1m <- promo1m
  panel$index <- index; panel$beta_it <- beta_it
  panel
}

# ---------------------------------------------------------------------------
# 4. Estimation: first stage, stocks, the six columns of Table 4, bootstrap
# ---------------------------------------------------------------------------

add_states <- function(df) {
  df$lf <- log1p(df$freq); df$lr <- log1p(df$rec); df$lf2 <- df$lf^2; df$lr2 <- df$lr^2
  df
}

controls_base <- c("post", "lf", "lf2", "lr", "lr2", "promo1m", "comp")
controls_traits <- c("z_view", "z_sports")

# A column specification of Table 4.
col_spec <- function(lambda_A = 0.7, cf = TRUE, lambda_CF = 0.95, interactions = TRUE, traits = TRUE,
                     trait_slopes = TRUE, random = FALSE, iv = "IV", d = "A", stock = "A_stock") {
  list(lambda_A = lambda_A, cf = cf, lambda_CF = lambda_CF, interactions = interactions, traits = traits,
       trait_slopes = trait_slopes, random = random, iv = iv, d = d, stock = stock)
}

second_formula <- function(spec) {
  rhs <- c("A_stock", if (spec$cf) "CF_stock", controls_base, if (spec$traits) controls_traits)
  if (spec$interactions) rhs <- c(rhs, paste0("A_stock:", c("post", "lf", "lf2", "lr", "lr2")))
  if (spec$traits && spec$trait_slopes) rhs <- c(rhs, paste0("A_stock:", controls_traits))
  if (spec$random) rhs <- c(rhs, "(1 + A_stock | hh)")
  stats::as.formula(paste("purchase ~", paste(rhs, collapse = " + ")))
}

# First stage with cf_residuals(); the control is accumulated within `group` (hh, or .draw in a bootstrap).
first_stage <- function(df, spec, group = "hh") {
  x <- c(controls_base, if (spec$traits) controls_traits)
  df <- cf_residuals(df, d = spec$d, z = spec$iv, x = x, name = "CF")
  df <- df[order(df[[group]], df$t), ]
  df$A_stock <- adstock_by(df[[spec$d]], df[[group]], spec$lambda_A)
  df$CF_stock <- adstock_by(df$CF, df[[group]], spec$lambda_CF)
  df
}

second_stage <- function(df, spec, nAGQ = 0) {
  f <- second_formula(spec)
  if (spec$random) {
    lme4::glmer(f, data = df, family = stats::binomial(link = "probit"), nAGQ = nAGQ)
  } else {
    stats::glm(f, data = df, family = stats::binomial(link = "probit"))
  }
}

fit_column <- function(df, spec, group = "hh", nAGQ = 0) {
  aug <- first_stage(df, spec, group)
  fit <- second_stage(aug, spec, nAGQ)
  list(fit = fit, data = aug, spec = spec)
}

coef_vector <- function(fit) {
  if (inherits(fit, "merMod")) {
    cf <- lme4::fixef(fit)
    vc <- lme4::VarCorr(fit)$hh
    sd <- attr(vc, "stddev"); cr <- attr(vc, "correlation")
    c(cf, sigma1 = unname(sd[1]), sigma2 = unname(sd[2]), rho = unname(cr[1, 2]))
  } else stats::coef(fit)
}

# Cluster bootstrap of the two-step estimator. Each draw relabels the sampled clusters
# (.draw), so a household drawn twice is accumulated twice, separately. unique_only = TRUE
# reproduces the authors' Step 3, which keeps each sampled household once.
cluster_bootstrap <- function(data, spec, n_boot = 50, seed = 1107, cluster = "hh", unique_only = FALSE,
                              workers = 1, nAGQ = 0, extra = NULL) {
  ids <- unique(data[[cluster]])
  rows_by <- split(seq_len(nrow(data)), data[[cluster]])
  set.seed(seed)
  draws <- lapply(seq_len(n_boot), function(b) sample(ids, size = length(ids), replace = TRUE))
  one <- function(pick) {
    if (unique_only) pick <- unique(pick)
    rb <- rows_by[as.character(pick)]
    df <- data[unlist(rb), , drop = FALSE]
    df$.draw <- rep(seq_along(pick), lengths(rb))
    if (spec$random) df$hh <- df$.draw
    fit <- fit_column(df, spec, group = ".draw", nAGQ = nAGQ)
    out <- coef_vector(fit$fit)
    if (!is.null(extra)) out <- c(out, extra(fit$fit, fit$data))
    out
  }
  res <- if (workers > 1) parallel::mclapply(draws, function(p) tryCatch(one(p), error = function(e) NULL), mc.cores = workers)
         else lapply(draws, function(p) tryCatch(one(p), error = function(e) NULL))
  ok <- !vapply(res, is.null, logical(1))
  B <- do.call(rbind, res[ok])
  list(draws = B, se = apply(B, 2, stats::sd), n_failed = sum(!ok))
}

# ---------------------------------------------------------------------------
# 5. Effects: average partial effects, elasticities, responsiveness curves
# ---------------------------------------------------------------------------

# Per-row derivative of the purchase probability with respect to the ad stock
# (cf_ape() contributions); same-day elasticity with respect to today's exposures,
# and the 30-day elasticity that accumulates the carryover.
elasticities <- function(fit, data, lambda_A = 0.7, horizon = 30) {
  m <- attr(cf_ape(fit, data, d = "A_stock"), "contributions")
  p <- as.numeric(stats::predict(fit, newdata = data, type = "response"))
  e0 <- mean(m * data$A) / mean(p)
  c(ape = mean(m), same_day = e0, long_run = e0 * sum(lambda_A^(0:(horizon - 1))))
}

# beta_it implied by the coefficients (Equation 13 and 15), for converted households.
beta_curve <- function(cf, freq, rec, post = 1) {
  g <- function(nm) if (nm %in% names(cf)) cf[[nm]] else 0
  lf <- log1p(freq); lr <- log1p(rec)
  g("A_stock") + g("A_stock:post") * post + g("A_stock:lf") * lf + g("A_stock:lf2") * lf^2 + g("A_stock:lr") * lr + g("A_stock:lr2") * lr^2
}

alpha_curve <- function(cf, freq, rec, post = 1) {
  g <- function(nm) if (nm %in% names(cf)) cf[[nm]] else 0
  lf <- log1p(freq); lr <- log1p(rec)
  g("(Intercept)") + g("post") * post + g("lf") * lf + g("lf2") * lf^2 + g("lr") * lr + g("lr2") * lr^2
}

# ---------------------------------------------------------------------------
# 6. Diagnostics: Table 2 moments, falsification checks, KS by network, decay grid
# ---------------------------------------------------------------------------

show_moments <- function(mom) {
  n <- mom$n; m <- mom$S1 / n
  V <- mom$S2 / n - tcrossprod(m)
  sd <- sqrt(diag(V)); R <- V / tcrossprod(sd)
  n2 <- mom$n_all; m2 <- mom$S1_all / n2; V2 <- mom$S2_all / n2 - tcrossprod(m2)
  list(mean = m, sd = sd, corr = R, n = n, corr_iv_target_all = V2[1, 2] / sqrt(V2[1, 1] * V2[2, 2]))
}

quantile_row <- function(x, label) {
  q <- stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
  data.frame(Variable = label, Mean = mean(x), SD = stats::sd(x), q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
}

ks_by_network <- function(ad_times) {
  do.call(rbind, lapply(sort(unique(ad_times$network)), function(k) {
    a <- ad_times$time[ad_times$network == k & ad_times$type == "Focal brand"]
    b <- ad_times$time[ad_times$network == k & ad_times$type == "Nonfocal brands"]
    ks <- suppressWarnings(stats::ks.test(a, b))
    data.frame(network = paste("Network", k), n_focal = length(a), n_nonfocal = length(b), D = unname(ks$statistic), p = ks$p.value)
  }))
}

# Out-of-sample log-likelihood over a grid of (lambda_A, lambda_CF): fit on days <= train_days,
# score the remaining days (footnote 17 of the paper).
decay_grid <- function(df, grid = c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95), train_days = 120, spec = col_spec(interactions = FALSE)) {
  x <- c(controls_base, if (spec$traits) controls_traits)
  df <- cf_residuals(df, d = spec$d, z = spec$iv, x = x, name = "CF")
  df <- df[order(df$hh, df$t), ]
  train <- df$t <= train_days
  res <- expand.grid(lambda_A = grid, lambda_CF = grid)
  res$loglik <- NA_real_
  for (i in seq_len(nrow(res))) {
    s <- spec; s$lambda_A <- res$lambda_A[i]; s$lambda_CF <- res$lambda_CF[i]
    df$A_stock <- adstock_by(df$A, df$hh, s$lambda_A); df$CF_stock <- adstock_by(df$CF, df$hh, s$lambda_CF)
    fit <- stats::glm(second_formula(s), data = df[train, ], family = stats::binomial(link = "probit"))
    p <- stats::predict(fit, newdata = df[!train, ], type = "response")
    y <- df$purchase[!train]
    res$loglik[i] <- sum(y * log(p) + (1 - y) * log(1 - p))
  }
  res
}

# ---------------------------------------------------------------------------
# 7. Monte Carlo under violations of the identifying assumptions
# ---------------------------------------------------------------------------

violation_mc <- function(n_sims = 20, n_hh = 2000, n_days = 90, type = c("none", "skip", "order"), pi = 0, seed = 1,
                         workers = 1, spec = col_spec(traits = FALSE, interactions = TRUE),
                         par = purchase_par(FALSE, FALSE, a0 = -3.0, delta = 0.015)) {
  type <- match.arg(type)
  one <- function(s) {
    sim <- simulate_exposure(n_hh, n_days, seed = seed * 1000 + s, n_nonfocal = 50000L, sample_frac = 0, par = par,
                             violation = if (type == "none") NULL else list(type = type, pi = pi), placebo = FALSE)
    pan <- add_states(simulate_purchases(sim$panel, par))
    fit_cf <- fit_column(pan, spec)
    naive <- stats::glm(second_formula(modifyList(spec, list(cf = FALSE))), data = fit_cf$data, family = stats::binomial(link = "probit"))
    mom <- show_moments(sim$moments)
    ks <- ks_by_network(sim$ad_times)
    tt <- stats::t.test(pan$IV)
    c(beta_cf = unname(stats::coef(fit_cf$fit)["A_stock"]), beta_naive = unname(stats::coef(naive)["A_stock"]),
      ape_cf = unname(cf_ape(fit_cf$fit, fit_cf$data, "A_stock")), ape_naive = unname(cf_ape(naive, fit_cf$data, "A_stock")),
      ape_true = mean(stats::dnorm(pan$index) * pan$beta_it),
      delta = unname(stats::coef(fit_cf$fit)["CF_stock"]), iv_mean = mean(pan$IV), iv_t = unname(tt$statistic),
      corr_iv_exp = mom$corr["IV", "A_exp"], corr_iv_phat = mom$corr["IV", "P_hat"], ks_p_min = min(ks$p),
      F = attr(fit_cf$data, "first_stage")$strength$F)
  }
  safe <- function(s) tryCatch(one(s), error = function(e) NULL)
  res <- if (workers > 1) parallel::mclapply(seq_len(n_sims), safe, mc.cores = workers) else lapply(seq_len(n_sims), safe)
  res <- res[!vapply(res, is.null, logical(1))]
  out <- as.data.frame(do.call(rbind, res))
  out$type <- type; out$pi <- pi
  out
}
