---
name: synthetic-control
description: Synthetic control with causalmetrics. Use when one or a few aggregate units (states, countries, firms) receive treatment, many untreated donors exist, and the pre-treatment window is long. Covers classic and demeaned synthetic control, ridge augmentation, synthetic DiD, placebo and conformal inference, the Ferman-Pinto specification test, and leave-one-out donor checks.
---

# Synthetic control with causalmetrics

## 1. Does this design apply?

The question: the effect of an intervention on one (or few) treated
aggregate units, estimated against a weighted combination of donors that
reproduces the treated unit's pre-treatment path. Identification leans on
a factor structure: shocks specific to the treated unit must be unrelated
to treatment; confounders that also move donors are absorbed by the
weights.

Checks to run first:

```r
length(unique(dat$id))                          # donors + treated
range(dat$time); t0 <- min(dat$time[dat$d == 1])  # adoption; T0 = pre-periods
table(tapply(dat$d, dat$id, max))               # how many treated units
```

- One row per unit-period, a 0/1 treatment column that switches on once.
- Long pre-period relative to the donor count (Ferman-Pinto request 3:
  donors few relative to T0, or good fit may be overfitting).
- At least two pre-periods per treated cohort (the estimator requires it).

Red flags pointing elsewhere: many treated units with staggered timing
and a usable comparison group -> `did-analysis` (or synthetic DiD here,
cohort by cohort); micro panels with thousands of treated units ->
`did-analysis`; donors themselves affected by the treatment (spillovers)
-> shrink the donor pool and say so.

## 2. Choose the variant

| Variant | Use when | Call |
|---|---|---|
| Classic SC | treated level inside the donors' convex hull | `synth_control()` |
| Demeaned SC | level differences; fit good after removing unit means | `synth_control(demean = TRUE)` (default recommendation) |
| Ridge-augmented | imperfect fit even demeaned | `synth_control(demean = TRUE, augment = "ridge")` |
| With covariates | predictors beyond outcome lags | `synth_control(x = , v = "mspe")` |
| Synthetic DiD | unit and time weights, DiD-style inference | `sdid_weights()` + `sdid_se()` |
| Staggered adoption | few cohorts of treated units | `synth_control(treated_units = "separate")` or SDID per cohort |

Tie-breakers: run classic and demeaned always; if the classic weights
cannot match the level (nonzero demeaned intercept), the demeaned
estimator is the honest one (Ferman-Pinto request 6). SDID when a DiD
audience wants unit-plus-time weighting and placebo/jackknife SEs.

## 3. Workflow

### Step 0: design argument and donor pool

Do: state why treatment timing is unrelated to treated-specific shocks;
list donors and exclude any plausibly affected by the treatment
(spillovers) or hit by their own contemporaneous policy.
Judge: every exclusion has a reason written down; the pool is not
selected on post-treatment outcomes.

### Step 1: fit

Do:
```r
fit <- synth_control(dat, id = "id", time = "time", y = "y", d = "d")
fit_dm <- synth_control(dat, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
fit$weights; fit$balance; fit$concentration
plot_synth(fit_dm, type = "path")
plot_synth(fit_dm, type = "demeaned_path")
```
Look: pre-treatment fit in levels and demeaned; weight concentration
(the l2 norm and effective donor count); the balance table when
covariates are used.
Judge: Ferman-Pinto requests 2-4: good fit after removing common trends,
donors few relative to T0, weights not concentrated on one or two donors
(diluted weights are the case with asymptotic-bias protection).
Fail: fit only in levels -> demeaned estimator; fit poor everywhere ->
augment with ridge or concede the design; weights on one donor -> the
"synthetic" control is that donor; step 5 decides if that is acceptable.

### Step 2: the estimate

Do:
```r
c(sc = fit$estimate, demeaned = fit_dm$estimate)
fit_aug <- synth_control(dat, id = "id", time = "time", y = "y", d = "d",
                         demean = TRUE, augment = "ridge")
```
Look: post-period gap path (`plot_synth(type = "gap")`), the average
post-treatment estimate.
Judge: estimates stable across classic/demeaned/augmented when the fit
is genuinely good; large spread is model dependence to report.

### Step 3: the specification test against DiD

Do:
```r
st <- synth_spec_test(fit_dm)
st$p_value; plot_synth(st)
```
Judge: the demeaned synthetic control and DiD should agree when either is
unbiased; rejection says both are suspect (Ferman-Pinto request 5),
usually selection on factor loadings.

### Step 4: inference with few treated units

Do:
```r
pl <- synth_placebo(fit_dm)                       # in-space: refit on every donor
pl; plot_synth(pl); plot_synth(pl, type = "ratio")
synth_placebo(fit_dm, type = "time", placebo_time = t0 - 5)$estimate
ci <- synth_conformal(fit_dm)                     # conformal p-values / interval
ci
```
Look: the treated unit's post/pre MSPE ratio against the placebo
distribution (the p-value is its rank); backdated placebo estimates near
zero; the conformal set.
Judge: rank-based p below 0.1 with, say, 20+ donors; in-time placebo
flat; conformal interval excluding zero for the horizons claimed.
Fail: treated ratio unexceptional -> the estimate is within the noise the
donor pool generates; report it as such.

### Step 5: fragility to donors

Do:
```r
synth_loo(fit_dm)      # refit dropping each positive-weight donor
```
Judge: the estimate survives dropping any single donor; if one donor's
removal kills it, the design is a two-unit comparative case study and
must be presented as one.

### Step 6: synthetic DiD (bridge to panel audiences)

Do:
```r
w <- sdid_weights(dat, id = "id", time = "time", y = "y", d = "d")
se <- sdid_se(w, method = "placebo", n_reps = 200, seed = 1)
plot_sdid(w)
```
Judge: SDID's unit weights plus time weights relax the level-matching
burden; placebo SEs need enough donors; with staggered cohorts run per
cohort and combine with treated-cell weights (each cohort needs >= 2
pre-periods).

## 4. Report

The design paragraph (treated-specific-shock argument, donor exclusions);
path and gap figures, plus the demeaned-fit figure; the weights table with
concentration and effective donors; the estimate across
classic/demeaned/augmented/SDID; the spec-test p; the placebo figure with
the rank p-value and the in-time placebo; leave-one-out; T0, donor count,
and the Ferman-Pinto checklist answered line by line.

## 5. Pitfalls

- A perfect pre-fit with many donors and short T0 is overfitting, not
  quality (request 3).
- Level-matching failure hides in path plots; always show the demeaned
  fit.
- In-space placebos need the same specification per donor; poor-fit
  donors can be excluded by an MSPE limit, but report how many.
- Cohorts with a single pre-period cannot be fit (regularization needs
  first differences); exclude and say so.
- `sdid_se` refits per replication; with very large donor pools this is
  slow - reduce `n_reps` and say so, or report the jackknife.
- Weights are the estimator: publish them.

## 6. Function reference

```r
synth_control(data, id, time, y, d, x = NULL, lags = "all",
              v = c("equal","mspe","regression"), demean = FALSE,
              constraints = c("simplex","nonnegative","none"), zeta = 0,
              augment = c("none","ridge"), lambda = NULL,
              treated_units = c("average","separate"), weights = NULL)
# $estimate, $weights, $balance, $concentration (l2, effective_donors),
# $intercept (demeaned), glance(): n_donors, t_pre; plot_synth(x, type =
# c("path","gap","demeaned_path","weights","balance"))
synth_placebo(x, type = c("space","time"), placebo_time = NULL, mspe_limit = Inf)
# $p_value (MSPE-ratio rank); plot_synth(pl), plot_synth(pl, type = "ratio")
synth_conformal(x, null = 0, q = 1, n_perm = 999, level = 0.95, per_period = FALSE)
synth_spec_test(x, q = 2)          # $p_value; demeaned SC vs DiD
synth_loo(x, min_weight = 0.001)
sdid_weights(data, id, time, y, d, estimator = c("sdid","sc","did","difp"))
sdid_se(x, method = c("placebo","bootstrap","jackknife"), n_reps = 200, seed = NULL)
sim_synth_panel(n_donors = 20, n_treated = 1, t_pre = 20, t_post = 5,
                n_factors = 2, effect = 0, seed = NULL)
# columns: id, time, y, d, y0, tau; treated unit(s) adopt at t_pre + 1
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_synth_panel(n_donors = 20, t_pre = 20, t_post = 5, effect = 2, seed = 8)
fit <- synth_control(d, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
stopifnot(abs(fit$estimate - 2) < 1)
pl <- synth_placebo(fit)
stopifnot(pl$p_value < 0.2)
lo <- synth_loo(fit)
st <- synth_spec_test(fit)
stopifnot(is.finite(st$p_value))
d0 <- sim_synth_panel(n_donors = 20, t_pre = 20, t_post = 5, effect = 0, seed = 8)
fit0 <- synth_control(d0, id = "id", time = "time", y = "y", d = "d", demean = TRUE)
stopifnot(abs(fit0$estimate) < 1)
cat("synthetic-control self-check passed\n")
```
