---
name: rd-analysis
description: Regression discontinuity with causalmetrics and the rdrobust family. Use when treatment (or its probability, or a policy slope) changes discontinuously at a known cutoff of a running variable - eligibility scores, vote shares, age or income thresholds. Covers sharp, fuzzy, kink, and local-randomization variants, manipulation and placebo diagnostics, covariate adjustment, weak-first-stage-robust fuzzy inference, and extrapolation.
---

# Regression discontinuity with causalmetrics

## 1. Does this design apply?

The question: the effect of a treatment assigned by a rule `D = 1{X >= c}`
(or whose probability or slope changes at `c`), identified by continuity
of the potential-outcome means at the cutoff. The estimand is local: the
effect at `X = c`, for compliers if fuzzy. Estimation runs on the
`rdrobust` family directly (the package does not wrap it); causalmetrics
adds the plotting, the diagnostics bundle, covariate adjustment,
weak-IV-robust fuzzy inference, kink reporting, and extrapolation.

Checks to run first:

```r
summary(dat$x); mean(dat$x >= 0)                      # cutoff normalized to 0
table(dat$d, dat$x >= 0)                              # sharp (diagonal) or fuzzy
length(unique(dat$x))                                 # mass points / discrete score
hist(dat$x, breaks = 100)                             # heaping near the cutoff?
```

- A known cutoff and a continuously measured (or finely discrete) running
  variable with mass on both sides.
- Knowledge of who could manipulate the score and how precisely; the
  design's credibility is exactly the impossibility of precise
  manipulation.

Red flags pointing elsewhere: no one is actually assigned by the rule
(first stage flat) -> there is no design; the score is coarse with few
mass points -> local randomization within a window (`RDHonest`,
`rdlocrand`) rather than local polynomials; treatment timing varies in a
panel around a dated policy -> `did-analysis`; the cutoff instrument is
to be used away from the threshold -> that is `iv-analysis` plus an
explicit extrapolation argument.

## 2. Choose the variant

| Variant | Use when | Main call |
|---|---|---|
| Sharp RD | D jumps 0 to 1 at c | `rdrobust::rdrobust(y, x, c = 0)` |
| Fuzzy RD | P(D=1) jumps at c | `rdrobust(y, x, fuzzy = d)`; weak case `rd_weak_iv()` |
| Kink | the policy's slope changes at c | `rd_kink()` |
| Local randomization | discrete score, small windows | `rdlocrand`/`RDHonest`; `rd_frame()` tidies |
| Covariate-adjusted | precision, or CIA-based extrapolation | `rd_adjust()`, `rd_extrapolate()` |
| Multiple cutoffs | different c by group | `rdmulti`; normalize-and-pool + per-cutoff |

Tie-breakers: sharp vs fuzzy is a fact of the data (the D-jump), not a
choice. Robust bias-corrected inference at the MSE-optimal bandwidth is
the default report. Honest (RDHonest) intervals accompany it when the
smoothness assumption is contested.

## 3. Workflow

### Step 0: the rule and the actors

Do: state the rule, who runs it, who knows their score, and how precisely
the score can be influenced.
Judge: manipulation must be imprecise for continuity to be credible.
Fail: precise manipulation (retaking, rounding by caseworkers) -> expect
the density test to fail; consider donut estimates and say what they
change.

### Step 1: see it - the RD plot

Do:
```r
rd_plot(dat, y = "y", x = "x", cutoff = 0, fit = "local", h = NULL, ci = TRUE)
```
Look: binned means, the local fits on each side, the jump; curvature and
outliers near the cutoff.
Judge: a jump visible against the binned noise; global polynomial
artifacts absent (the plot's fits are local, deliberately).
Fail: no visible jump does not kill the design (power), but a jump that
appears only under a global polynomial is a warning.

### Step 2: main estimate

Do:
```r
fit <- rdrobust::rdrobust(dat$y, dat$x, c = 0)   # sharp
summary(fit)
rd_frame(main = fit)                              # tidy conventional/bias-corrected/robust rows
# fuzzy:
ffit <- rdrobust::rdrobust(dat$y, dat$x, c = 0, fuzzy = dat$d)
```
Look: robust bias-corrected estimate and CI, MSE-optimal `h`, effective N
on each side.
Judge: report the robust row; conventional-at-optimal-h intervals
undercover by construction.
Fail: tiny effective N on one side -> widen the design discussion, not
the bandwidth by hand.

### Step 3: the diagnostics bundle (always)

Do:
```r
ck <- rd_checks(dat, y = "y", x = "x", cutoff = 0, covariates = zvars, d = NULL)
ck
plot_rd_checks(ck, "balance")
plot_rd_checks(ck, "sensitivity")
plot_rd_checks(ck, "placebo")
summary(rddensity::rddensity(dat$x, c = 0))       # manipulation test
dn <- rd_donut(dat, y = "y", x = "x", cutoff = 0)
plot_rd_checks(dn)
```
Look: density-test p; covariate jumps with the joint test; the estimate
across bandwidths; placebo cutoffs; donut radii.
Judge: density p above 0.05; no covariate jumps (joint test); estimate
stable across h and dead at placebo cutoffs; donut stable.
Fail: density rejects -> manipulation; lead with donut estimates and the
manipulation story, or abandon. A covariate jump -> the "cutoff" bundles
another change; identify and discuss it.

### Step 4: covariates - precision, not identification

Do:
```r
adj <- rd_adjust(dat, y = "y", x = "x", cutoff = 0, covariates = zvars, seed = 1)
rd_frame(raw = adj$fits$raw, adjusted = adj$fits$adjusted, methods = "robust")
```
Judge: the adjusted estimate should match the raw one with a smaller SE;
a moved point estimate signals covariate imbalance (back to step 3).

### Step 5: fuzzy specifics

Do:
```r
wk <- rd_weak_iv(dat, y = "y", d = "d", x = "x", cutoff = 0)
wk        # first-stage effective F at the cutoff + Anderson-Rubin interval
```
Judge: report the first stage and its effective F; when weak, the AR
interval replaces the Wald one. The estimand is the compliers-at-the-
cutoff effect; say so.

### Step 6: kink designs

Do:
```r
kk <- rd_kink(dat, y = "y", x = "x", cutoff = 0, d = "d", slope_change = 1, elasticity = TRUE)
plot_rd_kink(kk)
```
Judge: kinks estimate slope changes (first derivatives): noisier, more
bandwidth-sensitive; the four-panel plot (outcome, first stage, density,
covariates) is the standard evidence.

### Step 7: beyond the cutoff (only with an explicit assumption)

Do:
```r
ex <- rd_extrapolate(dat, y = "y", x = "x", cutoff = 0, covariates = zvars,
                     window = NULL, method = "aipw", seed = 1)
plot_rd_extrapolate(ex)
```
Judge: extrapolation swaps continuity for a CIA within a window; report
it as a different estimand under a stated assumption, next to (never
instead of) the cutoff effect.

## 4. Report

The rule paragraph (who, how precise); the RD plot with local fits plus
density and covariate plots; robust bias-corrected estimate at the
MSE-optimal h with effective Ns; the full checks bundle (density,
balance + joint test, placebo cutoffs, bandwidth sensitivity, donut);
fuzzy: first stage, effective F, AR interval; estimates with and without
covariates; the estimand sentence (at the cutoff; compliers if fuzzy;
window-average if local randomization; extrapolation assumption if any).

## 5. Pitfalls

- rdrobust's default `vce = "nn"`: match it when reproducing published
  numbers or comparisons will differ.
- Global high-order polynomials manufacture jumps; the plot's fit and the
  estimator are local by design.
- Discrete scores: cluster on score values, consider `RDHonest`;
  `rdlocrand::rdwinselect` can return NA windows with defaults - set
  `wmin`/`wstep`.
- Donut and density tests are complements: a passed density test with
  heaped data still warrants a donut row.
- Covariate adjustment that moves the point estimate is a red flag, not a
  fix.
- In fuzzy designs never instrument with the treatment itself elsewhere
  in the specification; the instrument is the threshold indicator.

## 6. Function reference

```r
# Estimation: rdrobust::rdrobust / rdbwselect / rdplot, RDHonest, rdlocrand, rdmulti (direct).
rd_plot(data, y, x, cutoff = 0, bins = c("qs","es"), p = 4,
        fit = c("polynomial","local","none"), h = NULL, ci = FALSE)
rd_bins(data, y, x, cutoff = 0, bins = c("qs","es"), n_bins = NULL)
rd_checks(data, y, x, cutoff = 0, covariates = NULL, d = NULL, cluster = NULL,
          h_grid = NULL, radius = NULL, cutoffs = NULL)   # + plot_rd_checks(x, type)
rd_balance(data, covariates, x, cutoff = 0)               # per-covariate + joint test
rd_placebo_cutoffs(data, y, x, cutoff = 0, cutoffs = NULL)
rd_donut(data, y, x, cutoff = 0, radius = NULL)
rd_sensitivity(data, y, x, cutoff = 0, h_grid = NULL)
rd_adjust(data, y, x, cutoff = 0, covariates, learner = NULL, folds = 5, seed = NULL)
# $fits$raw, $fits$adjusted (rdrobust objects)
rd_weak_iv(data, y, d, x, cutoff = 0, h = NULL)           # fuzzy: effective F + AR interval
rd_kink(data, y, x, cutoff = 0, d = NULL, slope_change = 1, elasticity = FALSE)
rd_extrapolate(data, y, x, cutoff = 0, covariates, method = c("regression","aipw"))
rd_frame(..., methods = c("conventional","bias_corrected","robust"))  # tidy comparison table
sim_rd(n, dgp = c("lee","fuzzy","covariates","cia","kink","discrete","manipulated","multi_cutoff"),
       tau = NULL, seed = NULL)
# lee: x, d, y (+ mu0, mu1), attr tau_true, cutoff; fuzzy adds z
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_rd(n = 3000, dgp = "lee", seed = 9)
tau <- attr(d, "tau_true")
fit <- rdrobust::rdrobust(d$y, d$x, c = 0)
est <- fit$coef["Robust", ]; se <- fit$se["Robust", ]
stopifnot(abs(est - tau) < 3 * se)
fr <- rd_frame(main = fit)
stopifnot(any(grepl("robust", fr$method, ignore.case = TRUE)))
dm <- sim_rd(n = 3000, dgp = "manipulated", seed = 9)
p_man <- rddensity::rddensity(dm$x, c = 0)$test$p_jk
p_ok <- rddensity::rddensity(d$x, c = 0)$test$p_jk
stopifnot(p_man < 0.05, p_ok > 0.05)
cat("rd-analysis self-check passed\n")
```
