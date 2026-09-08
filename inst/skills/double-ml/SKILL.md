---
name: double-ml
description: Double/debiased machine learning with causalmetrics. Use when identification is selection-on-observables (or its IV/panel analogues) but the nuisance functions need flexible ML - many covariates, unknown functional forms - or when a structured parametric outcome model with combinatorial treatments needs debiased inference (est_dml_structured). Covers plr vs irm, cross-fitting choices, learners, diagnostics, and sensitivity.
---

# Double machine learning with causalmetrics

## 1. Does this design apply?

The question: the same estimands as `selection-on-observables` (a
treatment coefficient, ATE, ATT) when the adjustment set is large or the
relationships are nonlinear enough that parametric nuisances are not
credible. DML adds no identification: conditional ignorability and
overlap still carry everything. What it adds is valid root-n inference
with ML nuisances, through Neyman-orthogonal scores and cross-fitting.

Checks: the same as `selection-on-observables` step 0-1 (adjustment-set
argument, overlap). Plus:

```r
length(xvars); nrow(dat)                 # p vs n motivates ML
sapply(dat[xvars], function(v) length(unique(v)))
```

Red flags: an instrument is available -> `iv-analysis` (est_dml pliv/iivm
live there); panel timing variation -> `did-analysis` (att_gt already
cross-fits); target is heterogeneity -> `hte-policy` (dr_scores reuses
these nuisances); the outcome model is a known parametric structure with
engineered treatments (demand-like settings) -> the structured branch in
Section 3, step S.

## 2. Choose the variant

| Variant | Model | Estimand | Use when |
|---|---|---|---|
| Partially linear (`model = "plr"`) | `Y = D*theta + g(X) + e` | overlap-weighted coefficient | effect plausibly constant, or a summary coefficient wanted |
| Interactive (`model = "irm"`) | fully heterogeneous | ATE or ATT | effects vary with X; binary D |
| Structured (`est_dml_structured`) | parametric outcome in engineered treatments | contrasts across treatment combinations | deep-learning first stage estimates theta(x) of a known link |

Tie-breakers: `irm` with `estimand = "ATE"` is the default causal
deliverable; `plr` when D is not binary or as the sensitivity vehicle
(`dml_sensitivity()` takes a plr fit). Under heterogeneity, the plr
coefficient is an overlap-weighted average, not the ATE; say which one
you report. `solve = "pooled"` (DML2) is the default; `"fold_average"`
(DML1) only for diagnostics. `n_rep >= 3` when n is small or learners are
unstable.

## 3. Workflow

### Step 0-1: inherit from selection-on-observables

The adjustment-set argument and the overlap check are identical; do them
first.

### Step 2: baseline with default learners

Do:
```r
fit_lm <- est_dml(dat, y = "y", d = "d", x = xvars, model = "irm",
                  estimand = "ATE", folds = 5, seed = 1)
tidy(fit_lm); glance(fit_lm)
```
Look: `glance()` reports the learners and nuisance RMSEs
(`rmse_p_hat`, `rmse_mu_hat`); `$diagnostics$fold_estimates` the
fold-level spread.
Judge: this is the parametric benchmark; keep it for the comparison
table.

### Step 3: ML nuisances, cross-fitted

Do:
```r
library(mlr3); library(mlr3learners)
fit_rf <- est_dml(dat, y = "y", d = "d", x = xvars, model = "irm", estimand = "ATE",
                  learner_p = lrn("classif.ranger", num.trees = 500, predict_type = "prob"),
                  learner_mu0 = lrn("regr.ranger", num.trees = 500),
                  learner_mu1 = lrn("regr.ranger", num.trees = 500),
                  folds = 5, n_rep = 3, seed = 1)
fit_rf$diagnostics$nuisance     # per-nuisance rmse and r2
```
Look: nuisance RMSE/R2 against the baseline learners; the propensity
range and clipped share; the spread across repetitions.
Judge: choose learners by nuisance fit, never by the theta they produce.
Remember the rate: orthogonality makes the bias second order, so
nuisances need only n^(-1/4) accuracy, but a forest whose propensity
RMSE sits at that boundary still leaves visible bias at moderate n (the
lecture's running teaching point) - more data or better learners, not a
different theta, is the fix.
Fail: propensity mass at the clip bounds -> trim and redescribe the
population, or switch to ATT.

### Step 4: external nuisances (optional)

Do: when the best predictions come from outside R (a Python stack), pass
out-of-fold predictions directly:
```r
fit_ext <- est_dml(dat, y = "y", d = "d", model = "irm", estimand = "ATE",
                   p_hat = dat$p_oof, mu0_hat = dat$mu0_oof, mu1_hat = dat$mu1_oof,
                   fold_id = dat$fold)
```
Judge: predictions must be out-of-fold with the fold map supplied;
in-sample predictions invalidate the inference (state in the report that
they are out-of-fold and how folds were built).

### Step 5: specification spread and solve method

Do: assemble the comparison: baseline learners, ML learners, plr vs irm,
pooled vs fold_average, and the repetition spread.
Judge: a wide learner-to-learner gap is model dependence; report the
range, not the favourite.

### Step 6: sensitivity (always)

Do:
```r
fit_plr <- est_dml(dat, y = "y", d = "d", x = xvars, model = "plr", seed = 1)
sens <- dml_sensitivity(fit_plr, r2_y = 0.05, r2_d = 0.05)
sens$robustness_value; sens$bounds
plot_dml_sensitivity(sens)
```
Look: the robustness value (the equal partial-R2 of an unobserved
confounder with outcome and treatment that drives theta to zero) and the
benchmark bounds.
Judge: compare the robustness value with the partial R2s of the strongest
observed covariates; below them, the headline must carry the warning.

### Step S: structured outcome models (the Farrell-Liang-Misra branch)

Use when the outcome follows a known parametric form in engineered
treatments (e.g. combinatorial feature bundles) and a flexible first
stage estimates the unit-level parameters theta(x). The first stage
(any engine; PyTorch script ships at `inst/python/dedl_nuisance.py`)
must be cross-fitted; the package then debiases contrasts:
```r
est <- est_dml_structured(theta_hat = theta_oof, y = dat$y, t_obs = t_obs,
                          targets = target_list, link = "gen_sigmoid",
                          fold_id = dat$fold, ridge = 0.01)
```
Judge: `Lambda(x)` near-singularity (saturated links) makes the
correction explode; a small `ridge` (0.01 in the Ye replication) is the
documented fix - report it. Contrasts to the best bundle come from
`contrast_best = TRUE`.

## 4. Report

Estimand and model named (plr coefficient vs ATE/ATT); learners per
nuisance with cross-fitted RMSE next to the estimate; propensity range,
clipped/trimmed shares, effective sample sizes; K, repetitions, solve
method; fold- and repetition-level spread; the learner comparison; the
sensitivity statement with benchmarks; the out-of-fold statement for any
external predictions. DML does not deliver identification - one placebo
or sensitivity analysis is part of the deliverable.

## 5. Pitfalls

- Choosing learners by the resulting theta is specification search;
  choose by nuisance fit.
- plr under heterogeneity is overlap-weighted; do not call it the ATE.
- Clipping (`p_clip`) is a numerical guard, trimming (`trim`) changes the
  population; report the latter.
- Structured DML without a ridge explodes when the link saturates.
- `dml_sensitivity()` needs a plr fit; run one alongside irm for this
  purpose.
- Repetition spread wider than the SE means the partition matters; raise
  `n_rep` and report the median-aggregated estimate.

## 6. Function reference

```r
est_dml(data, y, d, x = NULL, z = NULL, model = c("plr","irm","pliv","iivm"),
        estimand = c("ATE","ATT"),
        l_hat = NULL, m_hat = NULL, p_hat = NULL, mu0_hat = NULL, mu1_hat = NULL,
        fold_id = NULL, learner_l = NULL, learner_m = NULL, learner_p = NULL,
        learner_mu0 = NULL, learner_mu1 = NULL,
        folds = 5, n_rep = 1, solve = c("pooled","fold_average"),
        p_clip = c(0.01, 0.99), trim = NULL, seed = NULL)
# tidy(), glance() (learners, rmse_p_hat, rmse_mu_hat);
# $diagnostics: $nuisance (target, rmse, r2), $fold_estimates, $repetitions
dml_sensitivity(fit, r2_y = 0.05, r2_d = 0.05)
# $robustness_value, $bounds, $bounds_ci, $contour; plot_dml_sensitivity(x)
est_dml_structured(theta_hat, y, t_obs, targets, t0 = NULL,
                   link = c("gen_sigmoid","linear","logit","custom"),
                   G = NULL, G_grad = NULL, t_dist = NULL, fold_id = NULL,
                   ridge = 0, contrast_best = TRUE)
bind_scores(...)      # stack cm score objects across samples/arms
sim_hte(n, dgp = "smooth", seed)   # y, d, x1..x5, tau_true, p_true
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_hte(n = 2000, dgp = "smooth", seed = 5)
ate_true <- mean(d$tau_true)
irm <- est_dml(d, y = "y", d = "d", x = paste0("x", 1:5), model = "irm",
               estimand = "ATE", folds = 5, seed = 1)
td <- tidy(irm)
stopifnot(abs(td$estimate - ate_true) < 3 * td$std.error)
plr <- est_dml(d, y = "y", d = "d", x = paste0("x", 1:5), model = "plr", seed = 1)
sens <- dml_sensitivity(plr, r2_y = 0.05, r2_d = 0.05)
stopifnot(is.finite(sens$robustness_value), sens$robustness_value > 0)
stopifnot(all(c("rmse", "r2") %in% names(plr$diagnostics$nuisance)))
cat("double-ml self-check passed\n")
```
