---
name: selection-on-observables
description: Adjustment-based causal analysis (regression, IPW, matching, AIPW) with causalmetrics. Use when treatment was not randomized and identification rests on a defensible claim that observed pre-treatment covariates absorb all confounding. Covers estimand choice, overlap and trimming, est_aipw with internal or external nuisances, balance, estimator-family comparison, and sensitivity.
---

# Selection on observables with causalmetrics

## 1. Does this design apply?

The question: the ATE or ATT of a binary treatment when the only
identification claim is conditional ignorability, `Y(0), Y(1) independent
of D given X`, plus overlap `0 < P(D = 1 | X) < 1`. This claim cannot be
tested; the workflow's job is to make it explicit, check what is
checkable (overlap, balance), and quantify fragility (sensitivity).

Checks to run first:

```r
table(dat$d); mean(dat$d)
summary(glm(d ~ ., data = dat[, c("d", xvars)], family = binomial))  # crude selection scan
sapply(dat[xvars], function(v) sum(is.na(v)))
```

- A defensible adjustment set of pre-treatment covariates. The set is a
  design choice, argued variable by variable; it is not "everything in
  the file".
- Enough treated and control units across the covariate space.

Red flags pointing elsewhere: any usable design variation beats this —
randomization (`rct-analysis`), an instrument (`iv-analysis`), a cutoff
(`rd-analysis`), adoption timing in a panel (`did-analysis`), one treated
aggregate (`synthetic-control`). High-dimensional or clearly nonlinear
nuisances -> do this workflow but estimate with `double-ml`. Interest in
heterogeneity or targeting -> finish here, then `hte-policy`.

## 2. Choose the variant

| Variant | Use when | Main call |
|---|---|---|
| Regression adjustment | few covariates, outcome model trusted | `fixest::feols(y ~ d + x...)` as baseline only |
| IPW | propensity easier to model than outcome | `est_aipw(..., learner_mu0 = NULL)` comparison, or dr_scores type "ipw" |
| Matching | audience wants a design-like sample; ATT | `MatchIt::matchit()` as engine, then outcome model on matched data |
| AIPW (default) | both nuisances estimable; want double robustness | `est_aipw()` |
| AIPW with ML nuisances | flexible functional forms needed | `est_aipw(learner_* = mlr3 learners)` or external predictions |

Tie-breakers: AIPW is the default deliverable; regression and IPW are
comparison rows. ATT (`estimand = "ATT"`) when the treated population is
the policy question or when controls far from any treated unit make the
ATE ill-supported. Matching does not add identification over weighting;
use it when the trimmed, interpretable sample is itself wanted.

## 3. Workflow

### Step 0: the adjustment-set argument

Do: list candidate covariates; for each, state pre-treatment status and
the confounding path it blocks. Remove mediators, colliders, instruments
of D, and anything measured after treatment.
Judge: every kept variable has a one-line causal justification.
Fail: key confounder unmeasured and no proxy -> proceed only with the
sensitivity step promoted to headline status, or abandon the design.

### Step 1: overlap before anything

Do:
```r
ps <- glm(reformulate(xvars, "d"), data = dat, family = binomial)
dat$phat <- fitted(ps)
tapply(dat$phat, dat$d, summary)
hist(dat$phat[dat$d == 1]); hist(dat$phat[dat$d == 0])
```
Look: propensity ranges by arm; mass of controls where treated live and
vice versa.
Judge: no arm's support extends where the other has essentially no mass.
Fail: poor overlap -> trim (`trim = c(0.05, 0.95)` style bounds in the
estimator call) or switch to ATT; either changes the target population,
which must be reported (who was dropped, how the covariate means moved).

### Step 2: main estimate, AIPW

Do:
```r
fit <- est_aipw(dat, y = "y", d = "d", x = xvars, estimand = "ATE",
                folds = 5, seed = 1)
tidy(fit); glance(fit)
fit$diagnostics
# ATT version:
fit_att <- est_aipw(dat, y = "y", d = "d", x = xvars, estimand = "ATT", seed = 1)
```
Look: estimate and SE; `$diagnostics$propensity` (range, `n_clipped_*`,
common support), `$diagnostics$weights` (IPW tail, effective sample
sizes `ess_treated`/`ess_control`), `$diagnostics$fold_estimates`.
Judge: clipped share small (a few percent); nuisance models fit better
than the marginal mean.
Fail: large clipped share -> back to step 1; the estimand is drifting to
an overlap-weighted population.

### Step 3: nuisance flexibility

Do: refit with flexible learners and with external predictions when a
better stack exists outside R:
```r
library(mlr3); library(mlr3learners)
fit_rf <- est_aipw(dat, y = "y", d = "d", x = xvars,
                   learner_p = lrn("classif.ranger", predict_type = "prob"),
                   learner_mu0 = lrn("regr.ranger"), learner_mu1 = lrn("regr.ranger"),
                   folds = 5, seed = 1)
# external out-of-fold predictions (vignette usage_aipw_external_nuisances):
fit_ext <- est_aipw(dat, y = "y", d = "d", p_hat = dat$p_ext,
                    mu0_hat = dat$mu0_ext, mu1_hat = dat$mu1_ext)
```
Judge: parametric and flexible estimates within an SE of each other.
Fail: a large gap is model dependence; report both, prefer the flexible
one only with a stated reason, and consider `double-ml` (cross-fitting
discipline, sensitivity bounds).

### Step 4: balance after design

Do: for the weighted (or matched) sample, standardized differences of
each covariate; with matching:
```r
m <- MatchIt::matchit(reformulate(xvars, "d"), data = dat, method = "nearest", ratio = 1)
summary(m)             # balance table
md <- MatchIt::match.data(m, data = dat)
```
Judge: absolute standardized differences below ~0.1 after design.
Fail: rebuild the design (different distance, calipers, more flexible
propensity); balance failure is a design failure, not an estimation one.

### Step 5: estimator-family comparison

Do: one table with regression, IPW, matching+regression, AIPW
(parametric), AIPW (ML).
Judge: agreement across families says results are not artifacts of one
functional form; it does not certify ignorability, and the report must
say so.

### Step 6: placebo and sensitivity

Do: a placebo outcome (determined before treatment) through the same
pipeline should give zero. Then quantify hidden-confounding fragility;
with a DML fit available:
```r
dml <- est_dml(dat, y = "y", d = "d", x = xvars, model = "plr", seed = 1)
sens <- dml_sensitivity(dml, r2_y = 0.05, r2_d = 0.05)
plot_dml_sensitivity(sens)
```
Look: the confounder strength (partial R2 with treatment and outcome)
that drives the estimate to zero; benchmark against the observed
covariates' partial R2s.
Judge: a robustness value above the strongest observed covariate is
comfort; below it is a warning to print in the abstract.

## 4. Report

The adjustment-set paragraph (each variable's justification); overlap
figures and the trimming rule with the population change; the estimator
comparison table with nuisance choices named; the balance table; the
placebo result; the sensitivity statement with its benchmark. State
plainly: AIPW's double robustness covers nuisance misspecification, not
hidden confounding, bad controls, or missing support.

## 5. Pitfalls

- Bad controls do the most damage: one mediator or collider in `x`
  poisons every estimator equally; agreement across estimators does not
  detect it.
- Overlap and balance are different: overlap is support (checked before),
  balance is a post-design diagnostic (checked after).
- Trimming silently changes the estimand; always redescribe the analysis
  population.
- Extreme IPW weights: inspect the weight tail and effective sample size;
  AIPW tempers but does not remove the problem.
- Cross-fitting prevents overfit nuisances from leaking; it creates no
  identification.
- External nuisance predictions must be out-of-fold; in-sample
  predictions bias the score.

## 6. Function reference

```r
est_aipw(data, y, d, x = NULL, estimand = "ATE",           # or "ATT"
         p_hat = NULL, mu0_hat = NULL, mu1_hat = NULL,      # external out-of-fold preds
         learner_p = NULL, learner_mu0 = NULL, learner_mu1 = NULL,
         folds = 5, cross_fit = TRUE, p_clip = c(0.01, 0.99), trim = NULL,
         seed = NULL, conf_level = 0.95)
# tidy(): estimate/std.error/conf.*; glance(); $diagnostics: $sample, $propensity
# (raw_summary, n_clipped_low/high, common_support), $trimming, $weights
# (ipw_summary, ess_treated/ess_control), $nuisance, $fold_estimates
dr_scores(...)               # same nuisances, returns per-unit scores for hte-policy
est_dml(...) + dml_sensitivity(fit, r2_y, r2_d)   # sensitivity bounds (see double-ml skill)
sim_hte(n, dgp = c("smooth","unbalanced", ...), seed)  # y,d,x1..x5,tau_true,p_true
# MatchIt is the matching engine by design; the package does not wrap it.
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_hte(n = 2500, dgp = "smooth", seed = 3)   # assignment depends on x1, x2
ate_true <- mean(d$tau_true)
naive <- mean(d$y[d$d == 1]) - mean(d$y[d$d == 0])
fit <- est_aipw(d, y = "y", d = "d", x = paste0("x", 1:5), estimand = "ATE", seed = 1)
td <- tidy(fit)
stopifnot(abs(td$estimate - ate_true) < 3 * td$std.error)
stopifnot(abs(td$estimate - ate_true) < abs(naive - ate_true) + 0.05)  # beats naive
clip_n <- fit$diagnostics$propensity$n_clipped_low + fit$diagnostics$propensity$n_clipped_high
stopifnot(clip_n / nrow(d) < 0.05)
cat("selection-on-observables self-check passed\n")
```
