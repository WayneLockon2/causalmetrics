---
name: iv-analysis
description: Instrumental variables with causalmetrics. Use when an instrument moves treatment but plausibly affects the outcome only through it - lotteries, encouragements, judge/examiner assignment, shift-share exposure, policy discontinuities used as instruments. Covers 2SLS, weak-instrument-robust inference, complier description, DML-IV (PLIV/LATE), control functions, MTE, and exclusion sensitivity.
---

# Instrumental variables with causalmetrics

## 1. Does this design apply?

The question: the causal effect of an endogenous treatment `D` on `Y`,
using an instrument `Z` that (i) moves `D` (relevance, testable), (ii) is
as-good-as-random given controls `X` (independence), and (iii) affects
`Y` only through `D` (exclusion, not testable - it is the argument the
report stands on). With heterogeneous effects and a binary instrument,
the estimand is the LATE: the effect for compliers, whose treatment the
instrument changes.

Checks to run first:

```r
summary(fixest::feols(reformulate(c("z", xvars), "d"), data = dat, vcov = "hetero"))
cor(dat$z, dat[xvars])            # instrument balance scan
table(dat$z, dat$d)               # compliance pattern (binary case)
```

Red flags pointing elsewhere: the "instrument" is just another covariate
with no exclusion story -> `selection-on-observables`; the instrument is
a cutoff in a score -> `rd-analysis` (fuzzy RD is IV at the cutoff; that
skill calls back here for weak-IV tools); assignment itself randomized
and ITT is the deliverable -> `rct-analysis` first; the instrument moves
a mediator rather than the treatment -> `mediation-analysis`
(`mediate_iv`).

## 2. Choose the variant

| Variant | Use when | Main call |
|---|---|---|
| 2SLS (baseline) | low-dimensional linear controls | `fixest::feols(y ~ x1 | d ~ z)` (registry: fixest direct) |
| Partially linear IV | many/nonlinear controls, one endogenous D | `est_dml(model = "pliv")` |
| Doubly robust LATE | binary Z and D, heterogeneous effects | `est_dml(model = "iivm")` |
| Judge / leniency | random examiner assignment | `leniency_instrument()` then PLIV/2SLS |
| Shift-share | shares x aggregate shocks | `ssiv_rotemberg()`, `ssiv_shock_level()` |
| Control function | nonlinear second stage (probit/logit outcome or demand model) | `cf_residuals()` + `cf_bootstrap()` + `cf_ape()` |
| Marginal effects over p(Z) | policy question at other margins | `mte_curve()` |

Tie-breakers: always report 2SLS with the same controls and sample as the
first stage and reduced form. Move to `pliv` when the control set is the
threat, to `iivm` when the LATE interpretation is the point. Multiple
instruments: report per-instrument estimates and the 2SLS weights before
pooling.

## 3. Workflow

### Step 0: the exclusion argument

Do: write the structural equation, say what the error contains, and give
the mechanism for exclusion in one paragraph. Name what would violate it.
Judge: a named violation channel you can probe later (step 6) exists.
Fail: no articulable exclusion story -> stop; an instrument without a
mechanism is a specification search.

### Step 1: first stage, reduced form, IV - one sample, same controls

Do:
```r
fs <- iv_first_stage(dat, d = "d", z = "z", x = xvars, cluster = "cl")
fs                                       # F, robust F, effective F
rf <- fixest::feols(reformulate(c("z", xvars), "y"), data = dat, cluster = ~cl)
iv <- fixest::feols(y ~ x1 + x2 | d ~ z, data = dat, cluster = ~cl)
```
Look: the three F's; the reduced-form and first-stage signs (their ratio
is the IV estimate).
Judge: effective F above ~20-25 for reliable conventional inference; 10
is the old rule, not a safe harbor.
Fail: weak -> step 2 carries the inference regardless; do not drop the
design silently.

### Step 2: weak-instrument-robust inference (always report)

Do:
```r
ar <- iv_ar_confidence_set(dat, y = "y", d = "d", z = "z", x = xvars, cluster = "cl")
ar                                       # Anderson-Rubin set; may be unbounded
plot_ar_set(ar)
```
Judge: the AR set is valid at any strength; when it is unbounded or far
wider than the Wald interval, the Wald interval is not to be trusted.
Fail: unbounded set -> report it as the honest interval; the design
delivers sign/set information, not a point.

### Step 3: whose effect is it (binary Z, D)

Do:
```r
cp <- complier_profile(dat, d = "d", z = "z", covariates = xvars)
tidy(cp)                                 # complier shares and characteristics
w <- iv_late_weights(dat, y = "y", d = "d", z = "z", x = xvars)  # several instruments
tidy(w)
b <- iv_ate_bounds(dat, y = "y", d = "d", z = "z")               # from LATE toward ATE
```
Judge: state complier share and how compliers differ from the population;
with several instruments, report instrument-specific estimates and the
2SLS weights (negative weights are a warning).

### Step 4: orthogonal IV with cross-fitting

Do:
```r
library(mlr3); library(mlr3learners)
pliv <- est_dml(dat, y = "y", d = "d", z = "z", x = xvars, model = "pliv",
                learner_l = lrn("regr.ranger"), learner_m = lrn("regr.ranger"),
                learner_z = lrn("regr.ranger"), folds = 5, n_rep = 3, seed = 1)
tidy(pliv); glance(pliv)
late <- est_dml(dat, y = "y", d = "d", z = "z", x = xvars, model = "iivm",
                weak_iv = TRUE, seed = 1)
tidy(late)
```
Look: nuisance RMSEs; the residual first stage (instrument strength after
partialling out); with `weak_iv = TRUE`, the score-based AR set.
Judge: PLIV agrees with 2SLS when controls are simple; disagreement
means the control set mattered - report both.
Fail: residual first stage collapses after flexible partialling -> the
instrument works only through the controls; the design fails.

### Step 5: nonlinear outcomes - control function

Do:
```r
dat2 <- cf_residuals(dat, d = "d", z = "z", x = xvars, family = "gaussian")
second <- function(dd) glm(reformulate(c("d", xvars, "v_hat"), "y_bin"),
                           data = dd, family = binomial())
cb <- cf_bootstrap(dat2, first = function(dd) cf_residuals(dd, d = "d", z = "z",
                    x = xvars, family = "gaussian"), second = second,
                    n_boot = 199, seed = 1)
```
Judge: the control-function coefficient on `v_hat` is the endogeneity
test; report average partial effects (`cf_ape()`) with bootstrap SEs, not
raw coefficients.

### Step 6: exclusion sensitivity and constructed-instrument diagnostics

Do:
```r
pe <- iv_plausibly_exogenous(dat, y = "y", d = "d", z = "z", x = xvars,
                             gamma_grid = seq(-0.2, 0.2, 0.05))
tidy(pe)                        # bounds as direct effect gamma varies
# shift-share:
rot <- ssiv_rotemberg(dat, y = "y", d = "d", shares = share_cols, shocks = "shock")
tidy(rot)                       # which shocks carry the weight
sl <- ssiv_shock_level(dat, y = "y", d = "d", shares = share_cols, shocks = "shock")
# judges:
dat3 <- leniency_instrument(dat, examiner = "judge", d = "d", x = xvars)
```
Judge: how large a direct effect of Z on Y flips the conclusion; for
shift-share, report the Rotemberg weights (a few dominant shocks mean the
design is really those shocks' event studies); for judges, leave-one-out
leniency and balance across examiners.

### Step 7: MTE when the policy question is off the LATE margin

Do:
```r
mte <- mte_curve(dat, y = "y", d = "d", z = "z", x = xvars,
                 method = "polynomial", degree = 3, n_boot = 99, seed = 1)
plot_mte(mte); tidy(mte)
```
Judge: with a continuous-support propensity, the MTE curve says how the
effect varies with unobserved resistance; the LATE is one weighted
average of it. Flat curve -> LATE generalizes; steep -> say whose margin
the policy moves.

## 4. Report

The exclusion paragraph with its named threat; first stage, reduced form,
and IV in one table (same controls, same sample, clustered where the
instrument varies); the effective F and the AR interval regardless of F;
complier share and profile (or Rotemberg weights, or examiner
diagnostics); the DML-IV row with nuisance fits when used; the
plausibly-exogenous bounds; APEs for nonlinear outcomes.

## 5. Pitfalls

- Different controls in first stage and IV is the classic silent error;
  build all three regressions from one formula object.
- Effective F, not the homoskedastic F, is the relevant strength measure
  with robust/clustered errors.
- An unbounded AR set is a result, not a bug; report it.
- With several instruments, overidentification tests reject under
  heterogeneous LATEs even when every instrument is valid.
- `iv_ar_confidence_set()` takes controls linearly; with fixed effects,
  absorb them first or include dummies in `x`.
- Judge designs: leniency must exclude the own case (`leniency_instrument()`
  does leave-one-out) and balance must be shown across examiners.
- LATE weights can be negative with multiple instruments; report them.

## 6. Function reference

```r
iv_first_stage(data, d, z, x = NULL, cluster = NULL, weights = NULL)   # F, robust F, effective F
iv_ar_confidence_set(data, y, d, z, x = NULL, cluster = NULL, weights = NULL,
                     conf_level = 0.95, theta_grid = NULL)   # AR set; plot_ar_set()
iv_plausibly_exogenous(data, y, d, z, x = NULL, gamma_grid = seq(-0.5, 0.5, 0.1))
iv_ate_bounds(data, y, d, z, n_boot = 199)
complier_profile(data, d, z, covariates, x = NULL)           # tidy(): complier means
iv_late_weights(data, y, d, z, x = NULL)                     # per-instrument LATEs + weights
late_scores(fit); late_blp(scores, formula, newdata)         # complier heterogeneity (fit = iivm)
leniency_instrument(data, examiner, d, x = NULL, group = NULL, name = "leniency")
ssiv_rotemberg(data, y, d, shares, shocks, x = NULL)
ssiv_shock_level(data, y, d, shares, shocks, x = NULL, cluster = NULL)
mte_curve(data, y, d, z, x = NULL, method = c("polynomial","local"), degree = 3,
          grid = seq(0.05, 0.95, 0.05), n_boot = 99)
cf_residuals(data, d, z, x = NULL, family = c("gaussian","probit","logit"), name = "v_hat")
cf_bootstrap(data, first, second, n_boot = 499, cluster = NULL, ape = NULL)
cf_ape(fit, data, d, delta = NULL)
est_dml(..., model = "pliv" or "iivm", z = "z", weak_iv = TRUE)   # see double-ml skill
sim_iv(n, dgp = c("linear","weak","late","mte","probit_cf","shift_share","judge"),
       pi = 1, concentration = 10, seed = NULL)
# linear: y,d,z,x1..x5, attr theta = true coefficient; late: + type, tau_true
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_iv(n = 2000, dgp = "linear", seed = 11)
theta <- attr(d, "theta")
fs <- iv_first_stage(d, d = "d", z = "z", x = paste0("x", 1:5))
iv <- fixest::feols(y ~ x1 + x2 + x3 + x4 + x5 | d ~ z, data = d, vcov = "hetero")
ct <- fixest::coeftable(iv)["fit_d", ]
stopifnot(abs(ct["Estimate"] - theta) < 3 * ct["Std. Error"])
ar <- iv_ar_confidence_set(d, y = "y", d = "d", z = "z", x = paste0("x", 1:5))
dw <- sim_iv(n = 600, dgp = "weak", concentration = 2, seed = 2)
arw <- iv_ar_confidence_set(dw, y = "y", d = "d", z = "z", x = paste0("x", 1:5))
fsw <- iv_first_stage(dw, d = "d", z = "z", x = paste0("x", 1:5))
cat("strong-F ok, weak design flagged; iv-analysis self-check passed\n")
```
