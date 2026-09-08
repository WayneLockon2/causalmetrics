---
name: rct-analysis
description: Randomized experiment analysis. Use when treatment assignment was randomized (simple, stratified, clustered, or factorial) and the task is to estimate and report treatment effects with design-based inference, handle noncompliance, attrition, and multiple outcomes. Estimation runs on estimatr/fixest directly; causalmetrics supplies the judgment gates and downstream tools.
---

# Randomized experiment analysis

## 1. Does this design apply?

The question: the average effect of an intervention that was assigned by a
known random mechanism. Four variables must be distinguished before any
regression: assignment `Z`, take-up `D`, outcome `Y`, and whether the
outcome was observed `R`. The estimand ladder is: ITT (effect of `Z`,
always identified), then LATE/TOT (effect of `D` for compliers, needs the
instrument logic) only after the ITT is on the table.

Checks to run first:

```r
mean(dat$Z); table(dat$Z, dat$D, useNA = "ifany")   # assignment rate; compliance
tapply(is.na(dat$Y), dat$Z, mean)                    # missing outcomes by arm
table(dat$stratum, dat$Z)                            # if stratified: rate per stratum
```

- The randomization mechanism must be known (who was eligible, at what
  level, within which strata). If it cannot be stated, this is not an
  RCT analysis; load `selection-on-observables`.
- Clustered assignment (villages, cohorts, classrooms): the cluster count,
  not the individual count, drives inference.

Red flags pointing elsewhere: assignment probability depended on
covariates in unknown ways -> `selection-on-observables`; an encouragement
was randomized but the effect of take-up is the target -> do the ITT here,
then load `iv-analysis` for the LATE machinery; interest is in who
benefits or whom to treat -> finish the ITT, then load `hte-policy`.

## 2. Choose the variant

| Variant | Use when | Estimator |
|---|---|---|
| Difference in means / OLS on Z | simple randomization | `estimatr::lm_robust(Y ~ Z, se_type = "HC2")` |
| Covariate-adjusted (Lin) | pre-treatment covariates predict Y | `lm_robust(Y ~ Z * scale(x, scale = FALSE))` |
| Stratified | assignment within blocks | block fixed effects + `lm_robust`, or block-size-weighted |
| Clustered | assignment at group level | `fixest::feols(Y ~ Z, cluster = ~cluster_id)` |
| Randomization inference | small samples, exact mechanism known | `randomizr::declare_ra()` + `ri2::conduct_ri()` |
| Noncompliance (LATE) | D != Z for some units | `estimatr::iv_robust(Y ~ D | Z)`; diagnostics from `iv-analysis` |
| Multi-arm / dose | several treatments | one dummy per arm vs control; families defined upfront |

Tie-breakers: report the unadjusted ITT first, adjusted second; adjustment
uses pre-treatment covariates only and changes precision, not
identification. Cluster where assignment happened or where shocks are
shared, whichever is coarser.

## 3. Workflow

### Step 0: reconstruct the design

Do: write down N, arms, assignment level, strata, and intended
probabilities; verify against the data with the Step-1 checks.
Judge: realized assignment rates match the design (per stratum, per
cluster).
Fail: rates off in some strata -> assignment probabilities differ by
stratum; weight by inverse assignment probability or include stratum fixed
effects; if the mechanism cannot be reconstructed, treat as observational.

### Step 1: balance as an implementation check

Do:
```r
estimatr::lm_robust(x1 ~ Z, data = dat, se_type = "HC2")   # per covariate
# joint: F-test of Z on all covariates, or the omnibus regression x ~ Z
```
Look: standardized differences and the joint test.
Judge: balance checks assess implementation and precision; they do not
create or destroy identification. One significant covariate in ten is
expected.
Fail: systematic imbalance concentrated in variables tied to the
assignment process -> investigate the mechanism (broken randomization,
differential enrollment) before any effect estimate.

### Step 2: attrition before effects

Do:
```r
estimatr::lm_robust(is.na(Y) ~ Z, data = dat, se_type = "HC2")
tapply(is.na(dat$Y), dat$Z, mean)
```
Judge: missingness rates similar across arms and unrelated to Z.
Fail: differential attrition -> the simple ITT is biased for the original
population. Report Lee bounds (trim the arm with less attrition at the
attrition-rate gap) alongside the point estimate, and describe who is
missing.

### Step 3: the ITT

Do:
```r
fit0 <- estimatr::lm_robust(Y ~ Z, data = dat, se_type = "HC2")
fitx <- estimatr::lm_robust(Y ~ Z + x1 + x2, data = dat, se_type = "HC2")
# clustered assignment:
fitc <- fixest::feols(Y ~ Z, data = dat, cluster = ~cluster_id)
```
Look: point estimate, SE, and how much adjustment moved the SE (it should
shrink) versus the estimate (it should barely move).
Judge: adjusted and unadjusted estimates agree within a fraction of an SE.
Fail: adjustment moves the estimate a lot -> suspect bad controls
(post-treatment variables in `x`) or broken randomization; remove any
variable measured after assignment.

### Step 4: inference matched to the mechanism

Do: HC2 for individual randomization; CR2 or `cluster = ~` at the
assignment level; with < ~40 clusters, wild cluster bootstrap-t; and
randomization inference when the mechanism is exactly known:
```r
decl <- randomizr::declare_ra(N = nrow(dat), m = sum(dat$Z))
ri2::conduct_ri(Y ~ Z, declaration = decl, data = dat, sims = 999)
```
Judge: design-based and model-based p-values agree; when they differ,
report the randomization-inference one.

### Step 5: noncompliance

Do: report the first stage (share of compliers) and the LATE only after
the ITT:
```r
estimatr::lm_robust(D ~ Z, data = dat, se_type = "HC2")   # first stage
estimatr::iv_robust(Y ~ D | Z, data = dat, se_type = "HC2")
```
Judge: LATE = ITT / first stage; state that it is a complier effect. For
complier description, weak-assignment worries, or one-sided compliance
bounds, load `iv-analysis` (`complier_profile()`, `iv_ar_confidence_set()`,
`iv_ate_bounds()`).

### Step 6: multiple outcomes and subgroups

Do: define outcome families and subgroup lists before estimating; adjust
within family (Romano-Wolf or Westfall-Young via `ri2`/`wildwyoung`, or
report sharpened q-values); pre-registered subgroups only, or label the
rest exploratory.
Judge: a starred subgroup outside the pre-registered list is a hypothesis,
not a finding.
Fail: no pre-registration -> report the family-adjusted results and say
none of the subgroup analysis was prespecified (the Wheeler replication
shows the honest phrasing).

### Step 7: heterogeneity and targeting hand-off

Do: build doubly robust scores once the ITT stands, and hand off:
```r
sc <- dr_scores(dat, y = "Y", d = "Z", x = c("x1", "x2"), seed = 1)
```
Then load `hte-policy` for BLP/GATE, validation, and policy learning. In
an RCT the propensity is known; pass it via `p_hat = rep(mean(dat$Z), nrow(dat))`
for exactness.

## 4. Report

Assignment mechanism paragraph (level, strata, probabilities); a balance
table with the joint test; attrition by arm with bounds if differential;
ITT table unadjusted and adjusted with the SE type and cluster count
stated; first stage and LATE when compliance is imperfect; the
multiple-testing families; and the sentence saying which analyses were
prespecified.

## 5. Pitfalls

- Conditioning on any post-treatment variable (take-up, engagement,
  survey completion) in the main specification changes the estimand;
  never do it silently.
- "Controlling for" imbalance found in step 1 is fine for precision, but
  imbalanced variables tied to the mechanism signal a design problem
  adjustment cannot fix.
- Clustered assignment with individual-level SEs overstates precision by
  a factor that grows with cluster size; always cluster at assignment
  level.
- Randomization inference must replicate the actual mechanism (strata,
  clusters, arm sizes); a permutation that ignores strata tests the wrong
  null.
- LATE before ITT invites misreading; the ITT is the policy-relevant
  effect of offering the program.

## 6. Function reference

The estimators are `estimatr::lm_robust`, `estimatr::iv_robust`,
`fixest::feols`, `randomizr::declare_ra`, `ri2::conduct_ri` used directly
(the package deliberately does not wrap them). causalmetrics supplies:

```r
dr_scores(data, y, d, x = NULL, p_hat = NULL, learner_mu = NULL,
          type = c("dr","ipw","reg"), folds = 5, seed = NULL, p_clip = c(0.01, 0.99))
# cross-fitted doubly robust pseudo-outcomes; entry point to hte-policy
iv_ate_bounds(data, y, d, z, n_boot = 199)          # one-sided compliance bounds
complier_profile(data, d, z, covariates, x = NULL)  # who complies (kappa weights)
# table_task(), kable_notes(): reporting helpers
# Note: sim_hte() draws OBSERVATIONAL assignment (p_true depends on x);
# for an RCT self-check, randomize d by hand as below.
```

## 7. Self-check

```r
library(causalmetrics)
set.seed(7)
n <- 4000
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
d <- rbinom(n, 1, 0.5)                       # truly randomized
tau <- 1
y <- tau * d + x1 + 0.5 * x2^2 + rnorm(n)
dat <- data.frame(y, d, x1, x2, x3)
fit0 <- estimatr::lm_robust(y ~ d, data = dat, se_type = "HC2")
fitx <- estimatr::lm_robust(y ~ d + x1 + x2, data = dat, se_type = "HC2")
stopifnot(abs(fit0$coefficients["d"] - tau) < 3 * fit0$std.error["d"])
stopifnot(fitx$std.error["d"] < fit0$std.error["d"])          # adjustment buys precision
bal <- estimatr::lm_robust(x1 ~ d, data = dat, se_type = "HC2")
stopifnot(bal$p.value["d"] > 0.001)                            # randomization balanced
sc <- dr_scores(dat, y = "y", d = "d", x = c("x1", "x2", "x3"),
                p_hat = rep(0.5, n), seed = 1)
stopifnot(abs(mean(sc$score) - tau) < 0.15)
cat("rct-analysis self-check passed\n")
```
