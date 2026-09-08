---
name: did-analysis
description: Difference-in-differences with causalmetrics. Use when the outcome is observed for many units over time, some units switch into treatment at known dates and stay treated, and untreated or not-yet-treated units exist in every treated period. Covers 2x2, staggered Callaway-Sant'Anna, imputation, Sun-Abraham, synthetic DiD, and continuous doses, with diagnostics and sensitivity.
---

# Difference-in-differences analysis with causalmetrics

## 1. Does this design apply?

The question: what is the effect of an absorbing treatment on treated units,
`ATT(g, t)` for adoption cohort `g` at time `t`, and its aggregations (event
study, overall ATT, per-cohort ATT). Identification comes from parallel
trends: without treatment, treated and comparison outcomes would have moved
in parallel, possibly only conditional on covariates.

Data requirements, and the checks to run before anything else:

```r
table(table(dat$id))                      # panel balance: rows per unit
with(dat[!duplicated(dat$id), ], table(g))  # adoption cohorts; g = 0 or NA or Inf for never treated
with(dat, tapply(d, list(id_treated_once = ave(d, id, FUN = max)), mean))
any(with(dat, ave(d, id, FUN = function(v) any(diff(v) < 0))))  # reversals must be FALSE
```

- One row per unit and period (or repeated cross-sections; then
  `sampling = "rcs"`).
- A group column: the first treated period per unit, with never-treated
  coded `0`, `NA`, or `Inf` (the package convention).
- Untreated comparisons at every treated period; cohorts treated in the
  first observed period have no pre-period and must be dropped, not coded
  around.
- Cohorts with fewer than ~30 units make cell estimates noisy; flag them
  and consider merging or dropping (report the choice).

Red flags that point to another skill:

- Treatment switches on and off: the Callaway-Sant'Anna machinery here
  assumes absorbing treatment. De Chaisemartin-d'Haultfoeuille tools
  (`DIDmultiplegtDYN`, needs `library(polars)`) handle reversals; treat
  that as outside this skill's guarantees.
- One or two treated aggregate units with a long pre-period: load
  `synthetic-control`.
- Assignment by a cutoff in a score, everyone at once: load `rd-analysis`.
- No timing variation at all (one cross-section): load
  `selection-on-observables`.

## 2. Choose the variant

| Variant | Use when | Estimand | Main call |
|---|---|---|---|
| 2x2 / single-cohort event study | one adoption date | ATT, dynamic ATT | `att_gt()` (one group) or `fixest::feols` + `i()` |
| Callaway-Sant'Anna (default) | staggered adoption; covariates may be needed for parallel trends | ATT(g,t) + aggregations | `att_gt()` + `aggregate_att()` |
| Imputation (Borusyak-Jaravel-Spiess) | parallel trends trusted in all pre-periods; efficiency wanted | dynamic ATT | `did_imputation()` |
| Sun-Abraham | staggered, no covariates, fixest workflow | dynamic ATT | `fixest::sunab()` |
| Efficient (Roth-Sant'Anna) | treatment timing as-good-as random | simple/cohort/event | `staggered_efficient()` |
| Synthetic DiD | few treated units, or parallel trends shaky but pre-fit possible | ATT | `sdid_weights()` + `sdid_se()` |
| Continuous dose | treatment intensity varies across units | ATT by dose | `att_dose()` |

Tie-breakers: default to Callaway-Sant'Anna with `method = "dr"`. Use
`control_group = "notyet"` unless later adopters plausibly anticipate (then
`"never"`, and say why). TWFE with a single post dummy is a diagnostic
(step 2 below), never the headline estimate. If covariates are needed for
parallel trends, they must be time-invariant or pre-treatment values.

## 3. Workflow

### Step 0: audit and recode

Do:
```r
dat$g <- ifelse(is.na(dat$first_treat) | dat$first_treat > max(dat$time), 0, dat$first_treat)
dat <- dat[!(dat$g == min(dat$time) & dat$g != 0), ]  # drop no-pre-period cohorts
hist(dat$g[dat$g > 0 & !duplicated(dat$id)])
```
Look: cohort sizes, share never-treated, calendar span after each cohort.
Judge: every kept cohort has at least one pre-period and at least one
comparison cohort per post period.
Fail: merge adjacent tiny cohorts or drop them; report the sample change.

### Step 1: raw trajectories by cohort

Do:
```r
library(ggplot2)
agg <- aggregate(y ~ time + g, dat, mean)
ggplot(agg, aes(time, y, colour = factor(g))) + geom_line() +
  geom_vline(xintercept = sort(unique(dat$g[dat$g > 0])), linetype = "dashed")
```
Look: do cohorts ride rising trends into adoption? Do never-treated drift
away from everyone?
Judge: roughly parallel pre-adoption movement, cohort by cohort.
Fail: rising-into-adoption lines mean selection on trajectories. Plan
covariates in both nuisances, consider a matched design on recent outcomes
(MatchIt as engine, then re-estimate), and expect the honest bands of step
4 to carry the conclusion.

### Step 2: TWFE diagnostics (not estimates)

Do:
```r
bd <- bacon_decomp(dat, id = "id", time = "time", y = "y", d = "d")
bd$by_type
tw <- twfe_weights(dat, id = "id", time = "time", d = "d", y = "y")
c(tw$share_negative, tw$sum_negative)
```
Look: weight by comparison type; the share and sum of negative weights.
Judge: if later-vs-earlier comparisons or negative weights carry real
weight, any TWFE number in the literature you are comparing against is
suspect; either way proceed to step 3.
Fail: nothing fails here; this step only explains gaps you will see later.

### Step 3: main estimate

Do:
```r
cs <- att_gt(dat, id = "id", time = "time", group = "g", y = "y",
             x = c("x1", "x2"), method = "dr", control_group = "notyet",
             n_boot = 999, seed = 1)
dyn <- aggregate_att(cs, type = "dynamic", min_e = -4)
overall <- aggregate_att(cs, type = "simple")
overall$overall
plot_event_study(event_study_frame(CS = dyn))
```
Look: `dyn$by` (estimate, pointwise `conf.*`, uniform `band.*`),
`overall$overall`, `cs$pretest` (Wald test that all pre-treatment cells are
zero), propensity trimming messages.
Judge: pre-period cells jointly insignificant (pretest p > 0.10 is comfort,
not proof); post effects read off the uniform band, not pointwise
intervals; trimmed share small.
Fail: pretest rejects -> do not stop; go to step 4 with the violation
quantified, add covariates that plausibly restore parallel trends, or
switch comparison group; if a specific cohort drives it, inspect
`cs$att_gt` cell by cell.

### Step 4: pre-trend honesty (always, not only on failure)

Do:
```r
pp <- pretrend_power(dyn, seed = 1)
pp$detectable                         # slope detectable at 50/80% power
csu <- att_gt(dat, id = "id", time = "time", group = "g", y = "y",
              x = c("x1", "x2"), method = "dr", control_group = "notyet",
              base_period = "universal", n_boot = 0)
dynu <- aggregate_att(csu, type = "dynamic", min_e = -4, max_e = 4, n_boot = 0)
hi <- honestdid_inputs(dynu)
l1 <- HonestDiD::basisVector(1, hi$numPostPeriods)
orig <- HonestDiD::constructOriginalCS(betahat = hi$betahat, sigma = hi$sigma,
  numPrePeriods = hi$numPrePeriods, numPostPeriods = hi$numPostPeriods, l_vec = l1)
rm <- HonestDiD::createSensitivityResults_relativeMagnitudes(
  betahat = hi$betahat, sigma = hi$sigma, numPrePeriods = hi$numPrePeriods,
  numPostPeriods = hi$numPostPeriods, Mbarvec = seq(0.5, 2, 0.25), l_vec = l1,
  grid.lb = orig$lb - 2 * abs(orig$lb), grid.ub = orig$ub + 2 * abs(orig$ub))
```
Look: the detectable slope against the effect size; the smallest Mbar at
which zero enters the interval (breakdown value).
Judge: breakdown Mbar >= 1 means the effect survives violations as large
as the worst pre-period deviation; report the value either way. A tiny
detectable-slope bias relative to the effect strengthens a passing
pretest; a large one weakens it.
Fail: breakdown < 0.5 -> the conclusion rests on near-perfect parallel
trends; say so in the report and lean on design-based arguments, matching,
or a different comparison group. The universal base period is required
here; the default varying base gives wrong HonestDiD inputs.

### Step 5: composition and anticipation probes

Do:
```r
aggregate_att(cs, type = "group")$by       # per-cohort ATT
aggregate_att(cs, type = "calendar")$by
aggregate_att(cs, type = "dynamic", min_e = -4, balance_e = 4)$by  # fixed cohort mix
ant <- att_gt(dat, id = "id", time = "time", group = "g", y = "y",
              x = c("x1", "x2"), method = "dr", control_group = "notyet",
              anticipation = 1, n_boot = 499, seed = 1)
```
Look: do late cohorts (the ones with suspicious trajectories) carry the
effect? Does the balanced event study still attenuate/persist? Does the
impact effect move when the pre-adoption period is treated as affected?
Judge: stability across these cuts. A large anticipation shift means the
dating of treatment is itself a finding; report both datings.
Fail: composition explains the dynamic shape -> report the balanced
aggregation as the headline dynamic figure.

### Step 6: estimator triangulation

Do:
```r
imp <- did_imputation(dat, id = "id", time = "time", group = "g", y = "y", pre_window = 4)
dat$g_sa <- ifelse(dat$g == 0, 10000L, dat$g)
sa <- fixest::feols(y ~ sunab(g_sa, time) | id + time, data = dat, cluster = ~id)
fr <- event_study_frame(CS = dyn, Imputation = imp, SunAbraham = sa)
plot_event_study(fr)
```
Look: overlap of the three curves; imputation's pre-period residual
estimates are an extra pre-trend check.
Judge: agreement is a robustness statement, not proof; disagreement
localizes which assumption differs (imputation and Sun-Abraham lean on
pre-period trends and no covariates respectively).
Fail: none; report the comparison figure.

### Step 7: inference checks

Do: cluster at the treatment-assignment level (usually the unit; higher if
shocks are shared). With < ~40 clusters, add a wild cluster bootstrap-t
(write it in the note: refit under sign-flipped cluster residuals; see the
Wheeler replication for a 30-cluster implementation). When timing is
plausibly exchangeable across units:
```r
did_permutation_test(dat, id = "id", time = "time", group = "g", y = "y", n_perm = 999, seed = 1)
```
Judge: bootstrap/permutation p-values in the same region as analytic ones.
Fail: report the conservative one.

### Step 8: report (see Section 4)

## 4. Report

Minimum deliverable: (i) adoption histogram and cohort table; (ii) the
event-study figure with uniform bands and the estimator comparison
overlay; (iii) overall and per-cohort ATT table with standard errors and
the number of clusters; (iv) the pretest statistic, the detectable
pre-trend slope, and the HonestDiD breakdown value; (v) one paragraph
stating the comparison group, the covariates and why they restore
parallel trends, anticipation handling, and any sample drops or trimming.
Never headline a pooled TWFE coefficient.

## 5. Pitfalls

- Cohorts without a pre-period silently poison everything; drop them first.
- `honestdid_inputs()` requires `base_period = "universal"` in `att_gt()`.
- HonestDiD's relative-magnitudes default grid can miss the estimate; pass
  `grid.lb`/`grid.ub` explicitly.
- `event_study_frame()` reads fixest coefficients only in `name::value`
  form (use `i()` or `sunab()`, not hand-built dummies).
- Never-treated units for `fixest::sunab()` need a cohort value beyond the
  time range (e.g. 10000), not 0.
- Binned-endpoint TWFE event studies are a restriction; the saturated
  version equals universal-base `att_gt()` with one cohort exactly.
- Tiny cohorts can break the propensity step; `did`-style warnings about a
  group mean that cohort's cells are noisy, and `p_trim` may drop it.
- `aggregate_att(..., na.rm = TRUE)` is needed when subsetting makes some
  (g,t) cells empty.

## 6. Function reference

```r
att_gt(data, id, time, group, y, x = NULL, method = c("dr","reg","ipw"),
       control_group = c("notyet","never"), base_period = c("varying","universal"),
       anticipation = 0, sampling = c("panel","rcs"), weights = NULL, cluster = NULL,
       n_boot = 999, conf_level = 0.95, p_trim = 0.995, seed = NULL)
# returns: $att_gt (cell table), $pretest, $groups; feed to aggregate_att()

aggregate_att(x, type = c("dynamic","group","calendar","simple"),
              balance_e = NULL, min_e = -Inf, max_e = Inf, na.rm = FALSE)
# returns: $overall (estimate, std.error, conf.*), $by (per event/group/period,
# with band.low/band.high uniform bands), $inffunc

did_imputation(data, id, time, group, y, x = NULL, horizons = NULL, pre_window = 3,
               cluster = NULL, weights = NULL)      # $by_event, $pre
staggered_efficient(data, id, time, group, y,
                    estimand = c("simple","cohort","calendar","eventstudy"),
                    event_time = 0, control = c("notyet","last"))
bacon_decomp(data, id, time, y, d)                  # $decomposition, $by_type, $twfe
twfe_weights(data, id, time, d, y = NULL)           # $share_negative, $sum_negative
event_study_frame(..., ref = -1)                    # named objects -> tidy frame
plot_event_study(frame, band = TRUE)
pretrend_power(estimate, slope = NULL, target_power = c(0.5, 0.8), seed = NULL)
honestdid_inputs(x, ref = -1)                       # betahat, sigma, numPre/PostPeriods
did_permutation_test(data, id, time, group, y, n_perm = 999, seed = NULL)
sdid_weights(data, id, time, y, d, estimator = c("sdid","sc","did","difp"))
# block adoption only; run cohort by cohort for staggered; sdid_se(x) for SEs
att_dose(data, id, time, y, dose, dose_type = c("binned","spline"), n_bins = 4)
sim_did_panel(n_units = 500, n_periods = 10, groups = c(4, 7), never_share = 0.4, seed = NULL)
# columns: id, time, g, x1, x2, event_time, treated, tau (true effect), y
```

## 7. Self-check

Run before touching real data; every assertion should pass.

```r
library(causalmetrics)
d <- sim_did_panel(n_units = 400, n_periods = 10, groups = c(4, 7), seed = 42)
truth <- aggregate(tau ~ event_time, d[d$g > 0 & d$event_time >= 0, ], mean)
cs <- att_gt(d, id = "id", time = "time", group = "g", y = "y",
             x = c("x1", "x2"), method = "dr", control_group = "notyet",
             n_boot = 199, seed = 1)
dyn <- aggregate_att(cs, type = "dynamic", min_e = -3, seed = 1)
e0 <- dyn$by[dyn$by$event_time == 0, ]
stopifnot(abs(e0$estimate - truth$tau[truth$event_time == 0]) < 3 * e0$std.error)
pre <- dyn$by[dyn$by$event_time < 0, ]
stopifnot(all(abs(pre$estimate) < 4 * pre$std.error))
ov <- aggregate_att(cs, type = "simple", seed = 1)$overall
stopifnot(ov$estimate > 0, is.finite(ov$std.error))
bd <- bacon_decomp(d, id = "id", time = "time", y = "y", d = "treated")
stopifnot(abs(bd$check - bd$twfe) < 1e-6)
cat("did-analysis self-check passed\n")
```
