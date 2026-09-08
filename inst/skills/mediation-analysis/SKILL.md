---
name: mediation-analysis
description: Mediation and mechanisms with causalmetrics. Use when the question is how a treatment produces its effect - how much runs through a candidate mediator, what the direct effect is, or how a coefficient change decomposes. Covers regression and DML mediation under sequential ignorability, controlled direct effects with treatment-induced confounders, instrumented mediators, front-door adjustment, Gelbach decompositions, and mandatory sensitivity analysis.
---

# Mediation and mechanisms with causalmetrics

## 1. Does this design apply?

The question: not whether the treatment works, but how. With `M(d)` the
mediator under treatment `d` and `Y(d, m)` the outcome, the estimands:

- Total effect `TE = E[Y(1, M(1)) - Y(0, M(0))]` (the design skill's job).
- Controlled direct effect `CDE(m) = E[Y(1, m) - Y(0, m)]`.
- Natural direct/indirect `NDE = E[Y(1, M(0)) - Y(0, M(0))]`,
  `NIE = E[Y(1, M(1)) - Y(1, M(0))]`; `TE = NDE + NIE`.

The hard fact this skill enforces: randomizing D never identifies NDE/NIE
by itself, because `Y(1, M(0))` is cross-world. Identification needs
sequential ignorability (the mediator as-good-as-random given D and X - a
new, untestable assumption) or a design that moves the mediator
(instrument, front-door structure). Adding the mediator to a regression
and reading the D coefficient as "direct effect" is the central error:
conditioning on M opens the collider path through mediator-outcome
confounders.

Checks first: the total effect from the appropriate design skill; then

```r
summary(fixest::feols(reformulate(c("d", xvars), "m"), data = dat, vcov = "hetero"))  # D moves M?
```
No treatment effect on the mediator means no indirect effect to
decompose; stop here.

Red flags: the "mediator" is measured before treatment -> it is a
covariate or moderator, use the design skill; interest is effect
heterogeneity by M -> `hte-policy`; the mediator question is really an
IV question about M's own effect -> `iv-analysis`.

## 2. Choose the variant (the identification gate)

| Situation | Identified object | Function |
|---|---|---|
| Mediator as-good-as-random given D, X; parametric models fine | NDE, NIE, CDE, shares | `mediate_reg()` |
| Same assumption, flexible nuisances, root-n inference | NDE, NIE (+ CDE, binary M) | `mediate_dml()` |
| A post-treatment confounder of M and Y exists | CDE only (NDE/NIE not identified) | `mediate_cde()` |
| An instrument shifts the mediator | direct + indirect via IV system | `mediate_iv()` |
| D-Y confounding unobserved, M carries the whole effect and is clean | total effect of D | `front_door()` |
| Descriptive: how covariates absorb a coefficient | order-invariant accounting, not NDE/NIE | `gelbach_decomp()` |

Tie-breakers: state which row you are on before any estimate; the row is
the assumption. Several mediators that do not cause each other: pass a
vector to `mediate_reg()` (per-mediator NIEs sum to the joint one). If
one mediator causes another, per-mediator effects need path-specific
machinery outside this skill - report the joint NIE only.

## 3. Workflow

### Step 0: draw the graph, name the confounders

Do: write D -> M -> Y with the four confounding positions (D-Y, D-M,
M-Y, and treatment-induced M-Y). Name the candidate unobserved M-Y
confounder U and any treatment-induced L.
Judge: if L exists and matters, only the CDE row of the table is
available; say so up front.

### Step 1: effects on the mediator(s), no mediation assumption needed

Do:
```r
fixest::feols(reformulate(c("d", xvars), "m"), data = dat, vcov = "hetero")
```
Judge: report these first; they are identified by the design alone and
already constrain the story.

### Step 2: main estimate under sequential ignorability

Do:
```r
fit <- mediate_reg(dat, y = "y", d = "d", m = "m", x = xvars,
                   method = "delta", cluster = "cl")          # linear system
fit$effects; fit$shares; fit$potential
fit_int <- mediate_reg(dat, y = "y", d = "d", m = "m", x = xvars,
                       interaction = TRUE, n_sim = 1000, seed = 1)
plot_mediation(fit_int)
```
Look: total/NDE/NIE (+ both decompositions with interaction), the
proportion mediated with its SE, the four potential-outcome means.
Judge: total from the system matches the design skill's total; with
interaction, report both decompositions (pure/total); a proportion
mediated with a huge SE is reported as such, not rounded into a claim.
Fail: binary outcome or mediator -> `outcome = "logit"` /
`mediator = "logit"` (probability-scale effects by g-computation); wrong
functional forms suspected -> step 3.

### Step 3: flexible nuisances

Do:
```r
library(mlr3); library(mlr3learners)
fd <- mediate_dml(dat, y = "y", d = "d", m = "m", x = xvars,
                  learner_y = lrn("regr.ranger"),
                  learner_d = lrn("classif.ranger", predict_type = "prob"),
                  folds = 5, n_rep = 3, seed = 1)
fd$effects; fd$diagnostics
```
Look: the efficient-influence-function estimates; the two propensity
overlaps (`p(D|X)` and `p(D|M,X)`), trimmed share.
Judge: regression and DML agree -> functional form was not the issue;
disagree -> prefer DML with a stated reason and report both.

### Step 4: sensitivity (mandatory whenever natural effects are reported)

Do:
```r
sens <- mediate_sensitivity(fit, rho = seq(-0.8, 0.8, 0.05))
sens$rho_zero
plot_mediate_sensitivity(sens)
```
Look: the NIE as a function of the M-Y error correlation; `rho_zero`,
the correlation that erases it.
Judge: benchmark `rho_zero` against the residual-correlation shifts the
observed covariates produce; a small value goes in the abstract, not a
footnote. This rho analysis applies to the sequential-ignorability rows
(steps 2-3). On the CDE escape hatch there is no rho analogue in the
package; its honesty checks are the naive-vs-CDE contrast, the stated
no-D:M-interaction assumption of the demediation step, and, when doubted,
`interaction = TRUE` in `mediate_cde()`.

### Step 5: the escape hatches when sequential ignorability is untenable

Treatment-induced confounder L (observed):
```r
fc <- mediate_cde(dat, y = "y", d = "d", m = "m", x_pre = xvars,
                  x_post = "l", n_boot = 499, seed = 1)
fc$effects        # cde vs naive_direct; naive_direct is the D coefficient
                  # of the single regression y ~ d + m + x_pre + x_post
```
Instrument Z for the mediator:
```r
fi <- mediate_iv(dat, y = "y", d = "d", m = "m", z = "z", x = xvars,
                 fe = "unit_fe", cluster = "cl", weights = NULL,
                 homogeneity_by = "subgroup", n_boot = 499, seed = 1)
fi$effects; fi$first_stage_F; fi$homogeneity
```
Front-door structure:
```r
ff <- front_door(dat, y = "y", d = "d", m = "m", x = xvars,
                 method = "regression", n_boot = 499, seed = 1)
ff$effects        # front-door total vs the biased backdoor regression
```
Judge: each comes with its own gate - `mediate_cde` assumes no D:M
interaction in the demediation; `mediate_iv` needs instrument validity
given D, X plus a homogeneous mediator effect (the subgroup table is the
check; identity `total_z = direct + indirect` holds by construction);
`front_door` needs M to carry the whole effect and be unconfounded - the
exclusion argument is the report's core.

### Step 6: accounting decompositions, labelled as accounting

Do:
```r
g <- gelbach_decomp(dat, y = "y", d = "d", x_base = xvars, x_add = med_vars,
                    fe = "stratum", cluster = "cl", n_boot = 499, seed = 1)
g$coefficients; g$contributions; plot_gelbach(g)
```
Judge: contributions sum to the coefficient change exactly and are
order-invariant - but they are causal mediation only under the same
sequential ignorability as step 2; otherwise present them as "which
covariates absorb the coefficient", nothing more. Sequential addition of
mediators one at a time is never reported (order-dependent).

## 4. Report

The estimand row chosen and its assumption, stated before numbers; the
effect-on-mediator table; the effects table (both decompositions when
interaction is on) with the proportion mediated and its SE; the
sensitivity value with a benchmark; for escape-hatch rows, their specific
diagnostics (first-stage F and homogeneity table; naive-vs-CDE contrast;
backdoor-vs-front-door contrast); the graph with the assumed-absent
arrows named.

## 5. Pitfalls

- The Y ~ D + M regression's D coefficient is not the direct effect
  unless M is unconfounded given D, X; with an unobserved U on M and Y
  it is biased even when D is randomized.
- Proportion mediated explodes when the total effect is small or the
  channels offset; report the ratio with its SE, or not at all.
- Factor covariates are expanded by the g-computation internally; check
  levels present in both arms.
- The bootstrap drops resamples where a regression cannot be fit (rare
  binary covariates in cluster resamples) and warns; a large dropped
  share means the specification is too rich for the resample size.
- `mediate_iv` with weights carries them through every regression and
  the bootstrap; a weak mediator instrument makes the decomposition
  uninformative - report the F and, when marginal, an AR interval from
  `iv-analysis`.
- Parallel mediators must not cause each other; per-mediator NIEs are
  otherwise not separately meaningful.

## 6. Function reference

```r
mediate_reg(data, y, d, m, x = NULL, interaction = FALSE,
            outcome = c("linear","logit"), mediator = c("linear","logit"),
            method = c("simulation","bootstrap","delta"), m_ref = NULL,
            cluster = NULL, weights = NULL, n_boot = 499, n_sim = 1000, seed = NULL)
# $effects (total, nde, nie, nde_total, nie_pure, cde, nie_<m> per mediator),
# $shares (prop. mediated + Wheeler S1/S2 with suest-style SEs), $potential, $models
mediate_sensitivity(fit, rho = seq(-0.9, 0.9, 0.05))   # $curve, $rho_zero
mediate_dml(data, y, d, m, x, learner_y = NULL, learner_d = NULL, learner_nu = NULL,
            learner_m = NULL, m_ref = NULL, folds = 5, n_rep = 1, trim = 0.01, seed = NULL)
# $effects (total, nde, nie, nde_total, nie_pure, + cde with binary m_ref), $diagnostics, $psi
mediate_cde(data, y, d, m, x_pre = NULL, x_post = NULL, m_ref = 0,
            interaction = FALSE, cluster = NULL, n_boot = 499, seed = NULL)
# $effects: cde (sequential g-estimation), naive_direct (D coefficient of the
# one-step regression y ~ d + m + x_pre + x_post, biased by construction),
# delta (mediator effect from step 1), total (y ~ d + x_pre)
mediate_iv(data, y, d, m, z, x = NULL, fe = NULL, cluster = NULL, weights = NULL,
           homogeneity_by = NULL, n_boot = 499, seed = NULL)
# $effects: total, total_z, direct, indirect, first_stage, mediator_effect,
# naive_*; $first_stage_F; $homogeneity
front_door(data, y, d, m, x = NULL, method = c("regression","aipw"),
           folds = 5, n_boot = 499, seed = NULL)
# $effects: total, backdoor (+ d_on_m, m_on_y for regression); aipw needs binary m
gelbach_decomp(data, y, d, x_base = NULL, x_add, groups = NULL, fe = NULL,
               cluster = NULL, weights = NULL, n_boot = 499, seed = NULL)
# $coefficients (base/full/change), $contributions (Gamma, beta_full, share), $identity_gap
sim_mediation(n, dgp = c("linear","interaction","nonlinear","binary","correlated_errors",
              "post_treatment_confounder","iv_mediator","front_door","parallel"),
              rho = 0.5, seed = NULL)   # attr(d, "truth") holds every estimand
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_mediation(3000, dgp = "linear", seed = 1)
tr <- attr(d, "truth")
fit <- mediate_reg(d, y = "y", d = "d", m = "m", x = c("x1", "x2"), method = "delta")
e <- fit$effects
stopifnot(abs(e$estimate[e$term == "nie"] - tr$nie) < 3 * e$std.error[e$term == "nie"])
stopifnot(abs(e$estimate[e$term == "total"] - tr$total) < 3 * e$std.error[e$term == "total"])
dz <- sim_mediation(3000, dgp = "post_treatment_confounder", seed = 5)
fc <- mediate_cde(dz, y = "y", d = "d", m = "m", x_pre = c("x1", "x2"),
                  x_post = "z", n_boot = 49, seed = 1)
ec <- fc$effects
err_cde <- abs(ec$estimate[ec$term == "cde"] - attr(dz, "truth")$cde)
err_naive <- abs(ec$estimate[ec$term == "naive_direct"] - attr(dz, "truth")$cde)
stopifnot(err_cde < err_naive)
df <- sim_mediation(3000, dgp = "front_door", seed = 7)
ff <- front_door(df, y = "y", d = "d", m = "m", x = c("x1", "x2"), n_boot = 49, seed = 1)
ef <- ff$effects
stopifnot(abs(ef$estimate[ef$term == "total"] - attr(df, "truth")$total) <
          abs(ef$estimate[ef$term == "backdoor"] - attr(df, "truth")$total))
cat("mediation-analysis self-check passed\n")
```
