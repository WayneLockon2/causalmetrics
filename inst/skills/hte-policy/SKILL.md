---
name: hte-policy
description: Heterogeneous treatment effects and policy learning with causalmetrics. Use after a design (RCT, selection-on-observables, DML) delivers identification and the question becomes "for whom does it work" or "whom should we treat". Covers the doubly robust score, BLP/GATE inference, meta-learners, out-of-sample model selection, validation (calibration, TOC/QINI/AUTOC), and policy trees, budgets, and frontiers.
---

# Heterogeneous effects and policy learning with causalmetrics

## 1. Does this design apply?

Two questions, two objects: "for whom does it work" is the CATE
`tau(x) = E[Y(1) - Y(0) | X = x]`; "whom should we treat" is a policy
`pi(x)` judged by its value. Identification is inherited from the design
that produced the data (randomization or conditional ignorability); this
skill assumes it and organizes everything around one input: the
cross-fitted doubly robust pseudo-outcome from `dr_scores()`, whose
conditional mean is the CATE.

The goal-relaxation ladder (choose the lowest rung that answers the
question; each rung has honest inference, the higher ones do not):

1. Root-n summaries: best linear predictor of tau(X), group ATEs.
2. Flexible CATE models: no confidence intervals, only out-of-sample
   scores and validation.
3. Policies: regret bounds and held-out values, without needing the
   CATE's level.

Checks first: the design skill's overlap checks; then split the sample:

```r
set.seed(1)
idx <- sample(rep(1:3, length.out = nrow(dat)))   # train / select / test
```

Red flags: no credible design behind the data -> back to the design
skills; interest is one subgroup stated ex ante -> the design skill's
subgroup analysis suffices; "heterogeneity" in the outcome rather than
the effect -> descriptive statistics, not this skill.

## 2. Choose the variant

| Goal | Deliverable | Main calls |
|---|---|---|
| Inference on summaries | BLP coefficients, group ATEs with sup-t bands | `cate_blp()`, `cate_gate()` |
| Flexible tau(x) | predictions + out-of-sample comparison | `cate_learner()` (dr/r/x/t/s), `cate_score()`, `cate_ensemble()` |
| Is the model finding anything | calibration, TOC/QINI, AUTOC | `cate_validate()` |
| Whom to treat | tree/linear/budget policy + honest value | `policy_learn()`, `policy_value()` |
| Impact vs need trade-off | frontier across planner preferences | `policy_frontier()` |
| Several arms | per-arm scores and contrasts | `dr_scores(arms = )`, `contrast_scores()` |

Tie-breakers: DR-learner is the default meta-learner; R-learner when
propensities are extreme; T/S as baselines. Depth-2 trees are the default
policy class (implementable, auditable); `engine = "policytree"` for
exact trees when installed.

## 3. Workflow

### Step 0: the score object

Do:
```r
sc <- dr_scores(dat, y = "y", d = "d", x = xvars, folds = 5, seed = 1)
summary(sc$score); sc$ate; sc$diagnostics
# RCT: pass the known propensity for exactness
# sc <- dr_scores(dat, y = "y", d = "d", x = xvars, p_hat = rep(0.5, nrow(dat)))
```
Look: score mean (it is the ATE), score spread, propensity clipping.
Judge: the score mean matches the design skill's ATE; clipping small.
Fail: heavy clipping -> fix overlap upstream; everything downstream
inherits it.

### Step 1: is there heterogeneity at all (root-n rung)

Do:
```r
blp <- cate_blp(sc, formula = ~ x1 + x2, n_boot = 999, seed = 1)
tidy(blp); plot_cate_blp(blp)
gate <- cate_gate(sc, groups = "segment", n_boot = 999, seed = 1)
gate$table; plot_cate_gate(gate)
```
Look: BLP slope coefficients with uniform bands; GATE per group with
bands; the difference between largest and smallest group.
Judge: BLP slopes jointly zero -> report "no detectable heterogeneity
along these covariates" and stop the CATE modelling (the policy rung may
still be useless); bands, not pointwise stars, carry the claim.
Fail: nothing fails; this rung is always reported.

### Step 2: flexible models, trained and compared honestly

Do (train on fold 1, select on fold 2):
```r
library(mlr3); library(mlr3learners)
train <- dat[idx == 1, ]; select <- dat[idx == 2, ]; test <- dat[idx == 3, ]
sct <- dr_scores(train, y = "y", d = "d", x = xvars, seed = 1)
fits <- list(
  dr = cate_learner(scores = sct, x = xvars, method = "dr",
                    learner = lrn("regr.ranger", min.node.size = 20), seed = 1),
  r  = cate_learner(scores = sct, x = xvars, method = "r",
                    learner = lrn("regr.ranger", min.node.size = 20), seed = 1),
  t  = cate_learner(scores = sct, x = xvars, method = "t",
                    learner = lrn("regr.ranger"), seed = 1))
scs <- dr_scores(select, y = "y", d = "d", x = xvars, seed = 2)
tau_sel <- lapply(fits, predict, select)
comp <- cate_score(scs, dr = tau_sel$dr, r = tau_sel$r, t = tau_sel$t)
comp
ens <- cate_ensemble(scs, dr = tau_sel$dr, r = tau_sel$r, t = tau_sel$t, method = "q")
```
Look: DR-loss differences against the constant-ATE baseline, with
intervals; ensemble weights.
Judge: a model beats "constant" out of sample before its predictions are
used anywhere; ties -> Q-aggregation ensemble.
Fail: nothing beats constant -> heterogeneity is not learnable here at
this n; report rung-1 results only.

### Step 3: validation on the held-out third

Do:
```r
scv <- dr_scores(test, y = "y", d = "d", x = xvars, seed = 3)
tau_test <- predict(fits$dr, test)
val <- cate_validate(scv, tau_hat = tau_test, n_groups = 4, n_boot = 999, seed = 1)
val; plot_cate_validation(val)
```
Look: heterogeneity slope (target 1), calibration by quartile of
predicted tau, TOC/QINI curves with bands, AUTOC.
Judge: slope near 1 and monotone quartiles = the model ranks people
correctly; AUTOC's band excluding zero = targeting by this model beats
random.
Fail: slope near 0 with a significant rung-1 BLP -> the flexible model
overfit; fall back to the BLP as the deliverable.

### Step 4: policy learning (only now)

Do:
```r
pol <- policy_learn(sct, x = xvars, method = "tree", depth = 2, holdout = 0.5, seed = 1)
pol$policy; plot_policy_tree(pol)
pv <- policy_value(scv, policy = predict(pol, test), baseline = "none", n_boot = 999, seed = 1)
pv
# budgeted version:
tau_train <- predict(fits$dr, train)
polb <- policy_learn(sct, x = xvars, method = "budget", budget = 0.3,
                     tau_hat = tau_train, seed = 1)
```
Look: the tree; the held-out value against treat-everyone and
treat-no-one; the same with a per-unit `cost`.
Judge: the policy's held-out value beats both blanket policies by more
than its interval; report the value, not the in-sample gain.
Fail: no policy beats blanket treatment -> that is the finding; targeting
costs more than it buys here.

### Step 5: frontiers when the planner weighs need

Do:
```r
fr <- policy_frontier(scv, tau_hat = tau_test, y0_hat = predict_y0_test,
                      budget = 0.3, utility = "crra", curvature = c(0, 1, 3))
plot_policy_frontier(fr)
```
Judge: the frontier shows what impact is given up to reach the neediest;
the choice along it is the planner's, and the report says so.

### Step 6: several arms

Do:
```r
scm <- dr_scores(dat, y = "y", d = "arm", arms = TRUE, x = xvars, seed = 1)
c12 <- contrast_scores(scm, treat = "arm2", control = "arm1")
```
Then rungs 1-4 on each contrast; `policy_learn()` on multi-arm scores
picks the best arm per unit.

## 4. Report

The score construction (nuisances, folds, clipping); the sample split;
rung-1 BLP/GATE with uniform bands; the out-of-sample model comparison
table; validation on the untouched test set (slope, calibration, AUTOC
with band); the policy with its held-out value against both blanket
policies; the scale everything lives on (latent, probability, money) and
seed/split stability.

## 5. Pitfalls

- Never evaluate a learner on rows it trained on: DR/R labels are noisy
  and in-sample forests look brilliant while learning noise; larger leaf
  sizes and fresh samples are the defaults here.
- The `scores` argument partially matches a model named `s`/`sc` in
  `cate_score()`/`cate_ensemble()`; name models `s_learner` etc.
- AUTOC/QINI implementations differ in tie handling; this package matches
  grf exactly (rank-weighted integration), which is the comparison to
  trust.
- Binary outcomes: heterogeneity in the risk difference is not
  heterogeneity in the odds ratio; fix the scale before interpreting.
- Do not distill a policy from in-sample tau_hat; distill from held-out
  predictions or refit on the policy rung directly.
- Report the propensity design: with estimated propensities the bands are
  design-conditional.

## 6. Function reference

```r
dr_scores(data, y, d, x = NULL, p_hat = NULL, mu0_hat = NULL, mu1_hat = NULL,
          arms = NULL, learner_p = NULL, learner_mu = NULL, type = c("dr","ipw","reg"),
          folds = 5, seed = NULL, p_clip = c(0.01, 0.99), trim = NULL)
# $score (pseudo-outcome vector; mean = $ate), $data, $fold_id, $diagnostics
cate_blp(scores, formula = NULL, uniform = TRUE, n_boot = 999, seed = NULL)   # tidy(), glance()
cate_gate(scores, groups, n_boot = 999, seed = NULL)                          # $table + bands
cate_learner(scores = , x = , method = c("dr","r","x","t","s"), learner = ,
             learner_final = NULL, adapt = FALSE, seed = NULL)                # predict(fit, newdata)
cate_score(scores, ..., baseline = "constant")        # out-of-sample DR-loss comparison
cate_ensemble(scores, ..., method = c("q","convex","best","ols"))
cate_validate(scores, tau_hat, n_groups = 4, n_boot = 999, seed = NULL)
# slope, calibration groups, TOC/QINI, AUTOC; plot_cate_validation()
policy_learn(scores, x = NULL, method = c("tree","linear","classifier","budget"),
             depth = 2, cost = 0, budget = NULL, holdout = 0.5,
             engine = c("exhaustive","policytree"), seed = NULL)
policy_value(scores, policy, baseline = "none", cost = 0, n_boot = 0, seed = NULL)
policy_frontier(scores, tau_hat, y0_hat = NULL, budget = 0.3,
                utility = c("cara","crra"), curvature = c(0, 0.5, 1, 2, 5))
contrast_scores(x, treat, control); bind_scores(...)
sim_hte(n, dgp = c("smooth","simple_cate","complex_cate","unbalanced","binary_outcome","policy"))
# y, d, x1..x5 (smooth/policy), tau_true, p_true
```

## 7. Self-check

```r
library(causalmetrics)
d <- sim_hte(n = 4000, dgp = "smooth", seed = 21)
train <- d[1:2000, ]; test <- d[2001:4000, ]
sct <- dr_scores(train, y = "y", d = "d", x = paste0("x", 1:5), seed = 1)
stopifnot(abs(mean(sct$score) - mean(train$tau_true)) < 0.2)
blp <- cate_blp(sct, formula = ~ x1 + x2, n_boot = 199, seed = 1)
tb <- tidy(blp)
stopifnot(tb$estimate[tb$term == "x1"] > 2 * tb$std.error[tb$term == "x1"])  # true slope on x1
fit <- cate_learner(scores = sct, x = paste0("x", 1:5), method = "dr", seed = 1)
tau_hat <- predict(fit, test)
scv <- dr_scores(test, y = "y", d = "d", x = paste0("x", 1:5), seed = 2)
val <- cate_validate(scv, tau_hat = tau_hat, n_boot = 199, seed = 1)
stopifnot(cor(tau_hat, test$tau_true) > 0.3)
pol <- policy_learn(sct, x = paste0("x", 1:5), method = "tree", depth = 2, seed = 1)
pv <- policy_value(scv, policy = predict(pol, test))
stopifnot(is.finite(pv$value$estimate[1]))
cat("hte-policy self-check passed\n")
```
