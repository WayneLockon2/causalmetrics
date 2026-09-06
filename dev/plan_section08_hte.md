# Section 08 plan: Heterogeneous Treatment Effects and Policy Learning

Status: 2026-09-05 (late). Wayne approved the plan as written. Implemented the same day:
all package functions below (R/hte-internal.R, dr-scores.R, sim-hte.R, cate-learner.R,
cate-blp.R, cate-score.R, cate-validate.R, policy.R, tidiers), tests
(tests/testthat/test-hte-policy.R: BLP = lm + HC1; TOC/AUTOC/AUQC = grf's
rank_average_treatment_effect; exact policy trees = policytree; convex stacking = quadprog),
the usage vignette vignettes/usage_hte_policy.Rmd, helpers in
inst/vignette_helpers/08_hte/, and the
lecture inst/lectures/08_hte/lecture_08_heterogeneous_effects_and_policy.Rmd (9 sections +
appendix A.1-A.5 + references). Deviations from the plan: no `cate_forest` (grf unwrapped, as
planned); `policy_learn(method = "budget")` replaces a separate top-q helper; the
meta-learner Monte Carlo uses random-forest oracles and reports n = 500 and n = 2000, and a
second Monte Carlo (linear oracles, misspecified outcome model, correct propensity) isolates
double robustness. Replications remain Wayne's. Original proposal follows.

## Sources read (08_Heterogeneous_Treatment_Effects_and_Policy_Learning/materials + web)

| Source | What it contributes | Verdict |
|---|---|---|
| CausalML v0.1.2 ch. 14 (Statistical inference on HTE) | Why CATE is a function, not a parameter; DR pseudo-outcome `Y(eta)`; best linear predictor (BLP) of the CATE via OLS on the pseudo-outcome (Semenova-Chernozhukov 2021); GATEs as BLP with group dummies; pointwise vs uniform bands; honest/generalized random forests and the R-learner moment inside a forest; 401(k) and "welfare" survey examples | Core. Sections 2-4 and appendix |
| CausalML ch. 15 (Estimation and validation of HTE) | Meta-learners S/T/DR/R/X with the regression-oracle notation, error bounds (DR: product of nuisance errors; R: propensity to the fourth power), R-learner as the overlap-weighted projection, X-learner covariate shift and the domain-adapted X (DAX) and DRX variants; DR-loss scoring with confidence (Theorem 15.2.1), stacking, convex and Q-aggregation ensembles; validation: BLP heterogeneity test, calibration by CATE quartiles, TOC and QINI curves with one-sided uniform bands, AUTOC; policy learning as cost-sensitive classification with regret bound; welfare and Criteo applications | Core. Sections 3, 5, 6, 7 |
| MS&E 228 L16-17 (Syrgkanis) | The "ten goals" ladder from BLP to policy learning; multiplier bootstrap for uniform bands; GRF as a local moment; the model-selection and stacking slides in compact form | Framing device for Section 1; nothing new beyond the book |
| Childers, "Outcome and Effect Heterogeneity" | Distinguishes effects on distributions (QTE, distribution regression, counterfactual densities) from distributions of effects (Frechet-Hoeffding/Makarov bounds); quantile regression and its pitfalls; Breza-Kaur-Shamdasani example | Short section on "heterogeneity that is not the CATE" plus a warning box. No package code for QTE in v1 (quantreg and Counterfactual do it) |
| Facure ch. 18 | CATE as a derivative; interaction regressions `y ~ t*X` as the first CATE model; sensitivity partitions for prices | Motivation for Section 1 and for the continuous-treatment BLP (Section 4.4) |
| Facure ch. 21 | S/T/X learners with LightGBM; the two failure pictures (S shrinks the treatment away; T with an unbalanced design fits a linear model on the small arm) | The two pictures become Monte Carlo motivations in Section 3 |
| Facure ch. 23 | Binary outcomes: a covariate that only moves the baseline creates spurious CATE ranking because of the logistic curvature; continuous treatment: curvature of the response function | Section 8 "Traps" |
| jzhou.org meta-learner post | Q-learning vs A-learning taxonomy; multi-arm S/T/X/R; reference-free R-learner (Zhou-Zhang-Tu 2023); optimal-treatment recommendation from each learner; DTR pointer | Multi-arm extension paragraph in Section 3.6 and the "what a policy is" paragraph in Section 7. Do not implement multi-arm in v1 |
| Wayne's replication drafts | Athey-Wager 2021 (policy trees, GAIN data missing, simulation CSVs present); Haushofer et al. 2025 (targeting on impact vs deprivation with causal forests, CARA/CRRA planner); Yoganarasimhan-Barzegary-Pani 2023 (free-trial length as a multi-arm treatment, Python/R lasso, XGBoost, causal tree, causal forest, CART); Zhan-Ren-Athey-Zhou 2024 (policy learning with adaptively collected data, Python, AIPW with adaptive weights, policy trees) | Section 9 anchors. See "Replications" |

Installed locally: grf 2.6.1, DoubleML 1.0.2, mlr3 1.6.0, mlr3learners, ranger, glmnet, rpart, quantreg.
Not installed: policytree, rlearner, causalToolbox, xgboost, lightgbm, partykit, evalITR, SuperLearner.
Recommend installing `policytree` (CRAN) as the exact-tree reference for `policy_tree`-type
comparisons, and keeping grf as the causal-forest engine (Suggests).

## Scope decision (proposed, same spirit as Section 07)

The package does not re-implement random forests, boosting, or tree search. Learners come from
mlr3 (as in `est_dml()`), causal forests come from grf, exact policy trees from policytree.
The package supplies the causal layer those tools do not give in one place:

1. the cross-fitted DR pseudo-outcome (one object reused everywhere, built on `est_dml()` internals),
2. the meta-learner wrappers that turn "any regression oracle" into a CATE model with the right
   cross-fitting discipline (S, T, X, DR, R, plus the DRX/DAX variants),
3. inference on low-dimensional summaries (BLP, GATE, uniform bands via the existing multiplier
   bootstrap),
4. out-of-sample scoring, ensembling, and validation (DR score with CIs, Q-aggregation, BLP
   heterogeneity test, calibration, TOC/QINI/AUTOC with uniform bands),
5. policy evaluation and learning from the same scores (value with CI, EWM as weighted
   classification, depth-limited trees by exhaustive search for small depth, or policytree if
   installed), plus the welfare-curvature planner of Haushofer et al.

Everything sits on a single "scores" object so a lecture reader learns one idea (the DR
pseudo-outcome) and then sees it reused eight times. That is the pedagogical thread.

## Package functions

Names follow the existing verbs (`est_*` for estimators returning one number,
otherwise noun-ish). All accept `data` + column names like `att_gt()`; learners are mlr3
objects like `est_dml()`; every function that needs nuisances also accepts pre-computed
out-of-fold vectors (`p_hat`, `mu0_hat`, `mu1_hat`, `l_hat`, `m_hat`) so Python or grf nuisances
can be plugged in, as in `usage_dml_external_nuisances`.

### A. Scores (the foundation)

`dr_scores(data, y, d, x, p_hat = NULL, mu0_hat = NULL, mu1_hat = NULL, learner_p = NULL,
learner_mu = NULL, folds = 5L, fold_id = NULL, seed = NULL, p_clip = c(0.01, 0.99),
type = c("dr", "ipw", "reg"))`
-> class `cm_scores`: `score` (n vector, DR pseudo-outcome), `nuisance` (p, mu0, mu1),
`residuals` (`y - l(x)`, `d - m(x)` for the R-learner), `fold_id`, `d`, `y`, `x`, diagnostics.
Reuses `.cm_crossfit_predict()` and `.cm_score_irm()`. Multi-arm: not in v1 (document the
extension).

### B. CATE learners

`cate_learner(data, y, d, x, x_het = x, method = c("dr", "r", "x", "t", "s"),
learner = NULL, learner_final = NULL, scores = NULL, folds = 5L, seed = NULL,
final_holdout = 0)`
-> class `cm_cate`: `predict()` method on new data, `tau_hat` (in-sample fitted CATE),
`method`, `scores`, `learners`, the final-stage model.
- `dr`: regress the DR score on `x_het` with `learner_final` (this is Kennedy's DR-learner).
- `r`: weighted regression of `y_tilde / d_tilde` on `x_het` with weights `d_tilde^2`
  (Nie-Wager); `learner_final` must accept weights (ranger, glmnet, lm do).
- `x`: T-stage on each arm, imputed effects, propensity-weighted blend, then final regression;
  option `adapt = TRUE` gives the DAX weights `(1 - p)^2 / p` on the treated stage (CausalML 15.1.3).
- `t`, `s`: for the lecture's failure demonstrations. `s` optionally adds `d * x` interactions.
- All first stages cross-fitted; final stage on all rows (Remark 15.1.1).
`predict.cm_cate(object, newdata)`; `tidy()` returns per-row CATE.

`cate_forest(data, y, d, x, ..., engine = "grf")` is NOT proposed; users call
`grf::causal_forest()` directly and feed `predict()` output to the evaluation functions below
(registry rule: do not wrap). The usage vignette shows the two calls side by side.

### C. Inference on summaries

`cate_blp(scores, formula_or_x, vcov = "HC1", conf_level = 0.95, uniform = TRUE, n_boot = 999L)`
-> OLS of the DR score on a design (`~ income + I(income^2)`, or group dummies for GATEs);
pointwise and sup-t uniform bands using `.cm_multiplier_bootstrap()`; `predict()` on a grid;
plot method. GATE with `~ 0 + group` recovers the Section 14.2 boxes. Also accepts `cm_cate`
predictions as regressors (this is the heterogeneity test of Section 15.3 when the design is
`~ tau_hat_centered`).

`cate_gate(scores, groups, ...)` = thin convenience over `cate_blp` with group dummies, returns
the "GATE table" and a joint test of equality.

### D. Scoring, selection, ensembling

`cate_score(scores_test, ..., baseline = "constant")` where `...` are named `cm_cate` objects or
numeric prediction vectors -> DR-loss per model, difference to the constant-ATE model with the
CI of Theorem 15.2.1, normalized score `S(tau)` (15.2.16).

`cate_ensemble(scores_test, ..., method = c("q", "convex", "best", "ols"), lambda = NULL)`
-> weights on the simplex (Q-aggregation solved with a small projected-gradient or `quadprog`
for the convex case), intercept fixed at the training ATE (Remark 15.2.3), returns a `cm_cate`
whose `predict` combines the base models.

### E. Validation

`cate_validate(scores_test, tau_hat, n_groups = 4L, quantile_grid = seq(0.05, 1, by = 0.05),
n_boot = 999L, conf_level = 0.95)` -> class `cm_cate_val` with
- `blp`: the heterogeneity test (coef on centered CATE, ideally 1),
- `calibration`: DR GATE by CATE quantile group vs mean predicted CATE, CAL1/CAL2 scores,
  joint test that the groups differ,
- `toc`, `qini`: curves on the grid with pointwise and one-sided uniform bands (multiplier
  bootstrap on the influence functions of Theorem 15.3.1, including the estimated-ATE and
  estimated-share nuisance terms),
- `autoc`, `auqc`: point estimate, SE, one-sided CI,
- `group_diff`: covariate means top vs bottom group (Figure 15.16 table).
`plot(x, what = c("calibration", "toc", "qini"))`.
Tie-breaking rule of Remark 15.3.1 implemented (matters for tree-based CATE models).

### F. Policy evaluation and learning

`policy_value(scores, policy, baseline = c("none", "all", "control"), conf_level = 0.95)`
-> value of a fixed rule (`policy` = 0/1 vector, a function of the data, or a `cm_policy`)
with IF-based SE; difference between two policies with SE; cost argument (`cost` per treated
unit, or `budget` share for the constrained problem).

`policy_learn(scores, x, method = c("tree", "linear", "classifier"), depth = 2L,
budget = NULL, cost = 0, learner = NULL, honest_split = 0.5, seed = NULL)`
-> class `cm_policy`: `predict()`, `rule` (tree splits or coefficients), in-sample and
honest held-out value with SE, the weighted-classification reformulation (labels
`sign(score - cost)`, weights `|score - cost|`).
- `tree`: exact depth-1 and depth-2 search over all split pairs implemented in R on sorted
  covariates (O(p^2 n log n) for depth 2 with n <= ~1e4 features-by-thresholds; documented
  limits); if `policytree` is installed, `engine = "policytree"` delegates and we test equality.
- `linear`: EWM over `sign(x'b)` by a smoothed surrogate (logistic weighted classification),
  as Kitagawa-Tetenov's linear rules.
- `classifier`: any mlr3 classifier with weights (rpart, ranger) for the cost-sensitive view.
- `budget`: treat the top-q by CATE with the estimated quantile; returns the TOC-style value.

`policy_frontier(scores, tau_hat, y0_hat, utility = c("cara", "crra"), curvature, budget)`
-> the Haushofer et al. planner: rank units by marginal utility gain
`u(y0 + tau) - u(y0)` for a grid of curvature values, report overlap of the chosen set with
impact-only and deprivation-only targeting (their Figures 2-3, Tables 1-2). Small, but it is
exactly what the replication needs and no package offers it.

### G. Simulation and helpers

`sim_hte(n, dgp = c("simple_cate", "complex_baseline", "unbalanced", "binary_outcome",
"policy"), ...)` reproducing CausalML DGPs 1-3, Facure's binary-outcome trap, and a
policy DGP with a known optimal depth-2 tree, with the true CATE attached.

Plot helpers: `plot_cate_blp()`, `plot_cate_validation()`, `plot_policy_tree()` (text or
ggplot rendering of a depth-2 rule).

### Reuse map

| New function | Reuses |
|---|---|
| `dr_scores` | `.cm_crossfit_predict`, `.cm_score_irm`, `.cm_make_fold_sets`, `.cm_common_support`, `.cm_fit_quality` |
| `cate_learner` | `dr_scores`, `.cm_crossfit_predict` (arm subsets via `subset =`), mlr3 weights |
| `cate_blp`, `cate_validate` | `.cm_multiplier_bootstrap`, `.cm_if_analytic_se`, `.cm_if_wald` |
| `policy_value` | `.cm_if_analytic_se` |
| `policy_learn` (tree) | data.table sorted scans |

Tests: DR-learner vs a hand-rolled pipeline; R-learner vs `rlearner` formulas (not installed;
test against an lm-based R-learner closed form with known nuisances); BLP vs `lm()` + HC1;
policy tree vs `policytree` when available and vs brute force on tiny data; TOC/AUTOC vs
`grf::rank_average_treatment_effect()` (installed) for the point estimates; `policy_value` vs
the difference-in-means in an RCT with known propensity.

Performance: all score-level operations are vectorized; the multiplier bootstrap reuses the
chunked `crossprod` from DiD; depth-2 tree search uses cumulative sums over sorted thresholds,
no per-split refits.

## Teach with existing tools, no wrapper

- `grf::causal_forest()` (+ `variable_importance`, `best_linear_projection`,
  `rank_average_treatment_effect`) in the vignette and lecture; our functions accept its
  predictions and its `get_scores()` output.
- `DoubleML` is not needed (our `est_dml` covers it) but a footnote maps names.
- `quantreg::rq()` and distribution regression for QTE in the "not the CATE" section.

## Skip in v1 (mention only)

Multi-arm learners and reference-free R-learner; CFR-Net / TARNet; BART; conformal CATE
intervals; adaptive-data policy learning (Zhan et al.) beyond a paragraph; dynamic treatment
regimes; sensitivity of CATE to unobservables; continuous-treatment policy learning.

## Lecture skeleton: `inst/lectures/08_hte/lecture_08_heterogeneous_effects_and_policy.Rmd`

Style as lecture 07: short sentences, intuition -> simulation showing the failure -> the fix ->
code, proof skeletons in text, formal statements in the appendix, booktabs tables, footnotesize
code, every number computed in the document, Monte Carlos cached in
`inst/vignette_helpers/08_hte/`.

### 1. From "does it work" to "for whom" and "whom to treat"
1.1 The three questions: ATE, CATE tau(x) = E[Y(1) - Y(0) | X = x], policy pi(x) in {0,1}
1.2 Why CATE is different: a function, not a parameter; no n^-1/2 rate; no adaptive bands
    (Genovese-Wasserman); the ladder of goals (BLP -> tests -> forests -> MSE -> policies)
1.3 Identification: CATE = E[g(1,Z) - g(0,Z) | X] under conditional exogeneity; X can be a
    subset of Z (CATE on X vs the high-dimensional CATE delta(Z)); CATT = CATC = CATE on Z
1.4 The first CATE model: interaction regression `y ~ d * x` (Facure); when it is enough and
    what it hides. Code: an RCT with a known tau(x), two OLS fits, a sensitivity partition.
1.5 A roadmap and one idea that carries the whole lecture: the DR pseudo-outcome

### 2. One signal, many uses: the doubly robust pseudo-outcome
2.1 Y(eta) = H(mu)(Y - g(D,Z)) + g(1,Z) - g(0,Z), E[Y(eta0) | X] = tau(X); conditional
    Neyman orthogonality (skeleton; appendix A.1)
2.2 Cross-fitting the nuisances: why the final stage may use all rows (Remark 15.1.1)
2.3 Code: `dr_scores()`; showing that the score's conditional mean is the CATE on simulated
    data; the two-line bias formula bias(X; eta) = (H(mu0) - H(mu))(g0 - g)
2.4 Connection to lecture 04: the ATE is the mean of the score; everything below is a
    regression of the score on something

### 3. Meta-learners: turning a regression oracle into a CATE model
3.1 Oracle notation O_H({X, Y, W}); the meta-learning idea
3.2 S- and T-learners. Motivation experiment: S with lasso/boosting shrinks the treatment
    away (Facure/Chernozhukov picture); T with 5% treated fits a straight line on the treated
    arm (Kunzel picture). Code: `cate_learner(method = "s"/"t")`, MSE table
3.3 X-learner: CATT and CATC identification, imputed effects, propensity blend; when it
    wins (DGP 1) and when it loses (DGP 3); covariate shift and the DAX weights
3.4 DR-learner: regress Y(eta) on X; error bound r_n^2 + Err(g)^2 Err(H)^2 (skeleton;
    appendix A.2); variance from dividing by the propensity
3.5 R-learner: Robinson decomposition, weighted square loss with weights d_tilde^2,
    the overlap-weighted projection (15.1.8) proved in three lines; error bound with Err(mu)^4;
    stability under near-deterministic assignment
3.6 Guidelines table (S/T/X/DR/R rows; columns: needs outcome model, needs propensity,
    converges to, unstable when) and a Monte Carlo reproducing CausalML Figure 15.10
    (DGPs 1-3, 100 reps, cached). Multi-arm and reference-free R in one paragraph
3.7 Interpreting a black-box CATE: distillation tree on tau_hat, group differences

### 4. Inference on summaries of the CATE
4.1 Best linear predictor of the CATE: OLS of the score on p(X); beta_hat has the same
    asymptotics as with the true score (Semenova-Chernozhukov); HC sandwich
4.2 GATEs as BLP with group dummies; pointwise vs joint CIs; the "boxes" figure
4.3 Uniform bands for a curve x -> p(x)'beta: sup-t via multiplier bootstrap (reuse Section 07
    machinery); code: `cate_blp()` with a polynomial in income on the 401(k)-style simulation
4.4 Continuous treatment: BLP of the derivative E[y'(t) | X] via the partially linear
    interaction model (Facure ch. 18 with orthogonalization); code with `est_dml` residuals
4.5 Forests with confidence intervals: honesty, subsampling, balanced splits; GRF as a
    local moment; causal forest = R-learner moment inside a forest; DR forest = regression
    forest on the score; what the CI does and does not guarantee. Code: `grf::causal_forest`
    vs `cate_learner(method = "dr", learner_final = ranger)` on the same data

### 5. Choosing a CATE model out of sample
5.1 Why prediction metrics fail: the target is unobserved (Facure's "why prediction metrics
    are dangerous" point in two sentences)
5.2 The DR loss L(tau; eta) = E_n (Y(eta) - tau(X))^2 and why differences of losses are
    orthogonal and root-n normal (Theorem 15.2.1, skeleton; appendix A.3)
5.3 Comparing with confidence: the normalized score S(tau) and the "beat the constant model"
    CI; code: `cate_score()` on the meta-learners of Section 3
5.4 Ensembles: best-of, convex stacking, Q-aggregation and the log(M)/n rate; code:
    `cate_ensemble()`; a table showing the ensemble is never far from the best base learner
5.5 Stability: seeds and subsamples; what cannot be given (CIs on the ensemble)

### 6. Validating the chosen model
6.1 The three-way split (train / score / test) and what each set may see
6.2 Heterogeneity test: OLS of the score on (1, tau*(X) - mean); the coefficient is
    Cov(Y(1) - Y(0), tau*) / Var(tau*); ideal value 1
6.3 Calibration: DR GATE by CATE quartile vs mean prediction; CAL1, CAL2; the
    calibration-distortion decomposition (three-line proof)
6.4 Targeting curves: TOC(q) = GATE(q) - ATE, QINI(q) = TOC(q) P(top q), their covariance
    forms (appendix A.4), AUTOC and AUQC; one-sided uniform bands and the "heterogeneity
    statistic"; ties; the RATE view of Yadlowsky et al. (2021) as a remark
6.5 Code: `cate_validate()` end to end on the ensemble; a full validation figure panel

### 7. Policy learning
7.1 A policy is a map X -> {0,1}; value V(pi) = E[pi(X) Y(eta)]; value with a cost;
    the constrained (budget) version; inference on V(pi) for a fixed pi
7.2 Optimal unconstrained policy is the sign of the CATE; why sign is easier than magnitude
7.3 Empirical welfare maximization: argmax E_n[(2 pi(X) - 1) Y(eta)] equals weighted
    classification with labels sign(Y(eta)) and weights |Y(eta)| (three lines)
7.4 Restricted classes: depth-limited trees (Athey-Wager), linear rules (Kitagawa-Tetenov);
    regret bound sqrt(VC(Pi)/n) + product of nuisance errors (statement only; appendix A.5
    gives the reduction to a uniform law); why plug-in-and-threshold is not the same thing
7.5 Code: `policy_learn()` on the policy DGP; the true tree vs the learned tree; honest value
    on held-out data; comparison with `policytree` if installed; treating the top-q by CATE
    as the budget solution
7.6 Cautions: variance penalization, pessimism, distributionally robust versions, adaptive
    data (Zhan et al.), dynamic regimes; one paragraph each

### 8. Traps
8.1 Binary outcomes: a baseline-only covariate produces a CATE ranking through the logistic
    curvature (Facure ch. 23 reproduced with `sim_hte("binary_outcome")`); risk difference vs
    latent index; what to target
8.2 Extrapolation of parametric BLPs in heavy tails (CausalML Figure 15.34 message)
8.3 Specification search and small-sample noise (Childers); pre-registration of
    heterogeneity dimensions; the "surrogate" interpretation
8.4 Heterogeneity that is not the CATE: effects on distributions (QTE, distribution
    regression) vs distributions of effects (Frechet-Hoeffding/Makarov bounds); the DiD
    exception (Athey-Imbens 2006); code: `quantreg::rq` QTE on the RCT simulation

### 9. Practice
9.1 The replications (see below) and what each teaches
9.2 A reporting checklist for HTE and targeting papers (split design, nuisance learners,
    which score, which validation, policy class, value with CI, robustness across seeds)

### Appendix
A.1 Conditional orthogonality of the DR score and the bias identity
A.2 DR-learner and R-learner error bounds: proof sketch (Foster-Syrgkanis)
A.3 Theorem 15.2.1: normality of DR-loss differences; the variance lower bound (15.A)
A.4 TOC and QINI as covariances; influence functions with the estimated ATE and share
A.5 EWM regret: the reduction to a uniform law over Pi (Athey-Wager sketch)
A.6 References

### Usage vignette: `vignettes/usage_hte_policy.Rmd`
Simulated data only. Sections: build scores; five learners in five lines; BLP/GATE with
bands; grf side by side; score and ensemble; validate; learn and evaluate a policy; plug in
external nuisances (grf and Python) into the same pipeline.

## Replications (Wayne's, anchors in Section 9)

| Draft | Status in folder | Package hooks | Note |
|---|---|---|---|
| Athey-Wager 2021 | Simulation CSVs present; GAIN data private | `dr_scores` + `policy_learn(method = "tree", depth = 2)`; compare with policytree; Figure 2 from the saved CSV; Table III reconstruction | Install policytree; the IV variant of the score is Section 05 material, note it |
| Haushofer et al. 2025 | Full code + processed data; heavy causal-forest runs precomputed | `policy_frontier()` for the CARA/CRRA planner; grf predictions as input; overlap tables | Lightweight route only; document runtime |
| Yoganarasimhan et al. 2023 | Python + R, anonymized data, 7 trial-length arms | Multi-arm is out of v1: treat pairwise arms with `cate_learner` and compare lasso/forest/causal forest on the two-arm contrast; evaluate with `cate_validate` and `policy_value` | Good "model selection out of sample" case |
| Zhan et al. 2024 | Python, synthetic + OpenML | Teach only (adaptive weights) unless Wayne wants an R port of the synthetic experiment | Section 7.6 |

## Open decisions for Wayne

1. Function names: `dr_scores`, `cate_learner`, `cate_blp`, `cate_gate`, `cate_score`,
   `cate_ensemble`, `cate_validate`, `policy_value`, `policy_learn`, `policy_frontier`, `sim_hte`.
2. Implement the exact depth-2 tree in R (proposed) or require policytree (Suggests) for trees.
3. Install policytree now for testing.
4. Continuous-treatment CATE (Section 4.4): BLP only, or also an R-learner with continuous D?
5. Whether the QTE material (8.4) stays here or moves to a later "distributional effects" note.
