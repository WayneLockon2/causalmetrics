# Section 04 plan: Doubly Robust Estimation and Double Machine Learning

Status: APPROVED 2026-09-04. Package, tests, Python predictor, and usage vignette DONE
(see NEWS.md). Lecture 04 and the replications are being written by Wayne.

## Sources (notes folder, 04_Doubly_Robust_and_Double_ML/materials)

- CausalML book v0.1.2 (Chernozhukov, Hansen, Kallus, Spindler, Syrgkanis):
  chapter 9 = DML (PLM 9.2, IRM/ATE/ATT/GATE 9.3, generic DML 9.4, orthogonality
  calculations 9.B); chapter 4 = double lasso + Neyman orthogonality; chapter 3 =
  HD linear prediction; chapter 10 = feature engineering (NOT DML in this edition).
- MS&E 228 Lecture 13-14 (Syrgkanis): plug-in failure, AIPW as debiased formula,
  PLM/FWL generalization, cross-fitting variants.
- CMU "Semiparametric Functional Estimation" (Childers): influence functions,
  Neyman orthogonality, cross-fit DML theorem, PLM example.
- Chernozhukov et al. (2018, Econometrics Journal) replication bundle:
  401K (sipp1991.dta), Bonus (penn_jae.dat), AJR, Sim; Moment_Functions.R has the
  reference DML1/DML2 implementation (plinear, interactive, IV, LATE; median over splits).
- Ellickson, Kar, Reeder (2023, Marketing Science) public replication package
  (DigitalPromo_Replication.csv, EmailComponents.csv; NDA-limited, workflow only).

## Conceptual framing (to state in the lecture)

- Doubly robust = property of a SCORE (AIPW): bias is a product of the two nuisance
  errors, so consistency needs one correct nuisance and root-n inference needs the
  product of rates to vanish (each ~ n^-1/4). Not "switch to the other model".
- DML = recipe: (i) Neyman-orthogonal score, (ii) cross-fitting, (iii) ML nuisances
  with n^-1/4 rates. AIPW with cross-fitted ML nuisances IS DML for the ATE (IRM).
- PLM/Robinson = generalized FWL: residualize Y and D on X, regress residuals.
  Orthogonal but not doubly robust (partialling-out score); the IV-type score is DR.
- Other DML flavours: PLIV, interactive IV (LATE), GATE/ATT scores, DR-DiD, R/DR
  learners for CATE, DML1 vs DML2, repeated cross-fitting with median aggregation,
  learner selection by cross-fitted MSPE, double lasso as the linear special case.

## Deliverables

### 1. Package: `est_dml()` (R/est-dml.R) + shared internals

API sketch:
  est_dml(data, y, d, x = NULL,
          model = c("plr", "irm"), estimand = c("ATE", "ATT"),   # estimand: irm only
          l_hat, m_hat,                       # plr supplied nuisances E[Y|X], E[D|X]
          p_hat, mu0_hat, mu1_hat,            # irm supplied nuisances
          fold_id, learner_y, learner_d, learner_p, learner_mu0, learner_mu1,
          folds = 5, n_rep = 1, cross_fit = TRUE, seed = NULL,
          solve = c("pooled", "fold_average"),   # DML2 / DML1 in the paper
          p_clip = c(0.01, 0.99), trim = NULL,
          conf_level = 0.95, na_action = c("fail", "omit"))
Returns class cm_dml: estimate, std.error, conf.low/high, per-fold estimates,
per-repetition estimates and the median-aggregated result, score vector, residuals
(plr), nuisance predictions, diagnostics (cross-fitted RMSE / R2 of each nuisance,
var(D_tilde) as identification strength, propensity/weight/ESS for irm, fold
summary), call.

Internals: .cm_crossfit_predict() (generic K-fold out-of-fold prediction with
optional subset), reuse .cm_make_folds(), score functions .cm_score_plr(),
.cm_score_irm_ate(), .cm_score_irm_att(), linear-score inference helper
(theta = sum(psi_b)/sum(psi_a); SE via Jacobian), repeated-split aggregation per
Chernozhukov et al. (2018) sec 3.4 (median estimate; SE^2 = median(SE_s^2 + (theta_s - theta)^2)).
est_aipw() keeps its API and shares the internals; gains estimand = "ATT".
tidy()/glance() methods for cm_dml and cm_aipw (generics import) so modelsummary works.
Tests: FWL exactness (OLS nuisances, no cross-fit == lm coefficient), irm == est_aipw,
dml1 vs dml2, n_rep aggregation, supplied vs learned nuisances, ATT in a heterogeneous
DGP, validation errors, tidy/glance shapes.

### 2. Python contract extension

The predictor moved to inst/python/nuisance_predictor.py and gained
--model {irm,plr}: plr writes l_hat, m_hat, fold_id (regression for Y, regression or
classification for D). Same row-order contract. Single split only (n_rep = 1).

### 3. Vignettes

- vignettes/lecture_04_doubly_robust_and_double_ml.Rmd (Section 1 'The problem' written
  2026-09-05 by Claude on Wayne's request; Sections 2-8 and appendices pending). Helpers in
  inst/vignette_helpers/04_Doubly_Robust_and_Double_ML/dml_simulation_helpers.R; the Monte
  Carlo chunk is cached with cache.extra keyed on the helper file's md5. Agreed skeleton, building on lecture 03 (which already
  covers AIPW, double robustness, Robinson residualization pictures, the DR bias
  algebra, the IF standard error, and a cross-fitting recipe):
  1 The problem: ML plug-ins fail for causal parameters (regularization bias and
    own-observation/overfitting bias; opening simulation = book Fig 9.1/9.2 via
    est_dml(cross_fit = FALSE) and a naive plug-in written in the vignette).
  2 Just enough framework: target vs nuisance, moment M(theta, eta), Neyman
    orthogonality (Gateaux derivative), the n^-1/4 rate and product remainders,
    oracle equivalence. Double robustness as the special case with a product bias.
  3 Three views of DML: (a) partially linear model / generalized FWL (score,
    continuous D, overlap-weighted APE when PLM is wrong); (b) interactive model
    (AIPW = orthogonal score for ATE; ATT score; GATE); (c) high-dimensional linear
    special case = double lasso and why single selection fails.
  4 Cross-fitting: what it fixes, K-fold, why sample splitting restores the CLT
    without Donsker conditions (three-term decomposition, CMU notes).
  5 The generic recipe in practice: linear scores and the Jacobian, pooled vs
    fold-average solutions, repeated splits with median aggregation, learner choice
    by cross-fitted RMSE and ensembles, tuning inside folds, K, clipping/trimming.
  6 Diagnostics and failure modes: nuisance RMSE, residual treatment variance,
    fold spread, split dependence, learner sensitivity, weak overlap, bad controls
    and hidden confounding still bite.
  7 Practice: anchors (Chernozhukov 2018 Tables 1-2; Ellickson 2023), checklist.
  8 Extensions in one page each: PLIV/LATE (Sec 05), DR-DiD (Sec 07), CATE via
    R-/DR-learners (Sec 08), automatic DML / Riesz representers, TMLE.
  Appendices: A orthogonality calculations (PLM, AIPW; book 9.B); B proof sketch
  (influence, empirical-process, remainder terms); C variance, pooled/fold-average
  and median-aggregation formulas; D single vs double selection simulation.
- vignettes/usage_dml_external_nuisances.Rmd: est_dml with mlr3 learners, with the
  Python predictor, repeated cross-fitting, modelsummary table via tidy/glance.
- Update lecture 03 cross-references (04_double_ml -> new name).

### 4. Replications (inst/replications/04_Doubly_Robust_and_Double_ML/)

- chernozhukov2018double.Rmd: Table 1 (Penn bonus, PLR + IRM) and Table 2 (401(k)
  eligibility ATE, PLR + IRM) with est_dml(); Figures 1-2 simulation. Table 3 (LATE)
  and Table 4 (AJR PLIV) deferred until IV models exist (Section 05).
- ellickson2023estimating.Rmd: stage 1 orthogonalized scores per promotion pair via
  est_dml(model = "irm") ($score), stage 2 projection on email components via
  direct lm/fixest with clustered SEs. Workflow replication (public data is synthetic).
- data_raw/04_Doubly_Robust_and_Double_ML/{chernozhukov2018double, ellickson2023estimating}
  (git-ignored copies from the notes folder), read via params$data_path.

### 5. Verification

devtools::test(); smoke-test helper figures on strict pdf(); R CMD build + check --as-cran;
NEWS and DESCRIPTION updates (Suggests: glmnet for lasso learners).

## Order of work

1 package (est_dml, ATT, tidy/glance, tests, docs, check)  ->  2 Python + usage vignette
->  3 lecture 04 + helpers  ->  4 replications  ->  5 final check + commit points.

## Decisions taken

PLR + IRM now, IV models in Section 05; partialling-out score only (iv_type
deferred); n_rep default 1 with median aggregation available; anchors Chernozhukov
Tables 1-2 + Ellickson; double lasso as a conceptual section with glmnet learners;
Python script extended in place; separate usage vignette for est_dml().
