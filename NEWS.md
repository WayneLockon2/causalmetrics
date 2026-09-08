# causalmetrics (development version)

* Synthetic control: `synth_control()` (outcome lags and covariate
  predictors with equal, regression-based, or MSPE-optimized `V` weights;
  demeaned synthetic control of Ferman and Pinto; simplex, nonnegative, or
  unconstrained weights; ridge penalty; ridge augmentation of Ben-Michael,
  Feller, and Rothstein; one synthetic control per treated unit for
  staggered adoption; reuse of weights on another outcome), `synth_placebo()`
  (in-space and in-time placebos with MSPE-ratio p-values),
  `synth_conformal()` (Chernozhukov-Wuthrich-Zhu conformal p-values and
  confidence sets), `synth_spec_test()` (Ferman-Pinto demeaned-SC-versus-DiD
  test), `synth_loo()`, `plot_synth()`, `sim_synth_panel()` (linear factor
  model with selection on loadings or levels), with `tidy()`/`glance()`
  methods. Lecture 09 in `inst/lectures/09_sc/`.

* Regression discontinuity helpers around the `rdrobust` family (all in
  Suggests): `rd_adjust()` (flexible covariate adjustment with cross-fitted
  learners, Noack-Olma-Rothe 2024), `rd_extrapolate()` (effects away from the
  cutoff under conditional independence, Angrist-Rokkanen 2015),
  `rd_weak_iv()` (Anderson-Rubin sets for fuzzy RD on the kernel-weighted
  local sample), `rd_checks()` with `rd_balance()` (joint test from stacked
  influence functions), `rd_placebo_cutoffs()`, `rd_sensitivity()`,
  `rd_donut()`, and `plot_rd_checks()`, `rd_plot()`/`rd_bins()` (ggplot RD
  plots with IMSE-optimal bins), `rd_frame()` and tidiers for `RDHonest` and
  `rdhte` objects, `rd_kink()`/`plot_rd_kink()`, and `sim_rd()`.

* Mediation and mechanisms: `mediate_reg()` (product of coefficients,
  interaction and logit variants, g-computation with delta-method, bootstrap,
  or quasi-Bayesian inference; parallel mediators; the mediated shares `S1` and
  `S2` of Wheeler et al. 2022), `mediate_sensitivity()` (Imai-Keele-Yamamoto
  error-correlation curve), `mediate_dml()` (cross-fitted efficient influence
  function for natural direct and indirect effects and the controlled direct
  effect), `mediate_cde()` (sequential g-estimation with treatment-induced
  confounders), `gelbach_decomp()` (order-invariant decomposition of a
  coefficient change), `mediate_iv()` (direct and indirect effects with an
  instrumented mediator and a homogeneity check), `front_door()` (regression
  and influence-function estimators), plotting helpers, and `sim_mediation()`
  with nine designs and known truths. `mediate_reg()`, `mediate_iv()`, and
  `gelbach_decomp()` take a `weights` column. Vignette `usage_mediation`.
  Replication notes in `inst/replications/10_mediation/`: Wheeler et al.
  (2022, LinkedIn training; Tables 1-6, Figure 1 with the wild cluster
  bootstrap, Appendix Table C1), Zhang, Li, and Allenby (2024, text
  mediators; two-step version with the authors' MCMC as reference), and
  Cattan, Salvanes, and Tominey (2025, IV mediation; audit plus a design
  calibrated to their Table 4).

* Instrumental variables: `est_dml()` gains `model = "pliv"` (partially linear
  IV, one or several instruments) and `model = "iivm"` (the LATE as a ratio
  of two doubly robust scores) with a `weak_iv` Anderson-Rubin set; new
  `iv_first_stage()` (effective F), `iv_ar_confidence_set()`,
  `iv_plausibly_exogenous()`, `complier_profile()`, `late_scores()` and
  `late_blp()`, `iv_late_weights()`, `iv_ate_bounds()`,
  `leniency_instrument()`, `mte_curve()`, control-function helpers
  `cf_residuals()`, `cf_bootstrap()`, `cf_ape()`, shift-share diagnostics
  `ssiv_rotemberg()` and `ssiv_shock_level()`, `dml_sensitivity()`, and
  `sim_iv()`. Vignette `usage_iv` includes a Petrin-Train control function
  in a `mlogit` demand model.

* New `est_dml_structured()`: double machine learning for structured outcome
  models with combinatorial treatments (Farrell, Liang, and Misra 2020; Ye et
  al. 2025). Takes out-of-fold nuisance parameter functions from any engine,
  supports the generalized sigmoid, linear, logit, and custom links, and returns
  cross-fitted influence-function estimates, plug-in estimates, and gaps to the
  best combination. The PyTorch first stage lives in
  `inst/python/dedl_nuisance.py`. Replication of Ye et al. (2025) in
  `inst/replications/04_Doubly_Robust_and_Double_ML/ye2025deep.Rmd`.
* Heterogeneous treatment effects and policy learning, built on one object,
  the cross-fitted doubly robust pseudo-outcome of `dr_scores()`:
  `cate_learner()` (S, T, X, DR, and R meta-learners on `mlr3` learners, with
  prediction on new data), `cate_blp()` and `cate_gate()` (best linear
  predictors and group average effects with HC1 standard errors and
  simultaneous bands), `cate_score()` and `cate_ensemble()` (doubly robust
  loss with confidence intervals; best, convex, Q-aggregation, and least
  squares stacking), `cate_validate()` (heterogeneity test, calibration, TOC
  and QINI curves with one-sided simultaneous bands, AUTOC and AUQC matching
  `grf::rank_average_treatment_effect()`), `policy_value()`, `policy_learn()`
  (empirical welfare maximization over exact depth-1 and depth-2 trees,
  linear rules, weighted classifiers, and budgeted rules; matches
  `policytree` when both search all thresholds), `policy_frontier()`
  (targeting impact versus deprivation with CARA or CRRA planners), and
  `sim_hte()`. Plot helpers `plot_cate_blp()`, `plot_cate_gate()`,
  `plot_cate_validation()`, `plot_policy_tree()`, `plot_policy_frontier()`.
  `policy_value(n_boot = )` adds bootstrap percentile intervals and keeps the
  draws; `cm_policy` objects expose the fitted rule in `$model`.
  `policy_value(weights = )` and `policy_learn(weights = )` take unit weights
  (survey weights, or the time-decaying weights of policy learning with
  adaptively collected data) with the adaptively weighted standard error;
  `policy_learn(split_step = )` passes `split.step` to policytree.
* More than two arms: `dr_scores(arms = )` returns a `cm_scores_multi` object
  (one doubly robust score per arm, one-vs-rest propensities, all pairwise
  contrasts); `contrast_scores()` turns any pair into the binary object every
  other tool accepts; `policy_value()` prices arm-valued policies against a
  uniform baseline (`baseline = "<arm>"`, with `gain_pct` and arm shares);
  `policy_learn(method = "tree")` searches trees whose leaves are arms, exactly
  for depth 1 and 2 and through `policytree` for any depth. `bind_scores()`
  stacks score objects; `glance.cm_blp()`. Replications of Ellickson, Kar, and
  Reeder (2023) in `inst/replications/04_Doubly_Robust_and_Double_ML/` and
  Yoganarasimhan, Barzegary, and Pani (2023) in
  `inst/replications/08_Heterogeneous_Treatment_Effects_and_Policy_Learning/`
  (the latter on a population calibrated to the published moments, because
  the shipped data are randomized).

* New `est_dml()`: double/debiased machine learning for the partially linear
  model (partialling-out score, continuous or binary treatment) and the
  interactive model (doubly robust score, ATE and ATT), with cross-fitting,
  pooled or fold-averaged score solutions, repeated cross-fitting with median
  aggregation, nuisance fit and
  identification diagnostics, and supplied or `mlr3`-learned nuisances.
* `est_aipw()` gains `estimand = "ATT"` and shares its cross-fitting and score
  internals with `est_dml()`; `est_dml(model = "irm")` reproduces it exactly.
* `tidy()` and `glance()` methods for `cm_aipw` and `cm_dml` objects, so both
  estimators work in `modelsummary()` tables.
* The bundled Python predictor moved to `inst/python/nuisance_predictor.py` and
  gained `--model plr` (writes `l_hat`, `m_hat`) next to the default
  `--model irm`.
* New vignette `usage_dml_external_nuisances` on learners, external nuisance
  predictions, and repeated cross-fitting with `est_dml()`.

# causalmetrics 0.0.1

## Difference-in-differences toolkit (Section 07)

* `att_gt()` and `aggregate_att()`: Callaway-Sant'Anna group-time effects with
  never- or not-yet-treated comparisons, varying or universal base periods,
  anticipation, covariates by outcome regression, IPW, or doubly robust scores
  (reproducing `DRDID`/`did` to 1e-6), optional cross-fitted `mlr3` nuisances,
  sparse influence functions, clustered multiplier bootstrap, uniform bands.
* `did_imputation()`, `bacon_decomp()`, `twfe_weights()`, `event_study_frame()`,
  `plot_event_study()`, `pretrend_power()`, `honestdid_inputs()`,
  `did_permutation_test()`, `att_dose()`, `plot_att_dose()`.
* `sdid_weights()`, `sdid_se()`, `plot_sdid()`: synthetic DiD weights for a
  weighted `fixest` regression (also SC, DID, DIFP).
* `sim_did_panel()` simulator and the usage vignette `usage_did_variants`.

* `est_aipw()`: augmented inverse propensity weighting for the ATE with
  supplied or cross-fitted `mlr3` nuisance predictions, propensity clipping
  and trimming, and overlap, weight, and fold diagnostics.
* Paper-output helpers: `table_task()`, `wrap_latex_table()`, `kable_notes()`,
  `plot_task()`, and `wrap_latex_figure()`.
* Lecture-note vignettes for Section 2 (experiments and randomized trials)
  and Section 3 (selection on observables), plus a usage vignette on
  external (Python) nuisance predictions for `est_aipw()`.
* Paper replications in `inst/replications/`: Thornton (2008),
  Karaman (2026), and Imbens and Xu (2025).
* R Markdown templates for replication notes, lecture notes, homework,
  presentations, and working papers.
