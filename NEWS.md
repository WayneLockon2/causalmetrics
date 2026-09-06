# causalmetrics (development version)

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
