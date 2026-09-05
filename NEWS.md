# causalmetrics (development version)

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
