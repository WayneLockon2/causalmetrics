# causalmetrics

<!-- badges: start -->
<!-- badges: end -->

**causalmetrics** is one repository serving two purposes:

1. **A daily-use R toolkit for applied causal inference** in economics and
   marketing. The package implements the causal scores, diagnostics, and
   inference procedures that mature packages do not provide, and deliberately
   does *not* wrap what already works: `fixest`, `estimatr`, `rdrobust`,
   `MatchIt`, `mlr3`, and friends are used directly as engines.
2. **A self-contained course on causal inference**, written while learning it:
   one lecture note per identification strategy (knitted PDFs ship in
   `inst/lectures/`), each paired with replications of published papers
   (`inst/replications/`) that must run on this package, and an LLM-loadable
   *skill* (`inst/skills/`) that turns the lecture into an operational
   workflow.

The repository grows section by section; new estimators, lectures,
replications, and skills are added as I study them.

## Design principles

- The package owns the **estimand, the causal score, and the inference**;
  machine learning only produces nuisance functions (cross-fitted, or passed
  in as out-of-fold predictions from any external stack, including Python).
- **Engines are called directly.** No wrappers around `lm`/`fixest`/
  `estimatr`/`rdrobust`/`MatchIt`; where a function uses one internally (as
  `did_imputation()` or `mediate_iv()` use `fixest`), it is because the
  function owns an estimand those engines do not.
- **Estimates travel with their diagnostics.** Returned objects carry overlap,
  trimming, pre-tests, nuisance fits, and weight summaries; `tidy()` and
  `glance()` methods everywhere.
- **Everything is testable.** Estimators are tested against authoritative
  references (`did`/`DRDID`, `synthdid`, `policytree`, `grf`, `staggered`,
  hand-derived formulas) or against simulators with known truths
  (`sim_hte()`, `sim_iv()`, `sim_rd()`, `sim_did_panel()`,
  `sim_synth_panel()`, `sim_mediation()`).

## Installation

```r
# install.packages("pak")
pak::pak("WayneLockon2/causalmetrics")
```

Attaching the package also attaches the standard empirical toolkit
(tidyverse, data.table, fixest, estimatr, broom, modelsummary, kableExtra).

## Quick example

```r
library(causalmetrics)

# A staggered difference-in-differences with known truth
d  <- sim_did_panel(n_units = 400, n_periods = 10, groups = c(4, 7), seed = 42)
cs <- att_gt(d, id = "id", time = "time", group = "g", y = "y",
             x = c("x1", "x2"), method = "dr", control_group = "notyet")
aggregate_att(cs, type = "simple")$overall          # overall ATT
dyn <- aggregate_att(cs, type = "dynamic", min_e = -3)
plot_event_study(event_study_frame(CS = dyn))       # event study with uniform bands
```

## Methods in the package

Organized by the empirical logic of each design: how the effect is
estimated, how the identifying assumption is checked, how fragility is
quantified, and how the estimate is interpreted or decomposed.

### Randomized experiments (Section 02)

By design the package adds little here: estimation and randomization
inference stay in `estimatr`, `fixest`, `randomizr`, and `ri2`, used
directly. The package's contribution is downstream, once the ITT stands:
doubly robust scores with the known propensity for heterogeneity analysis
(`dr_scores(p_hat = )`), compliance bounds and complier description for
encouragement designs (`iv_ate_bounds()`, `complier_profile()`), and the
reporting helpers.

### Selection on observables (Section 03)

- **Estimation.** AIPW for the ATE and ATT (`est_aipw()`), with nuisances
  fit internally by `mlr3` learners or supplied as out-of-fold predictions
  from any external stack (including Python).
- **Design checks before estimates.** The returned object carries the
  propensity distribution, common support, clipped and trimmed counts, IPW
  weight tails, and effective sample sizes, so the overlap discussion and
  the change in target population under trimming are part of the output,
  not an afterthought.
- **Sensitivity.** Hidden-confounding bounds by partial R-squared
  benchmarking (`dml_sensitivity()` on a companion partially linear fit),
  reporting the confounder strength that would erase the estimate.

### Double machine learning (Section 04)

- **Estimation.** Neyman-orthogonal scores for the partially linear and
  interactive models (`est_dml()`, PLR/IRM, ATE/ATT), with K-fold
  cross-fitting, repeated partitions, and both DML1 and DML2 solutions.
- **Nuisance quality as a first-class output.** Per-nuisance cross-fitted
  RMSE and fold-level estimate spread are returned, so learners are chosen
  by fit, never by the estimate they produce.
- **Structured outcome models.** Debiased contrasts when a flexible first
  stage (for example a PyTorch network, script included) estimates the
  unit-level parameters of a known parametric outcome model
  (`est_dml_structured()`, with the ridge stabilization the saturated-link
  case requires).
- **Sensitivity and reuse.** `dml_sensitivity()` for robustness values;
  `bind_scores()` to carry scores into the heterogeneity toolkit.

### Instrumental variables (Section 05)

- **Estimation.** 2SLS itself stays in `fixest`; the package adds the
  orthogonal IV scores for many or nonlinear controls
  (`est_dml(model = "pliv")`) and the doubly robust LATE
  (`model = "iivm"`), plus control functions for nonlinear second stages
  (`cf_residuals()`, `cf_bootstrap()`, `cf_ape()`).
- **Instrument strength.** Conventional, robust, and effective first-stage
  F (`iv_first_stage()`), and weak-instrument-robust Anderson-Rubin
  confidence sets that remain valid at any strength
  (`iv_ar_confidence_set()`, or `weak_iv = TRUE` inside the DML scores).
- **Exclusion sensitivity.** Bounds under a direct instrument effect
  (`iv_plausibly_exogenous()`) and bounds from LATE toward the ATE under
  one-sided compliance (`iv_ate_bounds()`).
- **Whose effect it is.** Complier shares and characteristics via kappa
  weighting (`complier_profile()`), per-instrument LATEs with the 2SLS
  weights that combine them (`iv_late_weights()`), complier-level
  heterogeneity (`late_scores()`/`late_blp()`), and the marginal treatment
  effect curve when the policy question sits off the LATE margin
  (`mte_curve()`).
- **Constructed instruments.** Leave-one-out judge/examiner leniency
  (`leniency_instrument()`) and shift-share diagnostics: Rotemberg weights
  and shock-level equivalence (`ssiv_rotemberg()`, `ssiv_shock_level()`).

### Regression discontinuity (Section 06)

- **Estimation.** The `rdrobust` family is used directly; the package
  builds the surrounding workflow.
- **Seeing the design.** Binned RD plots with the estimator's own local
  fits rather than global polynomials (`rd_plot()`, `rd_bins()`).
- **Validity battery.** Manipulation (via `rddensity`), covariate balance
  with a joint test, placebo cutoffs, bandwidth sensitivity, and
  donut-hole estimates, bundled in one object (`rd_checks()`, or
  individually: `rd_balance()`, `rd_placebo_cutoffs()`,
  `rd_sensitivity()`, `rd_donut()`).
- **Precision and fuzzy designs.** Cross-fitted covariate adjustment that
  must move the standard error, not the estimate (`rd_adjust()`); for
  fuzzy designs, the first-stage effective F with an Anderson-Rubin
  interval when it is weak (`rd_weak_iv()`).
- **Beyond the sharp jump.** Regression kink designs with the standard
  four-panel evidence (`rd_kink()`), extrapolation away from the cutoff
  under an explicit conditional-independence assumption
  (`rd_extrapolate()`), and a tidy cross-package comparison table
  (`rd_frame()`).

### Difference-in-differences (Section 07)

- **Estimation.** Group-time effects `ATT(g,t)` with doubly robust,
  regression, or IPW scores, not-yet-treated or never-treated comparisons,
  covariates in both nuisances, panels or repeated cross-sections
  (`att_gt()`); aggregated to event studies, cohort, calendar, or overall
  effects with simultaneous confidence bands (`aggregate_att()`).
- **Why not TWFE.** The Goodman-Bacon decomposition of the pooled
  coefficient into its two-by-two comparisons (`bacon_decomp()`) and the
  de Chaisemartin-d'Haultfoeuille weights with their negative share
  (`twfe_weights()`) show what a naive regression would actually average.
- **Parallel trends, assessed rather than assumed.** Pre-treatment cells
  with a joint pre-test (inside `att_gt()`); the power of that pre-test
  against linear violations and the bias an undetected trend would leave
  (`pretrend_power()`); honest sensitivity intervals and breakdown values
  under bounded violations (`honestdid_inputs()` feeding the
  Rambachan-Roth `HonestDiD` machinery); anticipation windows as an
  explicit option (`anticipation = `).
- **Placebo and design-based inference.** Permutation tests of adoption
  timing (`did_permutation_test()`), multiplier-bootstrap bands
  throughout.
- **Alternative estimators, one comparison frame.** Imputation
  (`did_imputation()`), the efficient estimator under random timing
  (`staggered_efficient()`), synthetic DiD (`sdid_weights()`,
  `sdid_se()`), continuous doses (`att_dose()`), and Sun-Abraham via
  `fixest`, all overlaid with `event_study_frame()` and
  `plot_event_study()`.

### Heterogeneous effects and policy learning (Section 08)

- **One input for everything.** The cross-fitted doubly robust
  pseudo-outcome (`dr_scores()`), whose conditional mean is the CATE;
  multi-arm designs via per-arm scores and contrasts
  (`contrast_scores()`).
- **Honest summaries first.** Best linear predictor of the CATE and group
  average effects with simultaneous sup-t bands (`cate_blp()`,
  `cate_gate()`).
- **Flexible models, disciplined.** S/T/X/DR/R meta-learners on any
  `mlr3` learner (`cate_learner()`), compared out of sample by DR loss
  (`cate_score()`) and combined by Q-aggregation (`cate_ensemble()`).
- **Validation.** Calibration by predicted-effect groups, targeting
  curves, and rank-weighted average effects that match `grf` exactly
  (`cate_validate()`: TOC, QINI, AUTOC).
- **From effects to decisions.** Empirical-welfare policy trees, linear
  and budget rules with held-out values against blanket policies
  (`policy_learn()`, `policy_value()`), and impact-versus-need frontiers
  for planners with distributional preferences (`policy_frontier()`).

### Synthetic control (Section 09)

- **Estimation.** Donor weights from outcome lags and covariates with
  equal, regression, or MSPE-optimized predictor weighting; demeaned and
  ridge-augmented variants; separate fits per treated unit under
  staggered adoption (`synth_control()`).
- **Fit quality as identification evidence.** Pre-treatment balance,
  demeaned-path fit, and weight concentration (effective donors) reported
  by the object; the Ferman-Pinto specification test against DiD
  (`synth_spec_test()`).
- **Inference with one treated unit.** In-space and in-time placebos with
  MSPE-ratio rank p-values (`synth_placebo()`) and conformal p-values and
  intervals (`synth_conformal()`).
- **Fragility.** Leave-one-donor-out refits (`synth_loo()`); the
  synthetic DiD bridge for panel audiences (`sdid_weights()`,
  `sdid_se()`).

### Mediation and mechanisms (Section 10)

- **Estimands mapped to assumptions.** Natural direct and indirect
  effects under sequential ignorability, by parametric g-computation with
  delta-method, bootstrap, or quasi-Bayesian inference and parallel
  mediators (`mediate_reg()`, including the Wheeler-style mediated
  shares), or by the efficient influence function with cross-fitted ML
  nuisances (`mediate_dml()`).
- **Sensitivity is part of the estimate.** The indirect effect as a
  function of the mediator-outcome error correlation, with the value that
  erases it (`mediate_sensitivity()`).
- **When sequential ignorability fails.** Controlled direct effects by
  sequential g-estimation under treatment-induced confounding
  (`mediate_cde()`); direct/indirect decompositions with an instrumented
  mediator and its homogeneity check (`mediate_iv()`); front-door
  identification when treatment-outcome confounding is unobserved
  (`front_door()`).
- **Accounting, labelled as accounting.** Order-invariant Gelbach
  decompositions of coefficient changes (`gelbach_decomp()`), kept
  distinct from causal mediation.

### Simulators and reporting

Every section has a simulator with known truths stored in attributes
(`sim_hte()`, `sim_iv()`, `sim_rd()`, `sim_did_panel()`,
`sim_synth_panel()`, `sim_mediation()`), used by the tests, the lectures'
Monte Carlos, and the skills' self-checks. Reporting helpers
(`table_task()`, `kable_notes()`, `wrap_latex_table()`,
`wrap_latex_figure()`) produce the papers' table conventions.

Every function has a usage vignette; see `browseVignettes("causalmetrics")` or
[`vignettes/`](vignettes/).

## Lecture notes

Section-by-section notes (each PDF is knitted from the `.Rmd` beside it), with
proofs in appendices and Monte Carlo evidence computed in the document.

| Section | Topic | File |
|---|---|---|
| 02 | Experiments and Randomized Trials | [lecture\_02\_experiments\_and\_randomized\_trials.pdf](inst/lectures/02_rct/lecture_02_experiments_and_randomized_trials.pdf) |
| 03 | Selection on Observables | [lecture\_03\_selection\_on\_observables.pdf](inst/lectures/03_soo/lecture_03_selection_on_observables.pdf) |
| 04 | Doubly Robust Estimation and Double Machine Learning | [lecture\_04\_doubly\_robust\_and\_double\_ml.pdf](inst/lectures/04_dml/lecture_04_doubly_robust_and_double_ml.pdf) |
| 05 | Instrumental Variables and Unobserved Confounding | [lecture\_05\_instrumental\_variables.pdf](inst/lectures/05_iv/lecture_05_instrumental_variables.pdf) |
| 06 | Regression Discontinuity and Kink Designs | [lecture\_06\_regression\_discontinuity.pdf](inst/lectures/06_rd/lecture_06_regression_discontinuity.pdf) |
| 07 | Difference-in-Differences and Its Variants | [lecture\_07\_difference\_in\_differences.pdf](inst/lectures/07_did/lecture_07_difference_in_differences.pdf) |
| 08 | Heterogeneous Treatment Effects and Policy Learning | [lecture\_08\_heterogeneous\_effects\_and\_policy.pdf](inst/lectures/08_hte/lecture_08_heterogeneous_effects_and_policy.pdf) |
| 09 | Synthetic Control | [lecture\_09\_synthetic\_control.pdf](inst/lectures/09_sc/lecture_09_synthetic_control.pdf) |
| 10 | Mediation and Mechanisms | [lecture\_10\_mediation\_and\_mechanisms.pdf](inst/lectures/10_mediation/lecture_10_mediation_and_mechanisms.pdf) |

## Paper replications

Each replication note reproduces a published paper's tables and figures in
their published layout, using this package (plus the paper's own engines),
and records where the public materials fall short of the paper. Raw data live
in a git-ignored `data_raw/` folder (many are under NDA or too large); the
knitted PDFs are committed so the results are readable without the data.

| Section | Paper | Note |
|---|---|---|
| 02 | Thornton (2008, AER): The Demand for, and Impact of, Learning HIV Status | [thornton2008demand.pdf](inst/replications/02_Experiments_and_Randomized_Trials/thornton2008demand.pdf) |
| 02 | Karaman (2026): Review Solicitation and Future Spending | [karaman2026asymmetric.pdf](inst/replications/02_Experiments_and_Randomized_Trials/karaman2026asymmetric.pdf) |
| 03 | Imbens and Xu (2025): Comparing Experimental and Nonexperimental Methods Four Decades After LaLonde | [imbens2025comparing.pdf](inst/replications/03_Selection_on_Observables/imbens2025comparing.pdf) |
| 04 | Ellickson, Kar, and Reeder (2023, Marketing Science): Estimating Marketing Component Effects | [ellickson2023estimating.pdf](inst/replications/04_Doubly_Robust_and_Double_ML/ellickson2023estimating.pdf) |
| 04 | Ye et al. (2025, Management Science): Deep Learning-Based Causal Inference for Combinatorial Experiments | [ye2025deep.pdf](inst/replications/04_Doubly_Robust_and_Double_ML/ye2025deep.pdf) |
| 05 | Hsieh, Du, and Lu (2026, Marketing Science): Single-Source Data for TV Advertising, a Control-Function Strategy | [hsieh2026leveraging.pdf](inst/replications/05_IV_control_fun/hsieh2026leveraging.pdf) |
| 07 | Jiang, Uetake, and Yang (2026, Marketing Science): Premium Adoption in mHealth (Callaway–Sant'Anna workflow) | [jiang2026does.pdf](inst/replications/07_did/jiang2026does.pdf) |
| 08 | Yang, Eckles, Dhillon, and Aral (2024, Management Science): Targeting for Long-Term Outcomes | [yang2024targeting.pdf](inst/replications/08_Heterogeneous_Effects_and_Policy_Learning/yang2024targeting.pdf) |
| 08 | Yoganarasimhan, Barzegary, and Pani (2023, Management Science): Design and Evaluation of Optimal Free Trials | [yoganarasimhan2023design.pdf](inst/replications/08_Heterogeneous_Effects_and_Policy_Learning/yoganarasimhan2023design.pdf) |
| 08 | Zhan, Ren, Athey, and Zhou (2024, Management Science): Policy Learning with Adaptively Collected Data | [zhan2024policy.pdf](inst/replications/08_Heterogeneous_Effects_and_Policy_Learning/zhan2024policy.pdf) |
| 09 | Andersson (2019, AEJ: Economic Policy): Carbon Taxes and CO2 Emissions | [andersson2019carbon.pdf](inst/replications/09_Synthetic_Control/andersson2019carbon.pdf) |
| 10 | Wheeler, Garlick, Johnson, Shaw, and Gargano (2022, AEJ: Applied): LinkedIn(to) Job Opportunities | [wheeler2022linkedin.pdf](inst/replications/10_mediation/wheeler2022linkedin.pdf) |
| 10 | Zhang, Li, and Allenby (2024, Marketing Science): Text Analysis in Parallel Mediation | [zhang2024text.pdf](inst/replications/10_mediation/zhang2024text.pdf) |
| 10 | Cattan, Salvanes, and Tominey (2025, AER): First-Generation Elite (IV mediation audit + calibrated design) | [cattan2025elite.pdf](inst/replications/10_mediation/cattan2025elite.pdf) |

## Skills: the lectures as LLM-loadable workflows

[`inst/skills/`](inst/skills/) turns each lecture into a skill an LLM
assistant (e.g. Claude Code) can follow on real data: entry checks ("does
this design apply to these columns?"), a variant decision table, a workflow
whose every step carries *Do / Look / Judge / Fail* fields, reporting
minimums, pitfalls found in the replications, inline function signatures,
and a runnable self-check on the package's simulators.

- `causalmetrics-router` — start here: maps the question and the data shape
  to one of the nine design skills.
- `rct-analysis`, `selection-on-observables`, `double-ml`, `iv-analysis`,
  `rd-analysis`, `did-analysis`, `hte-policy`, `synthetic-control`,
  `mediation-analysis` — one per lecture.
- `package-registry` — development policy for contributing to the package
  itself.

To use them with Claude Code, symlink or copy the folders into
`~/.claude/skills/` (personal) or `<project>/.claude/skills/`; with any other
LLM, paste the router plus one design skill. `Rscript dev/check_skills.R`
re-runs every skill's self-check and audits every call in every skill against
the package's actual signatures, so the skills cannot drift from the API.

## Repository layout

```
R/, tests/            estimators, diagnostics, simulators, tests
vignettes/            usage vignettes, one per method family
inst/lectures/        lecture notes (Rmd + knitted PDF, per section)
inst/replications/    replication notes (Rmd + knitted PDF, per section)
inst/skills/          LLM-loadable analysis workflows
inst/python/          Python nuisance scripts (PyTorch first stages)
inst/vignette_helpers/ cached Monte Carlo helpers for the lectures
data_raw/             replication inputs (git-ignored; NDA/size)
dev/                  plans, skill checker, scenarios
```

## Roadmap

This is a living study project. Coming as I learn: more estimators and
diagnostics per section, new sections, additional replications of recent
economics and marketing papers, and continued refinement of the skills as
the workflows get exercised on real analyses. Suggestions and issues are
welcome.

## License

MIT © Weiheng Zhang
