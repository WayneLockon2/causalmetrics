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

| Area | Main functions |
|---|---|
| AIPW / doubly robust (Sec. 03) | `est_aipw()` (ATE/ATT; internal `mlr3` or external out-of-fold nuisances) |
| Double machine learning (Sec. 04) | `est_dml()` (PLR, IRM; repeated cross-fitting, DML1/DML2), `dml_sensitivity()`, `est_dml_structured()` (Farrell–Liang–Misra structured outcomes, e.g. deep-learning first stages), `bind_scores()` |
| Instrumental variables (Sec. 05) | `est_dml(model = "pliv"/"iivm")`, `iv_first_stage()` (effective F), `iv_ar_confidence_set()`, `iv_plausibly_exogenous()`, `iv_ate_bounds()`, `complier_profile()`, `late_scores()`/`late_blp()`, `iv_late_weights()`, `leniency_instrument()`, `mte_curve()`, control functions `cf_residuals()`/`cf_bootstrap()`/`cf_ape()`, shift-share `ssiv_rotemberg()`/`ssiv_shock_level()` |
| Regression discontinuity (Sec. 06) | built on the `rdrobust` family: `rd_plot()`, `rd_bins()`, `rd_checks()` (density, balance + joint test, placebo cutoffs, bandwidth sensitivity, donut), `rd_adjust()` (cross-fitted covariate adjustment), `rd_weak_iv()` (fuzzy + weak), `rd_extrapolate()`, `rd_kink()`, `rd_frame()` |
| Difference-in-differences (Sec. 07) | `att_gt()` (Callaway–Sant'Anna, DR/reg/IPW, panel/RCS), `aggregate_att()` (dynamic/group/calendar/simple, uniform bands), `did_imputation()`, `staggered_efficient()`, `bacon_decomp()`, `twfe_weights()`, `event_study_frame()`/`plot_event_study()`, `pretrend_power()`, `honestdid_inputs()` (Rambachan–Roth), `did_permutation_test()`, `att_dose()`, synthetic DiD `sdid_weights()`/`sdid_se()` |
| Heterogeneous effects & policy (Sec. 08) | `dr_scores()` (the shared DR pseudo-outcome, multi-arm via `arms`/`contrast_scores()`), `cate_blp()`/`cate_gate()` (sup-t bands), `cate_learner()` (S/T/X/DR/R), `cate_score()`/`cate_ensemble()` (Q-aggregation), `cate_validate()` (calibration, TOC/QINI/AUTOC matching `grf`), `policy_learn()`/`policy_value()`/`policy_frontier()` |
| Synthetic control (Sec. 09) | `synth_control()` (levels/demeaned/ridge-augmented, covariates, staggered `treated_units = "separate"`), `synth_placebo()` (space/time), `synth_conformal()`, `synth_spec_test()` (Ferman–Pinto vs DiD), `synth_loo()`, `plot_synth()` |
| Mediation & mechanisms (Sec. 10) | `mediate_reg()` (g-computation; delta/bootstrap/quasi-Bayesian; parallel mediators; Wheeler S1/S2 shares), `mediate_dml()` (efficient influence function), `mediate_cde()` (sequential g-estimation), `mediate_iv()`, `front_door()`, `gelbach_decomp()`, `mediate_sensitivity()` |
| Simulators & reporting | `sim_hte()`, `sim_iv()`, `sim_rd()`, `sim_did_panel()`, `sim_synth_panel()`, `sim_mediation()` (known truths in attributes); `table_task()`, `kable_notes()`, `wrap_latex_table()`/`wrap_latex_figure()` |

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
