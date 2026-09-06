# Section 07 plan: Difference-in-Differences and Panel Methods

Status: 2026-09-05 evening. Package functions, tests (matching did/DRDID to 1e-5 and synthdid to
1e-15), usage vignette `usage_did_variants`, and the full lecture 07 draft (12 sections + appendix,
34 pages) are written. Replications for 07 remain Wayne's. staggered_efficient() deferred.

## Sources read (notes folder, 07_Difference_in_Differences_and_Panel/materials)

- Pedro Sant'Anna, "Causal Inference using Difference-in-Differences" (Emory, Jan 2025),
  lectures 01-14, plus the NABE (Oct 2024) and USP (Nov 2022) condensed decks.
- CausalML book v0.1.2, chapter 16: DML for conditional DiD (ATT score on the
  differenced outcome), repeated-cross-section score in 16.A, minimum-wage example
  with ten learners (Tables 16.2-16.4).
- MS&E 228 (Syrgkanis) L18: DiD as ATT with outcome = Y2 - Y1, Riesz representer,
  CATT square-loss for heterogeneity. L19: dynamic treatment regimes, g-formula,
  dynamic DML, surrogates. L19 is not DiD; it belongs to the policy-learning /
  dynamic-treatment section, not here.
- Childers, "Difference in Differences" (panel FE basics, multi-valued treatment,
  feedback, extensions: interactive FE, matrix completion, synthetic control, ITS).
- Dehejia, "An Introduction to DiD" (2x2 intuition, Ashenfelter dip, DDD, compositional
  changes, fuzzy DiD).
- Wayne's eight replication drafts already in 07/replication:
  arkhangelsky2021 (SDID, manual QP implementation), de2020two (dCDH 2020),
  elenev2024 (staggered health policy, spillovers), freyaldenhoven2019 (pre-event
  trends), li2024 (forward DiD), nagengast2025 (staggered DiD in gravity/PPML),
  roth2022 (pre-test with caution), roth_santanna2023 (functional form).
- Installed locally: fixest, estimatr, sandwich, clubSandwich, glmnet. NOT installed:
  did, DRDID, HonestDiD, synthdid, didimputation, did2s, staggered, bacondecomp,
  panelView, fect, gsynth, augsynth, didFF, etwfe. The replications use only
  fixest/estimatr (+ quadprog for SDID).

## Verdict on each of Pedro's lectures

| Lecture | Include? | Where in the notes |
|---|---|---|
| 01 Introduction (popularity, potential outcomes Y_t(g), ATT(g,t), on/off & dose exercises) | Yes, short | Sec 1-2 |
| 02 Two-by-two (SUTVA, no anticipation, PT, DiD-by-hand = TWFE, influence functions, multiplier bootstrap) | Yes, core | Sec 2, 8, Appendix |
| 03 Clustering (few clusters: Donald-Lang, Conley-Taber, Ferman-Pinto, FRT, design-based) | Brief | Sec 8 |
| 04 Functional form (Roth-Sant'Anna 2023, PT of CDFs, mixture characterization, didFF test) | Yes, short | Sec 2 and 7 |
| 05 Covariates (TWFE-with-controls bias, OR/IPW/DR, Sant'Anna-Zhao Monte Carlo, panel vs RCS) | Yes, core | Sec 3 |
| 06 DML in DiD (lasso nuisances, orthogonality, cross-fitting) | Yes | Sec 3, ties to Section 04 |
| 07 Repeated cross-sections (stationarity, compositional changes, efficiency loss) | Short | Sec 3 and 8 |
| 08 Event studies (ATT(g,t) with two groups, long differences, simultaneous bands) | Yes, core | Sec 4 |
| 09 TWFE with multiple periods (static beta, all-leads-lags equivalence, BJS/Wooldridge baseline, Chen-Sant'Anna-Xie efficiency) | Yes | Sec 4 and 6 |
| 10 Pre-tests (Fadlon-Nielsen style simulation, anticipation, short vs long differences) | Yes | Sec 4 and 7 |
| 11 Staggered problems (Bacon decomposition, dCDH weights, Sun-Abraham contamination) | Yes, core | Sec 5 |
| 12 Callaway-Sant'Anna (identification, aggregation, estimation, inference) | Yes, core | Sec 6 |
| 13 On/off treatments (treatment sequences, dCDH 2024 aggregation, caveats) | Extension | Sec 10 |
| 14 Random timing (design-based, efficient estimator, FRT) | Extension | Sec 8 and 10 |

## Scope decision (Wayne, 2026-09-05)

No `est_did()`. The package does not wrap fixest or estimatr for any DiD variant.
The vignettes show how to run the standard pieces (2x2, TWFE, event-study dummies,
Sun-Abraham, clustering) with fixest/estimatr, and the package supplies only the
estimators and diagnostics those tools cannot produce. Goal: every advanced technique
in Sant'Anna's course can be performed with fixest/estimatr plus this package.

## Technique-to-tool map (Sant'Anna course)

| Course technique | Lecture | Tool |
|---|---|---|
| 2x2 DiD by hand, TWFE equivalence, unit-clustered SE | 02 | fixest/estimatr (vignette only) |
| Influence-function SE and multiplier bootstrap (unit or cluster level) | 02, 08, 12 | package internal `.cm_multiplier_bootstrap()`, exposed through `att_gt()` |
| Few clusters: wild cluster bootstrap, Fisher randomization test | 03 | teach; small `did_permutation_test()` (cluster-level FRT, studentized) |
| PT and functional form: PT of CDFs, implied-density monotonicity test (didFF) | 04 | teach; `did_ff_test()` later if time allows |
| Covariates: outcome regression, Hajek IPW, doubly robust (Sant'Anna-Zhao) | 05 | `att_gt(method = "reg"/"ipw"/"dr")`; a 2x2 is one (g,t) cell |
| DML nuisances, cross-fitting, lasso/forest first steps | 06 | same function; reuses `est_dml()` internals. Panel DR-DiD is also `est_dml(model = "irm", estimand = "ATT")` on the differenced outcome |
| Repeated cross-sections (stationary), efficiency loss vs panel | 07 | `att_gt(sampling = "rcs")` with the Sant'Anna-Zhao RCS score |
| Two-group event studies: ATT(g,t), long vs short differences, universal base period, limited anticipation | 08, 10 | `att_gt(base_period, anticipation)` + `aggregate_att()` |
| TWFE with all leads and lags; binned endpoints | 09 | fixest `i()` (vignette) |
| Imputation / Wooldridge baseline (BJS, Gardner) | 09, 12 | `did_imputation()` (first stage by fixest on untreated cells) |
| Efficient DiD under PT in all periods (Chen-Sant'Anna-Xie) | 09 | teach only (working paper) |
| Pre-test power and conditional bias (Roth 2022) | 10 | `pretrend_power()` |
| Sensitivity to PT violations (Rambachan-Roth) | 10 | HonestDiD in Suggests, fed by `att_gt()` influence functions |
| Bacon decomposition | 11 | `bacon_decomp()` |
| dCDH negative weights of static TWFE | 11 | `twfe_weights()` |
| Sun-Abraham interaction-weighted estimator | 11, 12 | fixest `sunab()` (vignette) |
| Callaway-Sant'Anna: never/not-yet comparison, aggregation (simple, group, calendar, dynamic, balanced), simultaneous bands | 12 | `att_gt()`, `aggregate_att()`, `plot_event_study()` |
| Treatment on and off: groups by first exposure, dCDH 2024 aggregation, per-unit normalization | 13 | `att_gt()` on G_start (mechanically identical, Pedro L13 slide 16) + FS normalization in `aggregate_att()`; teach the caveats |
| Random timing: Roth-Sant'Anna efficient estimator, FRT | 14 | `staggered_efficient()` (theta_0 - X'beta*, plug-in beta*, Neyman variance, studentized FRT); second priority |
| Comparing estimators on one plot | 11, 12 | `event_study_frame()` accepts fixest `i()`/`sunab`, `att_gt`, imputation |

## Package functions (all from scratch; regressions inside use fixest)

1. `att_gt()` -- Callaway-Sant'Anna ATT(g,t). Arguments: id, time, group (first
   treatment period, Inf/NA for never treated), outcome, covariates formula,
   `method = c("dr","reg","ipw")`, `control_group = c("notyet","never")`,
   `base_period = c("varying","universal")`, `anticipation = 0`,
   `sampling = c("panel","rcs")`, nuisance learners (glm/lm default, mlr3, external
   columns, Python) with optional cross-fitting, clipping, `cluster`, `n_boot`.
   Returns the (g,t) table, the n x (g,t) influence-function matrix, pointwise SE,
   sup-t simultaneous bands, pre-test of zero pre-treatment ATT(g,t)s.
   Covers 2x2 with covariates, two-group event studies, and staggered designs.
2. `aggregate_att()` -- simple, group, calendar, dynamic with `balance_e`,
   `min_e/max_e`, custom weights, optional dCDH per-unit-of-treatment normalization;
   SE via influence functions; uniform bands.
3. `did_imputation()` -- BJS/Gardner two-stage; conservative BJS SE + cluster bootstrap.
4. `bacon_decomp()`, `twfe_weights()` -- diagnostics for static TWFE.
5. `event_study_frame()`, `plot_event_study()`, `tidy()`/`glance()` methods.
6. `pretrend_power()` -- Roth (2022).
7. `did_permutation_test()` -- cluster-level FRT with a studentized statistic.
8. `staggered_efficient()` -- Roth-Sant'Anna (2023) design-based efficient estimator
   (second priority).
9. Synthetic DiD, split as Wayne asked:
   - `sdid_weights()` computes the synthetic part only: unit weights omega (simplex QP
     with ridge penalty zeta, quadprog or Frank-Wolfe), time weights lambda, and the
     regularization from the first differences of untreated units. Returns a tidy
     weight table plus the cell weights w_it = omega_i * lambda_t (treated cells 1)
     ready for `fixest::feols(y ~ d | unit + time, weights = w)`. Options
     `estimator = c("sdid","sc","did","difp")` give the corresponding weights
     (SC: no time weights, no intercept; DID: uniform; DIFP: intercept-free unit
     weights), so the same fixest call produces all four.
   - `sdid_se()` re-solves the weights under placebo, bootstrap, or jackknife
     resampling and re-runs the fixest step; returns SE and the resampled estimates.
   - `plot_sdid()` trajectories with lambda shading and unit-weight bars.
   - Optional `control_selection = "forward"` (Li 2024) inside `sdid_weights()`.
   - Promote Wayne's manual functions from the Arkhangelsky replication.

## Teach with existing tools, no wrapper

TWFE static/dynamic (fixest `i()`), Sun-Abraham (`sunab`), clustered SE and wild
cluster bootstrap (fixest/estimatr), HonestDiD (Suggests), didFF (later).

## Extensions: teach, no code

Continuous and multi-valued doses (Callaway-Goodman-Bacon-Sant'Anna 2024), RCS with
compositional changes (Sant'Anna-Xu 2023), nonlinear DiD (Wooldridge 2023,
Athey-Imbens CiC, Nagengast-Yotov PPML), efficient DiD (Chen-Sant'Anna-Xie 2024), CATT
and DR scores for heterogeneity (Callaway-Chen-Sant'Anna; DiD meta-learner 2025;
Section 08), spillovers/SUTVA (Elenev 2024), interactive FE and matrix completion
(Bai 2009; Athey et al. 2021), dynamic treatment regimes (Syrgkanis L19; Section 08/09).

## Skip

Model-based few-cluster derivations beyond one paragraph; Hagemann (2020); fuzzy DiD
details; surrogates (L19).

## Lecture skeleton (one long file on DiD and its variants)

### `lecture_07_difference_in_differences.Rmd` (subsection level, approved to write)

Style as lecture 04: each subsection = problem, intuition, code that shows the
failure and the fix, short sentences. Monte Carlos cached; helpers in
inst/vignette_helpers/07_Difference_in_Differences_and_Panel/.

1. The problem {#sec-problem}
   1.1 Same question, a new source of variation (before/after x treated/control; what PT buys over unconfoundedness; when DiD is the right tool)
   1.2 One experiment, the classic regression breaks three ways (4-cohort Monte Carlo: static TWFE, dynamic TWFE with bins, TWFE + covariates vs att_gt; table + figure)
   1.3 Reading the results (heterogeneous dynamics and forbidden comparisons; contamination of leads and lags; covariates)
   1.4 Roadmap
2. Just enough framework {#sec-framework}
   2.1 Potential outcomes indexed by adoption time (Y_t(g), G, never treated, ATT(g,t))
   2.2 Three assumptions (SUTVA, no anticipation, parallel trends: unconditional, never vs not-yet, conditional)
   2.3 Identification of the 2x2 ATT: skeleton (selection bias; PT imputes the counterfactual; DiD by hand = TWFE; fixest equivalence)
   2.4 From two periods to long differences (ATT(g,t) from Y_t - Y_{g-1}; short differences before treatment)
   2.5 What parallel trends really assumes (levels vs logs; PT of CDFs; the mixture characterization; code where the sign flips)
   2.6 Inference from influence functions (unit-level IF of the 2x2; clustering; multiplier bootstrap; IF standard error = fixest cluster SE)
3. Covariates: three faces of conditional parallel trends {#sec-covariates}
   3.1 Why covariates (covariate-specific trends; Ashenfelter dip)
   3.2 The tempting regression and why it fails (TWFE + X; moment restrictions; Kang-Schafer Monte Carlo)
   3.3 Outcome regression, IPW, doubly robust (formulas; att_gt(method =); four-DGP Monte Carlo)
   3.4 DR-DiD is the ATT score of Section 04 on the differenced outcome (est_dml on dY equals att_gt dr)
   3.5 Machine-learning nuisances and cross-fitting (learner_p, learner_or; forests example)
   3.6 Repeated cross-sections (stationarity; the two DR estimators; efficiency loss; sampling = "rcs")
4. Dynamics with two groups {#sec-dynamics}
   4.1 ATT(g,t) as an event study (event time; long vs short differences; varying vs universal base period)
   4.2 Anticipation (limited anticipation shifts the base period; example with one period of anticipation)
   4.3 TWFE event-study regressions (fixest i(); all leads and lags equivalence; what binning changes)
   4.4 Simultaneous confidence bands (why pointwise bands mislead; sup-t band)
5. Staggered adoption: what goes wrong {#sec-staggered-problems}
   5.1 Static TWFE as a weighted average (Bacon decomposition; bacon_decomp; forbidden comparisons)
   5.2 Negative weights (dCDH weights; twfe_weights)
   5.3 Dynamic TWFE and contamination (Sun-Abraham; Monte Carlo with all leads and lags)
   5.4 What is not the problem (heterogeneity is fine; the estimator is wrong)
6. Staggered adoption: reliable estimators {#sec-staggered-solutions}
   6.1 Identification, aggregation, inference as three steps
   6.2 The building block ATT(g,t) (never vs not-yet; att_gt; ATT(g,t) plot)
   6.3 Aggregations (simple, group, calendar, dynamic; balanced event window; aggregate_att; composition)
   6.4 Sun and Abraham (fixest sunab; equivalence with CS on the last-treated comparison)
   6.5 Imputation estimators (BJS/Gardner; did_imputation; efficiency vs assumptions; Wooldridge)
   6.6 Stacked DiD and local-projections DiD (fixest recipes; implicit weights)
   6.7 Comparing estimators (event_study_frame + plot_event_study; table: PT used, comparison group, covariates, efficiency)
   6.8 Covariates in staggered designs (att_gt dr with covariates; TWFE + X fails again)
7. Assessing parallel trends {#sec-assess}
   7.1 Pre-trends tests and what they cannot tell (power; pretrend_power)
   7.2 Conditioning on passing the pre-test (Roth 2022 bias; simulation)
   7.3 Sensitivity analysis (HonestDiD: relative magnitudes, smoothness, breakdown value; on aggregate_att output)
   7.4 Functional-form test (didFF idea; monotonicity)
   7.5 Placebos and design checks
8. Inference {#sec-inference}
   8.1 What level to cluster (sampling vs design view; unit at least; assignment level)
   8.2 Few clusters (why the CLT fails; wild bootstrap; permutation test; did_permutation_test)
   8.3 Design-based inference and random timing (Roth-Sant'Anna efficient estimator; FRT)
9. Beyond binary absorbing treatments {#sec-beyond}
   9.1 Continuous and multi-valued doses (level vs slope effects; strong PT; att_dose; two readings of TWFE)
   9.2 Treatments that turn on and off (sequences; first-exposure groups; dCDH 2024; caveats)
   9.3 Triple differences (DDD; when the naive regression fails with covariates)
   9.4 Nonlinear and distributional DiD (Poisson ratio-in-ratios; changes-in-changes; quantile DiD)
10. Synthetic difference-in-differences {#sec-sdid}
   10.1 When the comparison group is not credible (few aggregate units; Prop 99)
   10.2 Unit and time weights (SC, DID, SDID as weighted TWFE; regularization; sdid_weights + feols)
   10.3 Inference (placebo, bootstrap, jackknife; sdid_se)
   10.4 Diagnostics (weight concentration, pre-treatment fit; plot_sdid) and pointer to the SC note
11. Practice {#sec-practice}
   11.1 Replications and anchors
   11.2 A DiD checklist
12. Extensions {#sec-extensions} (efficient DiD, CATT and the DiD meta-learner, selection and PT, compositional changes, spillovers, interactive FE and matrix completion, dynamic regimes)
Appendix: A influence functions for the 2x2 (panel, RCS); B orthogonality and remainder of the DR-DiD score; C the Bacon decomposition; D identification of ATT(g,t) (never, not-yet); E multiplier bootstrap and simultaneous bands; References.

### Synthetic DiD: back in this note (Wayne, 2026-09-05)

SDID is a section of lecture 07 (after the staggered material), with the split agreed
earlier: `sdid_weights()` (unit and time weights, zeta, cell weights; `estimator =
sdid/sc/did/difp`), the DiD step in `fixest::feols(..., weights = w)`, `sdid_se()`
(placebo / bootstrap / jackknife re-solving the weights), `plot_sdid()`. Wayne's
separate synthetic control note (`causal_brave/05Synthetic_Control`) remains the deep
treatment of Abadie SC; lecture 07 cross-references it.

### Doubly robust DiD: no separate implementation

The panel DR-DiD score (Sant'Anna-Zhao 2020) is the AIPW ATT score of `est_dml()`
with the differenced outcome; `att_gt()` calls the same internals
(`.cm_crossfit_predict`, `.cm_score_irm(estimand = "ATT")`, `.cm_solve_linear_score`)
on the long difference Y_t - Y_{g-1} within each (g,t) cell and comparison set. Only the
repeated-cross-section score (four cells; CausalML 16.A, SZ2020 dr,rc) is new code.
The `did` package (Callaway) provides att_gt/aggte/ggdid and uses DRDID for the first
step; we do not depend on it. Install did + DRDID dev-only as numerical oracles for the
tests (not in Imports/Suggests); the `mpdta` county panel for the replication is
vendored into inst/extdata (public county data) or loaded from `did` in Suggests.
HonestDiD 0.2.8 installed 2026-09-05 (CRAN); ships LWdata_EventStudy and
BCdata_EventStudy; goes in Suggests.

### Continuous and multi-valued treatments (added on Wayne's request)

`att_dose()` -- Callaway, Goodman-Bacon, Sant'Anna (2024) for a dose D >= 0 with
untreated units at D = 0 (two periods first; staggered via (g,t) cells later).
Targets: level effects ATT(d|d) = E[Y_t(d) - Y_t(0) | D = d] and ATE(d); causal
responses ACRT(d|d) = dE[Y_t(l)|D=d]/dl at l = d and ACR(d). Identification: ATT(d|d)
under standard PT (dose-d units vs untreated); ACR needs "strong" PT (all dose groups
share the counterfactual trend at every dose). Estimation: `dose_type = c("binned",
"spline")`: binned dose groups run the 2x2 DiD per bin (reuse att_gt cell code);
the spline version regresses Delta Y on a B-spline of D among D > 0 and takes
ATT(d|d) = m(d) - E[Delta Y | D = 0] and ACRT(d|d) = m'(d); covariates via reg/ipw/dr
on the same cells; multiplier or nonparametric bootstrap uniform bands over d.
Diagnostics: the TWFE weights on ACR(l) and on ATE(l)/l (the two decompositions in
Pedro's USP slides 33-34) so the note can show the same beta with two readings.
Multi-valued discrete doses are the binned case with ACR(d_j) = ATE(d_j) - ATE(d_{j-1}).

### Performance design (Wayne's requirement: large panels)

- Data handling in data.table by reference; one long panel in memory; no per-cell
  copies (cells are index vectors).
- All regressions and logits through fixest (`feols`, `feglm`), which are C++ and
  handle millions of rows; nuisance fits per (g,t) cell are on subsets, so cost is
  O(#cells x n_cell), not O(#cells x n).
- Influence functions are sparse by construction: the IF of cell (g,t) is nonzero
  only for units in group g and its comparison set. Store them as a sparse matrix
  (Matrix package, dgCMatrix) or aggregate to the cluster level immediately
  (C x K instead of n x K). Option `keep_if = FALSE` to drop unit-level IFs.
- Multiplier bootstrap = one matrix product: (B x C) weight matrix times (C x K)
  cluster-summed IF matrix; BLAS handles B = 999, C = 1e5, K = 500 in seconds.
- Bacon decomposition: O(G^2) small regressions on collapsed group-by-time means, not
  on the microdata; equivalent for balanced panels.
- SDID weights: Frank-Wolfe with warm starts (as synthdid) for large donor pools;
  quadprog only for small N0. Placebo SEs re-solve weights B times; parallelize with
  future.apply.
- Rcpp: not for v1. The hot loops are already in fixest/BLAS; profile first on a
  synthetic 1e6 x 20 panel. Candidates if profiling demands it: the Frank-Wolfe loop
  in sdid_weights, the permutation test loop. Add Rcpp only with a benchmark in the
  usage vignette.
- Memory guardrails: `att_gt()` reports the (g,t) count and the IF matrix size before
  allocating; `allow_unbalanced_panel` with unit-level differencing done by keys.

### Cutting-edge DiD to consider

Implement (in addition to the map above): `att_dose()` (continuous DiD); LP-DiD of
Dube, Girardi, Jorda, Taylor (2023) as a vignette recipe with fixest (clean controls
per horizon), no function.
Teach in Sec 10 with one paragraph each: Chen-Sant'Anna-Xie (2024) efficient
DiD/event study; Callaway-Chen-Sant'Anna heterogeneity (CATT) and the DiD meta-learner
(2025, arXiv 2502.04699); Ghanem-Sant'Anna-Wuthrich (2022) selection and PT;
Marx-Tamer-Tang (2024) PT and dynamic choices; Sant'Anna-Xu (2023) compositional
changes; dCDH (2024) intertemporal effects and dCDH et al. (2024) continuous
treatments with stayers; Wooldridge (2023) nonlinear DiD and Lee-Wooldridge (2023);
Ortiz-Villavicencio-Sant'Anna (2025) triple differences; Goldsmith-Pinkham-Hull-
Kolesar (2024) contamination bias; Rambachan-Roth (2025) design-based uncertainty
and Athey-Imbens (2022); Arkhangelsky-Imbens (2024) panel survey and
Arkhangelsky-Imbens-Lei-Luo (2024) doubly robust panel identification; Imai-Kim-Wang
(2023) PanelMatch for on/off treatments; Liu-Wang-Xu (2024) counterfactual estimators
(fect) with placebo and equivalence tests; Freyaldenhoven et al. (2021) event-study
visualization; Bilinski-Hatfield (2018) and Dette-Schumann (2024) equivalence tests
for pre-trends; Butts (2023) spatial spillovers; Callaway-Drukker-Liu-Sant'Anna (2023)
DiD via ML.

### Additional DiD variants to cover (added 2026-09-05)

- Triple differences (DDD): Gruber 1994 style; Ortiz-Villavicencio and Sant'Anna 2025
  on why the naive DDD regression fails with covariates; short subsection in Sec 3/6.
- Stacked DiD (Cengiz et al. 2019; Baker-Larcker-Wang 2022): clean 2x2 stacks with
  fixest; one paragraph in Sec 6 with its implicit weights.
- dCDH DID_M (switchers vs stayers) and the on/off extension: Sec 5/10.
- Nonlinear and distributional DiD: log vs level (Roth-Sant'Anna 2023), Wooldridge
  2023 Poisson/exponential "ratio-in-ratios" (ties to the Nagengast-Yotov PPML
  replication), Athey-Imbens changes-in-changes, Callaway-Li 2019 quantile DiD: Sec 10.
- Fuzzy DiD (dCDH 2018) and continuous doses (CGS 2024): Sec 10.
- Spillovers and SUTVA (Elenev 2024; Butts 2021 spatial DiD): Sec 10.
- HonestDiD (Rambachan-Roth 2023) is a sensitivity/inference tool, not an estimator:
  Sec 7, applied to att_gt/aggregate_att output.

### Usage vignettes

`usage_did_staggered.Rmd` (att_gt, aggregate_att, plot against fixest sunab and
imputation), `usage_did_external_nuisances.Rmd` (DR-DiD with mlr3/Python nuisances;
extend inst/python/nuisance_predictor.py with `--model did`), no SDID usage vignette here.

## Replications

Existing drafts and where they anchor: de2020two (Sec 5), freyaldenhoven2019 and
roth2022 (Sec 7), roth_santanna2023 (Sec 2/7), elenev2024 (Sec 6, spillovers),
nagengast2025 (Sec 10, nonlinear), arkhangelsky2021 and li2024 (File 2).

Proposed additions, in priority order:
1. Callaway and Sant'Anna (2021) minimum wage and teen employment (county panel; data
   `mpdta` in `did`, also used in CausalML ch. 16). The canonical staggered anchor;
   validates att_gt and aggregate_att.
2. Sant'Anna and Zhao (2020): Monte Carlo Tables 1-8 (DR vs OR vs IPW vs TWFE) and the
   LaLonde DiD application; bridges Sections 03 and 04.
3. Goodman-Bacon (2021): decomposition on Stevenson-Wolfers unilateral divorce data
   (public via `bacondecomp`); validates bacon_decomp().
4. Rambachan and Roth (2023): HonestDiD on Lovenheim-Willen (2019) duty-to-bargain
   laws (data ships with HonestDiD).
5. Nice to have: Roth and Sant'Anna (2023) staggered rollout (Wood et al. police
   training, data in `staggered`); Borusyak-Jaravel-Spiess (2024) simulations;
   burtch2022peer (already in the 04 folder, a DiD design) could move here; one
   marketing paper with a public staggered-adoption package, Wayne to choose.

## Open decisions for Wayne

1. One lecture with two parts, or two files (DiD; synthetic control and panel)?
2. HonestDiD: depend on the CRAN package in Suggests, or re-implement?
3. Function names: att_gt, aggregate_att, did_imputation, bacon_decomp, twfe_weights,
   event_study_frame, plot_event_study, pretrend_power, did_permutation_test,
   staggered_efficient, sdid_weights, sdid_se, plot_sdid.
4. Install did, DRDID, synthdid, didimputation, bacondecomp as dev-only oracles for
   tests (not in Suggests). HonestDiD done.
5. Which marketing application to replicate.
