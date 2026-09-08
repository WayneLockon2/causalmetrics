---
name: causalmetrics-router
description: Entry point for causal analysis with the causalmetrics R package. Load first when the task is to estimate a causal effect from data and the design is not yet fixed - this skill maps the question and the data shape to the right design skill (rct, selection-on-observables, double-ml, iv, rd, did, hte-policy, synthetic-control, mediation).
---

# causalmetrics router: from question and data to a design skill

Answer two questions, load exactly one design skill (plus this file), and
follow its workflow. Each design skill is self-contained: entry checks,
variant choice, steps with Do/Look/Judge/Fail, reporting, signatures, and
a runnable self-check.

## Question 1: what creates the variation in treatment?

| You can say... | Load | Core check before committing |
|---|---|---|
| "Assignment was randomized (possibly within strata/clusters)" | `rct-analysis` | the mechanism can be written down; realized rates match it |
| "No design, but these pre-treatment covariates plausibly absorb confounding" | `selection-on-observables` | the adjustment set survives a variable-by-variable argument; overlap holds |
| "Same claim, but nuisances need ML (many covariates, unknown forms)" | `double-ml` | as above; learners chosen by nuisance fit |
| "An instrument/lottery/judge/shift-share moves treatment, excluded from the outcome" | `iv-analysis` | a written exclusion mechanism; first stage exists |
| "A rule assigns treatment at a cutoff of a running variable" | `rd-analysis` | mass on both sides; manipulation implausibly precise |
| "Units adopt at different times in a panel; untreated comparisons exist" | `did-analysis` | absorbing treatment; pre-periods for every kept cohort |
| "One or a few aggregate units treated; many donors; long pre-period" | `synthetic-control` | T0 large relative to donors; donor pool clean |

If several rows apply, prefer the stronger design: randomization beats
everything; a cutoff or instrument beats adjustment; panel timing beats a
single cross-section. Run the runner-up as a robustness row, not the
headline.

## Question 2: what is the question, beyond one average effect?

- "For whom does it work / whom should we treat?" -> finish the design
  skill, then load `hte-policy` (its input is `dr_scores()` built on the
  design's identification).
- "How does the effect come about (through which mediator)?" -> load
  `mediation-analysis`; bring the design skill's total effect with you.
- "How fragile is this to hidden confounding / trend violations?" -> stay
  in the design skill; each carries its sensitivity step
  (`dml_sensitivity`, HonestDiD, placebo distributions,
  `mediate_sensitivity`).

## Combination notes

- Fuzzy RD with a weak jump: `rd-analysis` step 5 uses the weak-IV tools
  (`rd_weak_iv`, Anderson-Rubin) from `iv-analysis`.
- RCT with noncompliance: ITT in `rct-analysis`, then the LATE machinery
  in `iv-analysis` (`complier_profile`, `iv_ate_bounds`).
- DiD with few treated units: `synthetic-control` (or SDID inside
  `did-analysis`).
- Any design's scores feed `hte-policy`; in an RCT pass the known
  propensity to `dr_scores(p_hat = )`.
- A mediator with its own instrument: `mediation-analysis` (`mediate_iv`),
  not plain `iv-analysis`.

## Shared minimums every skill assumes

- Seeds on every stochastic call; report them.
- Cluster standard errors at the level of assignment or shared shocks;
  with fewer than ~40 clusters, add a wild cluster bootstrap or
  permutation check.
- Estimates travel with their diagnostics; a table without its overlap /
  pretest / F / density row is incomplete.
- Sample changes (trimming, dropped cohorts, matched subsets) change the
  estimand; every skill's report section forces the statement.
- Simulators (`sim_hte`, `sim_iv`, `sim_rd`, `sim_did_panel`,
  `sim_synth_panel`, `sim_mediation`) carry known truths; run the skill's
  self-check before touching real data.

## What this package does not do

Structural/equilibrium models, forecasting, generic prediction, text or
image ML pipelines (engines like mlr3/MatchIt/rdrobust/fixest are used
directly where mature). For extending the package itself, load
`package-registry`, not an analysis skill.
