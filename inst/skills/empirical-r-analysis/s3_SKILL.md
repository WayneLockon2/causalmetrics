# Empirical R analysis rules

## General empirical rules

- Start from the paper question and estimand before writing estimator code.
- Separate paper logic, replication code, and package functions.
- Treat raw file reading as `data-raw/` or replication code, not exported package functionality.
- Treat variable construction, merges, reshapes, filters, and recodes as paper-specific unless they recur across designs at the diagnostic level.
- Do not promote a repeated R translation pattern to a package function unless it encodes reusable econometric judgment.
- Prefer direct calls to mature packages for generic modeling, matching, weighting, tables, and plotting.
- Report sample restrictions and target-population changes before reporting estimates.
- Keep estimator outputs tied to diagnostics and assumptions, not just point estimates.
- Be explicit when code reproduces a table, a diagnostic, a robustness check, or a data-preparation step.

## Selection-on-observables rules

- Define the estimand before choosing regression adjustment, matching, IPW, or AIPW.
- Treat the adjustment set as a design choice, not a prediction-variable dump.
- Use only pre-treatment covariates in the main adjustment set unless a paper gives a clear causal reason otherwise.
- Do not control for mediators, colliders, post-treatment variables, or treatment-only predictors as a default.
- State conditional ignorability, consistency, and overlap explicitly.
- Check overlap/common support before interpreting adjusted estimates.
- Distinguish overlap from balance: overlap is support; balance is a post-design diagnostic.
- Report trimming rules and how trimming changes the target population.
- Inspect propensity-score distributions, weight tails, and effective sample sizes.
- Compare estimator families to assess model dependence, but do not treat agreement as proof of ignorability.
- Use placebo outcomes, pre-treatment outcomes, validation samples, or sensitivity analysis whenever available.
- AIPW is doubly robust to nuisance-model misspecification, not to hidden confounding, bad controls, or lack of support.
- Cross-fitting reduces overfitting in nuisance prediction; it does not create identification.
- Keep matching, propensity-score estimation, and balance-table package wrappers out of `causalmetrics` unless a later registry decision explicitly promotes them.
- `est_aipw()` owns the causal score and inference; `mlr3` or Python should only produce nuisance predictions.
