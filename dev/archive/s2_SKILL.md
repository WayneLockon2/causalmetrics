## 12.1 General empirical rules

* Start from the paper’s causal question, estimand, and source of variation before reading implementation code.
* Separate raw ingestion, paper-specific cleaning, and design-level diagnostics.
* Treat original papers and original replication packages as authoritative; treat translated R notes as secondary.
* Do not infer package functions from repeated cleaning, recoding, reshaping, or table-formatting code.
* Prefer direct calls to established estimator, table, and plotting packages.
* Keep sample restrictions and variable construction close to the replication that requires them.
* Document whether a code pattern reflects paper logic or translation style.
* When unsure, put code in the vignette, not the exported package API.

## 12.2 Experiment-specific rules

* Always distinguish assignment `Z`, treatment take-up `D`, outcome `Y`, and outcome observation `R`.
* Report ITT before treatment-on-treated or LATE estimates when compliance is imperfect.
* Diagnose treatment assignment rates before estimating treatment effects.
* Balance checks assess implementation and precision; they do not create identification after randomization.
* Diagnose attrition and missing outcomes by assignment arm before interpreting treatment effects.
* Do not condition on post-treatment variables in the main specification.
* Use only pre-treatment covariates for regression adjustment.
* Cluster standard errors at the randomization level or the level of correlated shocks when appropriate.
* Match randomization inference to the actual assignment mechanism.
* Define outcome, subgroup, and treatment-arm families before applying multiple-testing adjustments.
* Treat spillovers and interference as design problems, not nuisance controls.
* Separate assignment effects from compliance-adjusted effects.
* Do not wrap standard estimators such as `lm`, `fixest`, `estimatr`, `ri2`, or `randomizr`.
* Do not promote R translation helpers to causalmetrics functions.
