---
name: package-registry
description: Development policy for the causalmetrics package itself. Load when adding, promoting, or reviewing package functions, vignettes, lectures, or replications - not when analyzing data. Encodes what belongs in the package versus the vignette or note, the no-wrapping rule, and promotion criteria.
---

# causalmetrics registry: what goes where

## The one rule

The package owns the causal score, its identification assumptions, and
its inference. Everything else - data cleaning, generic modelling,
matching, tables, plots of raw data - belongs to mature packages used
directly or to the analysis document.

## Never wrap

Do not wrap or re-export `lm`, `glm`, `fixest`, `estimatr`, `sandwich`,
`MatchIt`, `randomizr`, `ri2`, `rdrobust`/`rddensity`/`rdlocrand`/
`RDHonest`/`rdmulti`, `mlr3`, or table/plot infrastructure. Using them as
engines inside a function that owns an estimand and its inference is fine
(as `did_imputation()` and `mediate_iv()` use fixest, or `rd_adjust()`
feeds rdrobust); a thin convenience layer over their interfaces is not.

## Where code goes

- `R/` (exported): encodes reusable econometric judgment - a score, a
  diagnostic with a decision attached, an inference procedure. Takes
  `data` + column names, returns a `cm_*` object with `print()` and
  `tidy()`; simulators return known truths in attributes.
- Vignettes (`vignettes/usage_*.Rmd`): call patterns on simulated data;
  every exported function appears in one.
- Lectures (`inst/lectures/<NN>_<short>/`): the method, subsection by
  subsection, with cached Monte Carlos in
  `inst/vignette_helpers/<NN>_*/`; knitted PDFs in place; setup chunks
  `load_all()` the source tree. Folder names short (100-byte tarball
  path limit).
- Replications (`inst/replications/<NN>_<Topic>/<bibkey>.Rmd`): match
  the authors' scripts first, the paper second; published layouts; data
  in git-ignored `data_raw/<NN>_.../<bibkey>/` with `params$data_path`.
  One-off cleaning stays there, never in `R/`.
- Skills (`inst/skills/<name>/SKILL.md`): operational workflows per
  design; keep in sync with signatures (`dev/check_skills.R`).

## Promotion criteria

Promote replication or note code to `R/` only when all hold: it encodes
judgment reusable across papers (not translation style or repeated
cleaning); it owns an estimand/diagnostic/inference rather than
convenience; it can be tested against an authoritative reference
(another package to 1e-5-1e-15, a published table, or a simulator
truth); and a usage vignette section can show it on `sim_*()` data.
When unsure, put it in the vignette or note - demotion is harder than
promotion.

## Analysis-side conventions the code must support

Estimator outputs stay tied to their diagnostics and assumptions (the
object carries them; `print()` shows them). ML only produces nuisances;
external nuisance predictions must be out-of-fold with an explicit fold
map. Sample restrictions and target-population changes are reported by
the object where possible (trimming counts, dropped cohorts).

## Testing bar

Every exported estimator has a test against an external reference where
one exists (`did`/`DRDID`, `synthdid`, `policytree`, `grf`, `staggered`,
`fixest`, hand-computed formulas) and against simulator truths otherwise;
tolerance stated in the test. R CMD check clean apart from the known
long-path note.
