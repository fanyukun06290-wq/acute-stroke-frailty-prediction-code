# Acute stroke ICU prediction

Version 1.1.0 contains the source code for the revised BMC Medical Informatics and Decision Making manuscript. Start with [the portable run guide](PORTABLE_RUN_GUIDE.md) and [the reporting guide](REPORTING_RUN_GUIDE.md).

The study develops models in MIMIC-IV v3.1 and evaluates them in eICU-CRD v2.0. The target is recorded death after an ICU landmark, before hospital discharge and before day 30 from ICU admission. It is not complete 30-day all-cause mortality. The main transported models exclude retrospective stroke subtype and diagnosis-derived comorbidity information. No validated frailty instrument is included.

## Source layout

- `revision_2026-11-04/sql/03_mimic_cvst_v2.sql`: current MIMIC extraction definition.
- `revision_2026-11-04/sql/02_eicu_revision.sql`: current eICU extraction definition.
- `revision_2026-11-04/sql/01_mimic_revision.sql`: historical pre-CVST-label correction, retained for source comparison only.
- `revision_2026-11-04/R/`: fitting, uncertainty, sensitivity, checking and figure code. Historical `*subgroup_revision.R` files are not the current v2 subgroup route.
- `revision_2026-11-04/13_build_gcs_features.py` and `18_gcs_rule_sensitivity.py`: GCS transformations using locally obtained event exports.
- `revision_2026-11-04/37_prepare_cvst_v2_reporting.py`: v2 reporting builder; requires earlier computed outputs.
- `revision_2026-11-04/docs/`: original-analysis package versions. This is an observed environment record, not a dependency lockfile.
- `critical_care_code_release_2026-09-08/R/build_subgroup_threshold.R`: historical threshold-provenance input hashed by the revised DCA code; it is not part of the new primary fitting route.

## Reproduction status

Place the project directory under a working parent directory. R scripts use this retained relative layout. Install packages in your normal R library or set `STROKE_R_LIBRARY`; machine-specific library defaults were removed. Do not run analysis code against an existing frozen project directory: source hashes change after portability cleanup and saved caches must not be mixed between versions.

After obtaining authorized database access, the intended order is:

1. Obtain local cohort/analysis exports using the current SQL, with timing and source-quality summaries.
2. Run primary fitting via `R/start_revision_models.R`, then uncertainty and calibration routines.
3. Extract local GCS events and construct the neurological-reference/rule-sensitivity inputs before fitting the nested reference models.
4. Run missingness, endpoint, repeat-admission, complete-case, DCA, SHAP and corrected-subgroup evaluations.
5. Reconcile aggregate outputs with the exact manuscript/supplement version and archive the reviewed release.

The local replay used a fresh MIMIC extraction and hash-verified eICU and GCS event exports. Component checks covered primary and reference model fitting, uncertainty, neurological sensitivity, DCA, SHAP, corrected subgroups, complete cases, repeated records, and discharge-related endpoint analyses. The revision reports were checked against the manuscript and supplementary tables. Use `run_subgroup_fresh_v2.R` followed by `check_subgroup_fresh_v2.R` for a fresh subgroup analysis. The scripts depend on intermediate outputs documented in `revision_2026-11-04/runtime_stages.json`.

## Data and licensing

MIMIC/eICU require the applicable credentialing, training and data-use permissions. Keep exports, patient identifiers, predictions, fitted models, credentials and screenshots outside version control. The MIT license applies to code only and does not grant rights to source data. Do not publish generated output folders merely because code ran successfully.

## Citation status

Cite the archived record associated with the exact `v1.1.0` GitHub release once Zenodo issues the version-specific DOI. The historical [v1.0.0 archive](https://doi.org/10.5281/zenodo.22658666) documents the earlier analysis version.

## Aggregate reporting additions
See [REPORTING_RUN_GUIDE.md](REPORTING_RUN_GUIDE.md). Full clean-environment reproduction remains pending.

## GCS component replay
See `GCS_REPLAY_FIX.md` for the rule/directory prerequisites found during real local replay. Sixteen neurological reference models and sixteen missingness-sensitivity models were refitted and independently checked. Raw GCS event exports were reused after hash verification. This replay evaluates research models on the study cohorts.

## Fresh sensitivity checks
See FRESH_SENSITIVITY_CHECKS.md for the fresh subgroup checker and complete-case metadata prerequisites found in local replay. The guide lists the input and checker dependencies for these analyses.
