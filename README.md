# Acute stroke prognosis analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24128087.svg)](https://doi.org/10.5281/zenodo.24128087)

For the revised manuscript, use [GitHub release v1.1.0](https://github.com/fanyukun06290-wq/acute-stroke-frailty-prediction-code/releases/tag/v1.1.0) and cite the [version-specific Zenodo archive](https://doi.org/10.5281/zenodo.24128087). That release contains the revised extraction and analysis workflows and their reproduction instructions.

This repository contains the reproducible analysis code for the retrospective multicohort study of 30-day in-hospital mortality among critically ill patients with acute stroke. The study uses MIMIC-IV version 3.1 for model development and eICU-CRD version 2.0 for external validation.

## Original analysis on the default branch

The files and instructions below document the original analysis. Its source archive is [v1.0.0 on Zenodo](https://doi.org/10.5281/zenodo.22658666).

## Scope

The release contains read-only SQL extraction queries, R scripts for model development and validation, a reproducibility notebook, and the analysis plan. It does not contain MIMIC-IV or eICU-CRD data, patient-level exports, model objects, prediction files, database credentials, or manuscript files.

## Repository layout

- `sql/`: stroke phenotype and analysis-dataset extraction queries.
- `R/`: model development, validation, sensitivity, table, and figure scripts.
- `notebooks/`: portable R Markdown reproducibility notebook.
- `data/`: local input directory; intentionally empty in the public release.
- `results/`: local output directory; intentionally empty in the public release.
- `docs/`: analysis plan, code-availability wording, audit record, and publication steps.
- `LICENSE`: MIT license for the analysis code only; it does not license the MIMIC-IV or eICU-CRD data.

## Reproduce the analysis

1. Obtain credentialed access to both databases through the authorized PhysioNet procedures. Do not place database passwords or connection settings in this repository.
2. Run `sql/01_extract_mimic_iv_v3.1.sql` against the MIMIC-IV v3.1 schema and `sql/02_extract_eicu_crd_v2.0.sql` against the eICU-CRD v2.0 schema. Run the two files separately on their corresponding database connections.
3. Export the query results locally as `data/mimic_analysis.csv` and `data/eicu_analysis.csv`. These files are excluded by `.gitignore` and must never be committed.
4. Install the required R packages in a local R session:

```r
install.packages(c("tidymodels", "xgboost", "glmnet", "ranger"),
                 repos = "https://cloud.r-project.org")
```

5. From the repository root, run model development:

```text
Rscript R/run_stroke_models.R data/mimic_analysis.csv data/eicu_analysis.csv results/model_results
```

6. Run validation, figures, tables, and prespecified sensitivity analyses through `notebooks/acute_stroke_frailty_reproducibility.Rmd` or the individual scripts in `R/`. The default bootstrap setting is 1,000 patient-level resamples. Sensitivity analyses use `R/run_stroke_sensitivity.R`, followed by `R/run_sensitivity_validation.R` and `R/build_sensitivity_report.R`.

## Analysis safeguards

- MIMIC-IV is the development database; eICU-CRD is the external-validation database.
- The 24-hour model uses the prespecified landmark eligibility flag and the outcome `death_24h_to_30d_inhospital`.
- The 0–6-hour analysis is a sensitivity analysis using `death_30d_inhospital`.
- Patient-level identifiers are used only for grouped splitting and bootstrap resampling, not as model predictors.
- The models are research tools and are not validated for direct clinical decision-making.

## What must not be published

Never commit database exports, patient-level derived tables, prediction files, RDS model objects, `.Rhistory`, `.RData`, `.Renviron`, database passwords, Navicat screenshots, or local absolute paths. The repository's `.gitignore` blocks the main data and result extensions, but inspect `git status` before every commit.
