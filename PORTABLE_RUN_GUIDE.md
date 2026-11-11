# Portable extraction and fresh current-v2 route

Local unpublished candidate, 9 November 2026. No real source data were accessed to prepare this route. Adapter dry-runs and artificial-fixture checks do not establish full study reproduction.

Use a **new run directory**, not the original frozen analysis workspace. Retain the directory layout `revision_2026-11-04/` under the repository root. Commands below assume the repository root is the working directory and `python`, `Rscript` and, for direct extraction, `psql` are installed. Install the R packages listed in `revision_2026-11-04/docs/original_analysis_package_versions.csv` and Python pandas/numpy. The observed package list is not a complete environment lock.

## Check a query without data access

```
python revision_2026-11-04/portable_extract.py --source mimic --mode cohort --dry-run
python revision_2026-11-04/portable_extract.py --source mimic --mode extract --dry-run
python revision_2026-11-04/portable_extract.py --source eicu --mode extract --dry-run
python revision_2026-11-04/portable_gcs_extract.py --source mimic --dry-run
python revision_2026-11-04/register_navicat_export.py --source mimic --kind extract --dry-run
```

Dry-runs read code/SQL only and do not inspect patient files or start a database process. `gcs_*_fragment.sql` files are fragments used by the GCS adapter, not standalone queries to paste into Navicat.

## Existing Navicat route for cohort and analysis CSVs

Continue using the existing database connections. No Navicat installation, upgrade or connection modification is required by this code.

1. Run `sql/cohort_mimic_current.sql` and `sql/cohort_eicu_current.sql` on the corresponding authorized databases. Export UTF-8 CSVs with column headers to `private_data/mimic_cohort.csv` and `private_data/eicu_cohort.csv` inside the new project directory.
2. Run `sql/03_mimic_cvst_v2.sql` and `sql/02_eicu_revision.sql`. Export to `private_data/mimic_analysis_revision.csv` and `private_data/eicu_analysis_revision.csv`. Preserve missing values; do not replace them with zero.
3. Register each of the four exports locally. Example:

```
python revision_2026-11-04/register_navicat_export.py --source mimic --kind cohort --attest-current-sql
python revision_2026-11-04/register_navicat_export.py --source mimic --kind extract --attest-current-sql
python revision_2026-11-04/register_navicat_export.py --source eicu --kind cohort --attest-current-sql
python revision_2026-11-04/register_navicat_export.py --source eicu --kind extract --attest-current-sql
```

The attestation means the operator confirms the indicated SQL produced the export. The receipt checks file hashes, eligibility/outcome invariants and allowed predictor/missingness fields; it does **not** independently prove which query Navicat executed. The canonical SQL hash records the shared cohort definition; cohort-only queries are derived from that canonical CTE prefix.

## Alternative direct psql extraction

Keep credentials in an authorized private libpq service/password configuration or process environment. The adapters accept `PGSERVICE` or `PGDATABASE` and standard libpq settings. Optional `STROKE_MIMIC_SERVICE` and `STROKE_EICU_SERVICE` select source-specific service names. Do not put passwords in source files, command arguments or public configuration. Navicat settings are not automatically imported by libpq.

```
python revision_2026-11-04/portable_extract.py --source mimic --mode cohort
python revision_2026-11-04/portable_extract.py --source mimic --mode extract
```

Repeat for eICU after selecting its authorized connection. The adapter sets read-only sessions, timeouts and `ON_ERROR_STOP`. It refuses existing outputs or partials. Database stderr is captured but not echoed because it may contain query fragments or identifiers. A failed `.partial` remains private; inspect it locally and use a new run directory rather than overwriting the frozen evidence. SQL dry-run is not a PostgreSQL execution/plan check; `--mode explain` contacts the selected source when used without `--dry-run`.

## Primary models and GCS-dependent stages

After both extraction receipts exist:

```
Rscript revision_2026-11-04/R/start_revision_models.R revision_2026-11-04
Rscript revision_2026-11-04/R/run_revision_uncertainty.R revision_2026-11-04
```

The first command is computationally substantial. It changes run/source fingerprints after portability edits. Do not reuse old fitted/cache files in this fresh directory.

The GCS adapter requires the authorized source connection and the corresponding cohort CSV:

```
python revision_2026-11-04/portable_gcs_extract.py --source mimic
python revision_2026-11-04/portable_gcs_extract.py --source eicu
python revision_2026-11-04/13_build_gcs_features.py
python revision_2026-11-04/18_gcs_rule_sensitivity.py
```

The event fragments are byte-for-text identical to those in the original GCS extractor. The existing feature builder subsequently applies prediction-time availability and assessability rules. Event extraction alone is not evidence of harmonized neurological severity.

## Current CVST labels without the historical overlay

When the new current-v2 exports and new prediction files exist:

```
python revision_2026-11-04/prepare_fresh_v2.py --with-predictions
Rscript revision_2026-11-04/R/run_subgroup_fresh_v2.R revision_2026-11-04
```

The first command requires current-v2 extraction receipts and verifies that prediction labels already match the source labels. It creates byte-identical compatibility copies, edits no probabilities and does not apply the historical 18-record overlay. It rejects mismatched labels or historical source receipts rather than silently repairing them.

The second command uses the original subgroup bootstrap calculations with current paths and newly generated bootstrap samples for every group. It does not read historical subgroup caches. `run_subgroup_cvst_v2.R`, `check_cvst_v2_identity.R` and historical overlay/reuse checks remain provenance tools for the original revision, not this fresh-start entry point.

## Remaining reproduction work

The accompanying execution-order map identifies neurological-reference, missingness, endpoint, repeat-admission, complete-case, DCA/SHAP and reporting stages. A complete automated orchestration and reconciliation of every intermediate JSON/CSV/dictionary and every manuscript table/figure is still pending. The aggregate report builders and historical audit checks must not be assumed to form an executable clean-environment pipeline just because their files are present.

Before release: complete that dependency replay, compare aggregate outputs to the exact manuscript version, finalize CITATION/README metadata, and archive the verified version. The old DOI does not archive this portable candidate. This step did not publish a tag, DOI or release.
