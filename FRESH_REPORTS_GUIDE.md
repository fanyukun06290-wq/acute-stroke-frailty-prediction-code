# Four fresh-run reporting adapters

Run from the revision_2026-11-04 project directory:

    python build_fresh_reports.py phenotypes --plan
    python build_fresh_reports.py bounds --plan
    python build_fresh_reports.py time_flags --plan
    python build_fresh_reports.py index --plan

`--plan` prints the output path without reading data. Remove it to generate each report
from authorized LOCAL files. Index must run after the other three plus static metadata.
Each report refuses existing outputs/receipts and saves input/output hashes in private qc/.
No identifiers, patient rows, or probabilities are printed to stdout. Calculations necessarily
read local restricted records/predictions; no remote service is called by these scripts.

Phenotype counts use current *_cohort.csv files from canonical current SQL, with separate
coded-source, 6h and 24h populations. Source-level event counts remain undefined.
Bounds cover 18 measurement variables at two windows in two databases (72 rows), not age,
sex, full raw measurement capture, unit harmonization or clinical measurement validity.
SQL/dictionary bounds must agree; value bounds and 18 matching missingness flags are checked.
Time-flag reporting uses the internal MIMIC test partition and full eICU validation population,
requires matching prediction keys/truth, and checks actual hospital/ICU discharge differences.
Empty groups and no-event denominators are recorded as missing, never zero-risk claims.
Index is a fresh file-level SHA-256 register; it replaces the old narrative index for fresh runs.
It does not claim historical reviewer matrices exist or certify scientific correctness.

The two historical CVST correction/migration tables remain historical provenance.
All four adapters require actual local execution before real-data results are verified.
The 46-stage plan still omits some independent checkers/diagnostics; it is not a full verified DAG.
The old Zenodo DOI does not cover this unpublished candidate. Existing manuscript results
and response are not changed by preparing or testing these adapters.
