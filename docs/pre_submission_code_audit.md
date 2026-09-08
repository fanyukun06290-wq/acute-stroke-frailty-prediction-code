# Pre-submission code audit

Audit scope: the cleaned release tree for the acute-stroke prediction manuscript.

## Completed checks

- Removed legacy GPR exploratory scripts and unrelated manuscript/review files from the release tree.
- Excluded all CSV, RDS, RData, parquet, SQLite, database, spreadsheet, and archive files from the release tree.
- Removed machine-specific Windows paths from the R scripts and reproducibility notebook.
- Removed database host, port, username, and password information from the release tree.
- Converted the model, validation, sensitivity, table, and figure scripts to use command-line arguments or environment variables.
- Added `.gitignore`, `CITATION.cff`, input/output directory notices, code-availability wording, and publication instructions.
- Retained author names and affiliations for citation while removing direct contact email addresses from public code metadata.
- R syntax check passed for all 9 R scripts in the release tree.
- R Markdown code-chunk syntax check passed for the reproducibility notebook.

## Author checks still required

- Confirm that the final local SQL exports correspond to MIMIC-IV v3.1 and eICU-CRD v2.0 and are not copied into this repository.
- Review the staged Git file list before the first push.
- Select a software license.
- Replace the repository and DOI placeholders after GitHub and Zenodo publication.
- Re-run the complete analysis from the cleaned tree, if a new release is intended to support exact computational reproduction.
