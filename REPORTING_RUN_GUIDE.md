# Reporting adapters — candidate, not a released analysis

Run commands from `revision_2026-11-04` in a fresh authorized local workspace.
After model/validation aggregate outputs exist, run `python initialize_reporting_metadata.py` once.
Then use `python build_aggregate_reports.py --task TASK --dry-run` to inspect inputs and outputs.
TASK is primary, parameters, complete_cases, repeat, subgroups, endpoint, or units.
Remove `--dry-run` only when the listed aggregates are present. Existing reports cause refusal.
Run `python build_dca_presentation.py` after the DCA R analysis, and
`python build_shap_presentation.py` after the SHAP R analysis.
`python build_cohort_accounting.py` requires authorized local cohort files and split registries;
its original accounting rules are retained, but no reviewer attachment is needed.

The CSV dependency register covers all 88 supplement tables. A listed output is not proof
that its producer has been replayed. Unresolved inherited producers are explicitly marked.
The SHAP adapter retains frozen-study sample-size checks and is not a general cohort tool.
No new model fitting, real database extraction, or clinical validation is certified by these adapters.
Static metadata contain code patterns/measurement definitions, not observed patient records.
Historical CVST migration/scope tables are provenance records, not fresh-run numerical outputs.
The new candidate is not archived by the old Zenodo DOI. Review and a new release remain necessary.

Use `python run_reporting_stage.py TASK --plan` for a no-read plan; omit `--plan` to run the aggregate stage with input/output SHA-256 receipts. TASK also accepts dca or shap. A failed run has no success receipt; partial outputs are retained for local inspection.
