# Fresh runtime dependency changes

`runtime_stages.json` is an ordered dependency plan, not a claim of successful execution.
Commands with cwd=project_parent must run from the directory containing revision_2026-11-04.
Commands with cwd=project run inside revision_2026-11-04. Replace DATABASE with mimic/eicu
and run both extraction commands. No credentials belong in this plan.

Five R files use new runtime receipt paths instead of old audit manifests. Only path literals
were replaced; original candidates remain unchanged. This changes run hashes and must use
a fresh workspace without old fit/cache/output files. It is not a replay of the old cached run.

After each specified stage, `python build_runtime_receipt.py STAGE --plan` lists requirements;
omit --plan to hash the completed files locally. This reads local restricted files as bytes
to compute SHA-256; neither rows nor file contents are printed. Receipts are private QC.
Presence and hash receipts do NOT certify correctness or successful independent validation.
Statistical checker scripts must still run as appropriate before scientific acceptance.

The plan orders major computation/reporting stages, not every independent checker or diagnostic.
Four fresh producers remain unresolved: phenotype counts, new source index, input bounds,
and endpoint time-flag case mix. Two CVST correction tables remain historical provenance.
Full database execution, model refit, and numerical reproduction have not been performed.

Update: the four previously unresolved fresh report producers now have adapters; see FRESH_REPORTS_GUIDE.md. Numerical execution on real data remains pending.
