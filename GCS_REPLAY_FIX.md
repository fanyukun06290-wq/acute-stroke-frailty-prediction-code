# GCS replay portability fix

Actual local replay found two missing clean-workspace prerequisites: the original
10_第六步GCS构建规则.md and creation of gcs_results/. This candidate includes the original
rule and adds only out.mkdir(parents=True, exist_ok=True) to the GCS constructor.
Scoring, timestamp selection, missingness handling and model calculations are unchanged.
The local replay verified these prerequisites via its wrapper before reconstructing GCS.
Historical frozen exports and manuscript results are not overwritten.
Raw GCS event exports were reused after SHA-256 verification, not newly queried from databases.
The candidate contains no patient data, models, results or local failed-run logs.
This remains an unpublished candidate and is not covered by the old Zenodo DOI.
