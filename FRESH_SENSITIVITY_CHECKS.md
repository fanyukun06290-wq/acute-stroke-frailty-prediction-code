# Independent checks for fresh sensitivity runs

After R/run_subgroup_fresh_v2.R, run R/check_subgroup_fresh_v2.R from the repository parent
directory convention used by the other R scripts. It uses the fresh cache names; only source
entry/cache literals differ from the historical checker. Statistical checks are unchanged.
The historical CVST checker and cache-reuse route remain preserved as provenance.

Complete-case independent checks require model_results/encoded_predictors_24h.csv and
encoded_predictors_6h.csv as well as frozen recipe/model artifacts, raw inputs and split registries.
Current local component replays matched subgroup, complete-case, repeated-record and endpoint
outputs. The complete-case population reference uses the corrected CVST table, not its old version.
This is component-level numerical replay, not full fresh database extraction or clinical validation.
The old DOI does not archive this unpublished candidate. No private results are included.
