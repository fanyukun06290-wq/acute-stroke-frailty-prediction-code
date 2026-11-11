"""Fingerprint completed local stage files, not certify scientific correctness."""
from pathlib import Path
import argparse,json,hashlib
FILES={
 'model':['model_feature_contract.json','model_results/training_plan.json','model_results/selected_parameters_24h.json','model_results/selected_parameters_6h.json','private_splits/patients_24h.csv','private_splits/patients_6h.csv','private_predictions/predictions_24h.csv','private_predictions/predictions_6h.csv'],
 'gcs':['gcs_results/construction_plan.json','gcs_results/construction_complete.json','gcs_results/availability.csv','gcs_results/quality_audit_counts.csv','gcs_results/score_distributions.csv','private_data/mimic_gcs_features.csv','private_data/eicu_gcs_features.csv'],
 'clinical':['clinical_reference_results/metric_intervals.csv','clinical_reference_results/paired_differences.csv','private_predictions/clinical_ridge_predictions_24h.csv','private_predictions/clinical_ridge_predictions_6h.csv'],
 'neuro':['neuro_reference_results/fit_complete.json','neuro_reference_results/validation_complete.json','neuro_reference_results/metric_intervals.csv','neuro_reference_results/paired_differences.csv','private_predictions/neuro_predictions_24h.csv','private_predictions/neuro_predictions_6h.csv'],
 'missing':['missingness_sensitivity_results/fit_complete.json','missingness_sensitivity_results/validation_complete.json','missingness_sensitivity_results/metric_intervals.csv','missingness_sensitivity_results/paired_domain_differences.csv','private_predictions/missing_sensitivity_predictions_24h.csv','private_predictions/missing_sensitivity_predictions_6h.csv']}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def build(root,stage):
 target=root/'qc'/('runtime_'+stage+'_files.json')
 if target.exists():raise FileExistsError('Existing runtime receipt preserved')
 files=list(FILES[stage])
 if stage=='model':
  files += [f'private_models/{name}_{hour}h.{ext}' for hour in [24,6] for name,ext in [('recipe','rds'),('xgboost','ubj'),('random_forest','rds'),('elastic_net','rds')]]
 if any(not (root/n).is_file() or (root/n).stat().st_size==0 for n in files):raise FileNotFoundError('Required nonempty stage file absent')
 result=dict(stage=stage,scope='file_presence_and_sha256_only_not_analysis_validation',files=[dict(relative_path=n,sha256=sha(root/n)) for n in files])
 target.parent.mkdir(exist_ok=True)
 with target.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=FILES);p.add_argument('--plan',action='store_true');a=p.parse_args()
 if a.plan:print(json.dumps(dict(stage=a.stage,required_files=FILES[a.stage],files_read=False)))
 else:
  result=build(Path(__file__).resolve().parent,a.stage);print(json.dumps(dict(stage=a.stage,files=len(result['files']),scope=result['scope'])))
