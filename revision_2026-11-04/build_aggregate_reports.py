"""Format aggregate results only. Does not read patient records or fit models."""
from pathlib import Path
import argparse,json
import pandas as pd
METRICS=['auroc','auprc','brier_score','calibration_intercept','calibration_slope']
FORBIDDEN={'subject_id','stay_id','hadm_id','uniquepid','patient_group_id','group_id','probability','truth'}
TASKS={
 'primary':(['validation_results/primary_metric_intervals.csv'],['manuscript_revision_results/Table2_primary_models_24h.csv']),
 'parameters':(['model_results/selected_parameters_24h.json','model_results/selected_parameters_6h.json'],['reporting_audit_results/selected_parameter_index.csv']),
 'complete_cases':(['complete_case_results/complete_case_populations.csv','complete_case_results/complete_case_metric_intervals.csv','complete_case_results/missingness_by_event.csv'],['complete_case_results/complete_case_reporting.csv','complete_case_results/missingness_event_association.csv']),
 'repeat':(['patientflow_results/repeat_admission_populations.csv','patientflow_results/repeat_admission_metric_intervals.csv'],['patientflow_results/repeat_admission_reporting.csv']),
 'subgroups':(['phenotype_v2_results/subgroup_populations.csv','phenotype_v2_results/subgroup_metric_intervals.csv'],['phenotype_v2_results/subgroup_reporting_table.csv','phenotype_v2_results/external_phenotype_performance.csv','phenotype_v2_results/external_restricted_evaluation.csv','phenotype_v2_results/external_subgroup_threshold_010.csv','phenotype_v2_results/precision_warnings.csv']),
 'endpoint':(['endpoint_results/endpoint_sensitivity_intervals.csv','endpoint_results/competing_event_curves.csv','endpoint_results/fixed_risk_bin_cif.csv'],['endpoint_results/endpoint_sensitivity_reporting.csv','endpoint_results/day30_competing_event_summary.csv','endpoint_results/risk_bin_reporting.csv']),
 'units':(['qc/mimic_measurement_units.csv','qc/eicu_measurement_units.csv'],['phenotype_v2_results/observed_raw_units.csv'])}

def read(root,name):
 if any(p.startswith('private_') for p in Path(name).parts):raise ValueError('Private input forbidden in aggregate formatter')
 d=pd.read_csv(root/name)
 if FORBIDDEN.intersection(d.columns):raise ValueError('Identifier or row-level prediction field detected')
 return d

def display(row,style='to'):
 if pd.isna(row.estimate):return 'Not estimable'
 if pd.isna(row.lower) or pd.isna(row.upper):return f'{row.estimate:.3f} (CI unavailable)'
 if style=='bracket':return f'{row.estimate:.3f}\n[{row.lower:.3f}, {row.upper:.3f}]'
 return f'{row.estimate:.3f} ({row.lower:.3f} to {row.upper:.3f})'

def wide(pop,met,keys,metrics=METRICS):
 result=pop.copy()
 if result.duplicated(keys).any():raise ValueError('Duplicate population key')
 for metric in metrics:
  rows=met[met.metric.eq(metric)].copy()
  if rows.duplicated(keys).any():raise ValueError('Duplicate metric key')
  rows[metric]=[display(r) for r in rows.itertuples()]
  result=result.merge(rows[keys+[metric]],on=keys,how='left',validate='one_to_one')
 return result

def run(root,task):
 inputs,outputs=TASKS[task]
 if any((root/p).exists() for p in outputs):raise FileExistsError('Existing report outputs preserved')
 if any(not (root/p).is_file() for p in inputs):raise FileNotFoundError('Required aggregate input absent')
 tables=[]
 if task=='primary':
  met=read(root,inputs[0]);rows=[]
  for ds,label in [('MIMIC_internal_test','MIMIC internal'),('eICU_external_validation','eICU external')]:
   for model,name in [('xgboost','XGBoost'),('elastic_net','Elastic net'),('random_forest','Random forest')]:
    row={'Dataset':label,'Model':name}
    for metric,heading in zip(METRICS,['AUROC','AUPRC','Brier','Intercept','Slope']):
     q=met[met.window_hours.eq(24)&met.dataset.eq(ds)&met.model.eq(model)&met.metric.eq(metric)]
     if len(q)!=1:raise ValueError('Missing or duplicate primary performance cell')
     row[heading]=display(next(q.itertuples()),'bracket')
    rows.append(row)
  tables=[pd.DataFrame(rows)]
 elif task=='parameters':
  rows=[]
  for hour,path in zip([24,6],inputs):
   for model,objects in json.loads((root/path).read_text(encoding='utf-8')).items():
    for obj in objects:
     for key,value in obj.items():rows.append(dict(window_hours=hour,model=model,parameter=key,value=value))
  tables=[pd.DataFrame(rows)]
 elif task=='complete_cases':
  pop,met,by=[read(root,p) for p in inputs];keys=['window_hours','dataset','subset']
  pop=pop[pop.population.isin(['MIMIC_internal_test','eICU_external_validation'])].rename(columns={'population':'dataset'})
  # The original presentation keys/counts come from the metric table, as in the audited reporting helper.
  base=met[['window_hours','dataset','subset','records','events','clusters','precision_warning']].drop_duplicates()
  table=wide(base,met,keys)
  assoc=by.pivot(index=['window_hours','population','raw_input'],columns='event',values='missing_fraction').rename(columns={0:'nonevent_missing_fraction',1:'event_missing_fraction'}).reset_index()
  assoc['difference_event_minus_nonevent']=assoc.event_missing_fraction-assoc.nonevent_missing_fraction
  tables=[table,assoc]
 elif task=='repeat':
  pop,met=[read(root,p) for p in inputs];tables=[wide(pop,met,['window_hours','dataset','analysis'])]
 elif task=='subgroups':
  pop,met=[read(root,p) for p in inputs];keys=['window_hours','dataset','family','subgroup']
  table=wide(pop,met,keys);table['event_rate_pct']=100*table.event_rate
  threshold=wide(pop[pop.dataset.eq('eICU_external_validation')&pop.status.eq('Reported')],met[met.dataset.eq('eICU_external_validation')],keys,['sensitivity','specificity','ppv','npv','flagged_fraction'])
  tables=[table,table[table.dataset.eq('eICU_external_validation')&table.family.eq('Phenotype')],table[table.dataset.eq('eICU_external_validation')&table.family.eq('Restricted evaluation')],threshold,pop[pop.status.eq('Reported')&pop.precision_warning.ne('No predefined small-sample flag (not proof of adequate precision)')]]
 elif task=='endpoint':
  sens,curves,risk=[read(root,p) for p in inputs]
  base=sens[['window_hours','dataset','analysis','records','events']].drop_duplicates()
  bins=risk[['window_hours','dataset','risk_bin','records','events']].drop_duplicates()
  tables=[wide(base,sens,['window_hours','dataset','analysis']),curves[curves.icu_day.eq(30)].copy(),wide(bins,risk,['window_hours','dataset','risk_bin'],['death_cif','alive_discharge_cif','mean_predicted'])]
 elif task=='units':
  parts=[]
  for db,path in zip(['mimic','eicu'],inputs):
   d=read(root,path);d.insert(0,'source',db);parts.append(d)
  tables=[pd.concat(parts,ignore_index=True)]
 if len(tables)!=len(outputs):raise RuntimeError('Report output mismatch')
 for p,d in zip(outputs,tables):
  if FORBIDDEN.intersection(d.columns):raise ValueError('Row-level field in report output')
  dest=root/p;dest.parent.mkdir(parents=True,exist_ok=True)
  with dest.open('x',encoding='utf-8',newline='') as stream:d.to_csv(stream,index=False)
 return [{'file':p,'rows':len(d)} for p,d in zip(outputs,tables)]

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--task',choices=TASKS,required=True);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
 if a.dry_run:print(json.dumps({'task':a.task,'inputs':TASKS[a.task][0],'outputs':TASKS[a.task][1],'data_read':False}))
 else:print(json.dumps(run(a.root.resolve(),a.task)))
