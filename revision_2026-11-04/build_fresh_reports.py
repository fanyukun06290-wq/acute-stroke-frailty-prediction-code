"""Generate four fresh-run reports locally. Raw records never appear on stdout."""
from pathlib import Path
import argparse,json,hashlib,re
import numpy as np
import pandas as pd
OUTPUTS={'phenotypes':'phenotype_v2_results/phenotype_counts_cvst_v2.csv','bounds':'codebook_results/final_input_bounds_checks.csv','time_flags':'endpoint_results/time_flag_case_mix.csv','index':'phenotype_v2_results/authoritative_artifact_index.csv'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def binary(series):
 n=pd.to_numeric(series,errors='raise')
 if n.isna().any() or not n.isin([0,1]).all():raise ValueError('Invalid binary field')
 return n.astype(int)
def cohort(root,db):
 d=pd.read_csv(root/'private_data'/f'{db}_cohort.csv')
 if d.stay_id.isna().any() or not d.stay_id.is_unique or d.patient_group_id.isna().any():raise ValueError('Invalid cohort key')
 for h in [6,24]:
  eligible=binary(d[f'landmark_{h}h_eligible']);missing=binary(d.outcome_missing_flag)
  if (eligible.eq(1)&d.effective_icu_end_min.le(h*60)).any():raise ValueError('ICU landmark boundary violation')
  binary(d.loc[eligible.eq(1)&missing.eq(0),f'death_{h}h_to_30d_inhospital'])
 if (binary(d.landmark_24h_eligible)>binary(d.landmark_6h_eligible)).any():raise ValueError('Non-nested eligibility')
 return d
def phenotype_report(root):
 rows=[];inputs=[]
 for db in ['mimic','eicu']:
  inputs.append(f'private_data/{db}_cohort.csv');d=cohort(root,db)
  if d.audit_phenotype.isna().any():raise ValueError('Missing phenotype label')
  for h in [0,6,24]:
   a=d if h==0 else d[binary(d[f'landmark_{h}h_eligible']).eq(1)&binary(d.outcome_missing_flag).eq(0)]
   for label,g in a.groupby('audit_phenotype',sort=True):
    rows.append(dict(source=db,window_hours=h,phenotype=label,records=len(g),unique_patients=g.patient_group_id.nunique(),events=None if h==0 else int(binary(g[f'death_{h}h_to_30d_inhospital']).sum())))
 return pd.DataFrame(rows),inputs
def extract_bounds(text):
 block=text.split('bounds(variable,lower_bound,upper_bound) AS (VALUES ',1)[1].split('), lab_map',1)[0]
 pairs=re.findall(r"\('([^']+)',([\d.]+),([\d.]+)\)",block+')')
 result={v:(float(lo),float(hi)) for v,lo,hi in pairs}
 if len(result)!=18:raise ValueError('Expected 18 SQL measurement bounds')
 return result
def bounds_report(root):
 inputs=['sql/03_mimic_cvst_v2.sql','sql/02_eicu_revision.sql','phenotype_v2_results/measurement_dictionary.csv','model_feature_contract.json'];rows=[]
 maps=[extract_bounds((root/p).read_text(encoding='utf-8')) for p in inputs[:2]]
 if maps[0]!=maps[1]:raise ValueError('Cross-database bounds differ')
 dictionary=pd.read_csv(root/inputs[2]).set_index('variable')
 if not dictionary.index.is_unique or set(dictionary.index)!=set(maps[0]):raise ValueError('Measurement dictionary variable mismatch')
 for v,(lo,hi) in maps[0].items():
  if (float(dictionary.loc[v,'lower_inclusive']),float(dictionary.loc[v,'upper_inclusive']))!=(lo,hi):raise ValueError('SQL/dictionary bound mismatch')
 contract=json.loads((root/inputs[3]).read_text(encoding='utf-8'))
 for db in ['mimic','eicu']:
  path=f'private_data/{db}_analysis_revision.csv';inputs.append(path);d=pd.read_csv(root/path)
  if not d.stay_id.is_unique or d.stay_id.isna().any():raise ValueError('Invalid analysis keys')
  for h in [6,24]:
   a=d[binary(d[f'landmark_{h}h_eligible']).eq(1)&binary(d.outcome_missing_flag).eq(0)]
   if len(contract['windows'][str(h)]['predictors'])!=38:raise ValueError('Predictor contract changed')
   for v,(lo,hi) in maps[0].items():
    col=f'{v}_{"first_6h" if h==6 else "mean_24h"}';x=pd.to_numeric(a[col],errors='raise');z=x.dropna()
    bad=int(((z<lo-1e-9)|(z>hi+1e-9)|~np.isfinite(z)).sum())
    if bad or not x.isna().astype(int).eq(binary(a[f'{v}_missing_{h}h'])).all():raise ValueError('Measurement bounds or missingness flag violation')
    rows.append(dict(source=db,window_hours=h,variable=v,analysis_records=len(a),available_records=len(z),outside_bounds=bad,minimum=z.min(),maximum=z.max()))
 return pd.DataFrame(rows),inputs
def time_flag_report(root):
 rows=[];inputs=[]
 for h in [24,6]:
  path=f'private_predictions/predictions_{h}h.csv';inputs.append(path);pred=pd.read_csv(root/path)
  splitpath=f'private_splits/patients_{h}h.csv';inputs.append(splitpath);split=pd.read_csv(root/splitpath)
  if not split.group_id.is_unique:raise ValueError('Duplicate split patient')
  for db,ds in [('mimic','MIMIC_internal_test'),('eicu','eICU_external_validation')]:
   path=f'private_data/{db}_cohort.csv';inputs.append(path);raw=cohort(root,db)
   eligible=binary(raw[f'landmark_{h}h_eligible']).eq(1)&binary(raw.outcome_missing_flag).eq(0)
   if db=='mimic':eligible &=raw.patient_group_id.isin(split.loc[split.partition.eq('test'),'group_id'])
   p=pred[pred.model.eq('xgboost')&pred.dataset.eq(ds)]
   if p.stay_id.isna().any() or not p.stay_id.is_unique or set(p.stay_id)!=set(raw.loc[eligible,'stay_id']):raise ValueError('Prediction population mismatch')
   z=p[['stay_id','probability','truth']].merge(raw,on='stay_id',validate='one_to_one')
   probability=pd.to_numeric(z.probability,errors='raise').to_numpy()
   if not np.isfinite(probability).all() or ((probability<0)|(probability>1)).any():raise ValueError('Invalid prediction probability')
   if not z.truth.isin(['yes','no']).all():raise ValueError('Unknown prediction truth')
   yy=binary(z[f'death_{h}h_to_30d_inhospital']).to_numpy()
   expired=binary(z.hospital_expire_flag).eq(1).to_numpy() if db=='mimic' else z.hospital_status.eq('expired').to_numpy()
   observed=(expired&(z.death_offset_min.to_numpy()<43200)).astype(int)
   if not np.array_equal(yy,observed) or not np.array_equal(yy,z.truth.eq('yes').astype(int)):raise ValueError('Endpoint truth mismatch')
   terminal=np.where(expired,z.death_offset_min.to_numpy(),z.hospital_discharge_offset_min.to_numpy())
   if not np.isfinite(terminal).all() or (terminal<=h*60).any():raise ValueError('Terminal event at/before landmark')
   if db=='mimic':delay=(pd.to_datetime(z.icu_outtime)-pd.to_datetime(z.dischtime)).dt.total_seconds().to_numpy()/3600
   else:delay=(pd.to_numeric(z.unitdischargeoffset)-pd.to_numeric(z.hospitaldischargeoffset)).to_numpy()/60
   flag=binary(z.icu_end_after_hospital_flag).to_numpy()
   if not np.isfinite(delay).all() or not np.array_equal(flag,(delay>0).astype(int)):raise ValueError('Timing flag or discharge time mismatch')
   for value in [0,1]:
    mask=flag==value;values=delay[mask];n=int(mask.sum());events=int(yy[mask].sum())
    rows.append(dict(window_hours=h,dataset=ds,icu_end_after_hospital_flag=value,records=n,events=events,event_rate=float(yy[mask].mean()) if n else np.nan,mean_predicted=float(probability[mask].mean()) if n else np.nan,fraction_of_all_events=float(events/yy.sum()) if yy.sum() else np.nan,icu_end_minus_hospital_end_median_hours=float(np.median(values)) if n else np.nan,q25_hours=float(np.quantile(values,.25)) if n else np.nan,q75_hours=float(np.quantile(values,.75)) if n else np.nan))
 return pd.DataFrame(rows),list(dict.fromkeys(inputs))
def index_report(root):
 names=['sql/03_mimic_cvst_v2.sql','sql/02_eicu_revision.sql','model_feature_contract.json','phenotype_v2_results/phenotype_codes_cvst_v2.csv','phenotype_v2_results/measurement_dictionary.csv',OUTPUTS['phenotypes'],OUTPUTS['bounds'],OUTPUTS['time_flags']]
 if any(not(root/n).is_file() for n in names):raise FileNotFoundError('Fresh index prerequisite absent')
 return pd.DataFrame([dict(material=Path(n).stem,relative_path=n,version_instruction='Current local file fingerprint; not independent scientific validation',sha256=sha(root/n)) for n in names]),names
FUNCTIONS={'phenotypes':phenotype_report,'bounds':bounds_report,'time_flags':time_flag_report,'index':index_report}
def run(root,task):
 target=root/OUTPUTS[task];receipt=root/'qc'/f'fresh_report_{task}.json'
 if target.exists() or receipt.exists():raise FileExistsError('Existing output or receipt preserved')
 table,inputs=FUNCTIONS[task](root)
 target.parent.mkdir(parents=True,exist_ok=True)
 with target.open('x',encoding='utf-8',newline='') as f:table.to_csv(f,index=False)
 record=dict(task=task,scope='fresh_report_not_full_reproduction',inputs={n:sha(root/n) for n in inputs},output={OUTPUTS[task]:sha(target)},builder_sha256=sha(Path(__file__)))
 receipt.parent.mkdir(exist_ok=True)
 with receipt.open('x',encoding='utf-8') as f:json.dump(record,f,indent=2)
 return dict(task=task,rows=len(table),output=OUTPUTS[task])
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('task',choices=FUNCTIONS);p.add_argument('--plan',action='store_true');a=p.parse_args()
 if a.plan:print(json.dumps(dict(task=a.task,output=OUTPUTS[a.task],files_read=False)))
 else:print(json.dumps(run(Path(__file__).resolve().parent,a.task)))
