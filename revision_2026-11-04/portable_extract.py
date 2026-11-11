"""Read-only extraction adapter. Dry-run never reads patient files or starts psql."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,time
import pandas as pd

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def query_for(root,source,mode):
    name='03_mimic_cvst_v2.sql' if source=='mimic' else '02_eicu_revision.sql'
    path=root/'sql'/name
    text=path.read_text(encoding='utf-8').strip().rstrip(';')
    if mode=='explain':query='EXPLAIN (FORMAT JSON) '+text+';'
    elif mode=='cohort':
        marker=', feature_cohort AS ('
        if text.count(marker)!=1:raise ValueError('Unexpected cohort SQL structure')
        query='COPY ('+text.split(marker,1)[0]+' SELECT * FROM primary_cohort ORDER BY stay_id) TO STDOUT WITH CSV HEADER;'
    elif mode=='units':
        marker=', valid_long AS ('
        if text.count(marker)!=1:raise ValueError('Unexpected measurement SQL structure')
        tail=""" SELECT variable,event_source,coalesce(unit,'<missing>') AS unit,
 count(*) AS records,min(value) AS minimum,max(value) AS maximum,
 percentile_cont(0.5) WITHIN GROUP(ORDER BY value) AS median
 FROM raw_long GROUP BY variable,event_source,coalesce(unit,'<missing>')
 ORDER BY variable,event_source,unit"""
        query='COPY ('+text.split(marker,1)[0]+tail+') TO STDOUT WITH CSV HEADER;'
    else:query='COPY ('+text+') TO STDOUT WITH CSV HEADER;'
    return path,query

def validate_frame(df,contract,mode):
    if mode not in ('cohort','extract'):return {'source_records':len(df)}
    required={'stay_id','patient_group_id','landmark_6h_eligible','landmark_24h_eligible',
              'outcome_missing_flag','effective_icu_end_min','invalid_time_flag'}
    if not required.issubset(df.columns):raise ValueError('Required extraction fields absent')
    if not df.stay_id.is_unique or df.stay_id.isna().any():raise ValueError('Duplicate or missing stay key')
    if df.patient_group_id.isna().any():raise ValueError('Missing patient grouping key')
    flags=['landmark_6h_eligible','landmark_24h_eligible','outcome_missing_flag']
    if df[flags].isna().any().any() or not df[flags].isin([0,1]).all().all():raise ValueError('Invalid eligibility/outcome flags')
    if (df.landmark_24h_eligible>df.landmark_6h_eligible).any():raise ValueError('Landmarks are not nested')
    result={'source_records':len(df),'unique_patients':int(df.patient_group_id.nunique()),'windows':{}}
    for hour in [6,24]:
        elig=df[f'landmark_{hour}h_eligible'].eq(1);keep=elig & df.outcome_missing_flag.eq(0)
        y=df.loc[keep,f'death_{hour}h_to_30d_inhospital']
        ends=pd.to_numeric(df.loc[keep,'effective_icu_end_min'],errors='coerce')
        if ends.isna().any() or ends.le(hour*60).any():raise ValueError('ICU end does not exceed landmark')
        if y.isna().any() or not y.isin([0,1]).all():raise ValueError('Unknown/invalid target outcome')
        result['windows'][str(hour)]={'eligible_records':int(elig.sum()),'eligible_missing_outcomes':int((elig & ~keep).sum()),
          'analysis_records':int(keep.sum()),'unique_patients':int(df.loc[keep,'patient_group_id'].nunique()),'events':int(y.sum())}
    if mode=='extract':
        if {'frailty_comorbidity_count','stroke_group'}.intersection(df.columns):raise ValueError('Legacy primary predictor field present')
        for key,spec in contract['windows'].items():
            keep=df[spec['eligibility']].eq(1)&df.outcome_missing_flag.eq(0)
            if set(spec['predictors'])-set(df.columns):raise ValueError('Predictor allowlist field absent')
            for col in spec['predictors']:
                if '_missing_' not in col:continue
                value=col.replace('missing_6h','first_6h').replace('missing_24h','mean_24h')
                if not df.loc[keep,col].eq(df.loc[keep,value].isna().astype(int)).all():raise ValueError('Missingness indicator mismatch')
        result['feature_allowlist_checks']='PASS'
    result['invalid_time_records']=int(df.invalid_time_flag.astype(str).str.lower().isin(['1','t','true']).sum())
    return result

def connection_env(source,timeout_ms,planner_guard=False):
    env=os.environ.copy()
    service=env.get('STROKE_'+source.upper()+'_SERVICE')
    if service:env['PGSERVICE']=service
    if not env.get('PGSERVICE') and not env.get('PGDATABASE'):
        raise RuntimeError('Select an authorized database using PGSERVICE or PGDATABASE; keep credentials outside code')
    env.update(PGCONNECT_TIMEOUT='5',PGCLIENTENCODING='UTF8')
    # Deliberately replace, rather than inherit, PGOPTIONS so read-only cannot be disabled by an earlier setting.
    env['PGOPTIONS']=f'-c default_transaction_read_only=on -c statement_timeout={timeout_ms}'
    if planner_guard:env['PGOPTIONS']+=' -c enable_nestloop=off -c work_mem=64MB'
    return env

def run_psql(query,target,source,timeout_ms,psql='psql',planner_guard=False):
    part=target.with_suffix(target.suffix+'.partial')
    if target.exists() or part.exists():raise FileExistsError('Output/partial already exists; use a new run directory')
    env=connection_env(source,timeout_ms,planner_guard)
    target.parent.mkdir(parents=True,exist_ok=True)
    try:
        with part.open('xb') as stream:
            p=subprocess.run([psql,'-X','-w','-A','-t','-v','ON_ERROR_STOP=1'],input=query.encode('utf-8'),
              stdout=stream,stderr=subprocess.PIPE,env=env,timeout=timeout_ms/1000+30)
    except (subprocess.TimeoutExpired,OSError):
        raise RuntimeError('Extraction process failed or timed out; partial output retained privately') from None
    if p.returncode:
        # SQL errors may echo query fragments or selected IDs. Never print stderr or credentials.
        raise RuntimeError('Database extraction failed; partial output retained privately; no raw database error is printed')
    return part

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--source',choices=['mimic','eicu'],required=True)
    parser.add_argument('--mode',choices=['explain','cohort','extract','units'],required=True)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--planner-guard',action='store_true')
    parser.add_argument('--psql',default=os.environ.get('PSQL','psql'))
    args=parser.parse_args(argv);root=args.root.resolve()
    path,query=query_for(root,args.source,args.mode)
    if args.dry_run:
        print(json.dumps({'source':args.source,'mode':args.mode,'sql_file':path.name,'sql_sha256':sha(path),
          'query_sha256':hashlib.sha256(query.encode()).hexdigest(),'database_contacted':False,'patient_files_read':False}));return
    suffix={'explain':'explain.json','cohort':'cohort.csv','extract':'analysis_revision.csv','units':'measurement_units.csv'}[args.mode]
    target=root/('qc' if args.mode in ['explain','units'] else 'private_data')/f'{args.source}_{suffix}'
    report=root/'qc'/f'{args.source}_{args.mode}_summary.json'
    if report.exists():raise FileExistsError('Summary already exists; use a new run directory')
    start=time.monotonic();part=run_psql(query,target,args.source,900000,args.psql,args.planner_guard)
    if args.mode=='explain':json.loads(part.read_text(encoding='utf-8'));summary={'sql_parse_and_plan':'PASS'}
    else:
        frame=pd.read_csv(part)
        contract=json.loads((root/'model_feature_contract.json').read_text(encoding='utf-8'))
        summary=validate_frame(frame,contract,args.mode)
    summary.update(source=args.source,mode=args.mode,sql_sha256=sha(path),output_sha256=sha(part),
      output_file=str(target.relative_to(root)),elapsed_seconds=round(time.monotonic()-start,2),database_read_only=True)
    part.rename(target);report.parent.mkdir(parents=True,exist_ok=True)
    with report.open('x',encoding='utf-8') as stream:json.dump(summary,stream,indent=2)
    print(json.dumps(summary))

if __name__=='__main__':main()
