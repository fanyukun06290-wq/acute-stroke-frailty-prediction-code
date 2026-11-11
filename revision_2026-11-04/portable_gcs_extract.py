"""GCS adapter using unchanged source SQL fragments; dry-run is data-free."""
from pathlib import Path
import argparse,json,hashlib,os
import pandas as pd
from portable_extract import run_psql,sha

def query_for(root,source,selected_sql):
    body=(root/'sql'/f'gcs_{source}_fragment.sql').read_text(encoding='utf-8').rstrip().rstrip(';')
    return 'COPY (WITH selected(stay_id) AS ('+selected_sql+')'+body+') TO STDOUT WITH CSV HEADER;'

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--source',choices=['mimic','eicu'],required=True)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--psql',default=os.environ.get('PSQL','psql'))
    args=parser.parse_args(argv);root=args.root.resolve();source=args.source
    fragment=root/'sql'/f'gcs_{source}_fragment.sql'
    if args.dry_run:
        query=query_for(root,source,'SELECT NULL::bigint AS stay_id WHERE FALSE')
        print(json.dumps({'source':source,'fragment_sha256':sha(fragment),'empty_selection_query_sha256':hashlib.sha256(query.encode()).hexdigest(),
          'database_contacted':False,'patient_files_read':False}));return
    target=root/'private_data'/f'{source}_gcs_events.csv';report=root/'qc'/f'{source}_gcs_extract.json'
    if report.exists():raise FileExistsError('GCS summary already exists')
    cohort=pd.read_csv(root/'private_data'/f'{source}_cohort.csv')
    selected=cohort[cohort.landmark_6h_eligible.eq(1)&cohort.outcome_missing_flag.eq(0)]
    if selected.empty or selected.stay_id.isna().any() or not selected.stay_id.is_unique:raise ValueError('Invalid/empty eligible stay selection')
    numeric=pd.to_numeric(selected.stay_id,errors='coerce')
    if numeric.isna().any() or not numeric.eq(numeric.astype('int64')).all():raise ValueError('Non-integer stay key')
    values='VALUES '+','.join(f'({int(x)})' for x in numeric)
    query=query_for(root,source,values)
    part=run_psql(query,target,source,180000,args.psql,planner_guard=True)
    events=pd.read_csv(part)
    if not set(events.stay_id).issubset(set(numeric)):raise ValueError('GCS export contains an unselected stay')
    if not events.observation_offset.ge(0).all() or not events.observation_offset.lt(1440).all():raise ValueError('GCS observation outside extraction window')
    summary={'source':source,'rows':len(events),'stays':int(events.stay_id.nunique()),'output_sha256':sha(part),
      'fragment_sha256':sha(fragment),'database_read_only':True,'raw_includes_records_unavailable_at_landmark':True}
    part.rename(target);report.parent.mkdir(parents=True,exist_ok=True)
    with report.open('x',encoding='utf-8') as stream:json.dump(summary,stream,indent=2)
    print(json.dumps(summary))

if __name__=='__main__':main()
