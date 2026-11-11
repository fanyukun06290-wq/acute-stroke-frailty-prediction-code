"""Validate an existing local manual export; do not connect to a database."""
from pathlib import Path
import argparse,json
import pandas as pd
from portable_extract import validate_frame,sha

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--source',choices=['mimic','eicu'],required=True)
    p.add_argument('--kind',choices=['cohort','extract'],required=True)
    p.add_argument('--attest-current-sql',action='store_true',help='Operator attests that the selected current SQL produced this export; this is not independent verification')
    p.add_argument('--dry-run',action='store_true');args=p.parse_args(argv);root=args.root.resolve()
    sql=root/'sql'/('03_mimic_cvst_v2.sql' if args.source=='mimic' else '02_eicu_revision.sql')
    if args.dry_run:
        print(json.dumps({'route':'manual_export_registration','sql_file':sql.name,'sql_sha256':sha(sql),'patient_files_read':False,'database_contacted':False}));return
    if not args.attest_current_sql:raise ValueError('Manual export provenance requires an explicit operator attestation')
    name=f'{args.source}_cohort.csv' if args.kind=='cohort' else f'{args.source}_analysis_revision.csv'
    source=root/'private_data'/name;report=root/'qc'/f'{args.source}_{args.kind}_summary.json'
    if report.exists():raise FileExistsError('Existing receipt preserved; use a separate run directory')
    contract=json.loads((root/'model_feature_contract.json').read_text(encoding='utf-8'))
    summary=validate_frame(pd.read_csv(source),contract,args.kind)
    summary.update(source=args.source,mode=args.kind,sql_sha256=sha(sql),output_sha256=sha(source),output_file=str(source.relative_to(root)),
      database_contacted=False,provenance='Operator-attested manual export; query execution not independently verified',sql_execution_independently_verified=False)
    report.parent.mkdir(parents=True,exist_ok=True)
    with report.open('x',encoding='utf-8') as stream:json.dump(summary,stream,indent=2)
    print(json.dumps(summary))

if __name__=='__main__':main()
