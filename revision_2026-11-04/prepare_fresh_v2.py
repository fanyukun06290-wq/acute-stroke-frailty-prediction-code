"""Fresh-v2 compatibility copies, with no historical overlay or probability edit."""
from pathlib import Path
import argparse,json,shutil
import pandas as pd
from portable_extract import sha

def prepare(root,with_predictions=False):
    expected=sha(root/'sql/03_mimic_cvst_v2.sql')
    tasks=[]
    for stem,mode in [('mimic_cohort','cohort'),('mimic_analysis_revision','extract')]:
        source=root/'private_data'/f'{stem}.csv'
        report=json.loads((root/'qc'/f'mimic_{mode}_summary.json').read_text(encoding='utf-8'))
        if report.get('sql_sha256')!=expected or report.get('output_sha256')!=sha(source):
            raise ValueError('Fresh route requires verified current-v2 extraction hashes; historical exports are not accepted')
        target=root/'private_data'/('mimic_cohort_cvst_v2.csv' if mode=='cohort' else 'mimic_analysis_cvst_v2.csv')
        tasks.append((source,target))
    if with_predictions:
        source_frames={db:pd.read_csv(root/'private_data'/f'{db}_analysis_revision.csv',dtype=str,keep_default_na=False) for db in ['mimic','eicu']}
        for db,frame in source_frames.items():
            if not frame.stay_id.is_unique:raise ValueError('Duplicate source key')
        for hour in [6,24]:
            source=root/'private_predictions'/f'predictions_{hour}h.csv'
            predictions=pd.read_csv(source,dtype=str,keep_default_na=False)
            for db,dataset in [('mimic','MIMIC_internal_test'),('eicu','eICU_external_validation')]:
                frame=source_frames[db].set_index('stay_id')
                p=predictions[predictions.dataset.eq(dataset)]
                if p.empty:raise ValueError('Expected evaluation dataset absent')
                aligned=frame.reindex(p.stay_id)
                if aligned.audit_phenotype.isna().any() or not (p.phenotype.to_numpy()==aligned.audit_phenotype.to_numpy()).all():
                    raise ValueError('Prediction labels do not already match fresh-v2 extraction; no automatic relabeling is performed')
            tasks.append((source,root/'private_predictions'/f'predictions_{hour}h_cvst_v2.csv'))
    if any(target.exists() for _,target in tasks):raise FileExistsError('A compatibility output already exists; nothing overwritten')
    destination=root/'phenotype_v2_results/fresh_input_identity.json'
    if destination.exists():raise FileExistsError('Fresh-input receipt already exists')
    checks=[]
    for source,target in tasks:
        with source.open('rb') as a,target.open('xb') as b:shutil.copyfileobj(a,b)
        if sha(source)!=sha(target):raise RuntimeError('Compatibility copy hash mismatch')
        checks.append({'source':str(source.relative_to(root)),'output':str(target.relative_to(root)),'sha256':sha(target),'byte_identical':True})
    receipt={'route':'fresh_current_v2','historical_18_record_overlay_applied':False,'probabilities_modified':False,'files':checks}
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('x',encoding='utf-8') as stream:json.dump(receipt,stream,indent=2)
    return receipt

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--with-predictions',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args(argv)
    if args.dry_run:
        print(json.dumps({'route':'fresh_current_v2','requires_current_sql_hash':True,'historical_overlay':False,'database_contacted':False,'patient_files_read':False}));return
    print(json.dumps(prepare(args.root.resolve(),args.with_predictions)))

if __name__=='__main__':main()
