"""Data-free adapter and synthetic fresh-route checks. No real database is used."""
from pathlib import Path
import argparse,contextlib,copy,importlib.util,io,json,os,shutil,subprocess,sys,unittest
from unittest.mock import patch
import pandas as pd

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def run(repository,work,rscript=None):
    project=repository/'revision_2026-11-04';sys.path.insert(0,str(project))
    import portable_extract as ex
    import portable_gcs_extract as gx
    import prepare_fresh_v2 as fresh
    contract=json.loads((project/'model_feature_contract.json').read_text(encoding='utf-8'))
    work.mkdir(parents=True,exist_ok=True);results=[]
    def case(name,fn):
        fn();results.append({'check':name,'passed':True,'data_scope':'Artificial fixtures only; no real database'})
    def expect(kind,fn):
        try:fn()
        except kind:return
        raise AssertionError('Expected controlled rejection was not raised')
    def frame(base=1000,n=4):
        d=pd.DataFrame({'stay_id':range(base,base+n),'patient_group_id':range(base+100,base+100+n),
          'hospitalid':[1+(i%2) for i in range(n)],'landmark_6h_eligible':[1]*n,'landmark_24h_eligible':[1]*n,
          'outcome_missing_flag':[0]*n,'effective_icu_end_min':[2000]*n,'invalid_time_flag':[False]*n,
          'death_6h_to_30d_inhospital':[i%2 for i in range(n)],'death_24h_to_30d_inhospital':[i%2 for i in range(n)],
          'audit_phenotype':['AIS_STRICT']*n})
        for spec in contract['windows'].values():
            for key in spec['predictors']:
                if key=='sex':d[key]=['F' if i%2 else 'M' for i in range(n)]
                elif key=='age_years':d[key]=[60+i for i in range(n)]
                elif '_missing_' in key:d[key]=0
                else:d[key]=1.
        return d
    def dry_runs():
        with patch.object(ex.subprocess,'run',side_effect=AssertionError('No process is permitted')),patch.object(pd,'read_csv',side_effect=AssertionError('No patient file reads permitted')):
            for source in ['mimic','eicu']:
                for mode in ['explain','cohort','extract','units']:
                    out=io.StringIO()
                    with contextlib.redirect_stdout(out):ex.main(['--root',str(project),'--source',source,'--mode',mode,'--dry-run'])
                    report=json.loads(out.getvalue());assert not report['database_contacted'] and not report['patient_files_read']
                out=io.StringIO()
                with contextlib.redirect_stdout(out):gx.main(['--root',str(project),'--source',source,'--dry-run'])
                assert not json.loads(out.getvalue())['patient_files_read']
    case('Ten source/mode dry-runs cannot start a process or read a CSV',dry_runs)
    case('Both eligible landmarks and allowed inputs validate',lambda:ex.validate_frame(frame(),contract,'extract'))
    def boundary(hour):
        d=frame();d.loc[0,'effective_icu_end_min']=hour*60
        if hour==24:d.loc[0,'landmark_6h_eligible']=0;d.loc[0,'landmark_24h_eligible']=0;d.loc[1,'effective_icu_end_min']=1440
        expect(ValueError,lambda:ex.validate_frame(d,contract,'extract'))
    case('ICU end exactly at six-hour landmark is rejected',lambda:boundary(6))
    case('ICU end exactly at 24-hour landmark is rejected',lambda:boundary(24))
    def duplicate():
        d=frame();d.loc[1,'stay_id']=d.loc[0,'stay_id'];expect(ValueError,lambda:ex.validate_frame(d,contract,'extract'))
    case('Duplicate stay key is rejected',duplicate)
    def missingness():
        d=frame();key=next(k for k in contract['windows']['6']['predictors'] if '_missing_' in k);d.loc[0,key]=1
        expect(ValueError,lambda:ex.validate_frame(d,contract,'extract'))
    case('Inconsistent missingness indicator is rejected',missingness)
    def uncertain():
        d=frame();d.loc[0,'death_6h_to_30d_inhospital']=None
        expect(ValueError,lambda:ex.validate_frame(d,contract,'extract'))
        d.loc[0,'outcome_missing_flag']=1;ex.validate_frame(d,contract,'extract')
    case('Unknown target cannot enter known-outcome cohort',uncertain)
    def subprocess_contract():
        target=work/'mock_output.csv';seen={}
        def mock(args,**kwargs):
            seen.update(args=args,env=kwargs['env']);kwargs['stdout'].write(b'value\n1\n');return subprocess.CompletedProcess(args,0,b'',b'')
        with patch.dict(os.environ,{'PGSERVICE':'synthetic-service','PGOPTIONS':'-c default_transaction_read_only=off'},clear=True),patch.object(ex.subprocess,'run',side_effect=mock):
            part=ex.run_psql('SELECT 1;',target,'mimic',1000,psql='synthetic-psql')
        assert 'default_transaction_read_only=on' in seen['env']['PGOPTIONS'] and 'off' not in seen['env']['PGOPTIONS']
        assert '-w' in seen['args'] and '-h' not in seen['args'] and '-U' not in seen['args'] and part.exists()
        expect(FileExistsError,lambda:ex.run_psql('SELECT 1;',target,'mimic',1000))
    case('Read-only connection overrides inherited options and refuses existing partial',subprocess_contract)
    def no_error_echo():
        secret=b'SYNTHETIC_RESTRICTED_ERROR_SENTINEL'
        def mock(args,**kwargs):return subprocess.CompletedProcess(args,1,b'',secret)
        with patch.dict(os.environ,{'PGSERVICE':'synthetic-service'},clear=True),patch.object(ex.subprocess,'run',side_effect=mock):
            try:ex.run_psql('SELECT 1;',work/'error_output.csv','mimic',1000)
            except RuntimeError as e:assert secret.decode() not in str(e)
            else:raise AssertionError('Expected extraction failure')
    case('Raw database error text is not echoed',no_error_echo)
    gcs=load('synthetic_gcs_builder',project/'13_build_gcs_features.py')
    def gcs_rules():assert gcs.synthetic_checks()==15
    case('Fifteen inherited GCS timing and assessability scenarios',gcs_rules)
    fixture=work/'fresh_repository'/'revision_2026-11-04'
    shutil.copytree(project,fixture)
    for name in ['private_data','private_predictions','private_splits','private_cache','qc']:(fixture/name).mkdir(exist_ok=True)
    for db,base in [('mimic',1000),('eicu',2000)]:
        d=frame(base)
        for mode,name in [('cohort',f'{db}_cohort.csv'),('extract',f'{db}_analysis_revision.csv')]:
            p=fixture/'private_data'/name;d.to_csv(p,index=False)
            sql=fixture/'sql'/('03_mimic_cvst_v2.sql' if db=='mimic' else '02_eicu_revision.sql')
            (fixture/'qc'/f'{db}_{mode}_summary.json').write_text(json.dumps({'sql_sha256':ex.sha(sql),'output_sha256':ex.sha(p)}))
    for hour in [6,24]:
        parts=[]
        for db,dataset,base in [('mimic','MIMIC_internal_test',1000),('eicu','eICU_external_validation',2000)]:
            raw=frame(base)
            parts.append(pd.DataFrame({'model':['xgboost']*4,'dataset':[dataset]*4,'window_hours':[hour]*4,'stay_id':raw.stay_id,
             'group_id':raw.patient_group_id,'hospitalid':raw.hospitalid,'phenotype':raw.audit_phenotype,
             'truth':['yes' if i%2 else 'no' for i in range(4)],'probability':[.3,.6,.7,.4]}))
        pd.concat(parts,ignore_index=True).to_csv(fixture/'private_predictions'/f'predictions_{hour}h.csv',index=False)
        pd.DataFrame({'group_id':frame(1000).patient_group_id,'partition':['test']*4}).to_csv(fixture/'private_splits'/f'patients_{hour}h.csv',index=False)
    def stale_hash():
        path=fixture/'qc/mimic_extract_summary.json';saved=path.read_bytes();info=json.loads(saved);info['sql_sha256']='historical-version'
        path.write_text(json.dumps(info));expect(ValueError,lambda:fresh.prepare(fixture,True));path.write_bytes(saved)
    case('Fresh route rejects a historical SQL receipt',stale_hash)
    def alias():
        receipt=fresh.prepare(fixture,True);assert len(receipt['files'])==4 and all(x['byte_identical'] for x in receipt['files'])
        assert not receipt['historical_18_record_overlay_applied'] and not receipt['probabilities_modified']
        expect(FileExistsError,lambda:fresh.prepare(fixture,True))
    case('Fresh aliases preserve exact bytes and refuse rerun overwrite',alias)
    if rscript:
        def fresh_r():
            assert not list((fixture/'private_cache').iterdir())
            # Outputs are artificial and groups are deliberately below reporting size.
            process=subprocess.run([str(rscript),str(fixture/'R/run_subgroup_fresh_v2.R'),str(fixture)],cwd=fixture.parent,
              stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
            if process.returncode:
                raise RuntimeError('Synthetic R entry failed; inspect the test environment; no real patient data were used')
            files=list((fixture/'private_cache').glob('subgroup_fresh_v2_*.rds'));assert len(files)==4
            assert not list((fixture/'private_cache').glob('subgroup_6h_*.rds'))
            table=pd.read_csv(fixture/'phenotype_v2_results/subgroup_populations.csv')
            assert not table.empty and (fixture/'phenotype_v2_results/run_complete.json').exists()
        case('Fresh R subgroup entry starts without historical caches (small synthetic smoke test)',fresh_r)
    (work/'SYNTHETIC_ONLY.json').write_text(json.dumps({'artificial_data':True,'real_patient_data':False,'checks':len(results)}))
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repository',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--rscript',type=Path)
    args=p.parse_args();results=run(args.repository.resolve(),args.work.resolve(),args.rscript)
    for result in results:print(json.dumps(result))
