"""Independent subgroup counts/point estimates and publication tables, aggregate output only."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'phenotype_v2_results'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def collapse(x):
    allowed={'AIS_STRICT','ICH_STRICT','SAH','MIXED_PRIMARY_ACUTE','MIXED_WITH_AMBIGUOUS_SECONDARY','MIXED_WITH_CVST'}
    assert set(x)<=allowed
    return x.where(x.isin(['AIS_STRICT','ICH_STRICT','SAH']),'Mixed/ambiguous')
def auc(y,p):
    n1=y.sum();n0=len(y)-n1
    return (pd.Series(p).rank().to_numpy()[y==1].sum()-n1*(n1+1)/2)/(n1*n0) if n1*n0 else np.nan
def auprc(y,p):
    order=np.argsort(-p,kind='stable');yp=y[order];pp=p[order]
    last=np.r_[np.where(np.diff(pp)!=0)[0],len(pp)-1]
    tp=np.cumsum(yp)[last];fp=np.cumsum(1-yp)[last]
    r=np.r_[0,tp/y.sum()];precision=np.r_[1,tp/(tp+fp)]
    return np.sum(np.diff(r)*(precision[1:]+precision[:-1])/2)
def main():
    protected=json.loads((ROOT/'qc/shap_domain_artifact_manifest.json').read_text(encoding='utf-8'))
    for r in protected:assert sha(ROOT/r['relative_path'])==r['sha256']
    pop=pd.read_csv(OUT/'subgroup_populations.csv');met=pd.read_csv(OUT/'subgroup_metric_intervals.csv');mix=pd.read_csv(OUT/'phenotype_case_mix.csv')
    assert len(pop)==48 and len(mix)==32 and len(met)==748,(len(pop),len(mix),len(met))
    finite=met.lower.notna() & met.upper.notna()
    assert (met.loc[finite,'lower']<=met.loc[finite,'upper']).all()
    assert not met.duplicated(['window_hours','dataset','family','subgroup','metric']).any()
    assert len(pop[pop.status.ne('Reported')])==4
    checks=0
    for h in [24,6]:
        raw={db:pd.read_csv(ROOT/'private_data'/('mimic_analysis_cvst_v2.csv' if db=='mimic' else 'eicu_analysis_revision.csv'),usecols=['stay_id','patient_group_id','hospitalid','age_years','sex','audit_phenotype',f'landmark_{h}h_eligible','outcome_missing_flag',f'death_{h}h_to_30d_inhospital']) if db=='eicu' else
             pd.read_csv(ROOT/'private_data'/('mimic_analysis_cvst_v2.csv' if db=='mimic' else 'eicu_analysis_revision.csv'),usecols=['stay_id','patient_group_id','age_years','sex','audit_phenotype',f'landmark_{h}h_eligible','outcome_missing_flag',f'death_{h}h_to_30d_inhospital']) for db in ['mimic','eicu']}
        pred=pd.read_csv(ROOT/'private_predictions'/f'predictions_{h}h_cvst_v2.csv');pred=pred[pred.model.eq('xgboost')]
        for ds in ['MIMIC_internal_test','eICU_external_validation']:
            db='mimic' if ds=='MIMIC_internal_test' else 'eicu'
            r=raw[db];z=pred[pred.dataset.eq(ds)].merge(r,on='stay_id',validate='one_to_one',suffixes=('','_raw'))
            assert z[f'landmark_{h}h_eligible'].eq(1).all() and z.outcome_missing_flag.eq(0).all()
            assert z.group_id.eq(z.patient_group_id).all()
            ph=collapse(z.audit_phenotype);sex=z.sex.fillna('Unknown').where(z.sex.isin(['F','M']),'Unknown')
            masks={('Overall','All records'):np.ones(len(z),bool),('Age','<65'):z.age_years.lt(65),('Age','>=65'):z.age_years.ge(65),
                ('Sex','Female'):sex.eq('F'),('Sex','Male'):sex.eq('M'),('Sex','Unknown'):sex.eq('Unknown'),
                **{('Phenotype',k):ph.eq(k) for k in ['AIS_STRICT','ICH_STRICT','SAH','Mixed/ambiguous']},
                ('Restricted evaluation','Non-mixed AIS/ICH/SAH'):ph.ne('Mixed/ambiguous'),('Restricted evaluation','ICH/SAH combined'):ph.isin(['ICH_STRICT','SAH'])}
            for (family,label),mask in masks.items():
                zz=z[mask];row=pop[(pop.window_hours==h)&pop.dataset.eq(ds)&pop.family.eq(family)&pop.subgroup.eq(label)].iloc[0]
                if row.status=='Small cell suppressed':
                    assert 0<len(zz)<10 and pd.isna(row.records) and pd.isna(row.events)
                    continue
                assert row.records==len(zz) and row.unique_patients==zz.patient_group_id.nunique()
                y=zz.truth.eq('yes').to_numpy(dtype=int);p=zz.probability.to_numpy();e=int(y.sum());n=len(y)
                assert row.events==e and row.non_events==n-e
                assert np.array_equal(y,zz[f'death_{h}h_to_30d_inhospital'].to_numpy())
                if not n:continue
                positive=p>=.10;tp=int((y*positive).sum());fp=int(((1-y)*positive).sum());tn=int(((1-y)*~positive).sum());fn=int((y*~positive).sum())
                assert (row.tp,row.fp,row.tn,row.fn)==(tp,fp,tn,fn)
                m=met[(met.window_hours==h)&met.dataset.eq(ds)&met.family.eq(family)&met.subgroup.eq(label)].set_index('metric')
                expected={'auroc':auc(y,p),'auprc':auprc(y,p),'brier_score':np.mean((y-p)**2),'mean_predicted':p.mean(),'observed_rate':y.mean(),
                    'risk_bias':p.mean()-y.mean(),'sensitivity':tp/(tp+fn),'specificity':tn/(tn+fp),'ppv':tp/(tp+fp),'npv':tn/(tn+fn),'flagged_fraction':positive.mean()}
                for name,value in expected.items():assert np.isclose(m.loc[name,'estimate'],value,atol=1e-10), (h,ds,label,name)
                checks+=1
            # Age/sex/phenotype are separate exhaustive record partitions; patient counts can overlap.
            pp=pop[(pop.window_hours==h)&pop.dataset.eq(ds)]
            for family in ['Age','Sex','Phenotype']:
                g=pp[pp.family.eq(family)]
                if not g.status.eq('Small cell suppressed').any():
                    assert g.records.sum()==len(z) and g.events.sum()==z.truth.eq('yes').sum()
                else:
                    assert sum(np.asarray(mask,dtype=int) for (f,l),mask in masks.items() if f==family).tolist()==[1]*len(z)
    # Identical full-cohort weights must reproduce the previously frozen uncertainty results.
    baseline=pd.read_csv(ROOT/'validation_results/primary_metric_intervals.csv')
    baseline=baseline[baseline.model.eq('xgboost')]
    overall=met[met.family.eq('Overall')]
    c=baseline.merge(overall,on=['window_hours','dataset','metric'],suffixes=('_prior','_now'))
    assert len(c)==48
    for col in ['estimate','lower','upper']:
        assert np.allclose(c[col+'_prior'],c[col+'_now'],atol=1e-9,equal_nan=True),col
    # Compare aggregate composition to established extraction audit, with all mixed codes collapsed.
    for db in ['mimic','eicu']:
        path=ROOT/'private_data'/('mimic_analysis_cvst_v2.csv' if db=='mimic' else 'eicu_analysis_revision.csv')
        raw=pd.read_csv(path)
        for h in [6,24]:
            d=raw[(raw[f'landmark_{h}h_eligible']==1)&raw.outcome_missing_flag.eq(0)].copy()
            d['phenotype']=collapse(d.audit_phenotype)
            expected=d.groupby('phenotype').agg(records=('stay_id','size'),events=(f'death_{h}h_to_30d_inhospital','sum')).sort_index()
            actual=mix[mix.database.eq(db)&mix.window_hours.eq(h)&mix.partition.eq('All_analysis_records')].set_index('phenotype')[['records','events']].sort_index()
            pd.testing.assert_frame_equal(expected,actual,check_dtype=False)
    main_metrics=['auroc','auprc','brier_score','calibration_intercept','calibration_slope']
    report=pop.copy()
    for metric in main_metrics:
        x=met[met.metric.eq(metric)][['window_hours','dataset','family','subgroup','estimate','lower','upper']].copy()
        x[metric]=x.apply(lambda r:f'{r.estimate:.3f} ({r.lower:.3f} to {r.upper:.3f})' if pd.notna(r.lower) else ('Not estimable' if pd.isna(r.estimate) else f'{r.estimate:.3f} (CI unavailable)'),axis=1)
        report=report.merge(x[['window_hours','dataset','family','subgroup',metric]],on=['window_hours','dataset','family','subgroup'],how='left',validate='one_to_one')
    report['event_rate_pct']=100*report.event_rate
    report.to_csv(OUT/'subgroup_reporting_table.csv',index=False)
    report[report.dataset.eq('eICU_external_validation') & report.family.eq('Phenotype')].to_csv(OUT/'external_phenotype_performance.csv',index=False)
    report[report.dataset.eq('eICU_external_validation') & report.family.eq('Restricted evaluation')].to_csv(OUT/'external_restricted_evaluation.csv',index=False)
    threshold=pop[pop.dataset.eq('eICU_external_validation') & pop.status.eq('Reported')].copy()
    for metric in ['sensitivity','specificity','ppv','npv','flagged_fraction']:
        x=met[met.metric.eq(metric)][['window_hours','dataset','family','subgroup','estimate','lower','upper']].copy()
        x[metric]=x.apply(lambda r:f'{r.estimate:.3f} ({r.lower:.3f} to {r.upper:.3f})' if pd.notna(r.lower) else 'CI unavailable',axis=1)
        threshold=threshold.merge(x[['window_hours','dataset','family','subgroup',metric]],on=['window_hours','dataset','family','subgroup'],validate='one_to_one')
    threshold.to_csv(OUT/'external_subgroup_threshold_010.csv',index=False)
    flags=pop[(pop.status=='Reported')&pop.precision_warning.ne('No predefined small-sample flag (not proof of adequate precision)')]
    flags.to_csv(OUT/'precision_warnings.csv',index=False)
    (OUT/'aggregate_checks.json').write_text(json.dumps(dict(protected_shap_entries=len(protected),population_rows=len(pop),interval_rows=len(met),composition_rows=len(mix),point_estimate_groups_verified=checks,overall_metric_reproductions=len(c),precision_flagged_groups=len(flags),minimum_valid_replicates=int(met.valid_reps.min()),status='passed'),indent=2),encoding='utf-8')
    print('Independent aggregate checks passed. External phenotype/restriction summary:')
    print(report[report.dataset.eq('eICU_external_validation') & report.family.isin(['Phenotype','Restricted evaluation'])][['window_hours','subgroup','records','events','auroc','calibration_intercept','precision_warning']].to_string(index=False))
if __name__=='__main__':main()
