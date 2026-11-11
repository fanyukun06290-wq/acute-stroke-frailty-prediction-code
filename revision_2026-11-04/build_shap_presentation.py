"""Verify aggregate identities, preserve earlier artifacts, and prepare reporting tables."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'shap_domain_results'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    protected=[] # No historical review-manifest verification claimed.
    names=['domain_sizes.csv','domain_permutation_summary.csv','conceptual_variable_ranking.csv','rank_agreement.csv','external_domain_contributions.csv','external_top10_variables.csv','aggregate_checks.json']
    assert not any((OUT/n).exists() for n in names),'Existing report outputs preserved'
    s=pd.read_csv(OUT/'shap_summary_intervals.csv'); d=pd.read_csv(OUT/'encoded_domain_dictionary.csv')
    q=pd.read_csv(OUT/'reconstruction_checks.csv'); p=pd.read_csv(OUT/'domain_permutation_repetitions.csv')
    assert len(q)==4 and len(p)==4800 and (q.margin_additivity_error<2e-4).all() and (q.frozen_probability_error<1e-6).all()
    assert dict(zip(zip(q.window_hours,q.dataset),q.records))=={(24,'MIMIC_internal_test'):1303,(24,'eICU_external_validation'):3730,(6,'MIMIC_internal_test'):1490,(6,'eICU_external_validation'):4479}
    assert (s.valid_reps==1000).all() and np.isfinite(s[['estimate','lower','upper']]).all().all()
    assert (s.lower<=s.upper).all() and (s.estimate>=0).all()
    assert not d.feature.str.contains('frailty|comorbidity|stroke|phenotype',case=False).any()
    domain_sizes=[]
    for h,g in d.groupby('window_hours'):
        assert g.feature.is_unique and g.variable.nunique()==20
        for domain,k in g.groupby('domain'):
            domain_sizes.append(dict(window_hours=h,domain=domain,encoded_features=len(k),conceptual_variables=k.variable.nunique(),missingness_indicators=k.kind.eq('Missingness indicator').sum()))
    sizes=pd.DataFrame(domain_sizes);sizes.to_csv(OUT/'domain_sizes.csv',index=False)
    checks=0
    for (h,ds),g in s.groupby(['window_hours','dataset']):
        e=g[g.level.eq('Encoded feature')].set_index('name')
        dd=d[d.window_hours.eq(h)]
        for domain,k in dd.groupby('domain'):
            a=g[g.level.eq('Domain') & g.name.eq(domain)].set_index('metric').estimate
            total=e.loc[k.feature,'estimate'].sum()
            assert np.isclose(a['sum_absolute'],total)
            assert np.isclose(a['sum_absolute_per_encoded_feature'],total/len(k))
            assert np.isclose(a['sum_absolute_per_conceptual_variable'],total/k.variable.nunique())
            assert a['absolute_signed_sum']<=total+1e-12
            for v,vg in k.groupby('variable'):
                value=g[g.level.eq('Conceptual variable') & g.name.eq(v)].estimate.item()
                assert value<=e.loc[vg.feature,'estimate'].sum()+1e-12
            checks+=1
        assert np.isclose(g[g.level.eq('Encoding kind')].estimate.sum(),e.estimate.sum())
    agg=p.groupby(['window_hours','dataset','scope','domain','metric'],sort=False).performance_deterioration.agg(
        mean='mean',sd='std',q025=lambda x:x.quantile(.025),q975=lambda x:x.quantile(.975),repetitions='count').reset_index()
    assert len(agg)==96 and (agg.repetitions==50).all()
    agg['interval_type']='Monte Carlo shuffle spread, not sampling CI'
    agg.to_csv(OUT/'domain_permutation_summary.csv',index=False)
    rankings=[];agreement=[]
    for h in [24,6]:
        for level in ['Encoded feature','Conceptual variable']:
            k=s[(s.window_hours==h)&s.level.eq(level)]
            internal=k[k.dataset.eq('MIMIC_internal_test')].set_index('name').estimate
            external=k[k.dataset.eq('eICU_external_validation')].set_index('name').estimate
            agreement.append(dict(window_hours=h,level=level,spearman_rank_correlation=internal.rank().corr(external.rank()),top5_overlap=len(set(internal.nlargest(5).index)&set(external.nlargest(5).index)),features=len(internal)))
        for ds in q.dataset.unique():
            k=s[(s.window_hours==h)&s.dataset.eq(ds)&s.level.eq('Conceptual variable')].sort_values('estimate',ascending=False).copy()
            k['rank']=range(1,len(k)+1);rankings.append(k)
    pd.concat(rankings).to_csv(OUT/'conceptual_variable_ranking.csv',index=False)
    pd.DataFrame(agreement).to_csv(OUT/'rank_agreement.csv',index=False)
    ext=s[s.dataset.eq('eICU_external_validation') & s.level.eq('Domain')]
    ext.to_csv(OUT/'external_domain_contributions.csv',index=False)
    top=pd.concat(rankings);top=top[top.dataset.eq('eICU_external_validation')&top['rank'].le(10)]
    top.to_csv(OUT/'external_top10_variables.csv',index=False)
    result=dict(earlier_manifest_entries_verified=len(protected),reconstruction_sets_verified=len(q),aggregate_domain_identities_verified=checks,
        shap_summary_rows=len(s),permutation_rows=len(p),permutation_summary_rows=len(agg),valid_bootstrap_replicates=1000,
        status='aggregate_arithmetic_passed_not_manuscript_ready')
    (OUT/'aggregate_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))
    print('External top5 conceptual variables:')
    print(top[top['rank'].le(5)][['window_hours','rank','name','estimate']].to_string(index=False))
    print('External domain AUROC permutation summaries:')
    print(agg[agg.dataset.eq('eICU_external_validation')&agg.metric.eq('auroc')][['window_hours','scope','domain','mean','sd']].to_string(index=False))
if __name__=='__main__': main()
