from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'patientflow_results'
OUT.mkdir(exist_ok=True)
def boolean(x):
    values = x.astype(str).str.lower().str.strip()
    assert values.isin(['true', 'false', 't', 'f', '1', '0']).all()
    return values.isin(['true', 't', '1'])

def metadata(db):
    z = pd.read_csv(ROOT / 'private_data' / f'{db}_cohort.csv')
    assert z.stay_id.is_unique and z.patient_group_id.notna().all()
    admission = 'hadm_id' if db == 'mimic' else 'patienthealthsystemstayid'
    assert z[admission].is_unique
    assert z.groupby(admission).patient_group_id.nunique().max() == 1
    if db == 'mimic':
        assert ('MIMIC:' + z.subject_id.astype(str)).eq(z.patient_group_id.astype(str)).all()
    else:
        assert ('eICU:' + z.uniquepid.astype(str)).eq(z.patient_group_id.astype(str)).all()
    return z

def main():
    assert not any(((OUT / n).exists() for n in ['cohort_flow_counts.csv', 'repeat_admission_burden.csv', 'cross_window_overlap.csv', 'patient_partition_checks.csv'])), 'Existing report outputs preserved'
    flow = []
    repeat = []
    overlap = []
    splits = []
    for db in ['mimic', 'eicu']:
        raw = metadata(db)
        sets = {}
        patientsets = {}
        testsets = {}
        trainsets = {}
        for h in [24, 6]:
            invalid = boolean(raw.invalid_time_flag)
            valid = ~invalid
            eligible = raw[f'landmark_{h}h_eligible'].eq(1)
            assert np.array_equal(eligible, valid & raw.effective_icu_end_min.gt(h * 60))
            keep = eligible & raw.outcome_missing_flag.eq(0)
            sets[h] = set(raw.loc[keep, 'stay_id'])
            patientsets[h] = set(raw.loc[keep, 'patient_group_id'])
            masks = {'Source coded adult index stays': np.ones(len(raw), bool), 'Excluded invalid chronology': invalid, 'Valid chronology': valid, 'Excluded before/not reaching landmark': valid & ~eligible, 'Landmark eligible': eligible, 'Excluded unknown/conflicting outcome': eligible & ~raw.outcome_missing_flag.eq(0), 'Analysis cohort': keep}
            assert masks['Source coded adult index stays'].sum() == invalid.sum() + (valid & ~eligible).sum() + (eligible & ~raw.outcome_missing_flag.eq(0)).sum() + keep.sum()
            if db == 'mimic':
                registry = pd.read_csv(ROOT / 'private_splits' / f'patients_{h}h.csv')
                assert registry.group_id.is_unique and set(registry.group_id) == patientsets[h]
                mapped = raw.patient_group_id.map(registry.set_index('group_id').partition)
                train = keep & mapped.eq('train')
                test = keep & mapped.eq('test')
                trainsets[h] = set(raw.loc[train, 'patient_group_id'])
                testsets[h] = set(raw.loc[test, 'patient_group_id'])
                assert not trainsets[h] & testsets[h] and train.sum() + test.sum() == keep.sum()
                assert registry[registry.partition.eq('train')].fold.between(1, 5).all()
                assert registry[registry.partition.eq('test')].fold.isna().all()
                masks.update({'Training partition': train, 'Internal test partition': test})
                splits.append(dict(window_hours=h, training_patients=len(trainsets[h]), test_patients=len(testsets[h]), patient_overlap=0, training_records=int(train.sum()), test_records=int(test.sum()), training_folds=5, one_fold_per_training_patient=True))
            else:
                masks['External validation'] = keep
            for stage, m in masks.items():
                z = raw[m]
                known = z.outcome_missing_flag.eq(0)
                event = int(z[f'death_{h}h_to_30d_inhospital'].sum()) if stage in ['Analysis cohort', 'Training partition', 'Internal test partition', 'External validation'] else np.nan
                flow.append(dict(database=db, window_hours=h, stage=stage, records=len(z), unique_patients=z.patient_group_id.nunique(), events=event))
                if stage in ['Source coded adult index stays', 'Analysis cohort', 'Training partition', 'Internal test partition', 'External validation']:
                    counts = z.groupby('patient_group_id').size()
                    repeat.append(dict(database=db, window_hours=h, population=stage, records=len(z), unique_patients=len(counts), patients_with_multiple_records=int((counts > 1).sum()), excess_records=len(z) - len(counts), maximum_records_per_patient=int(counts.max()), fraction_records_from_repeat_patients=float(counts[counts > 1].sum() / len(z))))
        assert sets[24] <= sets[6] and patientsets[24] <= patientsets[6]
        r = dict(database=db, records_24h=len(sets[24]), records_6h=len(sets[6]), shared_records=len(sets[24] & sets[6]), records_6h_only=len(sets[6] - sets[24]), unique_patients_24h=len(patientsets[24]), unique_patients_6h=len(patientsets[6]), shared_patients=len(patientsets[24] & patientsets[6]))
        if db == 'mimic':
            r.update(internal_test_patients_24h=len(testsets[24]), internal_test_patients_6h=len(testsets[6]), internal_test_shared_patients=len(testsets[24] & testsets[6]), test24_train6_patients=len(testsets[24] & trainsets[6]), test6_train24_patients=len(testsets[6] & trainsets[24]))
        overlap.append(r)
    pd.DataFrame(flow).to_csv(OUT / 'cohort_flow_counts.csv', index=False)
    pd.DataFrame(repeat).to_csv(OUT / 'repeat_admission_burden.csv', index=False)
    pd.DataFrame(overlap).to_csv(OUT / 'cross_window_overlap.csv', index=False)
    pd.DataFrame(splits).to_csv(OUT / 'patient_partition_checks.csv', index=False)
if __name__=='__main__':main()
