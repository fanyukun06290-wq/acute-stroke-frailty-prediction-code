"""Outcome-blind score-rule variants; preserve original raw events and scores."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("frozen_gcs_builder",ROOT/"13_build_gcs_features.py")
builder=importlib.util.module_from_spec(spec); spec.loader.exec_module(builder)
VARIANTS=["strict_baseline","allow_reverse_recording","first_valid","last_valid",
          "motor_requires_complete","airway_60","airway_240"]

def available_events(events,eligible,hour):
    e=events[events.stay_id.isin(eligible.stay_id)].copy()
    ending=e.stay_id.map(eligible.set_index("stay_id").effective_icu_end_min)
    return e[(e.observation_offset>=0)&(e.observation_offset<hour*60)&(e.observation_offset<ending)&
             e.available_offset.notna()&(e.available_offset<hour*60)&(e.available_offset>=e.observation_offset)]

def replace_airway(exams,events,eligible,hour,tolerance):
    e=available_events(events,eligible,hour)
    e=e[e.kind=="airway"]
    base=exams.drop(columns=["airway_state","airway_offset"]).copy()
    if e.empty:
        base["airway_state"]="unknown"; base["airway_offset"]=np.nan
    else:
        a=e.groupby(["stay_id","observation_offset"]).agg(kind=("airway_type","first"),distinct=("airway_type","nunique")).reset_index()
        a["airway_state"]=np.where((a.distinct==1)&a.kind.eq("No Artificial Airway"),"no_artificial_airway",
            np.where((a.distinct==1)&a.kind.isin(["Oral ETT","Nasal ETT","Tracheostomy","Double-Lumen Tube","Cricothyrotomy"]),"artificial_airway","unknown"))
        a=a.rename(columns={"observation_offset":"airway_offset"})
        base=pd.merge_asof(base.sort_values("observation_offset"),a[["stay_id","airway_offset","airway_state"]].sort_values("airway_offset"),
            left_on="observation_offset",right_on="airway_offset",by="stay_id",direction="backward",tolerance=tolerance)
        base["airway_state"]=base.airway_state.fillna("unknown")
    block=base.airway_state.eq("artificial_airway")|(base.verbal.eq(1)&~base.airway_state.eq("no_artificial_airway"))
    base["airway_or_ambiguous_verbal_block"]=block
    full=base[[x+"_valid" for x in builder.COMPONENTS]].all(axis=1)
    base["gcs_valid"]=full&~base.quality_block&~base.direct_total_conflict&~block
    base["gcs_total"]=base.component_sum.where(base.gcs_valid)
    return base.sort_values(["stay_id","observation_offset"]).reset_index(drop=True)

def make_variant(events,cohort,hour,source,variant):
    adapted=events
    if variant=="allow_reverse_recording":
        adapted=events.copy()
        # This adapter passes effective information time to the frozen builder,
        # not a rewritten storage timestamp. Raw files remain byte-identical.
        adapted["actual_recording_offset"]=adapted.available_offset
        adapted["available_offset"]=np.maximum(adapted.available_offset,adapted.observation_offset)
    exams,eligible,_=builder.build_exams(adapted,cohort,hour,source)
    if variant.startswith("airway_") and source=="eicu":
        exams=replace_airway(exams,events,eligible,hour,int(variant.split("_")[1]))
    if variant=="motor_requires_complete":
        exams=exams.copy(); exams["motor_usable"]=exams.motor_usable&exams.gcs_valid
        exams["gcs_motor"]=exams.motor.where(exams.motor_usable)
    f=builder.summarize_features(exams,eligible,hour)
    if variant in ["first_valid","last_valid"]:
        selector="first" if variant=="first_valid" else "last"
        for variable in ["gcs_total","gcs_motor"]:
            valid=exams[exams[variable].notna()].sort_values("observation_offset")
            chosen=valid.drop_duplicates("stay_id",keep="first" if selector=="first" else "last").set_index("stay_id")
            f[variable+"_min"]=f.stay_id.map(chosen[variable])
            f[variable+"_selected_observation_offset"]=f.stay_id.map(chosen.observation_offset)
            f[variable+"_selected_available_offset"]=f.stay_id.map(chosen.available_offset)
    f["variant"]=variant
    f["selected_value_semantics"]="window_minimum" if variant not in ["first_valid","last_valid"] else variant
    f["availability_offset_semantics"]="max_observation_and_recording" if variant=="allow_reverse_recording" else "recorded_source_time"
    return f,exams

def synthetic_checks():
    cohort=pd.DataFrame({"stay_id":[1],"landmark_6h_eligible":[1],"outcome_missing_flag":[0],"effective_icu_end_min":[500]})
    def row(component,score,at=100,entry=None,invalid=False,kind="component",airway=""):
        return dict(stay_id=1,kind=kind,observation_offset=at,available_offset=at if entry is None else entry,
                    component=component,score=score,explicit_invalid=invalid,airway_type=airway)
    rows=[row("eye",4),row("verbal",5),row("motor",6,entry=90)]
    e=pd.DataFrame(rows)
    strict,_=make_variant(e,cohort,6,"mimic","strict_baseline")
    relaxed,_=make_variant(e,cohort,6,"mimic","allow_reverse_recording")
    assert pd.isna(strict.gcs_total_min.iloc[0]) and relaxed.gcs_total_min.iloc[0]==15
    late=pd.DataFrame([row("eye",4),row("verbal",5),row("motor",6,entry=360)])
    f,_=make_variant(late,cohort,6,"mimic","allow_reverse_recording")
    assert pd.isna(f.gcs_total_min.iloc[0])
    e=pd.DataFrame([row("eye",4),row("verbal",5),row("motor",2),row("eye",4,at=200),row("verbal",5,at=200),row("motor",6,at=200)])
    first,_=make_variant(e,cohort,6,"mimic","first_valid")
    last,_=make_variant(e,cohort,6,"mimic","last_valid")
    assert first.gcs_motor_min.iloc[0]==2 and last.gcs_motor_min.iloc[0]==6
    e=pd.DataFrame([row("eye",4),row("verbal",1,invalid=True),row("motor",2)])
    f,_=make_variant(e,cohort,6,"mimic","motor_requires_complete")
    assert pd.isna(f.gcs_motor_min.iloc[0])
    e=pd.DataFrame([row("eye",1,at=200),row("verbal",1,at=200),row("motor",1,at=200),
        row("airway",np.nan,at=100,kind="airway",airway="No Artificial Airway")])
    a,_=make_variant(e,cohort,6,"eicu","airway_60")
    b,_=make_variant(e,cohort,6,"eicu","airway_240")
    assert pd.isna(a.gcs_total_min.iloc[0]) and b.gcs_total_min.iloc[0]==3
    return 5

def main():
    out=ROOT/"gcs_rule_sensitivity_results"; out.mkdir(exist_ok=True)
    if (out/"score_construction_complete.json").exists():
        raise RuntimeError("Preserve completed score sensitivity outputs")
    checks=synthetic_checks(); comparison=[]; distribution=[]; manifests=[]
    for source in ["mimic","eicu"]:
        rawpath=ROOT/"private_data"/f"{source}_gcs_events.csv"
        originalpath=ROOT/"private_data"/f"{source}_gcs_features.csv"
        events=pd.read_csv(rawpath); cohort=pd.read_csv(ROOT/"private_data"/f"{source}_cohort.csv")
        original=pd.read_csv(originalpath); outputs=[]
        for h in [6,24]:
            expected=original[original.window_hours==h].set_index("stay_id")
            for variant in VARIANTS:
                f,exams=make_variant(events,cohort,h,source,variant)
                actual=f.set_index("stay_id").reindex(expected.index)
                assert set(f.stay_id)==set(expected.index) and f.stay_id.is_unique
                for variable,low,high in [("gcs_total",3,15),("gcs_motor",1,6)]:
                    a=actual[variable+"_min"]; b=expected[variable+"_min"]
                    both=a.notna()&b.notna(); changed=both&a.ne(b)
                    assert a.dropna().between(low,high).all() and a.dropna().eq(np.floor(a.dropna())).all()
                    assert (actual[variable+"_missing"]==a.isna().astype(int)).all()
                    for column in ["_selected_observation_offset","_selected_available_offset"]:
                        times=actual[variable+column].dropna(); assert times.ge(0).all() and times.lt(h*60).all()
                    delta=a[both]-b[both]
                    comparison.append(dict(database=source,window_hours=h,variant=variant,score_type=variable,records=len(a),
                        baseline_available=int(b.notna().sum()),variant_available=int(a.notna().sum()),newly_available=int((a.notna()&b.isna()).sum()),
                        newly_missing=int((a.isna()&b.notna()).sum()),changed_among_both_observed=int(changed.sum()),
                        unchanged_including_missing=int(((a.eq(b))|(a.isna()&b.isna())).sum()),
                        mean_difference_among_both=float(delta.mean()),mean_absolute_difference_among_both=float(delta.abs().mean())))
                    for value,count in delta.value_counts().sort_index().items():
                        distribution.append(dict(database=source,window_hours=h,variant=variant,score_type=variable,
                            score_difference=float(value),records=int(count)))
                    identity=variant=="strict_baseline" or (source=="eicu" and variant=="allow_reverse_recording") or \
                        (source=="mimic" and variant.startswith("airway_")) or (variable=="gcs_motor" and variant.startswith("airway_")) or \
                        (variable=="gcs_total" and variant=="motor_requires_complete")
                    if identity:
                        assert ((a.eq(b))|(a.isna()&b.isna())).all()
                outputs.append(f)
        target=ROOT/"private_data"/f"{source}_gcs_rule_variants.csv"
        pd.concat(outputs,ignore_index=True).to_csv(target,index=False)
        manifests.append(dict(database=source,raw_events_sha256=hashlib.sha256(rawpath.read_bytes()).hexdigest(),
            original_scores_sha256=hashlib.sha256(originalpath.read_bytes()).hexdigest(),
            variants_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    pd.DataFrame(comparison).to_csv(out/"score_changes.csv",index=False)
    pd.DataFrame(distribution).to_csv(out/"score_difference_distributions.csv",index=False)
    result=dict(status="score_variants_constructed_not_retrained_model_validation",posthoc=True,
        rules_sha256=hashlib.sha256((ROOT/"16_第九步GCS规则敏感性方案.md").read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sources=manifests,
        synthetic_case_groups_passed=checks,score_comparison_rows=len(comparison),outcomes_used_for_score_selection=False,
        frozen_builder_sha256=hashlib.sha256((ROOT/"13_build_gcs_features.py").read_bytes()).hexdigest())
    (out/"score_construction_complete.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(pd.DataFrame(comparison).query("variant=='allow_reverse_recording' or variant=='motor_requires_complete'")
          [["database","window_hours","variant","score_type","variant_available","changed_among_both_observed","newly_missing"]].to_string(index=False))
    print("GCS_RULE_VARIANTS_CONSTRUCTED")

if __name__=="__main__":
    main()
