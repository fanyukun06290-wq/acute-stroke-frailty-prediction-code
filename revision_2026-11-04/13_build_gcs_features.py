"""Construct time-valid documented GCS candidates; no outcomes used in selection."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
COMPONENTS = ["eye", "verbal", "motor"]
LIMITS = {"eye":4, "verbal":5, "motor":6, "direct_total":15}

def build_exams(events, cohort, hour, source):
    eligible = cohort[(cohort[f"landmark_{hour}h_eligible"] == 1) & (cohort.outcome_missing_flag == 0)]
    e = events[events.stay_id.isin(eligible.stay_id)].copy()
    e["explicit_invalid"] = e.explicit_invalid.astype(str).str.lower().isin(["true", "t", "1"])
    end = e.stay_id.map(eligible.set_index("stay_id").effective_icu_end_min)
    temporal = (e.observation_offset >= 0) & (e.observation_offset < hour*60) & (e.observation_offset < end)
    known = temporal & e.available_offset.notna() & (e.available_offset < hour*60) & (e.available_offset >= e.observation_offset)
    counts = {"raw_records_for_eligible_stays":len(e), "records_in_observation_window":int(temporal.sum()),
              "window_records_missing_or_late_availability":int((temporal & (e.available_offset.isna() | (e.available_offset >= hour*60))).sum()),
              "window_records_entry_before_observation":int((temporal & (e.available_offset < e.observation_offset)).sum())}
    e = e[known].copy()
    keys = ["stay_id", "observation_offset"]
    assessed = e[e.kind != "airway"].copy()
    if assessed.empty:
        raise ValueError("No exam records in candidate cohort")
    base = assessed.groupby(keys).agg(available_offset=("available_offset", "max"),
        explicit_quality_block=("explicit_invalid", "max")).reset_index()
    # A non-testable verbal component does not invalidate an observed motor component.
    quality = assessed[assessed.component == "quality"].groupby(keys).explicit_invalid.max()
    base["quality_block"] = pd.MultiIndex.from_frame(base[keys]).map(quality).fillna(False).astype(bool)
    c = assessed[assessed.component.isin(LIMITS)].copy()
    low = np.where(c.component == "direct_total", 3, 1)
    numeric_ok = c.score.notna() & c.score.ge(low) & c.score.le(c.component.map(LIMITS)) & c.score.eq(np.floor(c.score))
    c["invalid"] = c.explicit_invalid | ~numeric_ok
    grouped = c.groupby(keys+["component"]).agg(score=("score", "first"), distinct=("score", "nunique"), invalid=("invalid", "max"))
    for component in [*COMPONENTS, "direct_total"]:
        if component in grouped.index.get_level_values("component"):
            t = grouped.xs(component, level="component")
            index = pd.MultiIndex.from_frame(base[keys])
            base[component] = index.map(t.score)
            base[component+"_valid"] = (index.map(t.distinct).to_numpy() == 1) & ~index.map(t.invalid).fillna(True).to_numpy(dtype=bool)
            base[component+"_recorded"] = index.isin(t.index)
        else:
            base[component] = np.nan
            base[component+"_valid"] = False
            base[component+"_recorded"] = False
    base["airway_state"] = "unknown"
    base["airway_offset"] = np.nan
    if source == "eicu":
        airway = e[e.kind == "airway"].copy()
        if not airway.empty:
            a = airway.groupby(keys).agg(airway_type=("airway_type", "first"), distinct=("airway_type", "nunique")).reset_index()
            a["airway_state"] = np.where((a.distinct == 1) & a.airway_type.eq("No Artificial Airway"), "no_artificial_airway",
                np.where((a.distinct == 1) & a.airway_type.isin(["Oral ETT", "Nasal ETT", "Tracheostomy", "Double-Lumen Tube", "Cricothyrotomy"]), "artificial_airway", "unknown"))
            a = a.rename(columns={"observation_offset":"airway_offset"})
            base = pd.merge_asof(base.drop(columns=["airway_state", "airway_offset"]).sort_values("observation_offset"),
                a[["stay_id", "airway_offset", "airway_state"]].sort_values("airway_offset"),
                left_on="observation_offset", right_on="airway_offset", by="stay_id", direction="backward", tolerance=120)
            base["airway_state"] = base.airway_state.fillna("unknown")
    full = base[[x+"_valid" for x in COMPONENTS]].all(axis=1)
    base["component_sum"] = base[COMPONENTS].sum(axis=1, min_count=3)
    direct_conflict = base.direct_total_recorded & (~base.direct_total_valid | base.direct_total.ne(base.component_sum))
    base["direct_total_conflict"] = direct_conflict
    airway_block = pd.Series(False, index=base.index)
    if source == "eicu":
        airway_block = base.airway_state.eq("artificial_airway") | (base.verbal.eq(1) & ~base.airway_state.eq("no_artificial_airway"))
    base["airway_or_ambiguous_verbal_block"] = airway_block
    base["gcs_valid"] = full & ~base.quality_block & ~direct_conflict & ~airway_block
    base["motor_usable"] = base.motor_valid & ~base.quality_block
    base["gcs_total"] = base.component_sum.where(base.gcs_valid)
    base["gcs_motor"] = base.motor.where(base.motor_usable)
    base["window_hours"] = hour
    base = base.sort_values(keys).reset_index(drop=True)
    counts.update(exam_time_groups=len(base), complete_numeric_component_groups=int(full.sum()),
        explicit_quality_block_groups=int(base.quality_block.sum()),
        direct_total_conflict_groups=int((full & direct_conflict).sum()),
        airway_or_ambiguous_verbal_block_groups=int((full & airway_block).sum()),
        usable_complete_gcs_groups=int(base.gcs_valid.sum()), usable_motor_groups=int(base.motor_usable.sum()))
    return base, eligible, counts

def summarize_features(exams, eligible, hour):
    output = eligible[["stay_id"]].copy().reset_index(drop=True)
    output["window_hours"] = hour
    for variable in ["gcs_total", "gcs_motor"]:
        valid = exams[exams[variable].notna()]
        worst = valid.sort_values([variable, "observation_offset"], ascending=[True, False]).drop_duplicates("stay_id").set_index("stay_id")
        first = valid.sort_values("observation_offset").drop_duplicates("stay_id").set_index("stay_id")
        last = valid.sort_values("observation_offset").drop_duplicates("stay_id", keep="last").set_index("stay_id")
        output[variable+"_min"] = output.stay_id.map(worst[variable])
        output[variable+"_selected_observation_offset"] = output.stay_id.map(worst.observation_offset)
        output[variable+"_selected_available_offset"] = output.stay_id.map(worst.available_offset)
        output[variable+"_first"] = output.stay_id.map(first[variable])
        output[variable+"_last"] = output.stay_id.map(last[variable])
        output[variable+"_valid_exam_count"] = output.stay_id.map(valid.groupby("stay_id").size()).fillna(0).astype(int)
        output[variable+"_missing"] = output[variable+"_min"].isna().astype(int)
    return output

def synthetic_checks():
    cohort = pd.DataFrame({"stay_id":[1],"landmark_6h_eligible":[1],"outcome_missing_flag":[0],"effective_icu_end_min":[500]})
    def row(component, score, offset=100, available=None, invalid=False, kind="component", airway=""):
        return dict(stay_id=1,kind=kind,observation_offset=offset,available_offset=offset if available is None else available,
                    component=component,score=score,explicit_invalid=invalid,airway_type=airway)
    checks = 0
    def check(rows, source, total, motor):
        nonlocal checks
        e, c, _ = build_exams(pd.DataFrame(rows), cohort, 6, source)
        f = summarize_features(e,c,6).iloc[0]
        assert (pd.isna(f.gcs_total_min) if total is None else f.gcs_total_min == total)
        assert (pd.isna(f.gcs_motor_min) if motor is None else f.gcs_motor_min == motor)
        checks += 1
    normal = [row("eye",4),row("verbal",5),row("motor",6)]
    check(normal,"mimic",15,6)
    check([row("eye",4),row("verbal",5),row("motor",6,offset=101)],"mimic",None,6)
    check([row("eye",4),row("verbal",1,invalid=True),row("motor",6)],"mimic",None,6)
    check(normal+[row("motor",2)],"mimic",None,None)
    check([row("eye",4),row("verbal",5,available=360),row("motor",6)],"mimic",None,6)
    check(normal+[row("quality",np.nan,invalid=True,kind="quality")],"eicu",None,None)
    check(normal+[row("airway",np.nan,offset=90,kind="airway",airway="Oral ETT")],"eicu",None,6)
    low = [row("eye",1),row("verbal",1),row("motor",1)]
    check(low,"eicu",None,1)
    check(low+[row("airway",np.nan,offset=90,kind="airway",airway="No Artificial Airway")],"eicu",3,1)
    check(low+[row("airway",np.nan,offset=110,kind="airway",airway="No Artificial Airway")],"eicu",None,1)
    check(normal+[row("direct_total",14)],"eicu",None,6)
    check(normal+[row("direct_total",15)],"eicu",15,6)
    check(low,"mimic",3,1)
    check(normal+[row("eye",0,offset=360)],"mimic",15,6)
    check(normal+[row("eye",0,offset=200,available=190)],"mimic",15,6)
    return checks

def main():
    out = ROOT/"gcs_results"
    out.mkdir(parents=True, exist_ok=True)
    if (out/"construction_complete.json").exists():
        raise RuntimeError("Preserve completed GCS construction")
    synthetic = synthetic_checks()
    plan = {"version":"strict_documented_gcs_1.0","rules_sha256":hashlib.sha256((ROOT/"10_第六步GCS构建规则.md").read_bytes()).hexdigest(),
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"outcome_used_for_score_selection":False,
            "full_gcs_role":"restricted supplementary candidate, residual cross-database assessment differences",
            "motor_role":"documented motor component, not full GCS or proven unsedated examination"}
    coverage, reasons, distributions, audits = [], [], [], []
    for source in ["mimic","eicu"]:
        path = ROOT/"private_data"/f"{source}_gcs_events.csv"
        expected = json.loads((ROOT/"qc"/f"{source}_gcs_extract.json").read_text())["output_sha256"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        events = pd.read_csv(path)
        cohort = pd.read_csv(ROOT/"private_data"/f"{source}_cohort.csv")
        plan[source+"_events_sha256"] = expected
        all_exams, all_features = [], []
        for hour in [6,24]:
            exams, eligible, counts = build_exams(events,cohort,hour,source)
            f = summarize_features(exams,eligible,hour)
            assert len(f) == len(eligible) and f.stay_id.is_unique
            for variable, low, high in [("gcs_total",3,15),("gcs_motor",1,6)]:
                s = f[variable+"_min"].dropna()
                assert s.between(low,high).all() and s.eq(np.floor(s)).all()
                for field in ["selected_observation_offset","selected_available_offset"]:
                    values = f[variable+"_"+field].dropna()
                    assert values.ge(0).all() and values.lt(hour*60).all()
                assert (f[variable+"_missing"] == f[variable+"_min"].isna().astype(int)).all()
                for score,count in s.value_counts().sort_index().items():
                    distributions.append(dict(database=source,window_hours=hour,variable=variable,score=int(score),records=int(count)))
            assert not (f.gcs_total_min.notna() & f.gcs_motor_min.isna()).any()
            if source == "eicu":
                used = exams[exams.airway_offset.notna()]
                assert (used.airway_offset<=used.observation_offset).all()
                assert ((used.observation_offset-used.airway_offset)<=120).all()
            coverage.append(dict(database=source,window_hours=hour,records=len(f),
                full_gcs_available=int(f.gcs_total_min.notna().sum()),motor_available=int(f.gcs_motor_min.notna().sum()),
                full_gcs_missing_fraction=float(f.gcs_total_missing.mean()),motor_missing_fraction=float(f.gcs_motor_missing.mean())))
            reasons.append(dict(database=source,window_hours=hour,**counts))
            audits.append(dict(database=source,window_hours=hour,all_scores_in_range=True,selected_times_before_landmark=True,
                no_cross_time_component_fill=True,all_eligible_stays_retained=True))
            all_exams.append(exams); all_features.append(f)
        pd.concat(all_exams,ignore_index=True).to_csv(ROOT/"private_data"/f"{source}_gcs_exam_candidates.csv",index=False)
        pd.concat(all_features,ignore_index=True).to_csv(ROOT/"private_data"/f"{source}_gcs_features.csv",index=False)
    pd.DataFrame(coverage).to_csv(out/"availability.csv",index=False)
    pd.DataFrame(reasons).to_csv(out/"quality_audit_counts.csv",index=False)
    pd.DataFrame(distributions).to_csv(out/"score_distributions.csv",index=False)
    (out/"construction_plan.json").write_text(json.dumps(plan,indent=2),encoding="utf-8")
    complete = dict(status="documented_score_construction_complete_not_model_validation",synthetic_cases_passed=synthetic,
        cohort_checks=audits,model_fitting_performed=False,imputed_scores=False,
        limitation="Recorded-time grouping and unequal sedation/airway documentation do not establish fully equivalent unsedated GCS across databases")
    (out/"construction_complete.json").write_text(json.dumps(complete,indent=2),encoding="utf-8")
    print(pd.DataFrame(coverage).to_string(index=False))
    print(f"SYNTHETIC_TEMPORAL_AND_ASSESSABILITY_CASES_PASSED={synthetic}")

if __name__ == "__main__":
    main()
