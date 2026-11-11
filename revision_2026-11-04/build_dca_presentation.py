"""Aggregate presentation tables; never optimize thresholds or modify predictions."""
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"dca_threshold_results"

def contiguous_ranges(frame,criterion):
    f=frame.sort_values("threshold")
    keep=f.estimate.gt(1e-12) if criterion=="point_estimate_above_both_defaults" else f.lower.gt(1e-12)
    indices=(f.loc[keep,"threshold"]*100).round().astype(int).tolist()
    blocks=[]
    for value in indices:
        if not blocks or value>blocks[-1][-1]+1:
            blocks.append([value])
        else:
            blocks[-1].append(value)
    return "; ".join(f"{b[0]/100:.2f}" if len(b)==1 else f"{b[0]/100:.2f}–{b[-1]/100:.2f}" for b in blocks) or "none"

def main():
    names=['primary_operating_points.csv','primary_confusion_010.csv','primary_dca_010.csv','external_comparator_010.csv','descriptive_grid_ranges.csv']
    assert not any((OUT/n).exists() for n in names),'Existing report outputs preserved'
    points=pd.read_csv(OUT/"fixed_threshold_points.csv")
    intervals=pd.read_csv(OUT/"fixed_threshold_intervals.csv")
    curves=pd.read_csv(OUT/"decision_curve_intervals.csv")
    pairs=pd.read_csv(OUT/"paired_model_net_benefit.csv")
    primary=points[points.model=="xgboost"].copy()
    def display(r,metric):
        d=intervals[(intervals.window_hours==r.window_hours)&(intervals.dataset==r.dataset)&
            (intervals.model=="xgboost")&(intervals.threshold==r.threshold)&(intervals.metric==metric)].iloc[0]
        return f"{d.estimate:.3f} ({d.lower:.3f}–{d.upper:.3f})"
    for metric in ["sensitivity","specificity","ppv","npv","flagged_fraction"]:
        primary[metric+"_95ci"]=[display(r,metric) for r in primary.itertuples()]
    cols=["window_hours","dataset","threshold","records","events","sensitivity_95ci","specificity_95ci","ppv_95ci","npv_95ci","flagged_fraction_95ci"]
    primary[cols].to_csv(OUT/"primary_operating_points.csv",index=False)
    primary[primary.threshold==.10][["window_hours","dataset","records","events","tp","fp","tn","fn",
        "false_positives_per1000","false_negatives_per1000","residual_event_rate_negative"]].to_csv(OUT/"primary_confusion_010.csv",index=False)
    curves[(curves.model=="xgboost")&(curves.threshold==.10)&curves.metric.isin([
        "net_benefit","difference_vs_all","difference_vs_best_default"])].to_csv(OUT/"primary_dca_010.csv",index=False)
    comparison=pairs[(pairs.dataset=="eICU_external_validation")&(pairs.threshold==.10)].copy()
    comparison["comparator_role"]=comparison.comparator.map(lambda n:"primary_algorithm_comparator" if n in ["random_forest","elastic_net"] else
        "posthoc_missingness_sensitivity" if n.startswith("omit_") else "supplementary_clinical_neurological_reference")
    comparison.to_csv(OUT/"external_comparator_010.csv",index=False)
    rows=[]
    selected=curves[(curves.model=="xgboost")&(curves.metric=="difference_vs_best_default")]
    for (hour,dataset),d in selected.groupby(["window_hours","dataset"]):
        for criterion in ["point_estimate_above_both_defaults","pointwise_lower_delta_best_above_zero"]:
            rows.append(dict(window_hours=hour,dataset=dataset,criterion=criterion,
                observed_grid_ranges=contiguous_ranges(d,criterion),grid_step=.01,
                interpretation="Exploratory observed grid only; not simultaneous inference or a recommended clinical threshold range"))
    pd.DataFrame(rows).to_csv(OUT/"descriptive_grid_ranges.csv",index=False)
    print(pd.DataFrame(rows)[["window_hours","dataset","criterion","observed_grid_ranges"]].to_string(index=False))

if __name__=="__main__":
    main()
