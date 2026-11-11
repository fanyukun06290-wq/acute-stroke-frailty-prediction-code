"""Run a named aggregate report stage with input/output hash receipts."""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys
from build_aggregate_reports import TASKS
ROOT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('task',choices=list(TASKS)+['dca','shap']);p.add_argument('--plan',action='store_true');a=p.parse_args()
 if a.task in TASKS:
  ins,outs=TASKS[a.task];command=[sys.executable,str(ROOT/'build_aggregate_reports.py'),'--root',str(ROOT),'--task',a.task]
 else:
  folder='dca_threshold_results' if a.task=='dca' else 'shap_domain_results'
  ins=([n+'.csv' for n in ['fixed_threshold_points','fixed_threshold_intervals','decision_curve_intervals','paired_model_net_benefit']] if a.task=='dca' else [n+'.csv' for n in ['shap_summary_intervals','encoded_domain_dictionary','reconstruction_checks','domain_permutation_repetitions']])
  outs=([n+'.csv' for n in ['primary_operating_points','primary_confusion_010','primary_dca_010','external_comparator_010','descriptive_grid_ranges']] if a.task=='dca' else [n+'.csv' for n in ['domain_sizes','domain_permutation_summary','conceptual_variable_ranking','rank_agreement','external_domain_contributions','external_top10_variables']]+['aggregate_checks.json'])
  ins=[folder+'/'+n for n in ins];outs=[folder+'/'+n for n in outs];command=[sys.executable,str(ROOT/('build_'+a.task+'_presentation.py'))]
 if a.plan:
  print(json.dumps(dict(task=a.task,inputs=ins,outputs=outs,files_read=False)));return
 receipt=ROOT/'qc'/('reporting_'+a.task+'_receipt.json')
 if receipt.exists() or any((ROOT/n).exists() for n in outs):raise FileExistsError('Existing output or receipt preserved')
 if any(not (ROOT/n).is_file() for n in ins):raise FileNotFoundError('Required aggregate input absent')
 before={n:digest(ROOT/n) for n in ins};run=subprocess.run(command,cwd=ROOT,capture_output=True)
 if run.returncode:raise RuntimeError('Reporting stage failed; no successful receipt written. Inspect local inputs and adapter.')
 if any(digest(ROOT/n)!=h for n,h in before.items()):raise RuntimeError('Input changed during reporting')
 if any(not (ROOT/n).is_file() for n in outs):raise RuntimeError('Expected output absent')
 receipt.parent.mkdir(exist_ok=True)
 code={Path(command[1]).name:digest(Path(command[1]))}
 if a.task in TASKS:code['runner_imported_formatter']=digest(ROOT/'build_aggregate_reports.py')
 result=dict(task=a.task,status='aggregate_stage_completed_not_full_reproduction',inputs=before,outputs={n:digest(ROOT/n) for n in outs},code=code)
 with receipt.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
 print(json.dumps({'task':a.task,'outputs':len(outs),'receipt':str(receipt.relative_to(ROOT))}))
if __name__=='__main__':main()
