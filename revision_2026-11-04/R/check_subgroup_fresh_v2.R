if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
root<-file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_subgroup_fresh_v2.R"),encoding="UTF-8")
m<-read.csv(file.path(out,"subgroup_metric_intervals.csv"));ic<-rc<-0L
for(h in c(24L,6L))for(ds in c("MIMIC_internal_test","eICU_external_validation")){
 cache<-readRDS(file.path(root,"private_cache",sprintf("subgroup_fresh_v2_%sh_%s.rds",h,ds)))
 raw<-read_data(file.path(root,"private_data",if(ds=="MIMIC_internal_test")"mimic_analysis_cvst_v2.csv" else "eicu_analysis_revision.csv"),contract,h)
 pr<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh_cvst_v2.csv",h)))
 pr<-pr[pr$model=="xgboost" & pr$dataset==ds,];pr<-pr[order(pr$stay_id),]
 cl<-if(ds=="MIMIC_internal_test")pr$group_id else hospital_components(pr$group_id,pr$hospitalid)$row_component
 ci<-match(cl,unique(cl));ng<-max(ci)
 set.seed(cache$seed);draws<-lapply(1:3,function(b)sample.int(ng,ng,replace=TRUE))
 for(key in names(cache$groups)){
  g<-cache$groups[[key]];idx<-match(g$stay_id,pr$stay_id);stopifnot(!anyNA(idx))
  parts<-strsplit(key,"::",fixed=TRUE)[[1]]
  q<-m[m$window_hours==h & m$dataset==ds & m$family==parts[1] & m$subgroup==parts[2],]
  for(metric in names(g$point)){
   row<-q[q$metric==metric,];valid<-sum(is.finite(g$bootstrap[,metric]));v<-g$bootstrap[is.finite(g$bootstrap[,metric]),metric]
   stopifnot(row$valid_reps==valid,abs(row$estimate-g$point[metric])<1e-10)
   if(valid>=950)stopifnot(max(abs(c(row$lower,row$upper)-quantile(v,c(.025,.975),names=FALSE)))<1e-10)
   ic<-ic+1L
  }
  # Check optimized weight logic against explicit duplicated record rows, including calibration fits.
  for(b in 1:3){
   w<-tabulate(draws[[b]],nbins=ng)[ci[idx]];rep_idx<-rep(seq_along(g$truth),w)
   if(!length(rep_idx) || length(unique(g$truth[rep_idx]))<2)next
   direct<-point_scores(g$truth[rep_idx],g$p[rep_idx])
   metrics<-c("auroc","auprc","brier_score","log_loss","calibration_intercept","calibration_slope")
   stopifnot(max(abs(direct[metrics]-g$bootstrap[b,metrics]))<2e-6)
   op<-weighted_operating(g$truth[rep_idx],g$p[rep_idx],rep(1,length(rep_idx)))
   stopifnot(max(abs(op-g$bootstrap[b,names(op)]),na.rm=TRUE)<1e-10)
   rc<-rc+1L
  }
 }
}
write_json(list(status="passed",intervals_recomputed=ic,explicit_record_replication_checks=rc),file.path(out,"independent_r_checks.json"))
log_progress("Subgroup interval and explicit row-replication checks passed: ",ic," intervals; ",rc," bootstrap subgroup checks")
