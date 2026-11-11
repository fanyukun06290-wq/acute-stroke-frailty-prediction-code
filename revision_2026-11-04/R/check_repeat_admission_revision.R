if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
root<-file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_repeat_admission_revision.R"),encoding="UTF-8")
m<-read.csv(file.path(out,"repeat_admission_metric_intervals.csv"));ci_checks<-rep_checks<-0L
for(h in c(24L,6L))for(ds in c("MIMIC_internal_test","eICU_external_validation")){
 pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
 pred<-pred[pred$model=="xgboost" & pred$dataset==ds,];pred<-pred[order(pred$stay_id),]
 cache<-readRDS(file.path(root,"private_cache",sprintf("patient_repeat_%sh_%s.rds",h,ds)))
 stopifnot(identical(cache$stay_id,pred$stay_id))
 fixed<-select_patient_weights(pred,h,ds);y<-as.integer(pred$truth=="yes");p<-pred$probability
 cl<-if(ds=="MIMIC_internal_test")pred$group_id else hospital_components(pred$group_id,pred$hospitalid)$row_component
 idx<-match(cl,unique(cl));ng<-max(idx);set.seed(cache$seed)
 draws<-lapply(1:3,function(b)sample.int(ng,ng,replace=TRUE))
 for(a in names(cache$bootstrap)){
  q<-m[m$window_hours==h & m$dataset==ds & m$analysis==a,]
  for(metric in colnames(cache$bootstrap[[a]])){
   ci<-quantile_ci(cache$bootstrap[[a]][,metric]);r<-q[q$metric==metric,]
   stopifnot(abs(r$estimate-cache$points[[a]][metric])<1e-10,max(abs(c(r$lower,r$upper)-ci[c('lower','upper')]))<1e-10,
       r$valid_reps==ci['valid_reps'])
   ci_checks<-ci_checks+1L
  }
  # Explicit cluster replication preserves fixed within-patient weights and fixed selection.
  for(b in 1:3){
   w<-tabulate(draws[[b]],nbins=ng)[idx];rr<-rep(seq_along(y),w)
   yy<-y[rr];pp<-p[rr];ww<-fixed[[a]][rr]
   explicit<-weighted_metrics(prepare_probability(yy,pp),ww)
   stopifnot(max(abs(explicit-cache$bootstrap[[a]][b,]))<2e-6)
   rep_checks<-rep_checks+1L
  }
  # Calibration-in-the-large checked against a separate offset GLM implementation.
  keep<-fixed[[a]]>0;lp<-qlogis(pmin(pmax(p[keep],1e-6),1-1e-6))
  fit<-suppressWarnings(glm(y[keep]~1+offset(lp),weights=fixed[[a]][keep],family=binomial()))
  stopifnot(isTRUE(fit$converged),abs(coef(fit)[1]-cache$points[[a]]['calibration_intercept'])<2e-5)
 }
}
write_json(list(status='passed',intervals_recomputed=ci_checks,explicit_cluster_replication_checks=rep_checks,
  independent_offset_glm_checks=12),file.path(out,'independent_r_checks.json'))
log_progress('Repeat-admission interval and explicit weighted cluster replication checks passed')
