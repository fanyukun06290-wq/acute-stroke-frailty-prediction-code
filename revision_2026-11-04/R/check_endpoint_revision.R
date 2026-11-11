if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
root<-file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_endpoint_revision.R"),encoding="UTF-8")
cm<-read.csv(file.path(out,"competing_event_curves.csv"));sm<-read.csv(file.path(out,"endpoint_sensitivity_intervals.csv"))
pm<-read.csv(file.path(out,"endpoint_proxy_paired_differences.csv"));rm<-read.csv(file.path(out,"fixed_risk_bin_cif.csv"))
cc<-sc<-pc<-bc<-native_checks<-0L
for(h in c(24L,6L))for(ds in c("MIMIC_internal_test","eICU_external_validation")){
 db<-if(ds=="MIMIC_internal_test")"mimic" else "eicu"
 z<-read_endpoint(db,h,ds);pr<-z$pred
 cache<-readRDS(file.path(root,"private_cache",sprintf("endpoint_%sh_%s.rds",h,ds)))
 stopifnot(identical(cache$stay_id,pr$stay_id),identical(cache$status,z$status),max(abs(cache$time-z$time))<1e-12)
 cluster<-if(db=="mimic")pr$group_id else hospital_components(pr$group_id,pr$hospitalid)$row_component
 idx<-match(cluster,unique(cluster));ng<-max(idx);set.seed(cache$seed)
 weights<-matrix(0L,ng,1000);for(b in 1:1000)weights[,b]<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
 for(b in 0:3){
  w<-if(b==0)rep(1,nrow(pr)) else weights[idx,b]
  keep<-w>0;d<-data.frame(t=z$time[keep],status=factor(z$status[keep],levels=0:2,labels=c("censored","death","alive_discharge")),w=w[keep])
  fit<-survival::survfit(survival::Surv(t,status)~1,data=d,weights=w,se.fit=FALSE,timefix=FALSE)
  jj<-findInterval(cache$grid-h/24,fit$time)
  point<-if(b==0)cache$curve_point else cache$curve_bootstrap[b,,]
  for(metric in c("death_cif","alive_discharge_cif","remaining_hospitalized")){
   if(metric=="remaining_hospitalized")value<-c(1,fit$pstate[,1])[jj+1] else value<-c(0,fit$pstate[,if(metric=="death_cif")"death" else "alive_discharge"])[jj+1]
   stopifnot(max(abs(value-point[,metric]))<1e-9)
  }
  km<-survival::survfit(survival::Surv(t,status=="death")~1,data=d,weights=w,se.fit=FALSE,timefix=FALSE)
  kk<-findInterval(cache$grid-h/24,km$time)
  stopifnot(max(abs((1-c(1,km$surv)[kk+1])-point[,"one_minus_km_censor_discharge"]))<1e-9)
  native_checks<-native_checks+1L
 }
 q<-cm[cm$window_hours==h & cm$dataset==ds,]
 for(metric in colnames(cache$curve_point))for(t in seq_along(cache$grid)){
  row<-q[q$metric==metric & q$icu_day==cache$grid[t],]
  ci<-quantile_ci(cache$curve_bootstrap[,t,metric])
  stopifnot(max(abs(c(row$lower,row$upper)-ci[c("lower","upper")]))<1e-10,row$valid_reps==ci['valid_reps'])
  cc<-cc+1L
 }
 for(v in names(cache$sensitivity_bootstrap)){
  q<-sm[sm$window_hours==h & sm$dataset==ds & sm$analysis==v,]
  for(metric in colnames(cache$sensitivity_bootstrap[[v]])){
   row<-q[q$metric==metric,];ci<-quantile_ci(cache$sensitivity_bootstrap[[v]][,metric])
   stopifnot(max(abs(c(row$lower,row$upper)-ci[c("lower","upper")]))<1e-10,row$valid_reps==ci['valid_reps'])
   sc<-sc+1L
  }
  yy<-if(v=="hospital_death_discharge_proxy")z$yproxy else z$y
  mask<-if(v=="exclude_icu_end_after_hospital")z$data$icu_end_after_hospital_flag==0 else rep(TRUE,nrow(pr))
  for(b in 1:3){
   w<-weights[idx[mask],b];ix<-rep(seq_len(sum(mask)),w)
   direct<-point_scores(yy[mask][ix],pr$probability[mask][ix]);metrics<-c("auroc","auprc","brier_score","log_loss","calibration_intercept","calibration_slope")
   stopifnot(max(abs(direct[metrics]-cache$sensitivity_bootstrap[[v]][b,metrics]))<2e-6)
  }
 }
 delta<-cache$sensitivity_bootstrap$hospital_death_discharge_proxy-cache$sensitivity_bootstrap$primary
 for(metric in colnames(delta)){
  row<-pm[pm$window_hours==h & pm$dataset==ds & pm$metric==metric,];ci<-quantile_ci(delta[,metric])
  stopifnot(max(abs(c(row$lower,row$upper)-ci[c("lower","upper")]))<1e-10)
  pc<-pc+1L
 }
 bin<-cut(pr$probability,c(0,.05,.10,.20,1),include.lowest=TRUE,right=FALSE,labels=c("<0.05","0.05-<0.10","0.10-<0.20",">=0.20"))
 for(g in levels(bin)){
  m<-bin==g;if(sum(m)<10)next
  mat<-matrix(NA_real_,1000,3,dimnames=list(NULL,c("death_cif","alive_discharge_cif","mean_predicted")))
  for(b in 1:1000){w<-weights[idx[m],b];mat[b,]<-c(sum(w*(z$status[m]==1)),sum(w*(z$status[m]==2)),sum(w*pr$probability[m]))/sum(w)}
  for(metric in colnames(mat)){
   row<-rm[rm$window_hours==h & rm$dataset==ds & rm$risk_bin==g & rm$metric==metric,];ci<-quantile_ci(mat[,metric])
   stopifnot(max(abs(c(row$lower,row$upper)-ci[c("lower","upper")]))<1e-10)
   bc<-bc+1L
  }
 }
}
# A tied-event boundary example: discharge is a competing event, not a death or early unknown censor.
synthetic<-curve_engine(c(1,1,2,3),c(1L,2L,1L,0L),c(1,2,1,1),c(0,1,2,3))
stopifnot(max(abs(synthetic[,"death_cif"]-c(0,.2,.4,.4)))<1e-12,
 max(abs(rowSums(synthetic[,c("death_cif","alive_discharge_cif","remaining_hospitalized")])-1))<1e-12)
write_json(list(status="passed",curve_intervals_recomputed=cc,sensitivity_intervals_recomputed=sc,
 paired_intervals_recomputed=pc,risk_bin_intervals_recomputed=bc,native_survival_point_and_weighted_checks=native_checks,
 explicit_sensitivity_bootstrap_checks=36,tied_event_boundary_case=TRUE),file.path(out,"independent_r_checks.json"))
log_progress("Endpoint native survival, interval, risk-bin and explicit bootstrap checks passed")
