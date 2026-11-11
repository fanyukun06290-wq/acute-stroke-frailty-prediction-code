if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
root<-file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_complete_case_revision.R"),encoding="UTF-8")
met<-read.csv(file.path(out,"complete_case_metric_intervals.csv"));checked<-replicated<-0L;medians<-reconstruction<-list()
for(h in c(24L,6L)){
 raw<-read_data(file.path(root,"private_data/mimic_analysis_revision.csv"),contract,h)
 reg<-read.csv(file.path(root,"private_splits",sprintf("patients_%sh.csv",h)))
 training<-reg$partition[match(raw$metadata$group_id,reg$group_id)]=="train"
 recipe<-readRDS(file.path(root,"private_models",sprintf("recipe_%sh.rds",h)))
 step<-recipe$steps[[which(vapply(recipe$steps,function(s)inherits(s,"step_impute_median"),logical(1)))]]
 stopifnot(isTRUE(step$trained),length(step$medians)==37L)
 for(v in names(step$medians)){
  actual<-median(raw$x[[v]][training],na.rm=TRUE);saved<-unname(step$medians[[v]])
  stopifnot(abs(actual-saved)<1e-12)
  medians[[length(medians)+1L]]<-data.frame(window_hours=h,input=v,kind=if(grepl("_missing_",v))"Derived missingness flag" else "Raw numeric input",
   training_records=sum(training),observed_training_values=sum(!is.na(raw$x[[v]][training])),saved_training_median=saved,
   independent_training_median=actual,absolute_error=abs(actual-saved))
 }
 for(ds in c("MIMIC_internal_test","eICU_external_validation")){
  db<-if(ds=="MIMIC_internal_test")"mimic" else "eicu"
  data<-if(db=="mimic")raw else read_data(file.path(root,"private_data/eicu_analysis_revision.csv"),contract,h)
  p<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
  p<-p[p$model=="xgboost" & p$dataset==ds,];p<-p[order(p$stay_id),]
  ix<-match(p$stay_id,data$metadata$stay_id);cc<-completeness(data,h)$complete[ix]
  cache<-readRDS(file.path(root,"private_cache",sprintf("complete_case_%sh_%s.rds",h,ds)))
  stopifnot(identical(cache$stay_id,p$stay_id),identical(cache$complete,cc))
  # Confirm current raw inputs reproduce archived predictions under the frozen training recipe.
  x<-as.matrix(recipes::bake(recipe,new_data=data$x[ix,,drop=FALSE]));enc<-read.csv(file.path(root,"model_results",sprintf("encoded_predictors_%sh.csv",h)))[[1]]
  stopifnot(identical(colnames(x),enc),!anyNA(x))
  model<-xgboost::xgb.load(file.path(root,"private_models",sprintf("xgboost_%sh.ubj",h)))
  probability<-as.numeric(predict(model,xgboost::xgb.DMatrix(x,nthread=2)));error<-max(abs(probability-p$probability));stopifnot(error<1e-6)
  reconstruction[[length(reconstruction)+1L]]<-data.frame(window_hours=h,dataset=ds,records=nrow(p),maximum_frozen_probability_error=error)
  y<-as.integer(p$truth=="yes");cluster<-if(db=="mimic")p$group_id else hospital_components(p$group_id,p$hospitalid)$row_component
  idx<-match(cluster,unique(cluster));ng<-max(idx);set.seed(cache$seed);draws<-lapply(1:3,function(b)sample.int(ng,ng,replace=TRUE))
  for(a in names(cache$bootstrap)){
   q<-met[met$window_hours==h & met$dataset==ds & met$subset==a,]
   for(metric in colnames(cache$bootstrap[[a]])){
    ci<-quantile_ci(cache$bootstrap[[a]][,metric]);r<-q[q$metric==metric,]
    stopifnot(r$valid_reps==ci['valid_reps'],max(abs(c(r$raw_lower,r$raw_upper)-ci[c('lower','upper')]),na.rm=TRUE)<1e-10)
    if(is.finite(r$estimate) && r$valid_reps>=950)stopifnot(max(abs(c(r$lower,r$upper)-ci[c('lower','upper')]))<1e-10)
    checked<-checked+1L
   }
   m<-switch(a,"All records"=rep(TRUE,length(y)),"Complete cases"=cc,"Incomplete cases"=!cc)
   for(b in 1:3){
    w<-tabulate(draws[[b]],nbins=ng)[idx[m]];rr<-rep(seq_len(sum(m)),w)
    if(length(unique(y[m][rr]))<2)next
    direct<-point_scores(y[m][rr],p$probability[m][rr]);metrics<-c("auroc","auprc","brier_score","log_loss","calibration_intercept","calibration_slope")
    stopifnot(max(abs(direct[metrics]-cache$bootstrap[[a]][b,metrics]))<2e-6)
    replicated<-replicated+1L
   }
  }
 }
}
write.csv(do.call(rbind,medians),file.path(out,"training_imputation_medians.csv"),row.names=FALSE)
write.csv(do.call(rbind,reconstruction),file.path(out,"frozen_prediction_reconstruction.csv"),row.names=FALSE)
write_json(list(status="passed",intervals_recomputed=checked,explicit_bootstrap_checks=replicated,training_medians_verified=length(medians),
 frozen_prediction_sets=length(reconstruction),training_only_final_medians=TRUE,
 scope="Final recipe medians verified against original training partition; fold-only fitting supported by protected original source, not a new replay of all CV fits"),file.path(out,"independent_r_checks.json"))
log_progress("Complete-case checks passed; intervals=",checked,"; saved training medians=",length(medians),"; explicit bootstrap checks=",replicated)
