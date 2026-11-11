if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1L] else file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_revision_uncertainty.R"),encoding="UTF-8")
out<-file.path(root,"patientflow_results");dir.create(out,showWarnings=FALSE)
select_patient_weights<-function(pred,h,ds){
 priority<-vapply(pred$stay_id,function(id)digest::digest(paste("20261105",h,ds,id,sep="|"),algo="sha256",serialize=FALSE),character(1))
 ord<-order(as.character(pred$group_id),priority,pred$stay_id)
 selected<-rep(FALSE,nrow(pred));selected[ord[!duplicated(as.character(pred$group_id)[ord])]]<-TRUE
 count<-table(as.character(pred$group_id));equal<-1/as.numeric(count[as.character(pred$group_id)])
 stopifnot(sum(selected)==length(count),max(abs(tapply(equal,pred$group_id,sum)-1))<1e-12)
 list(primary=rep(1,nrow(pred)),one_record_per_patient=as.numeric(selected),equal_patient_weight=equal,priority=priority)
}
run_repeat<-function(){
 write_json(list(date="2026-11-05",protocol_sha256=sha_file(file.path(root,"26_第十四步重复住院与意见对照规则.md")),
   replicates=1000,selection="minimum UTF-8 SHA256(20261105|window_hours|dataset|stay_id) within patient",
   selection_repeated_in_bootstrap=FALSE,retraining=FALSE,calibration_updating=FALSE,
   target="fixed-model validation sensitivity, not first-admission model development"),file.path(out,"analysis_plan.json"))
 pops<-metrics<-list();np<-nm<-0L
 for(h in c(24L,6L))for(ds in c("MIMIC_internal_test","eICU_external_validation")){
  pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
  pred<-pred[pred$model=="xgboost" & pred$dataset==ds,];pred<-pred[order(pred$stay_id),]
  y<-as.integer(pred$truth=="yes");p<-pred$probability
  cluster<-if(ds=="MIMIC_internal_test")pred$group_id else hospital_components(pred$group_id,pred$hospitalid)$row_component
  ix<-match(cluster,unique(cluster));ng<-max(ix)
  seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(ds=="MIMIC_internal_test")1L else 4L)*1000L
  set.seed(seed);bw<-matrix(0L,ng,1000L);for(b in 1:1000)bw[,b]<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
  fixed<-select_patient_weights(pred,h,ds)
  write.csv(data.frame(stay_id=pred$stay_id,group_id=pred$group_id,priority=fixed$priority,
    selected=fixed$one_record_per_patient,patient_weight=fixed$equal_patient_weight),
    file.path(root,"private_splits",sprintf("patient_selection_%sh_%s.csv",h,ds)),row.names=FALSE)
  obj<-prepare_probability(y,p);cache<-list(seed=seed,stay_id=pred$stay_id,points=list(),bootstrap=list())
  for(a in c("primary","one_record_per_patient","equal_patient_weight")){
   w<-fixed[[a]];m<-w>0;pt<-weighted_metrics(obj,w)
   boot<-matrix(NA_real_,1000L,length(pt),dimnames=list(NULL,names(pt)))
   for(b in 1:1000)boot[b,]<-weighted_metrics(obj,w*bw[ix,b])
   np<-np+1L;pops[[np]]<-data.frame(window_hours=h,dataset=ds,analysis=a,records=sum(m),unique_patients=length(unique(pred$group_id[m])),
    raw_events=sum(y[m]),weight_sum=sum(w),weighted_event_rate=sum(w*y)/sum(w),
    clusters_represented=length(unique(cluster[m])),full_cohort_resampling_clusters=ng,
    records_removed=sum(!m),patient_weights_fixed=TRUE)
   for(metric in names(pt)){
    ci<-quantile_ci(boot[,metric]);show<-is.finite(pt[metric]) && ci['valid_reps']>=950;nm<-nm+1L
    metrics[[nm]]<-data.frame(window_hours=h,dataset=ds,analysis=a,metric=metric,estimate=unname(pt[metric]),
      lower=if(show)unname(ci['lower']) else NA_real_,upper=if(show)unname(ci['upper']) else NA_real_,valid_reps=unname(ci['valid_reps']))
   }
   cache$points[[a]]<-pt;cache$bootstrap[[a]]<-boot
   log_progress("Repeat-admission sensitivity: ",h,"h ",ds," / ",a,"; records=",sum(m),"; patients=",length(unique(pred$group_id[m])))
  }
  saveRDS(cache,file.path(root,"private_cache",sprintf("patient_repeat_%sh_%s.rds",h,ds)))
 }
 write.csv(do.call(rbind,pops),file.path(out,"repeat_admission_populations.csv"),row.names=FALSE)
 write.csv(do.call(rbind,metrics),file.path(out,"repeat_admission_metric_intervals.csv"),row.names=FALSE)
 write_json(list(status="computed_pending_checks",model_retraining=FALSE),file.path(out,"run_complete.json"))
}
if(sys.nframe()==0L)run_repeat()
