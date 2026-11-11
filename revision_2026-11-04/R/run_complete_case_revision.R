if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1L] else file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_revision_uncertainty.R"),encoding="UTF-8")
out<-file.path(root,"complete_case_results");dir.create(out,showWarnings=FALSE)
contract<-jsonlite::read_json(file.path(root,"model_feature_contract.json"),simplifyVector=TRUE)
completeness<-function(data,h){
 raw<-setdiff(names(data$x),grep("_missing_",names(data$x),value=TRUE))
 measurement<-setdiff(raw,c("age_years","sex"));stopifnot(length(raw)==20,length(measurement)==18)
 for(v in measurement){
  flag<-sub("_(mean_24h|first_6h)$",paste0("_missing_",h,"h"),v)
  stopifnot(flag %in% names(data$x),all(data$x[[flag]] %in% 0:1),all(as.integer(is.na(data$x[[v]]))==data$x[[flag]]),
    !any(is.infinite(data$x[[v]])))
 }
 stopifnot(!any(is.infinite(data$x$age_years)))
 miss<-is.na(as.data.frame(data$x[,raw,drop=FALSE]));nm<-rowSums(miss)
 list(raw=raw,miss=miss,n_missing=nm,complete=nm==0)
}
precision_flag<-function(n,e,cl){
 x<-c(if(n<100)"N<100",if(e<20)"events<20",if(n-e<20)"non-events<20",if(cl<20)"clusters<20")
 if(length(x))paste(x,collapse="; ") else "No predefined small-sample flag (not proof of adequate precision)"
}
run_cc<-function(){
 write_json(list(date="2026-11-05",protocol_sha256=sha_file(file.path(root,"28_第十五步完整病例与缺失数据规则.md")),
  raw_inputs_required=20,continuous_measurements=18,replicates=1000,minimum_valid_reps=950,
  fitting=FALSE,recalibration=FALSE,validation_subset_only=TRUE),file.path(out,"analysis_plan.json"))
 counts<-miss_rates<-by_event<-distributions<-intervals<-list();np<-nf<-ne<-nd<-ni<-0L
 for(h in c(24L,6L)){
  for(db in c("mimic","eicu")){
   all<-read_data(file.path(root,"private_data",paste0(db,"_analysis_revision.csv")),contract,h)
   masks<-list()
   if(db=="mimic"){
    registry<-read.csv(file.path(root,"private_splits",sprintf("patients_%sh.csv",h)))
    part<-registry$partition[match(all$metadata$group_id,registry$group_id)];stopifnot(!anyNA(part))
    masks<-list(MIMIC_all_analysis=rep(TRUE,length(part)),MIMIC_training=part=="train",MIMIC_internal_test=part=="test")
   }else masks<-list(eICU_external_validation=rep(TRUE,length(all$y)))
   cc<-completeness(all,h)
   for(population in names(masks)){
    mask<-masks[[population]];yy<-all$y[mask]
    for(a in c("All records","Complete cases","Incomplete cases")){
     m<-mask & switch(a,"All records"=rep(TRUE,length(mask)),"Complete cases"=cc$complete,"Incomplete cases"=!cc$complete)
     ph<-all$metadata$phenotype[m];n<-sum(m);ev<-sum(all$y[m]);age<-all$x$age_years[m]
     np<-np+1L;counts[[np]]<-data.frame(window_hours=h,population=population,subset=a,records=n,
      unique_patients=length(unique(all$metadata$group_id[m])),events=ev,event_rate=if(n)ev/n else NA_real_,
      fraction_of_parent=n/sum(mask),age_median=if(n)median(age) else NA_real_,
      male_fraction=if(n)mean(all$x$sex[m]=="M",na.rm=TRUE) else NA_real_,
      ais_records=sum(ph=="AIS_STRICT"),ich_records=sum(ph=="ICH_STRICT"),sah_records=sum(ph=="SAH"),mixed_records=sum(grepl("^MIXED_",ph)))
    }
    for(v in cc$raw){
     nf<-nf+1L;miss_rates[[nf]]<-data.frame(window_hours=h,population=population,raw_input=v,records=sum(mask),missing=sum(cc$miss[mask,v]),missing_fraction=mean(cc$miss[mask,v]))
     for(e in 0:1){
      m<-mask & all$y==e;ne<-ne+1L;by_event[[ne]]<-data.frame(window_hours=h,population=population,raw_input=v,event=e,records=sum(m),missing=sum(cc$miss[m,v]),missing_fraction=if(sum(m))mean(cc$miss[m,v]) else NA_real_)
     }
    }
    nb<-cut(cc$n_missing[mask],breaks=c(-1,0,3,6,12,20),labels=c("0","1-3","4-6","7-12","13-20"))
    for(label in levels(nb)){
     mm<-nb==label;nd<-nd+1L;distributions[[nd]]<-data.frame(window_hours=h,population=population,missing_input_count=label,
      records=sum(mm),events=sum(yy[mm]),fraction=mean(mm))
    }
   }
  }
  for(dataset in c("MIMIC_internal_test","eICU_external_validation")){
   db<-if(dataset=="MIMIC_internal_test")"mimic" else "eicu"
   all<-read_data(file.path(root,"private_data",paste0(db,"_analysis_revision.csv")),contract,h)
   pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
   pred<-pred[pred$model=="xgboost" & pred$dataset==dataset,];pred<-pred[order(pred$stay_id),]
   ix<-match(pred$stay_id,all$metadata$stay_id);stopifnot(!anyNA(ix),!anyDuplicated(ix),all((pred$truth=="yes")==all$y[ix]))
   stopifnot(identical(as.character(pred$group_id),as.character(all$metadata$group_id[ix])))
   cc<-completeness(all,h);complete<-cc$complete[ix];y<-all$y[ix];p<-pred$probability
   cluster<-if(db=="mimic")pred$group_id else hospital_components(pred$group_id,pred$hospitalid)$row_component
   idx<-match(cluster,unique(cluster));ng<-max(idx)
   seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(db=="mimic")1L else 4L)*1000L
   set.seed(seed);weights<-matrix(0L,ng,1000L);for(b in 1:1000)weights[,b]<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
   cache<-list(seed=seed,stay_id=pred$stay_id,complete=complete,points=list(),bootstrap=list())
   for(a in c("All records","Complete cases","Incomplete cases")){
    m<-switch(a,"All records"=rep(TRUE,length(y)),"Complete cases"=complete,"Incomplete cases"=!complete)
    n<-sum(m);e<-sum(y[m]);cl<-length(unique(cluster[m]));warn<-precision_flag(n,e,cl)
    if(n<10){
     for(metric in metric_names){
      ni<-ni+1L;intervals[[ni]]<-data.frame(window_hours=h,dataset=dataset,subset=a,records=n,events=e,clusters=cl,
       metric=metric,estimate=NA_real_,lower=NA_real_,upper=NA_real_,raw_lower=NA_real_,raw_upper=NA_real_,valid_reps=0,
       status=if(n==0)"Empty subset" else "N<10: performance not reported",precision_warning=warn)
     }
     next
    }
    obj<-prepare_probability(y[m],p[m]);pt<-weighted_metrics(obj,rep(1,n))
    boot<-matrix(NA_real_,1000L,length(pt),dimnames=list(NULL,names(pt)))
    for(b in 1:1000)boot[b,]<-weighted_metrics(obj,weights[idx[m],b])
    cache$points[[a]]<-pt;cache$bootstrap[[a]]<-boot
    for(metric in names(pt)){
     ci<-quantile_ci(boot[,metric]);show<-is.finite(pt[metric]) && ci['valid_reps']>=950
     ni<-ni+1L;intervals[[ni]]<-data.frame(window_hours=h,dataset=dataset,subset=a,records=n,events=e,clusters=cl,
      metric=metric,estimate=unname(pt[metric]),lower=if(show)unname(ci['lower']) else NA_real_,upper=if(show)unname(ci['upper']) else NA_real_,
      raw_lower=unname(ci['lower']),raw_upper=unname(ci['upper']),valid_reps=unname(ci['valid_reps']),
      status=if(!is.finite(pt[metric]))"Point estimate not estimable" else if(!show)"Insufficient valid replicates" else "Conditional percentile interval",precision_warning=warn)
    }
    log_progress("Complete-case sensitivity: ",h,"h ",dataset," / ",a,"; N=",n,"; events=",e)
   }
   saveRDS(cache,file.path(root,"private_cache",sprintf("complete_case_%sh_%s.rds",h,dataset)))
  }
 }
 for(pair in list(c("complete_case_populations.csv","counts"),c("raw_input_missingness.csv","miss_rates"),c("missingness_by_event.csv","by_event"),
   c("missing_input_distribution.csv","distributions"),c("complete_case_metric_intervals.csv","intervals")))
  write.csv(do.call(rbind,get(pair[2])),file.path(out,pair[1]),row.names=FALSE)
 write_json(list(status="computed_pending_checks",replicates=1000,model_retraining=FALSE),file.path(out,"run_complete.json"))
}
if(sys.nframe()==0L)run_cc()
