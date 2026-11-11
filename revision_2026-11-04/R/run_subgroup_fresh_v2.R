if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1L] else file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_revision_uncertainty.R"),encoding="UTF-8")
out<-file.path(root,"phenotype_v2_results");dir.create(out,showWarnings=FALSE)
contract<-jsonlite::read_json(file.path(root,"model_feature_contract.json"),simplifyVector=TRUE)
phenotype_collapse<-function(x){
 allowed<-c("AIS_STRICT","ICH_STRICT","SAH","MIXED_PRIMARY_ACUTE","MIXED_WITH_AMBIGUOUS_SECONDARY","MIXED_WITH_CVST")
 stopifnot(!anyNA(x),all(x %in% allowed))
 ifelse(x %in% c("AIS_STRICT","ICH_STRICT","SAH"),x,"Mixed/ambiguous")
}
get_groups<-function(age,sex,phenotype){
 sex<-trimws(as.character(sex));sex[is.na(sex)|!sex %in% c("F","M")]<-"Unknown"
 ph<-phenotype_collapse(phenotype)
 list(
  list(family="Overall",label="All records",mask=rep(TRUE,length(age))),
  list(family="Age",label="<65",mask=age<65),list(family="Age",label=">=65",mask=age>=65),
  list(family="Sex",label="Female",mask=sex=="F"),list(family="Sex",label="Male",mask=sex=="M"),list(family="Sex",label="Unknown",mask=sex=="Unknown"),
  list(family="Phenotype",label="AIS_STRICT",mask=ph=="AIS_STRICT"),list(family="Phenotype",label="ICH_STRICT",mask=ph=="ICH_STRICT"),
  list(family="Phenotype",label="SAH",mask=ph=="SAH"),list(family="Phenotype",label="Mixed/ambiguous",mask=ph=="Mixed/ambiguous"),
  list(family="Restricted evaluation",label="Non-mixed AIS/ICH/SAH",mask=ph!="Mixed/ambiguous"),
  list(family="Restricted evaluation",label="ICH/SAH combined",mask=ph %in% c("ICH_STRICT","SAH")))
}
operating_names<-c("sensitivity","specificity","ppv","npv","flagged_fraction")
weighted_operating<-function(y,p,w){
 r<-p>=.10;tp<-sum(w*y*r);fn<-sum(w*y*(!r));fp<-sum(w*(1-y)*r);tn<-sum(w*(1-y)*(!r))
 ratio<-function(a,b)if(b>0)a/b else NA_real_
 setNames(c(ratio(tp,tp+fn),ratio(tn,tn+fp),ratio(tp,tp+fp),ratio(tn,tn+fn),ratio(tp+fp,sum(w))),operating_names)
}
run_subgroups<-function(){
 write_json(list(date="2026-11-05",protocol_sha256=sha_file(file.path(root,"32_第十七步CVST标签修正规则.md")),
  model_retraining=FALSE,calibration_updating=FALSE,replicates=1000,threshold=.10,age_cutoff=65,
  minimum_valid_replicates=950,cluster_bootstrap="whole validation cohort, then restrict to subgroup",
  exploratory=TRUE),file.path(out,"analysis_plan.json"))
 populations<-intervals<-compositions<-list();ni<-np<-nc<-0L
 for(h in c(24L,6L)){
  # Composition includes training and full development-source cohort; not extra model evaluation.
  for(db in c("mimic","eicu")){
   raw<-read_data(file.path(root,"private_data",if(db=="mimic")"mimic_analysis_cvst_v2.csv" else "eicu_analysis_revision.csv"),contract,h)
   partitions<-list(All_analysis_records=seq_along(raw$y))
   if(db=="mimic"){
    split<-read.csv(file.path(root,"private_splits",sprintf("patients_%sh.csv",h)))
    pt<-split$partition[match(raw$metadata$group_id,split$group_id)];stopifnot(!anyNA(pt))
    partitions$Training<-which(pt=="train");partitions$Internal_test<-which(pt=="test")
   }
   for(part in names(partitions)){
    ii<-partitions[[part]];ph<-phenotype_collapse(raw$metadata$phenotype[ii])
    for(g in c("AIS_STRICT","ICH_STRICT","SAH","Mixed/ambiguous")){
     jj<-ii[ph==g];nc<-nc+1L
     compositions[[nc]]<-data.frame(window_hours=h,database=db,partition=part,phenotype=g,records=length(jj),
       unique_patients=length(unique(raw$metadata$group_id[jj])),events=sum(raw$y[jj]),denominator=length(ii),
       record_fraction=length(jj)/length(ii))
    }
   }
  }
  for(dataset in c("MIMIC_internal_test","eICU_external_validation")){
   db<-if(dataset=="MIMIC_internal_test")"mimic" else "eicu"
   raw<-read_data(file.path(root,"private_data",if(db=="mimic")"mimic_analysis_cvst_v2.csv" else "eicu_analysis_revision.csv"),contract,h)
   pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh_cvst_v2.csv",h)))
   pred<-pred[pred$model=="xgboost" & pred$dataset==dataset,];pred<-pred[order(pred$stay_id),]
   ii<-match(pred$stay_id,raw$metadata$stay_id);stopifnot(!anyNA(ii),!anyDuplicated(ii))
   y<-as.integer(pred$truth=="yes");p<-pred$probability
   stopifnot(all(y==raw$y[ii]),all(as.character(pred$group_id)==as.character(raw$metadata$group_id[ii])),all(pred$phenotype==raw$metadata$phenotype[ii]))
   groups<-get_groups(raw$x$age_years[ii],raw$x$sex[ii],pred$phenotype)
   cluster<-if(db=="mimic")pred$group_id else hospital_components(pred$group_id,pred$hospitalid)$row_component
   idx<-match(cluster,unique(cluster));ng<-max(idx)
   seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(db=="mimic")1L else 4L)*1000L
   set.seed(seed);weights<-matrix(0L,ng,1000L)
   for(b in 1:1000)weights[,b]<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
   cache<-list(seed=seed,groups=list(),records=nrow(pred),clusters=ng,scheme=if(db=="mimic")"patient" else "hospital_patient_component")
   for(j in seq_along(groups)){
    g<-groups[[j]];m<-g$mask;stopifnot(!anyNA(m));n<-sum(m);e<-sum(y[m]);ne<-n-e
    cluster_n<-length(unique(cluster[m]));status<-if(n==0)"Empty subgroup" else if(n<10)"Small cell suppressed" else "Reported"
    warning<-if(status!="Reported")status else paste(c(if(n<100)"N<100",if(e<20)"events<20",if(ne<20)"non-events<20",if(cluster_n<20)"clusters<20"),collapse="; ")
    if(!nchar(warning))warning<-"No predefined small-sample flag (not proof of adequate precision)"
    positive<-p[m]>=.10
    row<-data.frame(window_hours=h,dataset=dataset,family=g$family,subgroup=g$label,records=n,unique_patients=length(unique(pred$group_id[m])),
     events=e,non_events=ne,event_rate=if(n)e/n else NA_real_,hospitals=if(db=="eicu")length(unique(pred$hospitalid[m])) else NA_integer_,
     clusters=cluster_n,tp=sum(y[m]*positive),fp=sum((1-y[m])*positive),tn=sum((1-y[m])*(!positive)),fn=sum(y[m]*(!positive)),
     status=status,precision_warning=warning)
    # Suppress detailed public data for rare nonempty groups; private cache preserves audit.
    if(status=="Small cell suppressed")row[,c("records","unique_patients","events","non_events","event_rate","hospitals","clusters","tp","fp","tn","fn")]<-NA
    np<-np+1L;populations[[np]]<-row
    if(status!="Reported")next
    object<-prepare_probability(y[m],p[m]);points<-c(weighted_metrics(object,rep(1,n)),weighted_operating(y[m],p[m],rep(1,n)))
    # Single-class AUPRC is non-informative even where software would return 1.
    if(e==0 || ne==0)points[c("auroc","auprc","calibration_intercept","calibration_slope")]<-NA_real_
    boot<-matrix(NA_real_,1000,length(points),dimnames=list(NULL,names(points)))
    for(b in 1:1000){
     w<-weights[idx[m],b];boot[b,]<-c(weighted_metrics(object,w),weighted_operating(y[m],p[m],w))
    }
    for(metric in names(points)){
     ci<-quantile_ci(boot[,metric]);show<-ci['valid_reps']>=950 && is.finite(points[metric])
     ni<-ni+1L;intervals[[ni]]<-data.frame(window_hours=h,dataset=dataset,family=g$family,subgroup=g$label,metric=metric,
      estimate=unname(points[metric]),lower=if(show)unname(ci['lower']) else NA_real_,upper=if(show)unname(ci['upper']) else NA_real_,
      raw_lower=unname(ci['lower']),raw_upper=unname(ci['upper']),valid_reps=unname(ci['valid_reps']),
      interval_status=if(!is.finite(points[metric]))"Point estimate not estimable" else if(!show)"Insufficient valid replicates" else "Conditional percentile interval",
      precision_warning=warning)
    }
    cache$groups[[paste(g$family,g$label,sep="::")]]<-list(stay_id=pred$stay_id[m],truth=y[m],p=p[m],point=points,bootstrap=boot)
    log_progress(h,"h ",dataset," / ",g$family," ",g$label,"; N=",n,"; events=",e)
   }
   saveRDS(cache,file.path(root,"private_cache",sprintf("subgroup_fresh_v2_%sh_%s.rds",h,dataset)))
  }
 }
 write.csv(do.call(rbind,populations),file.path(out,"subgroup_populations.csv"),row.names=FALSE)
 write.csv(do.call(rbind,intervals),file.path(out,"subgroup_metric_intervals.csv"),row.names=FALSE)
 write.csv(do.call(rbind,compositions),file.path(out,"phenotype_case_mix.csv"),row.names=FALSE)
 write_json(list(status="computed_pending_checks",replicates=1000,model_retraining=FALSE),file.path(out,"run_complete.json"))
}
if(sys.nframe()==0L)run_subgroups()
