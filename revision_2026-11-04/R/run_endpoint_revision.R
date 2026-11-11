if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1] else file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_revision_uncertainty.R"),encoding="UTF-8")
out<-file.path(root,"endpoint_results");dir.create(out,showWarnings=FALSE)

read_endpoint<-function(db,h,dataset=NULL){
 raw<-read.csv(file.path(root,"private_data",paste0(db,"_cohort.csv")),stringsAsFactors=FALSE)
 elig<-raw[[paste0("landmark_",h,"h_eligible")]]==1
 keep<-elig & raw$outcome_missing_flag==0
 data<-raw[keep,,drop=FALSE]
 if(!is.null(dataset)){
  pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
  pred<-pred[pred$model=="xgboost" & pred$dataset==dataset,];pred<-pred[order(pred$stay_id),]
  idx<-match(pred$stay_id,data$stay_id);stopifnot(!anyNA(idx),!anyDuplicated(idx))
  data<-data[idx,,drop=FALSE]
  stopifnot(identical(as.character(data$patient_group_id),as.character(pred$group_id)))
 }else pred<-NULL
 expired<-if(db=="mimic")data$hospital_expire_flag==1 else data$hospital_status=="expired"
 stopifnot(!anyNA(expired),!anyNA(data$hospital_discharge_offset_min),all(data$hospital_discharge_offset_min>h*60),
  !anyNA(data$effective_icu_end_min),all(data$effective_icu_end_min>h*60),!anyNA(data$icu_end_after_hospital_flag))
 death<-data$death_offset_min;discharge<-data$hospital_discharge_offset_min
 stopifnot(all(is.finite(death[expired])),all(death[expired]>h*60),all(death[expired]<=discharge[expired]))
 y<-as.integer(expired & !is.na(death) & death>h*60 & death<43200)
 stored<-data[[sprintf("death_%sh_to_30d_inhospital",h)]]
 stopifnot(all(y==stored),!anyNA(stored))
 if(!is.null(pred))stopifnot(all(y==as.integer(pred$truth=="yes")))
 status<-ifelse(expired & death<43200,1L,ifelse(!expired & discharge<43200,2L,0L))
 stopifnot(!anyNA(status))
 # All nonterminal observations are administratively followed to day 30.
 exit<-ifelse(expired,death,discharge)
 time<-(pmin(exit,43200)-h*60)/1440
 stopifnot(all(time>0),all(time<=30-h/24),all(status==1L)==FALSE)
 yproxy<-as.integer(expired & discharge>h*60 & discharge<43200)
 list(data=data,pred=pred,y=y,yproxy=yproxy,expired=expired,time=time,status=status,
  source_summary=data.frame(window_hours=h,database=db,source_records=nrow(raw),landmark_eligible=sum(elig),
   eligible_outcome_unknown=sum(elig & raw$outcome_missing_flag==1),analysis_records=sum(keep),
   source_outcome_unknown=sum(raw$outcome_missing_flag==1)))
}

# General weighted AJ recursion and diagnostic KM, with administrative censoring at day 30.
curve_engine<-function(time,status,weights,grid){
 tt<-sort(unique(time));ti<-match(time,tt)
 counts<-rowsum(cbind(all=weights,death=weights*(status==1),alive=weights*(status==2)),ti,reorder=TRUE)
 stopifnot(nrow(counts)==length(tt))
 n<-sum(weights);risk<-n-c(0,head(cumsum(counts[,"all"]),-1L))
 ratio<-function(x)ifelse(risk>0,x/pmax(risk,1),0)
 s<-cumprod(pmax(0,1-ratio(counts[,"death"]+counts[,"alive"])))
 prev<-c(1,head(s,-1L));fd<-cumsum(prev*ratio(counts[,"death"]));fa<-cumsum(prev*ratio(counts[,"alive"]))
 km<-1-cumprod(pmax(0,1-ratio(counts[,"death"])))
 # With no earlier censoring AJ equals empirical cumulative incidence, including weighted bootstrap.
 stopifnot(max(abs(fd-cumsum(counts[,"death"])/n))<1e-9,max(abs(fa-cumsum(counts[,"alive"])/n))<1e-9)
 idx<-findInterval(grid,tt)
 cbind(death_cif=c(0,fd)[idx+1L],alive_discharge_cif=c(0,fa)[idx+1L],
  remaining_hospitalized=c(1,s)[idx+1L],one_minus_km_censor_discharge=c(0,km)[idx+1L])
}

run_endpoints<-function(){
 plan<-list(date="2026-11-05",protocol_sha256=sha_file(file.path(root,"24_第十三步结局与出院机制规则.md")),
  replicates=1000,landmark_hours=c(24,6),event_horizon="strictly before ICU-admission day 30",risk_cuts=c(0,.05,.10,.20,1),
  early_censoring="none among analyzed records; alive discharge is absorbing competing event",model_retraining=FALSE,
  sensitivities=c("hospital_death_discharge_proxy","exclude_icu_end_after_hospital"),survival_version=as.character(packageVersion("survival")))
 write_json(plan,file.path(out,"analysis_plan.json"))
 src<-pop<-curves<-risks<-sens<-paired<-list();ns<-np<-nc<-nr<-nm<-nd<-0L
 for(h in c(24L,6L)){
  for(db in c("mimic","eicu")){z<-read_endpoint(db,h);ns<-ns+1L;src[[ns]]<-z$source_summary}
  for(dataset in c("MIMIC_internal_test","eICU_external_validation")){
   db<-if(dataset=="MIMIC_internal_test")"mimic" else "eicu"
   z<-read_endpoint(db,h,dataset);d<-z$data;p<-z$pred$probability;n<-nrow(d)
   cluster<-if(db=="mimic")z$pred$group_id else hospital_components(z$pred$group_id,z$pred$hospitalid)$row_component
   idx<-match(cluster,unique(cluster));ng<-max(idx)
   seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(db=="mimic")1L else 4L)*1000L
   set.seed(seed);wmat<-matrix(0L,ng,1000)
   for(b in 1:1000)wmat[,b]<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
   grid<-seq(h/24,30,by=.25);relative_grid<-grid-h/24
   point<-curve_engine(z$time,z$status,rep(1,n),relative_grid)
   stopifnot(abs(tail(point[,"death_cif"],1)-mean(z$y))<1e-10)
   native_data<-data.frame(t=z$time,status=factor(z$status,levels=0:2,labels=c("censored","death","alive_discharge")))
   native<-survival::survfit(survival::Surv(t,status)~1,data=native_data,se.fit=FALSE,timefix=FALSE)
   ni<-findInterval(relative_grid,native$time)
   native_cif<-c(0,native$pstate[,"death"])[ni+1]
   err<-max(abs(point[,"death_cif"]-native_cif));stopifnot(err<1e-9)
   row<-data.frame(window_hours=h,dataset=dataset,records=n,unique_patients=length(unique(z$pred$group_id)),clusters=ng,
    death_before_day30=sum(z$status==1),alive_discharge_before_day30=sum(z$status==2),still_in_hospital_day30=sum(z$status==0),
    death_rate=mean(z$y),hospital_death_proxy_events=sum(z$yproxy),changed_endpoint_labels=sum(z$y!=z$yproxy),
    death_exact_day30=sum(z$expired & d$death_offset_min==43200,na.rm=TRUE),discharge_exact_day30=sum(d$hospital_discharge_offset_min==43200),
    icu_end_after_hospital=sum(d$icu_end_after_hospital_flag==1),native_aj_error=err,
    no_censoring_before_horizon=TRUE,death_time_source=if(db=="mimic")"Recorded deathtime" else "Unit/hospital death-discharge proxy",
    unit_death_proxy=if(db=="eicu")sum(d$death_time_source=="unit_death_discharge",na.rm=TRUE) else NA_integer_,
    hospital_death_proxy=if(db=="eicu")sum(d$death_time_source=="hospital_death_discharge",na.rm=TRUE) else NA_integer_,
    death_to_discharge_median_hours=median((d$hospital_discharge_offset_min-d$death_offset_min)[z$expired]/60),
    hospital_remaining_median_days=median(d$hospital_discharge_offset_min/1440-h/24))
   np<-np+1L;pop[[np]]<-row
   boot<-array(NA_real_,c(1000,nrow(point),ncol(point)),dimnames=list(NULL,NULL,colnames(point)))
   for(b in 1:1000)boot[b,,]<-curve_engine(z$time,z$status,wmat[idx,b],relative_grid)
   for(metric in colnames(point))for(t in seq_len(nrow(point))){
    ci<-quantile_ci(boot[,t,metric]);nc<-nc+1L
    curves[[nc]]<-data.frame(window_hours=h,dataset=dataset,icu_day=grid[t],metric=metric,estimate=point[t,metric],
      lower=unname(ci['lower']),upper=unname(ci['upper']),valid_reps=unname(ci['valid_reps']))
   }
   # Fixed risk strata and paired resampling, no external threshold optimization.
   bin<-cut(p,breaks=c(0,.05,.10,.20,1),include.lowest=TRUE,right=FALSE,labels=c("<0.05","0.05-<0.10","0.10-<0.20",">=0.20"))
   for(g in levels(bin)){
    mask<-bin==g;nn<-sum(mask);if(nn<10)next
    bp<-c(death_cif=mean(z$status[mask]==1),alive_discharge_cif=mean(z$status[mask]==2),mean_predicted=mean(p[mask]))
    br<-matrix(NA_real_,1000,3,dimnames=list(NULL,names(bp)))
    for(b in 1:1000){w<-wmat[idx[mask],b];den<-sum(w);if(den>0)br[b,]<-c(sum(w*(z$status[mask]==1)),sum(w*(z$status[mask]==2)),sum(w*p[mask]))/den}
    for(metric in names(bp)){
     ci<-quantile_ci(br[,metric]);nr<-nr+1L;risks[[nr]]<-data.frame(window_hours=h,dataset=dataset,risk_bin=g,records=nn,
      events=sum(z$y[mask]),metric=metric,estimate=bp[metric],lower=unname(ci['lower']),upper=unname(ci['upper']),valid_reps=unname(ci['valid_reps']))
    }
   }
   variants<-list(primary=list(y=z$y,mask=rep(TRUE,n)),hospital_death_discharge_proxy=list(y=z$yproxy,mask=rep(TRUE,n)),
    exclude_icu_end_after_hospital=list(y=z$y,mask=d$icu_end_after_hospital_flag==0))
   vb<-vp<-list()
   for(v in names(variants)){
    vv<-variants[[v]];m<-vv$mask;obj<-prepare_probability(vv$y[m],p[m]);pt<-weighted_metrics(obj,rep(1,sum(m)))
    bmat<-matrix(NA_real_,1000,length(pt),dimnames=list(NULL,names(pt)))
    for(b in 1:1000)bmat[b,]<-weighted_metrics(obj,wmat[idx[m],b])
    vb[[v]]<-bmat;vp[[v]]<-pt
    for(metric in names(pt)){
     ci<-quantile_ci(bmat[,metric]);show<-ci['valid_reps']>=950 && is.finite(pt[metric]);nm<-nm+1L
     sens[[nm]]<-data.frame(window_hours=h,dataset=dataset,analysis=v,records=sum(m),events=sum(vv$y[m]),metric=metric,estimate=pt[metric],
      lower=if(show)unname(ci['lower']) else NA_real_,upper=if(show)unname(ci['upper']) else NA_real_,valid_reps=unname(ci['valid_reps']))
    }
   }
   delta<-vb$hospital_death_discharge_proxy-vb$primary
   for(metric in colnames(delta)){
    ci<-quantile_ci(delta[,metric]);nd<-nd+1L
    paired[[nd]]<-data.frame(window_hours=h,dataset=dataset,comparison="hospital_proxy_minus_primary",metric=metric,
      difference=vp$hospital_death_discharge_proxy[metric]-vp$primary[metric],lower=unname(ci['lower']),upper=unname(ci['upper']),valid_reps=unname(ci['valid_reps']))
   }
   saveRDS(list(seed=seed,stay_id=z$pred$stay_id,time=z$time,status=z$status,grid=grid,curve_point=point,curve_bootstrap=boot,
    sensitivity_points=vp,sensitivity_bootstrap=vb),file.path(root,"private_cache",sprintf("endpoint_%sh_%s.rds",h,dataset)))
   log_progress("Endpoint analysis complete: ",h,"h ",dataset,"; N=",n,"; competing discharges=",sum(z$status==2))
  }
 }
 for(x in list(c("source_cohort_audit.csv","src"),c("endpoint_populations.csv","pop"),c("competing_event_curves.csv","curves"),
  c("fixed_risk_bin_cif.csv","risks"),c("endpoint_sensitivity_intervals.csv","sens"),c("endpoint_proxy_paired_differences.csv","paired")))
  write.csv(do.call(rbind,get(x[2])),file.path(out,x[1]),row.names=FALSE)
 write_json(list(status="computed_pending_checks",replicates=1000,model_retraining=FALSE),file.path(out,"run_complete.json"))
}
if(sys.nframe()==0L)run_endpoints()
