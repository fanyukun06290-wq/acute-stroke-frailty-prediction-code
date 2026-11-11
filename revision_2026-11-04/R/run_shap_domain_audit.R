# Local aggregate explanation audit. No fitting and no patient rows in console.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)) args[1] else file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_revision_uncertainty.R"),encoding="UTF-8")
out<-file.path(root,"shap_domain_results");dir.create(out,showWarnings=FALSE)
contract<-jsonlite::read_json(file.path(root,"model_feature_contract.json"),simplifyVector=TRUE)
domain_map<-function(features){
 v<-sub("_(mean_24h|first_6h|missing_24h|missing_6h)$","",features)
 v[grepl("^sex_",v)]<-"sex"
 d<-ifelse(v %in% c("age_years","sex"),"Demographic",
    ifelse(v %in% c("glucose","sodium","potassium","chloride","bicarbonate","bun","creatinine","calcium"),"Metabolic/renal",
    ifelse(v %in% c("wbc","platelets","hemoglobin","hematocrit"),"Inflammation-related hematologic",
    ifelse(v %in% c("heart_rate","respiratory_rate","spo2","sbp","dbp","map"),"Vital signs",NA))))
 stopifnot(!anyNA(d))
 data.frame(feature=features,variable=v,domain=d,kind=ifelse(grepl("_missing_",features),"Missingness indicator","Measured/demographic"))
}
load_validation<-function(h,dataset){
 db<-if(dataset=="MIMIC_internal_test") "mimic" else "eicu"
 dat<-read_data(file.path(root,"private_data",paste0(db,"_analysis_revision.csv")),contract,h)
 pred<-read.csv(file.path(root,"private_predictions",sprintf("predictions_%sh.csv",h)))
 pred<-pred[pred$model=="xgboost" & pred$dataset==dataset,]
 pred<-pred[order(pred$stay_id),]
 idx<-match(pred$stay_id,dat$metadata$stay_id);stopifnot(!anyNA(idx),!anyDuplicated(idx))
 stopifnot(identical(as.character(pred$group_id),as.character(dat$metadata$group_id[idx])),all((pred$truth=="yes")==dat$y[idx]))
 if(db=="mimic"){
  reg<-read.csv(file.path(root,"private_splits",sprintf("patients_%sh.csv",h)))
  expected<-reg$group_id[reg$partition=="test"]
  stopifnot(setequal(pred$stay_id,dat$metadata$stay_id[dat$metadata$group_id %in% expected]))
 }else stopifnot(nrow(pred)==nrow(dat$x))
 rec<-readRDS(file.path(root,"private_models",sprintf("recipe_%sh.rds",h)))
 x<-as.matrix(recipes::bake(rec,new_data=dat$x[idx,,drop=FALSE]));storage.mode(x)<-"double"
 enc<-read.csv(file.path(root,"model_results",sprintf("encoded_predictors_%sh.csv",h)))[[1]]
 stopifnot(identical(colnames(x),enc),!anyNA(x),all(is.finite(x)))
 list(x=x,y=dat$y[idx],pred=pred)
}
run_audit<-function(){
 plan<-list(date="2026-11-05",bootstrap_replicates=1000,permutation_replicates=50,
   native_exact_tree_shap=TRUE,unit="raw margin (log-odds)",background="native training tree path covers",
   fitting=FALSE,calibration_updating=FALSE,group_definition="sum of existing feature attributions, not recomputed coalition SHAP",
   protocol_sha256=sha_file(file.path(root,"20_第十一步SHAP与领域解释规则.md")))
 write_json(plan,file.path(out,"analysis_plan.json"))
 summaries<-maps<-checks<-permutations<-list();si<-mi<-ci<-pi<-0L
 for(h in c(24L,6L)){
  model<-xgboost::xgb.load(file.path(root,"private_models",sprintf("xgboost_%sh.ubj",h)))
  for(dataset in c("MIMIC_internal_test","eICU_external_validation")){
   z<-load_validation(h,dataset);dm<-xgboost::xgb.DMatrix(z$x,nthread=2L)
   phi_all<-predict(model,dm,predcontrib=TRUE,approxcontrib=FALSE)
   margin<-as.numeric(predict(model,dm,outputmargin=TRUE));p<-as.numeric(predict(model,dm))
   stopifnot(ncol(phi_all)==ncol(z$x)+1L)
   phi<-phi_all[,seq_len(ncol(z$x)),drop=FALSE];colnames(phi)<-colnames(z$x)
   bias<-phi_all[,ncol(phi_all)]
   e1<-max(abs(rowSums(phi)+bias-margin));e2<-max(abs(p-z$pred$probability))
   stopifnot(e1<2e-4,e2<1e-6,max(abs(plogis(margin)-p))<1e-6, diff(range(bias))<1e-6)
   map<-domain_map(colnames(z$x));map$window_hours<-h;mi<-mi+1L;maps[[mi]]<-map
   cluster<-if(dataset=="MIMIC_internal_test") z$pred$group_id else hospital_components(z$pred$group_id,z$pred$hospitalid)$row_component
   cluster_idx<-match(cluster,unique(cluster));ng<-max(cluster_idx)
   a<-abs(phi);labels<-data.frame(level="Encoded feature",name=map$feature,metric="mean_absolute",domain=map$domain)
   append_stat<-function(values,level,name,metric,domain){
    a<<-cbind(a,values)
    labels<<-rbind(labels,data.frame(level=level,name=name,metric=metric,domain=domain))
   }
   for(v in unique(map$variable)){
    cols<-which(map$variable==v)
    append_stat(abs(rowSums(phi[,cols,drop=FALSE])),"Conceptual variable",v,"absolute_signed_sum",map$domain[cols[1]])
   }
   for(d in unique(map$domain)){
    cols<-which(map$domain==d);nv<-length(unique(map$variable[cols]))
    append_stat(abs(rowSums(phi[,cols,drop=FALSE])),"Domain",d,"absolute_signed_sum",d)
    append_stat(rowSums(abs(phi[,cols,drop=FALSE])),"Domain",d,"sum_absolute",d)
    append_stat(rowSums(abs(phi[,cols,drop=FALSE]))/length(cols),"Domain",d,"sum_absolute_per_encoded_feature",d)
    append_stat(rowSums(abs(phi[,cols,drop=FALSE]))/nv,"Domain",d,"sum_absolute_per_conceptual_variable",d)
   }
   for(k in unique(map$kind))append_stat(rowSums(abs(phi[,map$kind==k,drop=FALSE])),"Encoding kind",k,"sum_absolute","All domains")
   cs<-rowsum(a,cluster_idx,reorder=FALSE);cn<-as.numeric(rowsum(rep(1,nrow(a)),cluster_idx,reorder=FALSE))
   # rowsum first-occurrence order matches cluster_idx=1..G.
   stopifnot(identical(as.integer(rownames(cs)),seq_len(ng)))
   seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(dataset=="MIMIC_internal_test")1L else 4L)*1000L
   set.seed(seed);boot<-matrix(NA_real_,1000,ncol(a))
   for(b in seq_len(1000)){
    w<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
    boot[b,]<-as.numeric(crossprod(w,cs))/sum(w*cn)
   }
   labels$estimate<-colMeans(a);labels$lower<-apply(boot,2,quantile,probs=.025,names=FALSE)
   labels$upper<-apply(boot,2,quantile,probs=.975,names=FALSE)
   labels$window_hours<-h;labels$dataset<-dataset;labels$records<-nrow(a);labels$clusters<-ng
   labels$valid_reps<-colSums(is.finite(boot));si<-si+1L;summaries[[si]]<-labels
   ci<-ci+1L;checks[[ci]]<-data.frame(window_hours=h,dataset=dataset,records=nrow(a),events=sum(z$y),encoded_features=ncol(phi),
      conceptual_variables=length(unique(map$variable)),clusters=ng,margin_additivity_error=e1,frozen_probability_error=e2,bias=mean(bias),
      total_mean_absolute=sum(colMeans(abs(phi))),missingness_share=sum(colMeans(abs(phi)) [map$kind=="Missingness indicator"])/sum(colMeans(abs(phi))))
   saveRDS(list(phi=phi,bias=bias,stay_id=z$pred$stay_id,bootstrap=boot,summary=labels,map=map),
      file.path(root,"private_cache",sprintf("shap_%sh_%s.rds",h,dataset)))
   log_progress("SHAP complete: ",h,"h ",dataset,"; records=",nrow(phi),"; additivity error=",signif(e1,3))
   base<-c(cv_scores(z$y,p),brier_score=mean((z$y-p)^2))
   scopes<-if(dataset=="MIMIC_internal_test")"Across records" else c("Across records","Within hospital")
   for(scope in scopes){
    set.seed(20261105L+h+if(scope=="Within hospital")1000L else 0L)
    for(b in seq_len(50)){
     idx<-if(scope=="Across records")sample.int(nrow(z$x)) else {
      ii<-seq_len(nrow(z$x));for(g in split(ii,z$pred$hospitalid))ii[g]<-g[sample.int(length(g))];ii
     }
     for(d in unique(map$domain)){
      xp<-z$x;cols<-which(map$domain==d);xp[,cols]<-z$x[idx,cols,drop=FALSE]
      pp<-as.numeric(predict(model,xgboost::xgb.DMatrix(xp,nthread=2L)))
      met<-c(cv_scores(z$y,pp),brier_score=mean((z$y-pp)^2))
      delta<-(base-met)*c(1,1,-1,-1)
      pi<-pi+1L;permutations[[pi]]<-data.frame(window_hours=h,dataset=dataset,scope=scope,repetition=b,domain=d,
         moved_fraction=mean(idx!=seq_along(idx)),metric=names(delta),performance_deterioration=as.numeric(delta))
     }
     if(b%%25==0)log_progress("Permutation: ",h,"h ",dataset," / ",scope," ",b,"/50")
    }
   }
  }
 }
 write.csv(do.call(rbind,summaries),file.path(out,"shap_summary_intervals.csv"),row.names=FALSE)
 write.csv(unique(do.call(rbind,maps)),file.path(out,"encoded_domain_dictionary.csv"),row.names=FALSE)
 write.csv(do.call(rbind,checks),file.path(out,"reconstruction_checks.csv"),row.names=FALSE)
 write.csv(do.call(rbind,permutations),file.path(out,"domain_permutation_repetitions.csv"),row.names=FALSE)
 write_json(list(status="computed_pending_independent_checks",R=R.version.string,xgboost=as.character(packageVersion("xgboost")),retraining=FALSE),file.path(out,"run_complete.json"))
}
if(sys.nframe()==0L)run_audit()
