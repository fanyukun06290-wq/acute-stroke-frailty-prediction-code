if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
root<-file.path(getwd(),"revision_2026-11-04")
source(file.path(root,"R/run_shap_domain_audit.R"),encoding="UTF-8")
summ<-read.csv(file.path(out,"shap_summary_intervals.csv"))
perm<-read.csv(file.path(out,"domain_permutation_repetitions.csv"))
nc<-pc<-0L
for(h in c(24L,6L))for(dataset in c("MIMIC_internal_test","eICU_external_validation")){
 z<-load_validation(h,dataset);s<-readRDS(file.path(root,"private_cache",sprintf("shap_%sh_%s.rds",h,dataset)))
 stopifnot(identical(s$stay_id,z$pred$stay_id),identical(colnames(s$phi),colnames(z$x)))
 model<-xgboost::xgb.load(file.path(root,"private_models",sprintf("xgboost_%sh.ubj",h)))
 margin<-as.numeric(predict(model,xgboost::xgb.DMatrix(z$x,nthread=2),outputmargin=TRUE))
 stopifnot(max(abs(rowSums(s$phi)+s$bias-margin))<2e-4)
 q<-summ[summ$window_hours==h & summ$dataset==dataset,]
 # Independent direct cluster-resampling of encoded features, compared against optimized bootstrap.
 cl<-if(dataset=="MIMIC_internal_test")z$pred$group_id else hospital_components(z$pred$group_id,z$pred$hospitalid)$row_component
 idx<-match(cl,unique(cl));ng<-max(idx)
 seed<-20261104L+(if(h==24L)0L else 4L)*1000L+(if(dataset=="MIMIC_internal_test")1L else 4L)*1000L
 set.seed(seed)
 for(b in 1:5){
  draw<-sample.int(ng,ng,replace=TRUE);weights<-tabulate(draw,nbins=ng)[idx]
  direct<-colSums(abs(s$phi)*weights)/sum(weights)
  stopifnot(max(abs(direct-s$bootstrap[b,seq_len(ncol(s$phi))]))<1e-12)
 }
 for(j in seq_len(ncol(s$phi))){
  row<-q[q$level=="Encoded feature" & q$name==colnames(s$phi)[j],]
  stopifnot(abs(mean(abs(s$phi[,j]))-row$estimate)<1e-12,
   max(abs(as.numeric(quantile(s$bootstrap[,j],c(.025,.975)))-c(row$lower,row$upper)))<1e-12)
  nc<-nc+1L
 }
 base<-c(cv_scores(z$y,z$pred$probability),brier_score=mean((z$y-z$pred$probability)^2))
 for(scope in if(dataset=="MIMIC_internal_test")"Across records" else c("Across records","Within hospital")){
  set.seed(20261105L+h+if(scope=="Within hospital")1000L else 0L)
  ii<-if(scope=="Across records")sample.int(nrow(z$x)) else {
   jj<-seq_len(nrow(z$x));for(g in split(jj,z$pred$hospitalid))jj[g]<-g[sample.int(length(g))];jj
  }
  if(scope=="Within hospital")stopifnot(all(z$pred$hospitalid[ii]==z$pred$hospitalid))
  for(d in unique(s$map$domain)){
   xx<-z$x;cols<-which(s$map$domain==d);xx[,cols]<-z$x[ii,cols,drop=FALSE]
   pp<-as.numeric(predict(model,xgboost::xgb.DMatrix(xx,nthread=2)))
   met<-c(cv_scores(z$y,pp),brier_score=mean((z$y-pp)^2))
   delta<-(base-met)*c(1,1,-1,-1)
   q<-perm[perm$window_hours==h & perm$dataset==dataset & perm$scope==scope & perm$domain==d & perm$repetition==1,]
   stopifnot(max(abs(delta-q$performance_deterioration[match(names(delta),q$metric)]))<1e-10)
   pc<-pc+1L
  }
 }
}
write_json(list(encoded_mean_and_interval_checks=nc,permutation_replay_checks=pc,bootstrap_weight_replay_sets=4,
  native_reconstruction_sets=4,status="passed"),file.path(out,"independent_r_checks.json"))
log_progress("Independent checks passed; encoded summaries=",nc,"; domain shuffle replays=",pc)
