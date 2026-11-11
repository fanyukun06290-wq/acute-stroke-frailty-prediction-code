Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root<-file.path(getwd(),'revision_2026-11-04');out<-file.path(root,'consistency_revision_results')
old<-read.csv(file.path(root,'complete_case_results/complete_case_populations.csv'),check.names=FALSE)
result<-old;contract<-jsonlite::read_json(file.path(root,'model_feature_contract.json'),simplifyVector=TRUE)
changes<-list();k<-0L;checked<-0L
for(h in c(24L,6L))for(db in c('mimic','eicu')){
 path<-file.path(root,'private_data',if(db=='mimic')'mimic_analysis_cvst_v2.csv' else 'eicu_analysis_revision.csv')
 d<-read.csv(path,check.names=FALSE)
 spec<-contract$windows[[as.character(h)]]
 keep<-!is.na(d[[spec$eligibility]])&d[[spec$eligibility]]==1&!is.na(d$outcome_missing_flag)&d$outcome_missing_flag==0
 d<-d[keep,];stopifnot(!anyDuplicated(d$stay_id))
 raw<-setdiff(spec$predictors,grep('_missing_',spec$predictors,value=TRUE));stopifnot(length(raw)==20L)
 d$sex<-trimws(as.character(d$sex));d$sex[!d$sex %in% c('F','M')]<-NA
 cc<-rowSums(is.na(d[,raw,drop=FALSE]))==0L
 if(db=='mimic'){
  reg<-read.csv(file.path(root,'private_splits',paste0('patients_',h,'h.csv')))
  part<-reg$partition[match(d$patient_group_id,reg$group_id)];stopifnot(!anyNA(part))
  masks<-list(MIMIC_all_analysis=rep(TRUE,nrow(d)),MIMIC_training=part=='train',MIMIC_internal_test=part=='test')
 }else masks<-list(eICU_external_validation=rep(TRUE,nrow(d)))
 ph<-d$audit_phenotype
 stopifnot(!anyNA(ph),all(ph %in% c('AIS_STRICT','ICH_STRICT','SAH')|grepl('^MIXED_',ph)))
 for(pop in names(masks))for(subset in c('All records','Complete cases','Incomplete cases')){
  mask<-masks[[pop]] & switch(subset,'All records'=rep(TRUE,nrow(d)),'Complete cases'=cc,'Incomplete cases'=!cc)
  ix<-which(result$window_hours==h&result$population==pop&result$subset==subset);stopifnot(length(ix)==1L)
  stopifnot(result$records[ix]==sum(mask),result$events[ix]==sum(d[[spec$outcome]][mask]),result$unique_patients[ix]==length(unique(d$patient_group_id[mask])))
  values<-c(ais_records=sum(ph[mask]=='AIS_STRICT'),ich_records=sum(ph[mask]=='ICH_STRICT'),sah_records=sum(ph[mask]=='SAH'),mixed_records=sum(grepl('^MIXED_',ph[mask])))
  stopifnot(sum(values)==sum(mask))
  for(col in names(values)){
   k<-k+1L;changes[[k]]<-data.frame(window_hours=h,population=pop,subset=subset,field=col,old_count=old[[col]][ix],v2_count=values[[col]],difference=values[[col]]-old[[col]][ix])
   result[[col]][ix]<-values[[col]]
  }
  checked<-checked+1L
 }
}
nonph<-setdiff(names(old),c('ais_records','ich_records','sah_records','mixed_records'))
stopifnot(identical(result[,nonph],old[,nonph]))
write.csv(result,file.path(out,'corrected_complete_case_populations_cvst_v2.csv'),row.names=FALSE)
write.csv(do.call(rbind,changes),file.path(out,'phenotype_count_changes.csv'),row.names=FALSE)
jsonlite::write_json(list(status='R_counts_recomputed',population_rows=checked,phenotype_counts=96L,nonphenotype_fields_unchanged=TRUE,model_refitting=FALSE),file.path(out,'R_count_checks.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Recomputed 24 population rows and 96 phenotype counts; nonphenotype fields unchanged\n')
