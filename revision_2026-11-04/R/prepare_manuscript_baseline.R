Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root<-file.path(getwd(),'revision_2026-11-04');out<-file.path(root,'manuscript_revision_results');dir.create(out,showWarnings=FALSE)
contract<-jsonlite::read_json(file.path(root,'model_feature_contract.json'),simplifyVector=TRUE)
map<-read.csv(file.path(root,'phenotype_v2_results/measurement_dictionary.csv'))
rows<-list();pop<-list();k<-j<-0L
for(h in c(24L,6L))for(db in c('mimic','eicu')){
 path<-file.path(root,'private_data',if(db=='mimic')'mimic_analysis_cvst_v2.csv' else 'eicu_analysis_revision.csv')
 d<-read.csv(path);d<-d[d[[paste0('landmark_',h,'h_eligible')]]==1&d$outcome_missing_flag==0,]
 parts<-list(All_analysis_records=seq_len(nrow(d)))
 if(db=='mimic'){
  s<-read.csv(file.path(root,'private_splits',paste0('patients_',h,'h.csv')))
  p<-s$partition[match(d$patient_group_id,s$group_id)];stopifnot(!anyNA(p))
  parts$Training<-which(p=='train');parts$Internal_test<-which(p=='test')
 }
 for(part in names(parts)){
  z<-d[parts[[part]],];n<-nrow(z);event<-sum(z[[paste0('death_',h,'h_to_30d_inhospital')]])
  j<-j+1L;pop[[j]]<-data.frame(window_hours=h,database=db,partition=part,records=n,patients=length(unique(z$patient_group_id)),events=event)
  add<-function(variable,label,unit,numerator=NA_real_,denominator=n,median=NA_real_,q1=NA_real_,q3=NA_real_,missing=0L,display){
   k<<-k+1L;rows[[k]]<<-data.frame(window_hours=h,database=db,partition=part,variable=variable,label=label,unit=unit,records=n,numerator=numerator,denominator=denominator,median=median,q1=q1,q3=q3,missing=missing,missing_pct=100*missing/n,display=display)
  }
  add('records','ICU records','count',display=format(n,big.mark=','))
  add('patients','Unique patients','count',display=format(length(unique(z$patient_group_id)),big.mark=','))
  add('events','Outcome events','n (%)',numerator=event,display=sprintf('%s (%.1f)',format(event,big.mark=','),100*event/n))
  variables<-c('age_years',map$variable)
  for(v in variables){
   col<-if(v=='age_years')v else paste0(v,if(h==24L)'_mean_24h' else '_first_6h')
   x<-z[[col]];stopifnot(all(is.na(x)|is.finite(x)))
   q<-quantile(x,c(.25,.5,.75),na.rm=TRUE,names=FALSE,type=7);miss<-sum(is.na(x))
   label<-if(v=='age_years')'Age' else map$description[match(v,map$variable)]
   unit<-if(v=='age_years')'years' else map$reporting_unit[match(v,map$variable)]
   add(v,label,unit,median=q[2],q1=q[1],q3=q[3],missing=miss,display=sprintf('%.1f [%.1f, %.1f]; %.1f%%',q[2],q[1],q[3],100*miss/n))
  }
  female<-sum(z$sex=='F',na.rm=TRUE);unknown<-sum(is.na(z$sex)|!z$sex %in% c('F','M'))
  add('sex_female','Female','n (%)',numerator=female,missing=unknown,display=sprintf('%s (%.1f)',format(female,big.mark=','),100*female/n))
  ph<-ifelse(z$audit_phenotype %in% c('AIS_STRICT','ICH_STRICT','SAH'),z$audit_phenotype,'Mixed/ambiguous')
  for(v in c('AIS_STRICT','ICH_STRICT','SAH','Mixed/ambiguous')){
   nn<-sum(ph==v);add(v,v,'n (%)',numerator=nn,display=sprintf('%s (%.1f)',format(nn,big.mark=','),100*nn/n))
  }
 }
}
data<-do.call(rbind,rows);write.csv(data,file.path(out,'baseline_long.csv'),row.names=FALSE)
write.csv(do.call(rbind,pop),file.path(out,'baseline_population_checks.csv'),row.names=FALSE)
for(h in c(24L,6L)){
 r<-data[data$window_hours==h,];r<-r[(r$database=='mimic'&r$partition %in% c('Training','Internal_test'))|(r$database=='eicu'&r$partition=='All_analysis_records'),]
 r$cohort<-ifelse(r$database=='eicu','eICU external',ifelse(r$partition=='Training','MIMIC training','MIMIC internal test'))
 tab<-unique(r[,c('variable','label','unit')])
 for(cohort in c('MIMIC training','MIMIC internal test','eICU external')){
  x<-r[r$cohort==cohort,];tab[[cohort]]<-x$display[match(tab$variable,x$variable)]
 }
 stopifnot(!anyNA(tab));write.csv(tab,file.path(out,paste0('Table1_',h,'h.csv')),row.names=FALSE)
}
cat('BASELINE_TABLES_CREATED_FROM_UNIMPUTED_V2_INPUTS; no data or models modified\n')
