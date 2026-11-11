# Independent aggregate check; no row identifiers are written or printed.
Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root<-file.path(getwd(),'revision_2026-11-04')
out<-file.path(root,'codebook_results')
bool<-function(x)tolower(as.character(x)) %in% c('1','t','true')
classify<-function(d){
 a<-bool(d$has_ais);i<-bool(d$has_ich);s<-bool(d$has_sah);v<-bool(d$has_cvst)
 u<-bool(d$has_unspecified_stroke);o<-bool(d$has_other_intracranial_hemorrhage);t<-bool(d$has_traumatic_intracranial_injury)
 n<-as.integer(a)+as.integer(i)+as.integer(s)
 ifelse(t&(n>0|v|u|o),'EXCLUDE_TRAUMATIC',
 ifelse(n>1,'MIXED_PRIMARY_ACUTE',ifelse(n==1&v,'MIXED_WITH_CVST',
 ifelse(n==1&(u|o),'MIXED_WITH_AMBIGUOUS_SECONDARY',ifelse(v,'CVST',
 ifelse(a,'AIS_STRICT',ifelse(i,'ICH_STRICT',ifelse(s,'SAH',ifelse(u,'UNSPECIFIED_STROKE',
 ifelse(o,'OTHER_INTRACRANIAL_HEMORRHAGE','NO_STROKE_GROUP'))))))))))
}
tt<-read.csv(file.path(out,'synthetic_flag_truth_table.csv'));stopifnot(nrow(tt)==128L,all(classify(tt)==tt$expected_label))
summary<-list()
for(src in c('mimic','eicu')){
 d<-read.csv(file.path(root,'private_data',paste0(src,'_cohort.csv')))
 stopifnot(all(classify(d)==d$audit_phenotype))
 summary[[src]]<-data.frame(source=src,records=nrow(d),mismatches=sum(classify(d)!=d$audit_phenotype))
}
u<-read.csv(file.path(out,'observed_raw_units.csv'));dict<-read.csv(file.path(out,'measurement_dictionary.csv'))
stopifnot(nrow(dict)==18L,nrow(u)==36L,all(dict$numeric_multiplier==1L),all(dict$lower_inclusive<dict$upper_inclusive))
b<-read.csv(file.path(out,'final_input_bounds_checks.csv'))
stopifnot(nrow(b)==72L,all(b$outside_bounds==0),sum(b$available_records)>0)
e<-read.csv(file.path(out,'eicu_phenotype_edge_aggregates.csv'));m<-read.csv(file.path(out,'mimic_phenotype_edge_aggregates.csv'))
getn<-function(d,n)d$records[match(n,d$diagnostic)]
stopifnot(getn(e,'original_source_cohort')==4673,getn(m,'original_source_cohort')==7594,
 getn(e,'newly_included_text_alternative')==0,getn(e,'icd10_concordant_only_original_guards_no_trauma')==0,
 getn(m,'new_cvst_flag_if_4376_I676_added')==18,getn(m,'24h_analysis_strict_label_to_mixed')==10,
 getn(m,'6h_analysis_strict_label_to_mixed')==12)
write.csv(do.call(rbind,summary),file.path(out,'independent_R_metadata_checks.csv'),row.names=FALSE)
jsonlite::write_json(list(status='bounded_checks_passed_with_open_CVST_label_issue',synthetic_cases=128L,
 source_metadata_rows=sum(vapply(summary,function(x)x$records,numeric(1))),measurements=18L,unit_rows=36L,bound_checks=72L,
 database_queries_read_only=TRUE,original_extraction_modified=FALSE,model_modified=FALSE),file.path(out,'independent_R_checks.json'),auto_unbox=TRUE,pretty=TRUE)
cat('CODEBOOK_INDEPENDENT_R_CHECKS_PASSED_WITH_OPEN_CVST_LABEL_ISSUE\n')
