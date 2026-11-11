Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
root<-file.path(getwd(),'revision_2026-11-04')
source(file.path(root,'R/run_revision_models.R'),encoding='UTF-8')
contract<-jsonlite::read_json(file.path(root,'model_feature_contract.json'),simplifyVector=TRUE)
rows<-list();k<-0L;input_checks<-list()
for(h in c(24L,6L)){
 a<-read_data(file.path(root,'private_data/mimic_analysis_revision.csv'),contract,h)
 b<-read_data(file.path(root,'private_data/mimic_analysis_cvst_v2.csv'),contract,h)
 stopifnot(identical(a$x,b$x),identical(a$y,b$y),identical(a$metadata$stay_id,b$metadata$stay_id),identical(a$metadata$group_id,b$metadata$group_id))
 input_checks[[as.character(h)]]<-data.frame(window_hours=h,analysis_records=length(a$y),predictor_columns=ncol(a$x),predictor_cells_changed=0L,outcomes_changed=0L)
 for(ds in c('MIMIC_internal_test','eICU_external_validation')){
  p<-readRDS(file.path(root,'private_cache',sprintf('subgroup_cvst_v2_%sh_%s.rds',h,ds)))
  stopifnot(p$reused_groups+p$recomputed_groups==11L)
  if(ds=='eICU_external_validation')stopifnot(p$reused_groups==11L,p$recomputed_groups==0L)
  k<-k+1L;rows[[k]]<-data.frame(window_hours=h,dataset=ds,whole_cohort_records=p$records,clusters=p$clusters,seed=p$seed,reused_groups=p$reused_groups,recomputed_groups=p$recomputed_groups)
 }
}
write.csv(do.call(rbind,rows),file.path(root,'phenotype_v2_results/bootstrap_reuse_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,input_checks),file.path(root,'phenotype_v2_results/R_predictor_identity.csv'),row.names=FALSE)
cat('R_NUMERIC_INPUT_AND_OUTCOME_IDENTITY_CONFIRMED; REUSE_REQUIRES_EXACT_RECORD_TRUTH_PROBABILITY_IDENTITY\n')
