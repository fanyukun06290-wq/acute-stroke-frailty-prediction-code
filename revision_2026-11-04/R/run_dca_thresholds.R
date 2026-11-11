# Frozen-prediction binary DCA and inherited fixed-threshold operating points.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "dca_threshold_results")
dir.create(out, showWarnings = FALSE)
thresholds <- (1:50)/100
fixed <- c(5,10,20,30,50)/100
model_names <- c("xgboost", "random_forest", "elastic_net", "clinical_ridge", "motor_B", "motor_BMI", "total_B", "total_BMI",
    "omit_neuro_flag_B", "omit_neuro_flag_BMI", "omit_all_flags_B", "omit_all_flags_BMI")
fixed_metrics <- c("sensitivity","specificity","ppv","npv","f1","flagged_fraction","residual_event_rate_negative",
    "false_positives_per1000","false_negatives_per1000")
dca_metrics <- c("net_benefit","difference_vs_all","difference_vs_none","difference_vs_best_default")
prepare_decision <- function(y,p) list(y=y,p=p,bin=findInterval(p,thresholds))
decision_counts <- function(o,w) {
    a <- rowsum(cbind(w*o$y,w*(1-o$y)),o$bin+1L,reorder=FALSE)
    b <- matrix(0,51L,2L); b[as.integer(rownames(a)),] <- a
    events <- sum(b[,1L]); nonevents <- sum(b[,2L]); n <- events+nonevents
    tp <- rev(cumsum(rev(b[,1L])))[-1L]; fp <- rev(cumsum(rev(b[,2L])))[-1L]
    list(n=n,events=events,tp=tp,fp=fp,tn=nonevents-fp,fn=events-tp)
}
safe_divide <- function(a,b) ifelse(b>0,a/b,NA_real_)
fixed_values <- function(c) {
    i <- match(fixed,thresholds); tp <- c$tp[i]; fp <- c$fp[i]; tn <- c$tn[i]; fn <- c$fn[i]
    z <- cbind(sensitivity=safe_divide(tp,tp+fn),specificity=safe_divide(tn,tn+fp),ppv=safe_divide(tp,tp+fp),
        npv=safe_divide(tn,tn+fn),f1=safe_divide(2*tp,2*tp+fp+fn),flagged_fraction=(tp+fp)/c$n,
        residual_event_rate_negative=safe_divide(fn,tn+fn),false_positives_per1000=1000*fp/c$n,false_negatives_per1000=1000*fn/c$n)
    z[,fixed_metrics,drop=FALSE]
}
net_benefits <- function(c) {
    odds <- thresholds/(1-thresholds)
    nb <- c$tp/c$n-c$fp/c$n*odds
    all <- c$events/c$n-(1-c$events/c$n)*odds
    cbind(net_benefit=nb,difference_vs_all=nb-all,difference_vs_none=nb,difference_vs_best_default=nb-pmax(0,all))
}
read_aligned_decisions <- function(hour,dataset) {
    readp <- function(name) read.csv(file.path(root,"private_predictions",paste0(name,"_",hour,"h.csv")))
    principal <- readp("predictions"); clinical <- readp("clinical_ridge_predictions")
    neuro <- readp("neuro_predictions"); sensitivity <- readp("missing_sensitivity_predictions")
    reference <- principal[principal$dataset==dataset & principal$model=="xgboost",]
    reference <- reference[order(reference$stay_id),]; stopifnot(!anyDuplicated(reference$stay_id))
    objects <- list()
    for (name in model_names) {
        raw <- if(name %in% c("xgboost","random_forest","elastic_net")) principal else if(name=="clinical_ridge") clinical else
            if(grepl("^omit_",name)) sensitivity else neuro
        d <- raw[raw$dataset==dataset & raw$model==name,]; stopifnot(nrow(d)==nrow(reference),!anyDuplicated(d$stay_id))
        d <- d[match(reference$stay_id,d$stay_id),]
        stopifnot(!anyNA(d$probability),identical(d$truth,reference$truth),identical(d$group_id,reference$group_id),
            identical(as.character(d$hospitalid),as.character(reference$hospitalid)))
        objects[[name]] <- prepare_decision(as.integer(d$truth=="yes"),d$probability)
    }
    list(reference=reference,objects=objects)
}
validate_counts_engine <- function() {
    y <- c(0,1,0,1,1,0,0,1); p <- c(.001,.05,.10,.10,.20,.30,.50,.99); w <- c(2,1,3,2,1,2,1,1)
    o <- prepare_decision(y,p); count <- decision_counts(o,w)
    for(j in seq_along(thresholds)) {
        flagged <- p>=thresholds[j]
        stopifnot(count$tp[j]==sum(w*y*flagged),count$fp[j]==sum(w*(1-y)*flagged),
            count$fn[j]==sum(w*y*!flagged),count$tn[j]==sum(w*(1-y)*!flagged))
    }
    all_negative <- fixed_values(decision_counts(prepare_decision(y,rep(0,length(y))),w))
    all_positive <- fixed_values(decision_counts(prepare_decision(y,rep(1,length(y))),w))
    stopifnot(all(is.na(all_negative[,"ppv"])),all(is.na(all_positive[,"npv"])))
    TRUE
}
run_decision_analysis <- function() {
    inputs <- c("qc/runtime_model_files.json","qc/runtime_clinical_files.json",
        "qc/runtime_neuro_files.json","qc/runtime_missing_files.json")
    for(hour in c(24L,6L)) for(prefix in c("predictions","clinical_ridge_predictions","neuro_predictions","missing_sensitivity_predictions"))
        inputs <- c(inputs,paste0("private_predictions/",prefix,"_",hour,"h.csv"))
    plan <- list(version="frozen_binary_dca_1.0",models=model_names,threshold_grid=thresholds,fixed_thresholds=fixed,primary_threshold=.10,
        original_threshold_script_sha256=sha_file(file.path(root,"../critical_care_code_release_2026-09-08/R/build_subgroup_threshold.R")),
        fixed_thresholds_inherited_not_newly_optimized=TRUE,dca_grid_added_during_revision=TRUE,
        rules_sha256=sha_file(file.path(root,"18_第十步DCA与固定阈值规则.md")),
        script_sha256=sha_file(file.path(root,"R/run_dca_thresholds.R")),
        inputs_sha256=setNames(lapply(inputs,function(p) sha_file(file.path(root,p))),inputs),
        endpoint="Recorded in-hospital death after landmark and by ICU admission day 30; not full 30-day all-cause mortality",
        intended_action="Hypothetical additional risk review, with no validated intervention or preference-elicited threshold",
        no_training_or_recalibration=TRUE,replicates=1000L,
        intervals="Conditional paired cluster percentile 95% intervals; pointwise grid, not simultaneous; exploratory multiple comparisons")
    hash <- digest::digest(plan,algo="sha256"); path <- file.path(out,"analysis_plan.json")
    if(file.exists(path)) stopifnot(identical(jsonlite::read_json(path)$run_hash,hash)) else write_json(list(run_hash=hash,plan=plan),path)
    stopifnot(validate_counts_engine())
    points <- metric_rows <- dca_rows <- pair_rows <- defaults <- populations <- list()
    replay <- 0L
    for(hour in c(24L,6L)) for(dataset in c("MIMIC_internal_test","eICU_external_validation")) {
        aligned <- read_aligned_decisions(hour,dataset); reference <- aligned$reference; objects <- aligned$objects
        cluster <- reference$group_id; scheme <- "patient"
        if(dataset=="eICU_external_validation") {cluster <- hospital_components(reference$group_id,reference$hospitalid)$row_component;scheme <- "hospital_patient_component"}
        group_index <- match(cluster,unique(cluster)); groups <- length(unique(cluster)); label <- paste(hour,dataset,sep="_")
        populations[[label]] <- data.frame(window_hours=hour,dataset=dataset,records=nrow(reference),patients=length(unique(reference$group_id)),
            events=sum(reference$truth=="yes"),scheme=scheme,resampling_blocks=groups)
        counts <- lapply(objects,decision_counts,w=rep(1,nrow(reference)))
        nb_points <- lapply(counts,net_benefits); fixed_points <- lapply(counts,fixed_values)
        old <- readRDS(file.path(root,"private_cache",paste0("uncertainty_",hour,"_",dataset,"_",scheme,".rds")))
        cache <- file.path(root,"private_cache",paste0("decision_validation_",label,".rds"))
        if(file.exists(cache)) {
            boot <- readRDS(cache);stopifnot(identical(boot$run_hash,hash),boot$completed==1000L)
        } else {
            boot <- list(run_hash=hash,completed=0L,replay_checks=0L,
                net_benefit=array(NA_real_,c(1000L,length(model_names),length(thresholds)),dimnames=list(NULL,model_names,NULL)),
                all_default=matrix(NA_real_,1000L,length(thresholds)),
                fixed=array(NA_real_,c(1000L,length(model_names),length(fixed),length(fixed_metrics)),dimnames=list(NULL,model_names,NULL,fixed_metrics)))
            set.seed(20261104L+(if(hour==24L)0L else 4L)*1000L+(if(dataset=="MIMIC_internal_test")1L else 4L)*1000L)
            for(b in 1:1000L) {
                w <- tabulate(sample.int(groups,groups,replace=TRUE),nbins=groups)[group_index]
                n <- sum(w); events <- sum(w*objects$xgboost$y)
                stopifnot(n==old$records[b],events==old$events[b])
                boot$all_default[b,] <- events/n-(1-events/n)*thresholds/(1-thresholds)
                if(b%%100L==0L) {
                    brier <- sum(w*(objects$xgboost$y-objects$xgboost$p)^2)/n
                    stopifnot(abs(brier-old$metrics[b,"xgboost","brier_score"])<1e-10);boot$replay_checks<-boot$replay_checks+1L
                }
                for(name in model_names) {
                    c <- decision_counts(objects[[name]],w)
                    boot$net_benefit[b,name,] <- net_benefits(c)[,"net_benefit"]
                    boot$fixed[b,name,,] <- fixed_values(c)
                }
                if(b%%250L==0L) log_progress("DCA/threshold ",label," ",b,"/1000")
            }
            boot$completed <- 1000L;saveRDS(boot,cache)
        }
        replay <- replay+boot$replay_checks
        for(j in seq_along(thresholds)) {
            all <- counts$xgboost$events/counts$xgboost$n-(1-counts$xgboost$events/counts$xgboost$n)*thresholds[j]/(1-thresholds[j])
            ci <- quantile_ci(boot$all_default[,j])
            defaults[[length(defaults)+1L]] <- data.frame(window_hours=hour,dataset=dataset,threshold=thresholds[j],all_flagged=all,
                all_lower=ci[["lower"]],all_upper=ci[["upper"]],none_flagged=0)
        }
        for(name in model_names) {
            for(j in seq_along(fixed)) {
                k <- match(fixed[j],thresholds);c <- counts[[name]]
                points[[length(points)+1L]] <- data.frame(window_hours=hour,dataset=dataset,model=name,threshold=fixed[j],
                    primary_threshold=fixed[j]==.10,records=c$n,events=c$events,tp=c$tp[k],fp=c$fp[k],tn=c$tn[k],fn=c$fn[k],
                    as.list(fixed_points[[name]][j,]))
                for(metric in fixed_metrics) {
                    ci <- quantile_ci(boot$fixed[,name,j,metric])
                    metric_rows[[length(metric_rows)+1L]] <- data.frame(window_hours=hour,dataset=dataset,scheme=scheme,model=name,
                        threshold=fixed[j],metric=metric,estimate=fixed_points[[name]][j,metric],lower=ci[["lower"]],upper=ci[["upper"]],
                        valid_reps=ci[["valid_reps"]],status=if(ci[["valid_reps"]]==0)"undefined_no_valid_replicates" else
                            if(ci[["valid_reps"]]<950)"limited_valid_replicates" else "computed")
                }
            }
            for(j in seq_along(thresholds)) for(metric in dca_metrics) {
                v <- switch(metric,net_benefit=boot$net_benefit[,name,j],difference_vs_none=boot$net_benefit[,name,j],
                    difference_vs_all=boot$net_benefit[,name,j]-boot$all_default[,j],
                    difference_vs_best_default=boot$net_benefit[,name,j]-pmax(0,boot$all_default[,j]))
                ci <- quantile_ci(v);stopifnot(ci[["valid_reps"]]>=950L)
                dca_rows[[length(dca_rows)+1L]] <- data.frame(window_hours=hour,dataset=dataset,scheme=scheme,model=name,
                    threshold=thresholds[j],metric=metric,estimate=nb_points[[name]][j,metric],lower=ci[["lower"]],upper=ci[["upper"]],valid_reps=ci[["valid_reps"]])
            }
            if(name!="xgboost") for(j in seq_along(thresholds)) {
                ci <- quantile_ci(boot$net_benefit[,"xgboost",j]-boot$net_benefit[,name,j])
                pair_rows[[length(pair_rows)+1L]] <- data.frame(window_hours=hour,dataset=dataset,comparator=name,threshold=thresholds[j],
                    difference_xgboost_minus_comparator=nb_points$xgboost[j,"net_benefit"]-nb_points[[name]][j,"net_benefit"],
                    lower=ci[["lower"]],upper=ci[["upper"]],valid_reps=ci[["valid_reps"]])
            }
        }
    }
    write.csv(do.call(rbind,points),file.path(out,"fixed_threshold_points.csv"),row.names=FALSE)
    write.csv(do.call(rbind,metric_rows),file.path(out,"fixed_threshold_intervals.csv"),row.names=FALSE)
    write.csv(do.call(rbind,dca_rows),file.path(out,"decision_curve_intervals.csv"),row.names=FALSE)
    write.csv(do.call(rbind,pair_rows),file.path(out,"paired_model_net_benefit.csv"),row.names=FALSE)
    write.csv(do.call(rbind,defaults),file.path(out,"default_strategy_intervals.csv"),row.names=FALSE)
    write.csv(do.call(rbind,populations),file.path(out,"evaluation_populations.csv"),row.names=FALSE)
    write_json(list(status="exploratory_frozen_prediction_dca_complete_not_clinical_impact_validation",run_hash=hash,
        counts_engine_edge_cases_passed=TRUE,original_bootstrap_replay_checks=replay,models_refitted=FALSE,
        intervals="pointwise conditional percentile; no simultaneous inference",new_threshold_optimization=FALSE),file.path(out,"complete.json"))
}
if(sys.nframe()==0L) run_decision_analysis()
