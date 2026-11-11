Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root <- file.path(getwd(),"revision_2026-11-04"); out <- file.path(root,"dca_threshold_results")
plan <- jsonlite::read_json(file.path(out,"analysis_plan.json"))
sha <- function(p) digest::digest(file=p,algo="sha256")
stopifnot(identical(sha(file.path(root,"R/run_dca_thresholds.R")),plan$plan$script_sha256),
    identical(sha(file.path(root,"18_第十步DCA与固定阈值规则.md")),plan$plan$rules_sha256))
for(p in names(plan$plan$inputs_sha256)) stopifnot(identical(sha(file.path(root,p)),plan$plan$inputs_sha256[[p]]))
fixed <- read.csv(file.path(out,"fixed_threshold_intervals.csv"))
curves <- read.csv(file.path(out,"decision_curve_intervals.csv"))
pairs <- read.csv(file.path(out,"paired_model_net_benefit.csv"))
defaults <- read.csv(file.path(out,"default_strategy_intervals.csv"))
pop <- read.csv(file.path(out,"evaluation_populations.csv"))
stopifnot(nrow(fixed)==2160L,nrow(curves)==9600L,nrow(pairs)==2200L,nrow(defaults)==200L,nrow(pop)==4L)
ci <- function(v) as.numeric(quantile(v[is.finite(v)],c(.025,.975),type=7))
counts <- c(fixed=0L,curves=0L,pairs=0L,defaults=0L)
for(i in seq_len(nrow(pop))) {
    r <- pop[i,]; label <- paste(r$window_hours,r$dataset,sep="_")
    boot <- readRDS(file.path(root,"private_cache",paste0("decision_validation_",label,".rds")))
    stopifnot(identical(boot$run_hash,plan$run_hash),boot$completed==1000L,boot$replay_checks==10L)
    select <- function(x) x[x$window_hours==r$window_hours & x$dataset==r$dataset,]
    f <- select(fixed); c <- select(curves); p <- select(pairs); d <- select(defaults)
    for(k in seq_len(nrow(f))) {
        j <- match(f$threshold[k],c(.05,.10,.20,.30,.50))
        v <- boot$fixed[,f$model[k],j,f$metric[k]]
        if(f$valid_reps[k]>0) stopifnot(max(abs(ci(v)-c(f$lower[k],f$upper[k])))<1e-10)
        stopifnot(sum(is.finite(v))==f$valid_reps[k]);counts["fixed"]<-counts["fixed"]+1L
    }
    for(k in seq_len(nrow(c))) {
        j <- round(c$threshold[k]*100)
        v <- switch(c$metric[k],net_benefit=boot$net_benefit[,c$model[k],j],difference_vs_none=boot$net_benefit[,c$model[k],j],
            difference_vs_all=boot$net_benefit[,c$model[k],j]-boot$all_default[,j],
            difference_vs_best_default=boot$net_benefit[,c$model[k],j]-pmax(0,boot$all_default[,j]))
        stopifnot(max(abs(ci(v)-c(c$lower[k],c$upper[k])))<1e-10,sum(is.finite(v))==c$valid_reps[k]);counts["curves"]<-counts["curves"]+1L
    }
    for(k in seq_len(nrow(p))) {
        j <- round(p$threshold[k]*100);v <- boot$net_benefit[,"xgboost",j]-boot$net_benefit[,p$comparator[k],j]
        stopifnot(max(abs(ci(v)-c(p$lower[k],p$upper[k])))<1e-10);counts["pairs"]<-counts["pairs"]+1L
    }
    for(k in seq_len(nrow(d))) {
        j<-round(d$threshold[k]*100)
        stopifnot(max(abs(ci(boot$all_default[,j])-c(d$all_lower[k],d$all_upper[k])))<1e-10,d$none_flagged[k]==0)
        counts["defaults"]<-counts["defaults"]+1L
    }
}
jsonlite::write_json(list(status="passed_all_cached_dca_and_threshold_interval_checks",intervals_verified=as.list(counts),
    configurations=4L,replicates_each=1000L,pointwise_not_simultaneous=TRUE,models_refitted=FALSE),
    file.path(root,"qc/dca_threshold_interval_checks.json"),auto_unbox=TRUE,pretty=TRUE)
cat("DCA_THRESHOLD_INTERVAL_CHECKS_PASSED\n")
