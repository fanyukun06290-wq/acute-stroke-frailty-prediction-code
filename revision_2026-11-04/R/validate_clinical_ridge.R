# Frozen supplementary model; replay the original paired cluster weights exactly.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "clinical_reference_results")
fit_plan <- jsonlite::read_json(file.path(out, "clinical_ridge_plan.json"))
old_plan <- jsonlite::read_json(file.path(root, "validation_results/uncertainty_plan.json"))
plan <- list(replicates = 1000L, original_uncertainty_hash = old_plan$run_hash,
    new_model_hash = fit_plan$run_hash, script_sha256 = sha_file(file.path(root, "R/validate_clinical_ridge.R")),
    target = "Conditional on frozen models; exploratory paired evaluation; not a validated score comparison",
    rule = "Replay original patient/hospital/connected-block seeds in original sorted-stay order")
plan$prediction_hashes <- lapply(c(24L, 6L), function(h) sha_file(file.path(root, "private_predictions", paste0("clinical_ridge_predictions_", h, "h.csv"))))
hash <- digest::digest(plan, algo = "sha256")
planfile <- file.path(out, "validation_plan.json")
if (file.exists(planfile)) stopifnot(identical(jsonlite::read_json(planfile)$run_hash, hash)) else
    write_json(list(run_hash = hash, plan = plan), planfile)
metric_rows <- paired_rows <- bin_rows <- list(); counter <- 0L
for (hour in c(24L, 6L)) for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
    aligned <- align_predictions(file.path(root, "private_predictions", paste0("predictions_", hour, "h.csv")), dataset)
    reference <- aligned$reference
    new <- read.csv(file.path(root, "private_predictions", paste0("clinical_ridge_predictions_", hour, "h.csv")))
    new <- new[new$dataset == dataset, ]
    stopifnot(!anyDuplicated(new$stay_id), nrow(new) == nrow(reference))
    new <- new[match(reference$stay_id, new$stay_id), ]
    stopifnot(!anyNA(new$probability), identical(new$truth, reference$truth), identical(new$group_id, reference$group_id),
        identical(as.character(new$hospitalid), as.character(reference$hospitalid)))
    object <- prepare_probability(as.integer(new$truth == "yes"), new$probability)
    schemes <- list(patient = reference$group_id)
    if (dataset == "eICU_external_validation") {
        schemes$hospital <- reference$hospitalid
        schemes$hospital_patient_component <- hospital_components(reference$group_id, reference$hospitalid)$row_component
    }
    for (scheme in names(schemes)) {
        counter <- counter + 1L
        label <- paste(hour, dataset, scheme, sep = "_")
        old <- readRDS(file.path(root, "private_cache", paste0("uncertainty_", label, ".rds")))
        stopifnot(identical(old$run_hash, old_plan$run_hash), old$completed == 1000L)
        cluster <- schemes[[scheme]]; group_index <- match(cluster, unique(cluster)); groups <- length(unique(cluster))
        stopifnot(groups == old$clusters)
        cache <- file.path(root, "private_cache", paste0("clinical_ridge_uncertainty_", label, ".rds"))
        if (file.exists(cache)) {
            boot <- readRDS(cache); stopifnot(identical(boot$run_hash, hash), boot$completed == 1000L)
        } else {
            boot <- list(run_hash = hash, completed = 0L, metrics = matrix(NA_real_, 1000L, length(metric_names), dimnames = list(NULL, metric_names)),
                bin_rates = matrix(NA_real_, 1000L, max(object$bins)), replay_checks = 0L)
            set.seed(20261104L + counter * 1000L)
            for (b in 1:1000) {
                weights <- tabulate(sample.int(groups, groups, replace = TRUE), nbins = groups)[group_index]
                stopifnot(sum(weights) == old$records[b], sum(weights*object$y) == old$events[b])
                if (b %% 100L == 0L) {
                    check <- weighted_metrics(aligned$objects$xgboost, weights)
                    stopifnot(max(abs(check[c("auroc", "brier_score")] - old$metrics[b, "xgboost", c("auroc", "brier_score")])) < 1e-10)
                    boot$replay_checks <- boot$replay_checks + 1L
                }
                boot$metrics[b, ] <- weighted_metrics(object, weights)
                by_bin <- rowsum(cbind(weights*object$y, weights), object$bins, reorder = TRUE)
                boot$bin_rates[b, ] <- by_bin[, 1L]/by_bin[, 2L]
                if (b %% 250L == 0L) log_progress("clinical ridge ", label, " bootstrap ", b, "/1000")
            }
            boot$completed <- 1000L
            saveRDS(boot, cache)
        }
        point <- weighted_metrics(object, rep(1, nrow(reference)))
        for (metric in metric_names) {
            ci <- quantile_ci(boot$metrics[, metric])
            stopifnot(ci[["valid_reps"]] >= 950L)
            metric_rows[[length(metric_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset,
                scheme = scheme, model = "clinical_ridge", metric = metric, estimate = point[metric],
                lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
        }
        for (metric in c("auroc", "auprc", "brier_score", "log_loss")) {
            difference <- old$metrics[, "xgboost", metric] - boot$metrics[, metric]
            ci <- quantile_ci(difference)
            paired_rows[[length(paired_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset,
                scheme = scheme, comparator = "clinical_ridge", metric = metric,
                difference_xgboost_minus_comparator = old$point_metrics$xgboost[metric] - point[metric],
                lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
        }
        for (bin in seq_len(max(object$bins))) {
            keep <- object$bins == bin; ci <- quantile_ci(boot$bin_rates[, bin])
            bin_rows[[length(bin_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset, scheme = scheme,
                model = "clinical_ridge", bin = bin, records = sum(keep), events = sum(object$y[keep]),
                mean_predicted = mean(object$p[keep]), observed_rate = mean(object$y[keep]),
                lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
        }
    }
}
write.csv(do.call(rbind, metric_rows), file.path(out, "metric_intervals.csv"), row.names = FALSE)
write.csv(do.call(rbind, paired_rows), file.path(out, "paired_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, bin_rows), file.path(out, "calibration_bins.csv"), row.names = FALSE)
write_json(list(status = "conditional_validation_complete_not_validated_score_requirement_resolved", run_hash = hash,
    schemes = counter, replicates_per_scheme = 1000L, paired_weights_replayed = TRUE, main_model_changed = FALSE),
    file.path(out, "validation_complete.json"))
