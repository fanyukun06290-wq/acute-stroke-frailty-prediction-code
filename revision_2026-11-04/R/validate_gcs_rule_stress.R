# Paired discrimination/loss uncertainty for a fixed-model input stress test.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "gcs_rule_sensitivity_results")
models <- c("omit_all_flags_B", "omit_all_flags_BMI", "total_B", "total_BMI")
variants <- c("strict_baseline", "allow_reverse_recording", "first_valid", "last_valid", "motor_requires_complete", "airway_60", "airway_240")
four_metrics <- c("auroc", "auprc", "brier_score", "log_loss")
fast_metrics <- function(o, weights) {
    n <- sum(weights); events <- sum(weights*o$y); nonevents <- n-events
    pp <- pmin(pmax(o$p, 1e-6), 1-1e-6)
    brier <- sum(weights*(o$y-o$p)^2)/n
    loss <- -sum(weights*(o$y*log(pp)+(1-o$y)*log(1-pp)))/n
    auc <- pr <- NA_real_
    if (events > 0 && nonevents > 0) {
        i <- o$order
        grouped <- rowsum(cbind(weights[i]*o$y[i], weights[i]*(1-o$y[i])), o$score_group, reorder = FALSE)
        grouped <- grouped[rowSums(grouped) > 0, , drop = FALSE]
        tp <- cumsum(grouped[, 1L]); fp <- cumsum(grouped[, 2L])
        recall <- tp/events; fpr <- fp/nonevents; precision <- tp/(tp+fp)
        auc <- sum(diff(c(0, fpr))*(c(0, head(recall, -1L))+recall)/2)
        pr <- sum(diff(c(0, recall))*(c(1, head(precision, -1L))+precision)/2)
    }
    setNames(c(auc, pr, brier, loss), four_metrics)
}
evaluation <- jsonlite::read_json(file.path(out, "evaluation_plan.json"))
missing_plan <- jsonlite::read_json(file.path(root, "missingness_sensitivity_results/validation_plan.json"))
neuro_plan <- jsonlite::read_json(file.path(root, "neuro_reference_results/validation_plan.json"))
plan <- list(version = "frozen_rule_stress_intervals_1.0", posthoc = TRUE, models_refitted = FALSE,
    evaluation_hash = evaluation$run_hash, original_missingness_hash = missing_plan$run_hash, original_neuro_hash = neuro_plan$run_hash,
    script_sha256 = sha_file(file.path(root, "R/validate_gcs_rule_stress.R")), replicates = 1000L, models = models,
    prediction_hashes = setNames(lapply(c(24L, 6L), function(h)
        sha_file(file.path(root, "private_predictions", paste0("gcs_rule_stress_", h, "h.csv")))), c("24", "6")),
    inference = "fixed-model paired input perturbation, not a refitted model or new independent validation; no multiplicity adjustment",
    uncertainty_metrics = four_metrics, calibration_change_intervals = FALSE)
hash <- digest::digest(plan, algo = "sha256")
path <- file.path(out, "validation_plan.json")
if (file.exists(path)) stopifnot(identical(jsonlite::read_json(path)$run_hash, hash)) else
    write_json(list(run_hash = hash, plan = plan), path)
paired <- domain <- populations <- list(); engine_checks <- 0L; replay_checks <- 0L
for (hour in c(24L, 6L)) {
    raw <- read.csv(file.path(root, "private_predictions", paste0("gcs_rule_stress_", hour, "h.csv")))
    original_missing <- read.csv(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv")))
    original_neuro <- read.csv(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")))
    for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
        reference <- raw[raw$dataset == dataset & raw$model == models[1L] & raw$variant == "strict_baseline", ]
        reference <- reference[order(reference$stay_id), ]; stopifnot(!anyDuplicated(reference$stay_id))
        objects <- base_objects <- list(); identities <- list()
        for (model in models) {
            d <- if (grepl("^omit_", model)) original_missing else original_neuro
            d <- d[d$dataset == dataset & d$model == model, ]
            d <- d[match(reference$stay_id, d$stay_id), ]; stopifnot(!anyNA(d$probability))
            base_objects[[model]] <- prepare_probability(as.integer(d$truth == "yes"), d$probability)
            for (variant in variants) {
                d <- raw[raw$dataset == dataset & raw$model == model & raw$variant == variant, ]
                d <- d[match(reference$stay_id, d$stay_id), ]
                stopifnot(!anyNA(d$probability), identical(d$truth, reference$truth), identical(d$group_id, reference$group_id),
                    identical(as.character(d$hospitalid), as.character(reference$hospitalid)))
                key <- paste(variant, model, sep = "/")
                identities[[key]] <- max(abs(d$probability-base_objects[[model]]$p)) < 1e-14
                objects[[key]] <- if (identities[[key]]) base_objects[[model]] else
                    prepare_probability(as.integer(d$truth == "yes"), d$probability)
            }
        }
        if (engine_checks == 0L) for (seed in c(20261105L, 20261106L)) {
            set.seed(seed); n <- nrow(reference); w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
            duplicated <- rep(seq_len(n), w)
            for (o in objects) {
                fast <- fast_metrics(o, w); expected <- cv_scores(o$y[duplicated], o$p[duplicated])
                stopifnot(max(abs(fast[c("auroc", "auprc", "log_loss")]-expected)) < 1e-9,
                    abs(fast["brier_score"]-mean((o$y[duplicated]-o$p[duplicated])^2)) < 1e-12)
            }
            engine_checks <- engine_checks+length(objects)
        }
        cluster <- reference$group_id; scheme <- "patient"
        if (dataset == "eICU_external_validation") {
            cluster <- hospital_components(reference$group_id, reference$hospitalid)$row_component
            scheme <- "hospital_patient_component"
        }
        group_index <- match(cluster, unique(cluster)); groups <- length(unique(cluster))
        label <- paste(hour, dataset, sep = "_")
        missing_cache <- readRDS(file.path(root, "private_cache", paste0("missing_validation_", hour, "_", dataset, "_all_eligible.rds")))
        neuro_cache <- readRDS(file.path(root, "private_cache", paste0("neuro_validation_", hour, "_", dataset, "_all_eligible_motor_total.rds")))
        stopifnot(identical(missing_cache$run_hash, missing_plan$run_hash), identical(neuro_cache$run_hash, neuro_plan$run_hash))
        baseline_array <- array(NA_real_, c(1000L, length(models), 4L), dimnames = list(NULL, models, four_metrics))
        for (model in models) baseline_array[, model, ] <-
            if (grepl("^omit_", model)) missing_cache$metrics[, model, four_metrics] else neuro_cache$metrics[, model, four_metrics]
        cache <- file.path(root, "private_cache", paste0("gcs_rules_validation_", label, ".rds"))
        if (file.exists(cache)) {
            boot <- readRDS(cache); stopifnot(identical(boot$run_hash, hash), boot$completed == 1000L)
        } else {
            boot <- list(run_hash = hash, completed = 0L, baseline_vectors_checked = 0L,
                metrics = array(NA_real_, c(1000L, length(objects), 4L), dimnames = list(NULL, names(objects), four_metrics)))
            seed <- 20261104L+(if (hour == 24L) 0L else 4L)*1000L+(if (dataset == "MIMIC_internal_test") 1L else 4L)*1000L
            set.seed(seed)
            for (b in 1:1000L) {
                w <- tabulate(sample.int(groups, groups, replace = TRUE), nbins = groups)[group_index]
                if (b %% 100L == 0L) for (model in models) {
                    stopifnot(max(abs(fast_metrics(base_objects[[model]], w)-baseline_array[b, model, ])) < 1e-9)
                    boot$baseline_vectors_checked <- boot$baseline_vectors_checked+1L
                }
                for (model in models) for (variant in variants) {
                    key <- paste(variant, model, sep = "/")
                    boot$metrics[b, key, ] <- if (identities[[key]]) baseline_array[b, model, ] else fast_metrics(objects[[key]], w)
                }
                if (b %% 250L == 0L) log_progress("GCS rule stress intervals ", label, " ", b, "/1000")
            }
            boot$completed <- 1000L; saveRDS(boot, cache)
        }
        replay_checks <- replay_checks+boot$baseline_vectors_checked
        points <- lapply(objects, fast_metrics, weights = rep(1, nrow(reference)))
        base_points <- lapply(base_objects, fast_metrics, weights = rep(1, nrow(reference)))
        for (model in models) for (variant in setdiff(variants, "strict_baseline")) for (metric in four_metrics) {
            key <- paste(variant, model, sep = "/")
            diff <- boot$metrics[, key, metric]-baseline_array[, model, metric]
            ci <- quantile_ci(diff); stopifnot(ci[["valid_reps"]] >= 950L)
            paired[[length(paired)+1L]] <- data.frame(window_hours = hour, dataset = dataset, scheme = scheme,
                model = model, variant = variant, metric = metric, difference = points[[key]][metric]-base_points[[model]][metric],
                lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]],
                unchanged_predictions = identities[[key]])
        }
        base_increment <- baseline_array[, "omit_all_flags_BMI", ]-baseline_array[, "omit_all_flags_B", ]
        base_increment_point <- base_points$omit_all_flags_BMI-base_points$omit_all_flags_B
        for (variant in variants) for (metric in four_metrics) {
            a <- paste(variant, "omit_all_flags_BMI", sep = "/"); b <- paste(variant, "omit_all_flags_B", sep = "/")
            inc <- boot$metrics[, a, metric]-boot$metrics[, b, metric]
            delta <- inc-base_increment[, metric]
            ci <- quantile_ci(inc); change <- quantile_ci(delta)
            stopifnot(ci[["valid_reps"]] >= 950L, change[["valid_reps"]] >= 950L)
            point <- points[[a]][metric]-points[[b]][metric]
            domain[[length(domain)+1L]] <- data.frame(window_hours = hour, dataset = dataset, scheme = scheme, variant = variant,
                metric = metric, domain_increment = point, lower = ci[["lower"]], upper = ci[["upper"]],
                increment_change_vs_strict = point-base_increment_point[metric], change_lower = change[["lower"]],
                change_upper = change[["upper"]], valid_reps = ci[["valid_reps"]])
        }
        populations[[label]] <- data.frame(window_hours = hour, dataset = dataset, scheme = scheme, records = nrow(reference),
            patients = length(unique(reference$group_id)), events = sum(reference$truth == "yes"), resampling_blocks = groups)
    }
}
write.csv(do.call(rbind, paired), file.path(out, "paired_input_rule_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, domain), file.path(out, "domain_increment_stress.csv"), row.names = FALSE)
write.csv(do.call(rbind, populations), file.path(out, "evaluation_populations.csv"), row.names = FALSE)
write_json(list(status = "fixed_model_rule_stress_uncertainty_complete_not_refitted_validation", run_hash = hash,
    fast_metric_engine_duplicated_record_checks = engine_checks, baseline_bootstrap_vectors_reproduced = replay_checks,
    models_refitted = FALSE, calibration_change_intervals = FALSE), file.path(out, "validation_complete.json"))
