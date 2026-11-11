# Replay prior validation weights for post-hoc missingness sensitivity comparisons.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "missingness_sensitivity_results")
strategies <- c("omit_neuro_flag", "omit_all_flags")
contrasts <- list(c("BM", "B"), c("BI", "B"), c("BMI", "B"), c("BMI", "BM"), c("BMI", "BI"))
fit_plan <- jsonlite::read_json(file.path(out, "training_plan.json"))
prior_plan <- jsonlite::read_json(file.path(root, "neuro_reference_results/validation_plan.json"))
plan <- list(version = "posthoc_missing_flags_intervals_1.0", fit_hash = fit_plan$run_hash,
    prior_validation_hash = prior_plan$run_hash, posthoc = TRUE, replicates = 1000L,
    script_sha256 = sha_file(file.path(root, "R/validate_missingness_sensitivity.R")),
    prediction_hashes = setNames(lapply(c(24L, 6L), function(h)
        sha_file(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", h, "h.csv")))), c("24", "6")),
    prior_prediction_hashes = setNames(lapply(c(24L, 6L), function(h)
        sha_file(file.path(root, "private_predictions", paste0("neuro_predictions_", h, "h.csv")))), c("24", "6")),
    comparison = "same weights within nested models and paired to original motor_B and motor_BMI",
    inference = "conditional frozen-model percentile 95% intervals; posthoc exploratory multiple comparisons")
hash <- digest::digest(plan, algo = "sha256")
path <- file.path(out, "validation_plan.json")
if (file.exists(path)) stopifnot(identical(jsonlite::read_json(path)$run_hash, hash)) else
    write_json(list(run_hash = hash, plan = plan), path)
metric_rows <- domain_rows <- original_rows <- populations <- list()
for (hour in c(24L, 6L)) {
    raw <- read.csv(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv")))
    prior_raw <- read.csv(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")))
    for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
        for (population in if (dataset == "MIMIC_internal_test") "all_eligible" else c("all_eligible", "neuro_observed")) {
            chosen <- raw[raw$dataset == dataset, ]
            reference <- chosen[chosen$model == "omit_neuro_flag_B", ]
            if (population == "neuro_observed") reference <- reference[reference$neuro_observed, ]
            reference <- reference[order(reference$stay_id), ]
            stopifnot(!anyDuplicated(reference$stay_id))
            models <- unique(chosen$model); objects <- old_objects <- list()
            aligned_object <- function(d) {
                d <- d[match(reference$stay_id, d$stay_id), ]
                stopifnot(!anyNA(d$probability), identical(d$truth, reference$truth), identical(d$group_id, reference$group_id),
                    identical(as.character(d$hospitalid), as.character(reference$hospitalid)))
                prepare_probability(as.integer(d$truth == "yes"), d$probability)
            }
            for (model in models) objects[[model]] <- aligned_object(chosen[chosen$model == model, ])
            for (domain in c("B", "BMI")) old_objects[[domain]] <- aligned_object(
                prior_raw[prior_raw$dataset == dataset & prior_raw$model == paste0("motor_", domain), ])
            cluster <- reference$group_id; scheme <- "patient"; hospitals <- NA_integer_
            if (dataset == "eICU_external_validation") {
                cluster <- hospital_components(reference$group_id, reference$hospitalid)$row_component
                scheme <- "hospital_patient_component"; hospitals <- length(unique(reference$hospitalid))
            }
            group_index <- match(cluster, unique(cluster)); groups <- length(unique(cluster))
            label <- paste(hour, dataset, population, sep = "_")
            old_label <- paste(label, if (population == "all_eligible") "motor_total" else "motor", sep = "_")
            old <- readRDS(file.path(root, "private_cache", paste0("neuro_validation_", old_label, ".rds")))
            stopifnot(identical(old$run_hash, prior_plan$run_hash), old$completed == 1000L)
            populations[[label]] <- data.frame(window_hours = hour, dataset = dataset, population = population,
                records = nrow(reference), patients = length(unique(reference$group_id)),
                events = sum(reference$truth == "yes"), hospitals = hospitals, resampling_blocks = groups)
            seed <- if (population == "all_eligible") 20261104L + (if (hour == 24L) 0L else 4L)*1000L +
                (if (dataset == "MIMIC_internal_test") 1L else 4L)*1000L else 20261104L+30000L+hour*100L+1L
            cache <- file.path(root, "private_cache", paste0("missing_validation_", label, ".rds"))
            checkpoint <- sub("\\.rds$", "_checkpoint.rds", cache)
            if (file.exists(cache)) {
                boot <- readRDS(cache); stopifnot(identical(boot$run_hash, hash), boot$completed == 1000L)
            } else {
                if (file.exists(checkpoint)) {
                    boot <- readRDS(checkpoint); stopifnot(identical(boot$run_hash, hash))
                    assign(".Random.seed", boot$rng, envir = .GlobalEnv); start <- boot$completed+1L
                } else {
                    set.seed(seed); start <- 1L
                    boot <- list(run_hash = hash, completed = 0L, replay_checks = 0L,
                        metrics = array(NA_real_, c(1000L, length(models), length(metric_names)), dimnames = list(NULL, models, metric_names)))
                }
                if (start <= 1000L) for (b in start:1000L) {
                    weights <- tabulate(sample.int(groups, groups, replace = TRUE), nbins = groups)[group_index]
                    if (b %% 100L == 0L) {
                        check <- weighted_metrics(old_objects$B, weights)
                        stopifnot(max(abs(check[c("auroc", "brier_score")]-old$metrics[b, "motor_B", c("auroc", "brier_score")])) < 1e-10)
                        boot$replay_checks <- boot$replay_checks+1L
                    }
                    for (model in models) boot$metrics[b, model, ] <- weighted_metrics(objects[[model]], weights)
                    if (b %% 250L == 0L) {
                        boot$completed <- b; boot$rng <- .Random.seed; saveRDS(boot, checkpoint)
                        log_progress("Missingness sensitivity validation ", label, " ", b, "/1000")
                    }
                }
                saveRDS(boot, cache)
            }
            stopifnot(boot$completed == 1000L, boot$replay_checks == 10L)
            points <- lapply(objects, weighted_metrics, weights = rep(1, nrow(reference)))
            old_points <- lapply(old_objects, weighted_metrics, weights = rep(1, nrow(reference)))
            for (model in models) for (metric in metric_names) {
                ci <- quantile_ci(boot$metrics[, model, metric]); stopifnot(ci[["valid_reps"]] >= 950L)
                metric_rows[[length(metric_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset, population = population,
                    scheme = scheme, model = model, metric = metric, estimate = points[[model]][metric],
                    lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
            }
            for (strategy in strategies) {
                for (contrast in contrasts) for (metric in c("auroc", "auprc", "brier_score", "log_loss")) {
                    a <- paste(strategy, contrast[1L], sep = "_"); b <- paste(strategy, contrast[2L], sep = "_")
                    ci <- quantile_ci(boot$metrics[, a, metric]-boot$metrics[, b, metric])
                    domain_rows[[length(domain_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset, population = population,
                        strategy = strategy, augmented = a, reference = b, contrast = paste(contrast, collapse = "_minus_"), metric = metric,
                        difference = points[[a]][metric]-points[[b]][metric], lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
                }
                for (domain in c("B", "BMI")) for (metric in c("auroc", "auprc", "brier_score", "log_loss")) {
                    a <- paste(strategy, domain, sep = "_"); b <- paste0("motor_", domain)
                    ci <- quantile_ci(boot$metrics[, a, metric]-old$metrics[, b, metric])
                    original_rows[[length(original_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset, population = population,
                        strategy = strategy, model = a, original_model = b, metric = metric,
                        difference_new_minus_original = points[[a]][metric]-old_points[[domain]][metric],
                        lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
                }
            }
        }
    }
}
write.csv(do.call(rbind, metric_rows), file.path(out, "metric_intervals.csv"), row.names = FALSE)
write.csv(do.call(rbind, domain_rows), file.path(out, "paired_domain_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, original_rows), file.path(out, "paired_original_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, populations), file.path(out, "evaluation_populations.csv"), row.names = FALSE)
write_json(list(status = "posthoc_conditional_sensitivity_complete_not_new_independent_validation", run_hash = hash,
    resampling_configurations = length(populations), replicates_each = 1000L, original_models_changed = FALSE),
    file.path(out, "validation_complete.json"))
