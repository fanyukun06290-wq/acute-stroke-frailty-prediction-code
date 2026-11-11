# Conditional paired validation of the frozen same-learner nested models.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "neuro_reference_results")
fit_plan <- jsonlite::read_json(file.path(out, "training_plan.json"))
old_plan <- jsonlite::read_json(file.path(root, "validation_results/uncertainty_plan.json"))
contrasts <- list(c("BM", "B"), c("BI", "B"), c("BMI", "B"), c("BMI", "BM"), c("BMI", "BI"))
plan <- list(version = "nested_neuro_intervals_1.0", fit_hash = fit_plan$run_hash, replicates = 1000L,
    script_sha256 = sha_file(file.path(root, "R/validate_neuro_nested.R")),
    prediction_hashes = setNames(lapply(c(24L, 6L), function(h)
        sha_file(file.path(root, "private_predictions", paste0("neuro_predictions_", h, "h.csv")))), c("24", "6")),
    primary = "patient bootstrap internally; connected hospital/shared-patient block bootstrap externally",
    sensitivity = "same frozen predictions restricted to observed neurological score in eICU; no refitting",
    contrasts = contrasts, interval = "95% percentile, type 7, conditional on frozen models; exploratory multiple comparisons",
    primary_seeds = c(20262104L, 20265104L, 20266104L, 20269104L))
hash <- digest::digest(plan, algo = "sha256")
path <- file.path(out, "validation_plan.json")
if (file.exists(path)) stopifnot(identical(jsonlite::read_json(path)$run_hash, hash)) else
    write_json(list(run_hash = hash, plan = plan), path)
metric_rows <- pair_rows <- populations <- list()
for (hour in c(24L, 6L)) {
    raw <- read.csv(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")))
    for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
        configurations <- list(list(population = "all_eligible", neuro = c("motor", "total")))
        if (dataset == "eICU_external_validation") configurations <- c(configurations,
            list(list(population = "neuro_observed", neuro = "motor"), list(population = "neuro_observed", neuro = "total")))
        for (configuration in configurations) {
            chosen <- raw[raw$dataset == dataset & raw$neuro %in% configuration$neuro, ]
            reference <- chosen[chosen$model == paste0(configuration$neuro[1L], "_B"), ]
            if (configuration$population == "neuro_observed") reference <- reference[reference$neuro_observed, ]
            reference <- reference[order(reference$stay_id), ]
            stopifnot(!anyDuplicated(reference$stay_id))
            models <- unique(chosen$model); objects <- list()
            for (model in models) {
                d <- chosen[chosen$model == model, ]; d <- d[match(reference$stay_id, d$stay_id), ]
                stopifnot(!anyNA(d$probability), identical(d$truth, reference$truth), identical(d$group_id, reference$group_id),
                    identical(as.character(d$hospitalid), as.character(reference$hospitalid)))
                objects[[model]] <- prepare_probability(as.integer(d$truth == "yes"), d$probability)
            }
            cluster <- reference$group_id
            scheme <- "patient"; hospitals <- cross_patients <- NA_integer_
            if (dataset == "eICU_external_validation") {
                comp <- hospital_components(reference$group_id, reference$hospitalid)
                cluster <- comp$row_component; scheme <- "hospital_patient_component"
                hospitals <- length(unique(reference$hospitalid)); cross_patients <- comp$cross_hospital_patients
            }
            group_index <- match(cluster, unique(cluster)); groups <- length(unique(cluster))
            label <- paste(hour, dataset, configuration$population, paste(configuration$neuro, collapse = "_"), sep = "_")
            populations[[label]] <- data.frame(window_hours = hour, dataset = dataset, population = configuration$population,
                neuro_types = paste(configuration$neuro, collapse = ","), records = nrow(reference),
                patients = length(unique(reference$group_id)), events = sum(reference$truth == "yes"),
                hospitals = hospitals, patients_crossing_hospitals = cross_patients, resampling_blocks = groups)
            primary <- configuration$population == "all_eligible"
            seed <- if (primary) 20261104L + (if (hour == 24L) 0L else 4L)*1000L +
                (if (dataset == "MIMIC_internal_test") 1L else 4L)*1000L else
                20261104L + 30000L + hour*100L + if (identical(configuration$neuro, "motor")) 1L else 2L
            old <- original <- NULL
            if (primary) {
                old <- readRDS(file.path(root, "private_cache", paste0("uncertainty_", hour, "_", dataset, "_", scheme, ".rds")))
                stopifnot(identical(old$run_hash, old_plan$run_hash), old$clusters == groups)
                original <- align_predictions(file.path(root, "private_predictions", paste0("predictions_", hour, "h.csv")), dataset)
                stopifnot(identical(original$reference$stay_id, reference$stay_id))
            }
            cache <- file.path(root, "private_cache", paste0("neuro_validation_", label, ".rds"))
            checkpoint <- sub("\\.rds$", "_checkpoint.rds", cache)
            if (file.exists(cache)) {
                boot <- readRDS(cache); stopifnot(identical(boot$run_hash, hash), boot$completed == 1000L)
            } else {
                if (file.exists(checkpoint)) {
                    boot <- readRDS(checkpoint); stopifnot(identical(boot$run_hash, hash))
                    assign(".Random.seed", boot$rng, envir = .GlobalEnv); start <- boot$completed+1L
                } else {
                    set.seed(seed); start <- 1L
                    boot <- list(run_hash = hash, completed = 0L, metrics = array(NA_real_, c(1000L, length(models), length(metric_names)),
                        dimnames = list(NULL, models, metric_names)), weight_replay_checks = 0L)
                }
                if (start <= 1000L) for (b in start:1000L) {
                    weights <- tabulate(sample.int(groups, groups, replace = TRUE), nbins = groups)[group_index]
                    if (primary) {
                        stopifnot(sum(weights) == old$records[b], sum(weights*(reference$truth == "yes")) == old$events[b])
                        if (b %% 100L == 0L) {
                            check <- weighted_metrics(original$objects$xgboost, weights)
                            stopifnot(max(abs(check[c("auroc", "brier_score")] - old$metrics[b, "xgboost", c("auroc", "brier_score")])) < 1e-10)
                            boot$weight_replay_checks <- boot$weight_replay_checks+1L
                        }
                    }
                    for (model in models) boot$metrics[b, model, ] <- weighted_metrics(objects[[model]], weights)
                    if (b %% 250L == 0L) {
                        boot$completed <- b; boot$rng <- .Random.seed
                        saveRDS(boot, checkpoint)
                        log_progress("Nested validation ", label, " ", b, "/1000")
                    }
                }
                saveRDS(boot, cache)
            }
            stopifnot(boot$completed == 1000L)
            points <- lapply(objects, weighted_metrics, weights = rep(1, nrow(reference)))
            for (model in models) for (metric in metric_names) {
                ci <- quantile_ci(boot$metrics[, model, metric])
                stopifnot(ci[["valid_reps"]] >= 950L)
                metric_rows[[length(metric_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset,
                    population = configuration$population, scheme = scheme, model = model, metric = metric,
                    estimate = points[[model]][metric], lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
            }
            for (neuro in configuration$neuro) for (contrast in contrasts) for (metric in c("auroc", "auprc", "brier_score", "log_loss")) {
                a <- paste(neuro, contrast[1L], sep = "_"); b <- paste(neuro, contrast[2L], sep = "_")
                ci <- quantile_ci(boot$metrics[, a, metric] - boot$metrics[, b, metric])
                stopifnot(ci[["valid_reps"]] >= 950L)
                pair_rows[[length(pair_rows)+1L]] <- data.frame(window_hours = hour, dataset = dataset,
                    population = configuration$population, scheme = scheme, neuro = neuro,
                    augmented = a, reference = b, contrast = paste(contrast, collapse = "_minus_"), metric = metric,
                    difference = points[[a]][metric]-points[[b]][metric], lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]])
            }
        }
    }
}
write.csv(do.call(rbind, metric_rows), file.path(out, "metric_intervals.csv"), row.names = FALSE)
write.csv(do.call(rbind, pair_rows), file.path(out, "paired_differences.csv"), row.names = FALSE)
write.csv(do.call(rbind, populations), file.path(out, "evaluation_populations.csv"), row.names = FALSE)
write_json(list(status = "conditional_nested_validation_complete_not_final_paper_ready", run_hash = hash,
    resampling_configurations = length(populations), replicates_each = 1000L, same_weights_within_comparisons = TRUE,
    main_model_changed = FALSE), file.path(out, "validation_complete.json"))
