# Fixed-model measurement substitution. No fitting, refitting, or recalibration.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/fit_neuro_nested.R"), encoding = "UTF-8")
out <- file.path(root, "gcs_rule_sensitivity_results")
variants <- c("strict_baseline", "allow_reverse_recording", "first_valid", "last_valid", "motor_requires_complete", "airway_60", "airway_240")
stress_models <- c("motor_B", "motor_BMI", "total_B", "total_BMI", "omit_neuro_flag_B", "omit_neuro_flag_BMI", "omit_all_flags_B", "omit_all_flags_BMI")
evaluate_rule_stress <- function() {
    inputs <- c("qc/runtime_neuro_files.json", "qc/runtime_missing_files.json",
        "gcs_rule_sensitivity_results/score_construction_complete.json", "private_data/mimic_gcs_rule_variants.csv", "private_data/eicu_gcs_rule_variants.csv")
    plan <- list(version = "frozen_gcs_rule_stress_1.0", posthoc = TRUE, models_refitted = FALSE, variants = variants,
        frozen_models = stress_models, script_sha256 = sha_file(file.path(root, "R/evaluate_gcs_rule_stress.R")),
        rules_sha256 = sha_file(file.path(root, "16_第九步GCS规则敏感性方案.md")),
        inputs_sha256 = setNames(lapply(inputs, function(p) sha_file(file.path(root, p))), inputs),
        interpretation = "Input-rule stress test of fixed models, not performance of models trained using new score definitions")
    hash <- digest::digest(plan, algo = "sha256")
    path <- file.path(out, "evaluation_plan.json")
    if (file.exists(path)) stop("Preserve existing stress evaluation")
    write_json(list(run_hash = hash, plan = plan), path)
    contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
    points <- manifest <- list(); maximum_baseline_difference <- 0
    for (hour in c(24L, 6L)) {
        mimic <- load_neuro("mimic", hour, contract); external <- load_neuro("eicu", hour, contract)
        registry <- read.csv(file.path(root, "private_splits", paste0("patients_", hour, "h.csv")))
        test <- which(registry$partition[match(mimic$metadata$group_id, registry$group_id)] == "test")
        sets <- list(MIMIC_internal_test = list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ], source = "mimic"),
            eICU_external_validation = list(x = external$x, y = external$y, meta = external$metadata, source = "eicu"))
        original <- read.csv(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")))
        missing <- read.csv(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv")))
        score_files <- list(mimic = read.csv(file.path(root, "private_data/mimic_gcs_rule_variants.csv")),
            eicu = read.csv(file.path(root, "private_data/eicu_gcs_rule_variants.csv")))
        predictions <- list()
        for (name in stress_models) {
            prefix <- if (grepl("^omit_", name)) "missing_sensitivity_" else "neuro_"
            modelpath <- file.path(root, "private_models", paste0(prefix, name, "_", hour, "h.rds"))
            model <- readRDS(modelpath)
            manifest[[paste(name, hour)]] <- list(path = modelpath, sha256 = sha_file(modelpath))
            for (dataset in names(sets)) {
                d <- sets[[dataset]]
                for (variant in variants) {
                    f <- score_files[[d$source]]
                    f <- f[f$window_hours == hour & f$variant == variant, ]
                    stopifnot(!anyDuplicated(f$stay_id))
                    idx <- match(d$meta$stay_id, f$stay_id); stopifnot(!anyNA(idx))
                    f <- f[idx, ]; x <- d$x
                    for (variable in c("gcs_motor", "gcs_total")) {
                        x[[paste0(variable, "_min")]] <- f[[paste0(variable, "_min")]]
                        x[[paste0(variable, "_missing")]] <- f[[paste0(variable, "_missing")]]
                    }
                    m <- as.matrix(recipes::bake(model$recipe, new_data = x[, model$raw_columns, drop = FALSE]))
                    stopifnot(identical(colnames(m), model$encoded_columns), !anyNA(m))
                    p <- as.numeric(predict(model$fit, newx = m, s = model$lambda, type = "response"))
                    stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1))
                    if (variant == "strict_baseline") {
                        old <- if (grepl("^omit_", name)) missing else original
                        old <- old[old$dataset == dataset & old$model == name, ]
                        old <- old[match(d$meta$stay_id, old$stay_id), ]
                        stopifnot(!anyNA(old$probability), identical(old$truth, ifelse(d$y == 1, "yes", "no")))
                        delta <- max(abs(p-old$probability)); stopifnot(delta < 1e-12)
                        maximum_baseline_difference <- max(maximum_baseline_difference, delta)
                    }
                    points[[length(points)+1L]] <- data.frame(window_hours = hour, dataset = dataset, model = name,
                        variant = variant, records = length(p), patients = length(unique(d$meta$group_id)), events = sum(d$y),
                        as.list(point_scores(d$y, p)))
                    predictions[[length(predictions)+1L]] <- data.frame(window_hours = hour, dataset = dataset, model = name,
                        variant = variant, d$meta[, c("stay_id", "group_id", "hospitalid")],
                        truth = ifelse(d$y == 1, "yes", "no"), probability = p)
                }
            }
        }
        write.csv(do.call(rbind, predictions), file.path(root, "private_predictions", paste0("gcs_rule_stress_", hour, "h.csv")), row.names = FALSE)
        log_progress("GCS fixed-model stress ", hour, "h complete")
    }
    points <- do.call(rbind, points)
    baseline <- points[points$variant == "strict_baseline", ]
    for (metric in c("auroc", "auprc", "brier_score", "log_loss", "calibration_intercept", "calibration_slope")) {
        key <- paste(points$window_hours, points$dataset, points$model)
        base_key <- paste(baseline$window_hours, baseline$dataset, baseline$model)
        points[[paste0("difference_", metric)]] <- points[[metric]]-baseline[[metric]][match(key, base_key)]
    }
    write.csv(points, file.path(out, "performance_points.csv"), row.names = FALSE)
    write_json(manifest, file.path(out, "frozen_model_manifest.json"))
    write_json(list(status = "fixed_model_score_perturbation_evaluated_not_retrained_validation", run_hash = hash,
        prediction_groups = nrow(points), frozen_models_reloaded = 16L, maximum_strict_baseline_prediction_difference = maximum_baseline_difference,
        model_training_performed = FALSE, score_imputation_updated = FALSE), file.path(out, "evaluation_complete.json"))
}
if (sys.nframe() == 0L) evaluate_rule_stress()
