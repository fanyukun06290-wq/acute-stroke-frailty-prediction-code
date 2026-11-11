Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/fit_missingness_sensitivity.R"), encoding = "UTF-8")
plan <- jsonlite::read_json(file.path(out, "training_plan.json"))
stopifnot(identical(sha_file(file.path(root, "R/fit_missingness_sensitivity.R")), plan$plan$script_sha256),
    identical(sha_file(file.path(root, "14_第八步缺失指示稳健性规则.md")), plan$plan$rules_sha256))
for (p in names(plan$plan$inputs_sha256)) stopifnot(identical(sha_file(file.path(root, p)), plan$plan$inputs_sha256[[p]]))
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
selected <- read.csv(file.path(out, "selected_parameters.csv"))
max_error <- 0; median_count <- 0L
for (hour in c(24L, 6L)) {
    mimic <- load_neuro("mimic", hour, contract); external <- load_neuro("eicu", hour, contract)
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", hour, "h.csv")))
    mapped <- match(mimic$metadata$group_id, registry$group_id)
    train <- which(registry$partition[mapped] == "train"); test <- which(registry$partition[mapped] == "test")
    saved <- read.csv(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv")))
    stopifnot(!length(intersect(mimic$metadata$group_id[train], mimic$metadata$group_id[test])))
    for (strategy in strategies) for (domain in domains) {
        name <- paste(strategy, domain, sep = "_"); tag <- paste0(name, "_", hour, "h")
        model <- readRDS(file.path(root, "private_models", paste0("missing_sensitivity_", tag, ".rds")))
        stopifnot(identical(model$run_hash, plan$run_hash), identical(model$raw_columns, sensitivity_columns(mimic$x, strategy, domain)),
            !"gcs_motor_missing" %in% model$raw_columns,
            !any(grepl("(group_id|stay_id|phenotype|comorbidity|outcome|death|offset)", model$raw_columns)))
        if (strategy == "omit_all_flags") stopifnot(!any(grepl("_missing($|_)", model$raw_columns)))
        grid <- read.csv(file.path(out, paste0("oof_grid_", tag, ".csv")))
        fold <- read.csv(file.path(out, paste0("fold_grid_", tag, ".csv")))
        stopifnot(nrow(grid) == 37L, nrow(fold) == 185L, all(is.finite(grid$log_loss)),
            abs(model$lambda-grid$lambda[order(grid$log_loss, -grid$lambda)[1L]]) < 1e-14)
        oof <- read.csv(file.path(root, "private_predictions", paste0("missing_oof_", tag, ".csv")))
        idx <- match(oof$group_id, registry$group_id)
        stopifnot(all(registry$partition[idx] == "train"), identical(oof$fold, registry$fold[idx]),
            !anyDuplicated(unique(oof[, c("group_id", "fold")])$group_id))
        r <- selected[selected$window_hours == hour & selected$model == name, ]
        s <- cv_scores(as.integer(oof$truth == "yes"), oof$probability)
        stopifnot(max(abs(s-unlist(r[, names(s)]))) < 1e-9)
        step <- which(vapply(model$recipe$steps, inherits, logical(1), what = "step_impute_median"))
        medians <- unlist(model$recipe$steps[[step]]$medians, use.names = TRUE)
        expected <- vapply(mimic$x[train, names(medians), drop = FALSE], median, numeric(1), na.rm = TRUE)
        stopifnot(max(abs(medians-expected)) < 1e-12); median_count <- median_count+length(medians)
        for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
            d <- if (dataset == "MIMIC_internal_test")
                list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ]) else external
            if (dataset == "eICU_external_validation") d$meta <- d$metadata
            x <- as.matrix(recipes::bake(model$recipe, new_data = d$x[, model$raw_columns, drop = FALSE]))
            p <- as.numeric(predict(model$fit, newx = x, s = model$lambda, type = "response"))
            stored <- saved[saved$model == name & saved$dataset == dataset, ]
            stored <- stored[match(d$meta$stay_id, stored$stay_id), ]
            stopifnot(!anyNA(stored$probability), identical(stored$truth, ifelse(d$y == 1, "yes", "no")),
                identical(stored$group_id, d$meta$group_id))
            delta <- max(abs(stored$probability-p)); stopifnot(delta < 1e-12); max_error <- max(max_error, delta)
        }
    }
}
write_json(list(status = "passed_bounded_training_and_prediction_checks", models_reloaded = 16L,
    full_prediction_sets_reproduced = 32L, maximum_probability_difference = max_error,
    training_medians_verified = median_count, candidate_fold_scores = 2960L,
    excluded_missingness_predictors_verified = TRUE, posthoc = TRUE, original_models_changed = FALSE),
    file.path(root, "qc/missing_sensitivity_fit_checks.json"))
cat("MISSINGNESS_SENSITIVITY_FIT_CHECKS_PASSED\n")
