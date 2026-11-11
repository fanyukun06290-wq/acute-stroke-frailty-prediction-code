# Verify saved model transport, patient partitions, and training-only medians.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R", "run_revision_models.R"), encoding = "UTF-8")
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
plan <- jsonlite::read_json(file.path(root, "model_results", "training_plan.json"), simplifyVector = FALSE)
stopifnot(identical(plan$fingerprints$contract_sha256, sha_file(file.path(root, "model_feature_contract.json"))))
versions <- read.csv(file.path(root, "model_results", "package_versions.csv"))
stopifnot(all(versions$version == vapply(versions$package, function(x) as.character(packageVersion(x)), character(1))))
for (source in c("mimic", "eicu")) stopifnot(identical(
    plan$fingerprints[[paste0(source, "_sha256")]],
    sha_file(file.path(root, "private_data", paste0(source, "_analysis_revision.csv")))))
checks <- list()
for (hour in c(24, 6)) {
    tag <- paste0(hour, "h")
    complete <- jsonlite::read_json(file.path(root, "model_results", paste0("complete_", tag, ".json")))
    stopifnot(identical(complete$run_hash, plan$run_hash),
        identical(plan$fingerprints$script_sha256, sha_file(file.path(root, "R", "run_revision_models.R"))))
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", tag, ".csv")))
    stopifnot(!anyDuplicated(registry$group_id),
        !length(intersect(registry$group_id[registry$partition == "train"],
                          registry$group_id[registry$partition == "test"])))
    mimic <- read_data(file.path(root, "private_data", "mimic_analysis_revision.csv"), contract, hour)
    eicu <- read_data(file.path(root, "private_data", "eicu_analysis_revision.csv"), contract, hour)
    mapped <- match(mimic$metadata$group_id, registry$group_id)
    stopifnot(!anyNA(mapped))
    train <- which(registry$partition[mapped] == "train")
    test <- which(registry$partition[mapped] == "test")
    rec <- readRDS(file.path(root, "private_models", paste0("recipe_", tag, ".rds")))
    stopifnot(setequal(rec$var_info$variable[rec$var_info$role == "predictor"], names(mimic$x)))
    median_step <- which(vapply(rec$steps, inherits, logical(1), what = "step_impute_median"))
    learned <- unlist(rec$steps[[median_step]]$medians, use.names = TRUE)
    expected <- vapply(mimic$x[train, names(learned), drop = FALSE], median, numeric(1), na.rm = TRUE)
    stopifnot(max(abs(learned-expected)) < 1e-12)
    parameters <- jsonlite::read_json(file.path(root, "model_results", paste0("selected_parameters_", tag, ".json")), simplifyVector = TRUE)
    preds <- read.csv(file.path(root, "private_predictions", paste0("predictions_", tag, ".csv")))
    errors <- list()
    for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
        item <- if (dataset == "MIMIC_internal_test")
            list(x = mimic$x[test, , drop = FALSE], y = mimic$y[test], meta = mimic$metadata[test, ])
        else list(x = eicu$x, y = eicu$y, meta = eicu$metadata)
        matrix <- as.matrix(bake(rec, new_data = item$x)); storage.mode(matrix) <- "double"
        p <- list(
            xgboost = as.numeric(predict(xgboost::xgb.load(file.path(root, "private_models", paste0("xgboost_", tag, ".ubj"))),
                xgboost::xgb.DMatrix(matrix))),
            random_forest = predict(readRDS(file.path(root, "private_models", paste0("random_forest_", tag, ".rds"))),
                as.data.frame(matrix), num.threads = 2L)$predictions[, "yes"],
            elastic_net = as.numeric(predict(readRDS(file.path(root, "private_models", paste0("elastic_net_", tag, ".rds"))),
                newx = matrix, s = parameters$elastic_net$lambda, type = "response")),
            clinical_reference = predict_clinical(readRDS(file.path(root, "private_models", paste0("clinical_reference_", tag, ".rds"))), matrix),
            constant_training_prevalence = rep(readRDS(file.path(root, "private_models", paste0("constant_training_prevalence_", tag, ".rds"))), length(item$y)))
        for (name in names(p)) {
            saved <- preds[preds$dataset == dataset & preds$model == name, ]
            stopifnot(identical(as.character(saved$stay_id), as.character(item$meta$stay_id)),
                      identical(as.character(saved$group_id), as.character(item$meta$group_id)),
                      identical(saved$truth, ifelse(item$y == 1, "yes", "no")))
            difference <- max(abs(saved$probability-p[[name]]))
            stopifnot(difference < 1e-10)
            errors[[paste(dataset, name, sep = "/")]] <- difference
        }
    }
    cv <- read.csv(file.path(root, "model_results", paste0("cv_fold_scores_", tag, ".csv")))
    stopifnot(nrow(cv) == (20L+9L+125L+2L)*5L, all(is.finite(cv$auroc)))
    oof <- read.csv(file.path(root, "private_predictions", paste0("selected_oof_", tag, ".csv")))
    stopifnot(all(oof$group_id %in% registry$group_id[registry$partition == "train"]))
    oof_group_folds <- unique(oof[, c("group_id", "fold")])
    stopifnot(!anyDuplicated(oof_group_folds$group_id))
    checks[[tag]] <- list(patient_partition_overlap = 0, patient_fold_overlap = 0,
        exact_training_medians_verified = TRUE,
        saved_models_reproduce_all_internal_external_probabilities = TRUE,
        maximum_probability_difference = max(unlist(errors)),
        candidate_fold_results = nrow(cv), predictor_allowlist_verified = TRUE)
}
write_json(list(scope = "Training and point-prediction artifact checks only; not final paper integrity certification",
                checks = checks), file.path(root, "qc", "model_artifact_checks.json"))
cat("MODEL_ARTIFACT_CHECKS_PASS\n")
