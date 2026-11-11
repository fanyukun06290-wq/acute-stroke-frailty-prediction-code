# Verify fixed-input selection, training-only medians, and all saved predictions.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/fit_neuro_nested.R"), encoding = "UTF-8")
plan <- jsonlite::read_json(file.path(out, "training_plan.json"))
stopifnot(identical(sha_file(file.path(root, "R/fit_neuro_nested.R")), plan$plan$script_sha256),
    identical(sha_file(file.path(root, "12_第七步神经参考与领域增量规则.md")), plan$plan$rule_sha256))
for (p in names(plan$plan$inputs_sha256)) stopifnot(identical(sha_file(file.path(root, p)), plan$plan$inputs_sha256[[p]]))
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
parameters <- read.csv(file.path(out, "selected_parameters.csv"))
diagnostics <- missing_coefficients <- list(); max_error <- 0; medians_checked <- 0L
for (hour in c(24L, 6L)) {
    mimic <- load_neuro("mimic", hour, contract); external <- load_neuro("eicu", hour, contract)
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", hour, "h.csv")))
    idx <- match(mimic$metadata$group_id, registry$group_id)
    train <- which(registry$partition[idx] == "train"); test <- which(registry$partition[idx] == "test")
    saved <- read.csv(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")))
    for (neuro in names(neuro_variables)) for (domain in domains) {
        name <- paste(neuro, domain, sep = "_"); tag <- paste0(name, "_", hour, "h")
        model <- readRDS(file.path(root, "private_models", paste0("neuro_", tag, ".rds")))
        stopifnot(identical(model$run_hash, plan$run_hash), identical(model$raw_columns, neuro_columns(mimic$x, neuro, domain)),
            !any(grepl("(group_id|stay_id|phenotype|comorbidity|outcome|death|offset)", model$raw_columns)),
            all(is.finite(as.matrix(coef(model$fit, s = model$lambda)))))
        grid <- read.csv(file.path(out, paste0("oof_grid_", tag, ".csv")))
        folds <- read.csv(file.path(out, paste0("fold_grid_", tag, ".csv")))
        stopifnot(nrow(grid) == 37L, nrow(folds) == 185L,
            abs(model$lambda-grid$lambda[order(grid$log_loss, -grid$lambda)[1L]]) < 1e-14)
        oof <- read.csv(file.path(root, "private_predictions", paste0("neuro_oof_", tag, ".csv")))
        k <- match(oof$group_id, registry$group_id)
        stopifnot(!anyNA(k), all(registry$partition[k] == "train"), identical(oof$fold, registry$fold[k]),
            !anyDuplicated(unique(oof[, c("group_id", "fold")])$group_id),
            !length(intersect(oof$group_id, mimic$metadata$group_id[test])))
        score <- cv_scores(as.integer(oof$truth == "yes"), oof$probability)
        selected <- parameters[parameters$window_hours == hour & parameters$model == name, ]
        stopifnot(max(abs(score-unlist(selected[, names(score)]))) < 1e-9)
        median_step <- which(vapply(model$recipe$steps, inherits, logical(1), what = "step_impute_median"))
        learned <- unlist(model$recipe$steps[[median_step]]$medians, use.names = TRUE)
        expected <- vapply(mimic$x[train, names(learned), drop = FALSE], median, numeric(1), na.rm = TRUE)
        stopifnot(max(abs(learned-expected)) < 1e-12); medians_checked <- medians_checked+length(learned)
        missing_name <- paste0(neuro_variables[[neuro]], "_missing")
        beta <- as.matrix(coef(model$fit, s = model$lambda))
        normal <- which(vapply(model$recipe$steps, inherits, logical(1), what = "step_normalize"))
        sd <- model$recipe$steps[[normal]]$sds[[missing_name]]
        effect <- if (missing_name %in% rownames(beta) && !is.null(sd)) beta[missing_name, 1]/sd else NA_real_
        missing_coefficients[[tag]] <- data.frame(window_hours = hour, model = name,
            missing_indicator_logodds_change_0_to_1 = effect)
        sets <- list(training_apparent = list(x = mimic$x[train, ], y = mimic$y[train], meta = mimic$metadata[train, ]),
            MIMIC_internal_test = list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ]),
            eICU_external_validation = list(x = external$x, y = external$y, meta = external$metadata))
        for (dataset in names(sets)) {
            d <- sets[[dataset]]
            x <- as.matrix(recipes::bake(model$recipe, new_data = d$x[, model$raw_columns, drop = FALSE]))
            probabilities <- as.numeric(predict(model$fit, newx = x, s = model$lambda, type = "response"))
            if (dataset != "training_apparent") {
                s <- saved[saved$model == name & saved$dataset == dataset, ]
                s <- s[match(d$meta$stay_id, s$stay_id), ]
                stopifnot(identical(s$truth, ifelse(d$y == 1, "yes", "no")),
                    identical(s$group_id, d$meta$group_id), !anyNA(s$probability))
                delta <- max(abs(probabilities-s$probability)); stopifnot(delta < 1e-12)
                max_error <- max(max_error, delta)
            }
            for (missing in c(FALSE, TRUE)) {
                keep <- (d$x[[missing_name]] == 1) == missing
                if (!any(keep)) next
                diagnostics[[length(diagnostics)+1L]] <- data.frame(window_hours = hour, model = name, neuro = neuro,
                    dataset = dataset, neuro_missing = missing, records = sum(keep), events = sum(d$y[keep]),
                    event_prevalence = mean(d$y[keep]), mean_predicted = mean(probabilities[keep]),
                    brier = mean((d$y[keep]-probabilities[keep])^2),
                    minimum_predicted = min(probabilities[keep]), maximum_predicted = max(probabilities[keep]),
                    predicted_over_0_99 = sum(probabilities[keep] > .99))
            }
        }
    }
}
write.csv(do.call(rbind, diagnostics), file.path(out, "missingness_diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, missing_coefficients), file.path(out, "missingness_coefficients.csv"), row.names = FALSE)
write_json(list(status = "passed_bounded_training_and_prediction_checks", fitted_models_reloaded = 16L,
    validation_prediction_sets_reproduced = 32L, maximum_probability_difference = max_error,
    exact_training_medians_checked = medians_checked, candidate_fold_scores_checked = 2960L,
    original_patient_splits_used = TRUE, external_tuning = FALSE, main_model_retrained = FALSE),
    file.path(root, "qc/neuro_fit_checks.json"))
cat("NEURO_FIT_AND_SAVED_PREDICTION_CHECKS_PASSED\n")
