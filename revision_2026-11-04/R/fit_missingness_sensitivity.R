# Post-hoc sensitivity after external missingness-shift diagnosis; retain all originals.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/fit_neuro_nested.R"), encoding = "UTF-8")
out <- file.path(root, "missingness_sensitivity_results")
dir.create(out, showWarnings = FALSE)
strategies <- c("omit_neuro_flag", "omit_all_flags")
sensitivity_columns <- function(x, strategy, domain) {
    columns <- neuro_columns(x, "motor", domain)
    if (strategy == "omit_neuro_flag") columns <- setdiff(columns, "gcs_motor_missing") else
        columns <- columns[!grepl("_missing($|_)", columns)]
    expected <- if (strategy == "omit_neuro_flag") switch(domain, B = 15L, BM = 31L, BI = 23L, BMI = 39L) else
        switch(domain, B = 9L, BM = 17L, BI = 13L, BMI = 21L)
    stopifnot(length(columns) == expected, !"gcs_motor_missing" %in% columns)
    columns
}
run_sensitivity_fit <- function() {
    inputs <- c("model_feature_contract.json", "qc/runtime_neuro_files.json", "private_data/mimic_analysis_revision.csv",
        "private_data/eicu_analysis_revision.csv", "private_data/mimic_gcs_features.csv", "private_data/eicu_gcs_features.csv",
        "private_splits/patients_24h.csv", "private_splits/patients_6h.csv", "R/fit_neuro_nested.R")
    plan <- list(version = "posthoc_missing_flags_1.0", posthoc_external_diagnosis_triggered = TRUE,
        strategies = strategies, domains = domains, alpha = 0, lambda_grid = lambda_grid,
        selection = "same grouped-training OOF log loss rule; no new external selection",
        rules_sha256 = sha_file(file.path(root, "14_第八步缺失指示稳健性规则.md")),
        script_sha256 = sha_file(file.path(root, "R/fit_missingness_sensitivity.R")),
        inputs_sha256 = setNames(lapply(inputs, function(p) sha_file(file.path(root, p))), inputs),
        same_training_partitions = TRUE, original_models_changed = FALSE)
    hash <- digest::digest(plan, algo = "sha256")
    path <- file.path(out, "training_plan.json")
    if (file.exists(path)) stopifnot(identical(jsonlite::read_json(path)$run_hash, hash)) else
        write_json(list(run_hash = hash, plan = plan), path)
    contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
    points <- parameters <- predictors <- diagnostics <- list()
    for (hour in c(24L, 6L)) {
        mimic <- load_neuro("mimic", hour, contract)
        registry <- read.csv(file.path(root, "private_splits", paste0("patients_", hour, "h.csv")))
        mapped <- match(mimic$metadata$group_id, registry$group_id)
        stopifnot(!anyNA(mapped))
        train <- which(registry$partition[mapped] == "train"); test <- which(registry$partition[mapped] == "test")
        folds <- registry$fold[mapped][train]; predictions <- list()
        for (strategy in strategies) for (domain in domains) {
            name <- paste(strategy, domain, sep = "_"); tag <- paste0(name, "_", hour, "h")
            cols <- sensitivity_columns(mimic$x, strategy, domain); predictors[[tag]] <- cols
            modelpath <- file.path(root, "private_models", paste0("missing_sensitivity_", tag, ".rds"))
            cachepath <- file.path(root, "private_cache", paste0("missing_fit_", tag, ".rds"))
            if (file.exists(cachepath)) {
                item <- readRDS(cachepath); stopifnot(identical(item$run_hash, hash))
                selected <- item$selected; model <- readRDS(modelpath)
            } else {
                oof <- matrix(NA_real_, length(train), length(lambda_grid)); scores <- list()
                for (k in 1:5) {
                    a <- which(folds != k); b <- which(folds == k)
                    stopifnot(!length(intersect(mimic$metadata$group_id[train[a]], mimic$metadata$group_id[train[b]])))
                    prep <- prepare_matrix(mimic$x[train[a], cols, drop = FALSE], list(mimic$x[train[b], cols, drop = FALSE]))
                    fit <- fit_neuro_path(prep$train, mimic$y[train[a]])
                    p <- as.matrix(predict(fit, newx = prep$new[[1L]], s = lambda_grid, type = "response"))
                    oof[b, ] <- p
                    scores[[k]] <- do.call(rbind, lapply(seq_along(lambda_grid), function(j)
                        data.frame(fold = k, lambda = lambda_grid[j], as.list(cv_scores(mimic$y[train[b]], p[, j])))))
                }
                stopifnot(!anyNA(oof))
                grid <- do.call(rbind, lapply(seq_along(lambda_grid), function(j)
                    data.frame(lambda = lambda_grid[j], as.list(cv_scores(mimic$y[train], oof[, j])))))
                selected <- grid[order(grid$log_loss, -grid$lambda), ][1L, ]
                write.csv(grid, file.path(out, paste0("oof_grid_", tag, ".csv")), row.names = FALSE)
                write.csv(do.call(rbind, scores), file.path(out, paste0("fold_grid_", tag, ".csv")), row.names = FALSE)
                prep <- prepare_matrix(mimic$x[train, cols, drop = FALSE])
                model <- list(run_hash = hash, raw_columns = cols, recipe = prep$recipe, encoded_columns = colnames(prep$train),
                    fit = fit_neuro_path(prep$train, mimic$y[train]), lambda = selected$lambda,
                    strategy = strategy, domain = domain, hour = hour)
                saveRDS(model, modelpath); model <- readRDS(modelpath)
                write.csv(data.frame(model = name, mimic$metadata[train, ], fold = folds,
                    truth = ifelse(mimic$y[train] == 1, "yes", "no"), probability = oof[, match(selected$lambda, lambda_grid)]),
                    file.path(root, "private_predictions", paste0("missing_oof_", tag, ".csv")), row.names = FALSE)
                saveRDS(list(run_hash = hash, selected = selected), cachepath)
            }
            parameters[[tag]] <- data.frame(window_hours = hour, model = name, strategy = strategy, domain = domain, selected)
            external <- load_neuro("eicu", hour, contract)
            sets <- list(training_apparent = list(x = mimic$x[train, ], y = mimic$y[train], meta = mimic$metadata[train, ]),
                MIMIC_internal_test = list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ]),
                eICU_external_validation = list(x = external$x, y = external$y, meta = external$metadata))
            for (dataset in names(sets)) {
                d <- sets[[dataset]]
                x <- as.matrix(recipes::bake(model$recipe, new_data = d$x[, model$raw_columns, drop = FALSE]))
                stopifnot(identical(colnames(x), model$encoded_columns), !anyNA(x))
                p <- as.numeric(predict(model$fit, newx = x, s = model$lambda, type = "response"))
                stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1))
                observed <- d$x$gcs_motor_missing == 0
                for (missing in c(FALSE, TRUE)) {
                    keep <- !observed == missing
                    if (!any(keep)) next
                    diagnostics[[length(diagnostics)+1L]] <- data.frame(window_hours = hour, model = name, strategy = strategy,
                        domain = domain, dataset = dataset, neuro_missing = missing, records = sum(keep), events = sum(d$y[keep]),
                        observed_rate = mean(d$y[keep]), mean_predicted = mean(p[keep]), brier = mean((d$y[keep]-p[keep])^2),
                        minimum_predicted = min(p[keep]), maximum_predicted = max(p[keep]))
                }
                if (dataset == "training_apparent") next
                predictions[[paste(name, dataset)]] <- data.frame(model = name, strategy = strategy, domain = domain, dataset = dataset,
                    window_hours = hour, d$meta, neuro_observed = observed, truth = ifelse(d$y == 1, "yes", "no"), probability = p)
                populations <- list(all_eligible = seq_along(p))
                if (dataset == "eICU_external_validation") populations$neuro_observed <- which(observed)
                for (population in names(populations)) {
                    keep <- populations[[population]]
                    points[[length(points)+1L]] <- data.frame(window_hours = hour, model = name, strategy = strategy, domain = domain,
                        dataset = dataset, population = population, records = length(keep), patients = length(unique(d$meta$group_id[keep])),
                        events = sum(d$y[keep]), as.list(point_scores(d$y[keep], p[keep])))
                }
            }
            log_progress("Missingness sensitivity ", tag, " fit/evaluated; lambda=", signif(model$lambda, 5))
        }
        write.csv(do.call(rbind, predictions), file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv")), row.names = FALSE)
    }
    write.csv(do.call(rbind, points), file.path(out, "performance_points.csv"), row.names = FALSE)
    write.csv(do.call(rbind, parameters), file.path(out, "selected_parameters.csv"), row.names = FALSE)
    write.csv(do.call(rbind, diagnostics), file.path(out, "missingness_diagnostics.csv"), row.names = FALSE)
    write_json(predictors, file.path(out, "raw_predictor_sets.json"))
    write_json(list(status = "posthoc_sensitivity_fits_complete_not_final_validation", run_hash = hash,
        fitted_models = 16L, old_models_changed = FALSE), file.path(out, "fit_complete.json"))
}
if (sys.nframe() == 0L) run_sensitivity_fit()
