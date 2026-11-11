# Supplementary clinical-variable reference, fixed before its fit/evaluation.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R", "run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "clinical_reference_results")
dir.create(out, showWarnings = FALSE)
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
grid <- 10^seq(2, -4, length.out = 37)
plan <- list(version = "clinical_ridge_1.0", alpha = 0, lambda = grid,
    selection = "minimum pooled patient-grouped training OOF log loss; ties larger lambda",
    role = "revision supplementary clinical-variable comparator, not validated severity score",
    predictors = "age sex and six vital signs with missingness; frozen training-only preprocessing",
    main_model_refit = FALSE, class_weighting = FALSE, recalibration = FALSE,
    script_sha256 = sha_file(file.path(root, "R/run_clinical_ridge.R")),
    rules_sha256 = sha_file(file.path(root, "08_第五步临床对照规则.md")),
    original_model_manifest_sha256 = sha_file(file.path(root, "qc/runtime_model_files.json")))
run_hash <- digest::digest(plan, algo = "sha256")
planfile <- file.path(out, "clinical_ridge_plan.json")
if (file.exists(planfile)) stop("Preserve existing clinical-ridge results; do not overwrite")
write_json(list(run_hash = run_hash, plan = plan), planfile)
all_scores <- all_predictions <- diagnostics <- list()
for (hour in c(24L, 6L)) {
    tag <- paste0(hour, "h")
    mimic <- read_data(file.path(root, "private_data/mimic_analysis_revision.csv"), contract, hour)
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", tag, ".csv")))
    mapped <- match(mimic$metadata$group_id, registry$group_id)
    stopifnot(!anyNA(mapped))
    train <- which(registry$partition[mapped] == "train")
    test <- which(registry$partition[mapped] == "test")
    folds <- registry$fold[mapped][train]
    oof <- matrix(NA_real_, length(train), length(grid))
    fit_path <- function(x, y) glmnet::glmnet(x, y, family = "binomial", alpha = 0,
        lambda = grid, standardize = FALSE, intercept = TRUE, thresh = 1e-9, maxit = 100000L)
    for (k in 1:5) {
        a <- which(folds != k); b <- which(folds == k)
        stopifnot(!length(intersect(mimic$metadata$group_id[train[a]], mimic$metadata$group_id[train[b]])))
        prep <- prepare_matrix(mimic$x[train[a], , drop = FALSE], list(mimic$x[train[b], , drop = FALSE]))
        cols <- clinical_columns(prep$train)
        fit <- fit_path(prep$train[, cols, drop = FALSE], mimic$y[train[a]])
        oof[b, ] <- as.matrix(predict(fit, newx = prep$new[[1L]][, cols, drop = FALSE], s = grid, type = "response"))
    }
    stopifnot(!anyNA(oof), all(oof > 0 & oof < 1))
    scores <- do.call(rbind, lapply(seq_along(grid), function(j) data.frame(lambda = grid[j],
        as.list(cv_scores(mimic$y[train], oof[, j])))))
    selected <- scores[order(scores$log_loss, -scores$lambda), ][1L, ]
    write.csv(scores, file.path(out, paste0("training_oof_grid_", tag, ".csv")), row.names = FALSE)
    write_json(as.list(selected), file.path(out, paste0("selected_lambda_", tag, ".json")))
    recipe <- readRDS(file.path(root, "private_models", paste0("recipe_", tag, ".rds")))
    transform <- function(x) {
        m <- as.matrix(recipes::bake(recipe, new_data = x)); storage.mode(m) <- "double"; m
    }
    xt <- transform(mimic$x[train, , drop = FALSE]); cols <- clinical_columns(xt)
    model <- list(fit = fit_path(xt[, cols, drop = FALSE], mimic$y[train]),
        lambda = selected$lambda, columns = colnames(xt)[cols], run_hash = run_hash,
        recipe_sha256 = sha_file(file.path(root, "private_models", paste0("recipe_", tag, ".rds"))))
    path <- file.path(root, "private_models", paste0("clinical_ridge_", tag, ".rds"))
    saveRDS(model, path)
    model <- readRDS(path)
    write.csv(data.frame(model = "clinical_ridge", mimic$metadata[train, ], fold = folds,
        truth = ifelse(mimic$y[train] == 1, "yes", "no"), probability = oof[, match(selected$lambda, grid)]),
        file.path(root, "private_predictions", paste0("clinical_ridge_oof_", tag, ".csv")), row.names = FALSE)
    # External data are first loaded only after tuning and final fit are saved.
    external <- read_data(file.path(root, "private_data/eicu_analysis_revision.csv"), contract, hour)
    datasets <- list(MIMIC_internal_test = list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ]),
        eICU_external_validation = list(x = external$x, y = external$y, meta = external$metadata))
    predictions <- list()
    for (dataset in names(datasets)) {
        d <- datasets[[dataset]]; xv <- transform(d$x)[, model$columns, drop = FALSE]
        p <- as.numeric(predict(model$fit, newx = xv, s = model$lambda, type = "response"))
        stopifnot(all(is.finite(p)), all(p > 0 & p < 1))
        score <- point_scores(d$y, p)
        all_scores[[length(all_scores)+1L]] <- data.frame(window_hours = hour, model = "clinical_ridge", dataset = dataset,
            records = length(p), patients = length(unique(d$meta$group_id)), events = sum(d$y), as.list(score))
        predictions[[dataset]] <- data.frame(model = "clinical_ridge", dataset = dataset, window_hours = hour,
            d$meta, truth = ifelse(d$y == 1, "yes", "no"), probability = p)
        diagnostics[[paste0(tag, "_", dataset)]] <- list(minimum_probability = min(p), maximum_probability = max(p),
            predictions_below_1e_minus6 = sum(p < 1e-6), predictions_above_1_minus_1e_minus6 = sum(p > 1-1e-6),
            finite_coefficients = all(is.finite(as.matrix(coef(model$fit, s = model$lambda)))))
    }
    pred <- do.call(rbind, predictions)
    write.csv(pred, file.path(root, "private_predictions", paste0("clinical_ridge_predictions_", tag, ".csv")), row.names = FALSE)
    all_predictions[[tag]] <- pred
    log_progress("Clinical ridge ", tag, " complete; lambda=", signif(selected$lambda, 5))
}
write.csv(do.call(rbind, all_scores), file.path(out, "performance_points.csv"), row.names = FALSE)
write_json(diagnostics, file.path(out, "probability_diagnostics.json"))
write_json(list(status = "fitted_and_point_evaluated_not_final_reference_requirement_resolved", run_hash = run_hash,
    main_model_changed = FALSE, validated_score = FALSE), file.path(out, "fit_complete.json"))
