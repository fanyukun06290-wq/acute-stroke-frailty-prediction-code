# Reviewer 2: same-learner neurological reference and laboratory-domain increments.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "neuro_reference_results")
dir.create(out, showWarnings = FALSE)
neuro_variables <- c(motor = "gcs_motor", total = "gcs_total")
domains <- c("B", "BM", "BI", "BMI")
lambda_grid <- 10^seq(2, -4, length.out = 37)
neuro_columns <- function(x, neuro, domain) {
    variables <- c("age_years", "sex", grep("^(heart_rate|respiratory_rate|spo2|sbp|dbp|map)_", names(x), value = TRUE),
        paste0(neuro_variables[[neuro]], c("_min", "_missing")))
    if (domain %in% c("BM", "BMI")) variables <- c(variables,
        grep("^(glucose|sodium|potassium|chloride|bicarbonate|bun|creatinine|calcium)_", names(x), value = TRUE))
    if (domain %in% c("BI", "BMI")) variables <- c(variables,
        grep("^(wbc|platelets|hemoglobin|hematocrit)_", names(x), value = TRUE))
    stopifnot(!anyDuplicated(variables), all(variables %in% names(x)),
        length(variables) == switch(domain, B = 16L, BM = 32L, BI = 24L, BMI = 40L))
    variables
}
load_neuro <- function(source_name, hour, contract) {
    d <- read_data(file.path(root, "private_data", paste0(source_name, "_analysis_revision.csv")), contract, hour)
    g <- read.csv(file.path(root, "private_data", paste0(source_name, "_gcs_features.csv")))
    g <- g[g$window_hours == hour, ]
    stopifnot(!anyDuplicated(g$stay_id), setequal(g$stay_id, d$metadata$stay_id))
    g <- g[match(d$metadata$stay_id, g$stay_id), ]
    for (variable in unname(neuro_variables)) {
        value <- g[[paste0(variable, "_min")]]
        missing <- g[[paste0(variable, "_missing")]]
        stopifnot(identical(as.integer(is.na(value)), as.integer(missing)))
        d$x[[paste0(variable, "_min")]] <- as.numeric(value)
        d$x[[paste0(variable, "_missing")]] <- as.numeric(missing)
    }
    d
}
fit_neuro_path <- function(x, y) glmnet::glmnet(x, y, family = "binomial", alpha = 0,
    lambda = lambda_grid, standardize = FALSE, intercept = TRUE, control = list(thresh = 1e-9, maxit = 100000L))
run_neuro_fit <- function() {
    contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
    inputs <- c("model_feature_contract.json", "qc/runtime_gcs_files.json", "qc/runtime_model_files.json",
        "private_data/mimic_analysis_revision.csv", "private_data/eicu_analysis_revision.csv",
        "private_data/mimic_gcs_features.csv", "private_data/eicu_gcs_features.csv",
        "private_splits/patients_24h.csv", "private_splits/patients_6h.csv")
    plan <- list(version = "nested_neuro_ridge_1.0", alpha = 0, lambda_grid = lambda_grid,
        neuro_types = neuro_variables, domains = domains,
        selection = "pooled grouped-training OOF log loss minimum; ties larger lambda",
        preprocessing = "separate training-only recipe per feature set per fold; final refit on original 80% training partition",
        rule_sha256 = sha_file(file.path(root, "12_第七步神经参考与领域增量规则.md")),
        script_sha256 = sha_file(file.path(root, "R/fit_neuro_nested.R")),
        inputs_sha256 = setNames(lapply(inputs, function(p) sha_file(file.path(root, p))), inputs),
        external_tuning = FALSE, main_model_changed = FALSE, supplementary_revision_analysis = TRUE)
    hash <- digest::digest(plan, algo = "sha256")
    planfile <- file.path(out, "training_plan.json")
    if (file.exists(planfile)) stopifnot(identical(jsonlite::read_json(planfile)$run_hash, hash)) else
        write_json(list(run_hash = hash, plan = plan), planfile)
    points <- all_parameters <- predictor_lists <- counts <- list()
    for (hour in c(24L, 6L)) {
        mimic <- load_neuro("mimic", hour, contract)
        registry <- read.csv(file.path(root, "private_splits", paste0("patients_", hour, "h.csv")))
        idx <- match(mimic$metadata$group_id, registry$group_id)
        stopifnot(!anyNA(idx))
        train <- which(registry$partition[idx] == "train"); test <- which(registry$partition[idx] == "test")
        folds <- registry$fold[idx][train]
        stopifnot(!length(intersect(mimic$metadata$group_id[train], mimic$metadata$group_id[test])))
        new_predictions <- list()
        for (neuro in names(neuro_variables)) for (domain in domains) {
            name <- paste(neuro, domain, sep = "_"); tag <- paste0(name, "_", hour, "h")
            cols <- neuro_columns(mimic$x, neuro, domain)
            predictor_lists[[tag]] <- cols
            modelpath <- file.path(root, "private_models", paste0("neuro_", tag, ".rds"))
            cachepath <- file.path(root, "private_cache", paste0("neuro_fit_", tag, ".rds"))
            if (file.exists(cachepath)) {
                fit_cache <- readRDS(cachepath); stopifnot(identical(fit_cache$run_hash, hash))
                model <- readRDS(modelpath)
                selected <- fit_cache$selected
            } else {
                oof <- matrix(NA_real_, length(train), length(lambda_grid)); fold_scores <- list()
                for (k in 1:5) {
                    a <- which(folds != k); b <- which(folds == k)
                    stopifnot(!length(intersect(mimic$metadata$group_id[train[a]], mimic$metadata$group_id[train[b]])))
                    prepared <- prepare_matrix(mimic$x[train[a], cols, drop = FALSE], list(mimic$x[train[b], cols, drop = FALSE]))
                    fit <- fit_neuro_path(prepared$train, mimic$y[train[a]])
                    p <- as.matrix(predict(fit, newx = prepared$new[[1L]], s = lambda_grid, type = "response"))
                    oof[b, ] <- p
                    fold_scores[[k]] <- do.call(rbind, lapply(seq_along(lambda_grid), function(j)
                        data.frame(fold = k, lambda = lambda_grid[j], as.list(cv_scores(mimic$y[train[b]], p[, j])))))
                }
                stopifnot(!anyNA(oof), all(is.finite(oof)), all(oof >= 0 & oof <= 1))
                pooled <- do.call(rbind, lapply(seq_along(lambda_grid), function(j)
                    data.frame(lambda = lambda_grid[j], as.list(cv_scores(mimic$y[train], oof[, j])))))
                selected <- pooled[order(pooled$log_loss, -pooled$lambda), ][1L, ]
                write.csv(pooled, file.path(out, paste0("oof_grid_", tag, ".csv")), row.names = FALSE)
                write.csv(do.call(rbind, fold_scores), file.path(out, paste0("fold_grid_", tag, ".csv")), row.names = FALSE)
                prepared <- prepare_matrix(mimic$x[train, cols, drop = FALSE])
                model <- list(run_hash = hash, recipe = prepared$recipe, raw_columns = cols,
                    encoded_columns = colnames(prepared$train), fit = fit_neuro_path(prepared$train, mimic$y[train]),
                    lambda = selected$lambda, hour = hour, neuro = neuro, domain = domain)
                saveRDS(model, modelpath)
                model <- readRDS(modelpath)
                write.csv(data.frame(model = name, mimic$metadata[train, ], fold = folds,
                    truth = ifelse(mimic$y[train] == 1, "yes", "no"), probability = oof[, match(selected$lambda, lambda_grid)]),
                    file.path(root, "private_predictions", paste0("neuro_oof_", tag, ".csv")), row.names = FALSE)
                saveRDS(list(run_hash = hash, selected = selected), cachepath)
            }
            all_parameters[[tag]] <- data.frame(window_hours = hour, model = name, neuro = neuro, domain = domain, selected)
            # Selection and training fit are already saved before validation predictions.
            external <- load_neuro("eicu", hour, contract)
            sets <- list(MIMIC_internal_test = list(x = mimic$x[test, ], y = mimic$y[test], meta = mimic$metadata[test, ]),
                eICU_external_validation = list(x = external$x, y = external$y, meta = external$metadata))
            for (dataset in names(sets)) {
                d <- sets[[dataset]]
                x <- as.matrix(recipes::bake(model$recipe, new_data = d$x[, model$raw_columns, drop = FALSE]))
                stopifnot(identical(colnames(x), model$encoded_columns), !anyNA(x))
                p <- as.numeric(predict(model$fit, newx = x, s = model$lambda, type = "response"))
                stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1))
                observed <- d$x[[paste0(neuro_variables[[neuro]], "_missing")]] == 0
                new_predictions[[paste(name, dataset)]] <- data.frame(model = name, neuro = neuro, domain = domain,
                    dataset = dataset, window_hours = hour, d$meta, neuro_observed = observed,
                    truth = ifelse(d$y == 1, "yes", "no"), probability = p)
                populations <- list(all_eligible = seq_along(p))
                if (dataset == "eICU_external_validation") populations$neuro_observed <- which(observed)
                for (population in names(populations)) {
                    keep <- populations[[population]]
                    points[[length(points)+1L]] <- data.frame(window_hours = hour, model = name, neuro = neuro, domain = domain,
                        dataset = dataset, population = population, records = length(keep),
                        patients = length(unique(d$meta$group_id[keep])), events = sum(d$y[keep]),
                        as.list(point_scores(d$y[keep], p[keep])))
                }
            }
            for (part in c("train", "test")) {
                keep <- if (part == "train") train else test
                counts[[paste(tag, part)]] <- data.frame(window_hours = hour, model = name, partition = part,
                    records = length(keep), observed_neuro = sum(mimic$x[[paste0(neuro_variables[[neuro]], "_missing")]][keep] == 0))
            }
            log_progress("Nested neurological model ", tag, " fit/evaluated; lambda=", signif(selected$lambda, 5))
        }
        write.csv(do.call(rbind, new_predictions), file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv")), row.names = FALSE)
    }
    write.csv(do.call(rbind, points), file.path(out, "performance_points.csv"), row.names = FALSE)
    write.csv(do.call(rbind, all_parameters), file.path(out, "selected_parameters.csv"), row.names = FALSE)
    write.csv(do.call(rbind, counts), file.path(out, "training_neuro_availability.csv"), row.names = FALSE)
    write_json(predictor_lists, file.path(out, "raw_predictor_sets.json"))
    write_json(list(status = "nested_model_fits_and_points_complete_not_final_paper_ready", run_hash = hash,
        fitted_models = 16L, main_model_changed = FALSE, external_tuning = FALSE), file.path(out, "fit_complete.json"))
}
if (sys.nframe() == 0L) run_neuro_fit()
