# Reviewer revision: training-only preprocessing and one frozen model per window.
# Run: Rscript R/run_revision_models.R /absolute/path/to/revision_2026-11-04
# Uses only private local CSV exports. Never connects to a database.

options(stringsAsFactors = FALSE)
local_library <- Sys.getenv("STROKE_R_LIBRARY", "")
if (dir.exists(local_library)) .libPaths(unique(c(local_library, .libPaths())))
required <- c("recipes", "yardstick", "xgboost", "glmnet", "ranger", "jsonlite", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(recipes))

log_progress <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""), "\n")
    flush.console()
}

sha_file <- function(path) digest::digest(file = path, algo = "sha256")
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
    pretty = TRUE, dataframe = "rows", na = "null", digits = 16)

make_plan <- function() {
    set.seed(20261104L)
    xgb <- data.frame(candidate_id = sprintf("XGB_%02d", 1:20),
        trees = sample(c(400L, 800L, 1200L), 20, replace = TRUE),
        max_depth = sample(2:5, 20, replace = TRUE),
        eta = exp(runif(20, log(0.02), log(0.10))),
        min_child_weight = sample(c(1, 3, 5, 10, 15), 20, replace = TRUE),
        gamma = runif(20, 0, 1), subsample = runif(20, 0.65, 1),
        colsample_bytree = sample(c(0.8, 1), 20, replace = TRUE))
    rf <- expand.grid(mtry = c(6L, 12L, 24L), min_node_size = c(5L, 20L, 50L))
    rf$candidate_id <- sprintf("RF_%02d", seq_len(nrow(rf)))
    rf$trees <- 500L
    en <- expand.grid(alpha = c(0, 0.25, 0.5, 0.75, 1), lambda = 10^seq(0, -5, length.out = 25))
    en$candidate_id <- sprintf("EN_%03d", seq_len(nrow(en)))
    list(version = "revision_models_1.0", seed = 20261104L, threads = 2L,
        primary_model = "xgboost", windows = c(24L, 6L),
        training_fraction_by_patient = 0.8, grouped_folds = 5L,
        tuning_selection = "highest mean fold AUROC; ties: lower mean log loss, then candidate ID",
        final_fit = "Refit on all 80% training patients only; use same fitted model and recipe for both validation sets",
        class_imbalance = "No oversampling, undersampling, or class weighting; unweighted probability estimation",
        preprocessing = "Fold-only median imputation, unknown/novel nominal levels, dummy encoding, zero-variance and exact-linear-combination removal, normalization",
        clinical_reference = "Untuned logistic regression: age, sex, 6 vital signs and their missingness indicators; no laboratory inputs; not a validated severity score",
        pending_reference_requirement = "Reviewer 2 text after clinically meaningful reference model is incomplete; clinical-score comparator requirement remains unresolved",
        xgboost = xgb, random_forest = rf, elastic_net = en)
}

read_data <- function(path, contract, hour) {
    raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    spec <- contract$windows[[as.character(hour)]]
    keep <- !is.na(raw[[spec$eligibility]]) & raw[[spec$eligibility]] == 1 &
        !is.na(raw$outcome_missing_flag) & raw$outcome_missing_flag == 0
    raw <- raw[keep, , drop = FALSE]
    y <- as.integer(raw[[spec$outcome]])
    features <- unlist(spec$predictors, use.names = FALSE)
    forbidden <- unlist(contract$excluded_predictors, use.names = FALSE)
    stopifnot(length(features) == 38L, !any(features %in% forbidden),
        all(features %in% names(raw)), !anyNA(y), all(y %in% c(0L, 1L)),
        !anyNA(raw$patient_group_id), !anyDuplicated(raw$stay_id))
    x <- raw[, features, drop = FALSE]
    x$sex <- factor(trimws(as.character(x$sex)), levels = c("F", "M"))
    for (name in setdiff(names(x), "sex")) x[[name]] <- as.numeric(x[[name]])
    metadata <- data.frame(stay_id = raw$stay_id, group_id = raw$patient_group_id,
        hospitalid = if ("hospitalid" %in% names(raw)) raw$hospitalid else NA_integer_,
        phenotype = raw$audit_phenotype)
    list(x = x, y = y, metadata = metadata)
}

patient_partition <- function(data, seed) {
    groups <- aggregate(data$y, list(group_id = data$metadata$group_id), max)
    names(groups)[2] <- "event_stratum"
    groups <- groups[order(groups$group_id), , drop = FALSE]
    groups$partition <- "test"
    set.seed(seed)
    for (stratum in 0:1) {
        idx <- which(groups$event_stratum == stratum)
        selected <- sample(idx, floor(0.8 * length(idx)))
        groups$partition[selected] <- "train"
    }
    groups$fold <- NA_integer_
    set.seed(seed + 1L)
    for (stratum in 0:1) {
        idx <- which(groups$event_stratum == stratum & groups$partition == "train")
        groups$fold[idx] <- sample(rep(1:5, length.out = length(idx)))
    }
    mapped <- match(data$metadata$group_id, groups$group_id)
    train <- which(groups$partition[mapped] == "train")
    test <- which(groups$partition[mapped] == "test")
    stopifnot(!length(intersect(data$metadata$group_id[train], data$metadata$group_id[test])))
    list(registry = groups, train = train, test = test, fold = groups$fold[mapped][train])
}

prepare_matrix <- function(train_x, new_x = list()) {
    rec <- recipe(~ ., data = train_x) |>
        step_unknown(all_nominal_predictors(), new_level = "unknown") |>
        step_novel(all_nominal_predictors(), new_level = "novel") |>
        step_impute_median(all_numeric_predictors()) |>
        step_dummy(all_nominal_predictors(), one_hot = FALSE) |>
        step_zv(all_predictors()) |>
        step_lincomb(all_numeric_predictors()) |>
        step_normalize(all_numeric_predictors())
    rec <- prep(rec, training = train_x, retain = FALSE)
    bake_matrix <- function(x) {
        m <- as.matrix(bake(rec, new_data = x))
        storage.mode(m) <- "double"
        if (anyNA(m) || any(!is.finite(m))) stop("Non-finite values after preprocessing")
        m
    }
    list(recipe = rec, train = bake_matrix(train_x), new = lapply(new_x, bake_matrix))
}

cv_scores <- function(y, p) {
    if (anyNA(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid probabilities")
    truth <- factor(y, levels = c(0, 1), labels = c("no", "yes"))
    pp <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    c(auroc = yardstick::roc_auc_vec(truth, p, event_level = "second"),
      auprc = yardstick::pr_auc_vec(truth, p, event_level = "second"),
      log_loss = -mean(y * log(pp) + (1-y) * log(1-pp)))
}

clinical_columns <- function(x) {
    vital <- "^(heart_rate|respiratory_rate|spo2|sbp|dbp|map)_"
    which(colnames(x) == "age_years" | grepl("^sex_", colnames(x)) | grepl(vital, colnames(x)))
}

fit_clinical <- function(x, y) {
    columns <- clinical_columns(x)
    fit <- glm.fit(cbind(`(Intercept)` = 1, x[, columns, drop = FALSE]), y,
        family = binomial(), control = glm.control(maxit = 100))
    if (!isTRUE(fit$converged) || anyNA(fit$coefficients)) stop("Clinical reference failed to converge/full rank")
    list(coefficients = fit$coefficients, columns = colnames(x)[columns])
}

predict_clinical <- function(model, x) as.numeric(plogis(
    cbind(1, x[, model$columns, drop = FALSE]) %*% model$coefficients))

xgb_parameters <- function(row, seed, threads) list(objective = "binary:logistic",
    eval_metric = "logloss", tree_method = "hist", nthread = threads, seed = seed,
    max_depth = as.integer(row$max_depth), eta = row$eta,
    min_child_weight = row$min_child_weight, gamma = row$gamma,
    subsample = row$subsample, colsample_bytree = row$colsample_bytree,
    scale_pos_weight = 1, verbosity = 0)

fit_forest <- function(row, x, y, seed, threads) ranger::ranger(
    x = as.data.frame(x), y = factor(y, levels = c(0, 1), labels = c("no", "yes")),
    probability = TRUE, num.trees = as.integer(row$trees), mtry = min(as.integer(row$mtry), ncol(x)),
    min.node.size = as.integer(row$min_node_size), splitrule = "gini",
    num.threads = threads, seed = seed, importance = "none", oob.error = FALSE)

fit_enet_path <- function(alpha, x, y, plan) glmnet::glmnet(x, y, family = "binomial",
    alpha = alpha, lambda = sort(unique(plan$elastic_net$lambda), decreasing = TRUE),
    standardize = FALSE, intercept = TRUE, control = list(thresh = 1e-7, maxit = 100000L))

run_fold <- function(x, y, training_idx, validation_idx, plan, seed) {
    transformed <- prepare_matrix(x[training_idx, , drop = FALSE], list(x[validation_idx, , drop = FALSE]))
    xt <- transformed$train; xv <- transformed$new[[1L]]
    yt <- y[training_idx]; yv <- y[validation_idx]
    dtrain <- xgboost::xgb.DMatrix(xt, label = yt, nthread = plan$threads)
    dvalid <- xgboost::xgb.DMatrix(xv, nthread = plan$threads)
    pred <- list(xgboost = matrix(NA_real_, nrow(xv), nrow(plan$xgboost)),
                 random_forest = matrix(NA_real_, nrow(xv), nrow(plan$random_forest)),
                 elastic_net = matrix(NA_real_, nrow(xv), nrow(plan$elastic_net)))
    for (j in seq_len(nrow(plan$xgboost))) {
        row <- plan$xgboost[j, ]
        model <- xgboost::xgb.train(xgb_parameters(row, seed+j, plan$threads), dtrain,
            nrounds = as.integer(row$trees), verbose = 0)
        pred$xgboost[, j] <- as.numeric(predict(model, dvalid))
        if (j %% 5 == 0) log_progress("  XGBoost candidates ", j, "/", nrow(plan$xgboost))
    }
    for (j in seq_len(nrow(plan$random_forest))) {
        model <- fit_forest(plan$random_forest[j, ], xt, yt, seed+1000L+j, plan$threads)
        pred$random_forest[, j] <- predict(model, as.data.frame(xv), num.threads = plan$threads)$predictions[, "yes"]
    }
    for (alpha in unique(plan$elastic_net$alpha)) {
        idx <- which(plan$elastic_net$alpha == alpha)
        model <- fit_enet_path(alpha, xt, yt, plan)
        pred$elastic_net[, idx] <- as.matrix(predict(model, newx = xv,
            s = plan$elastic_net$lambda[idx], type = "response"))
    }
    pred$clinical_reference <- matrix(predict_clinical(fit_clinical(xt, yt), xv), ncol = 1L)
    pred$constant_training_prevalence <- matrix(mean(yt), nrow(xv), 1L)
    rows <- list()
    for (name in names(pred)) {
        ids <- if (name %in% c("xgboost", "random_forest", "elastic_net")) plan[[name]]$candidate_id else name
        for (j in seq_len(ncol(pred[[name]]))) {
            score <- cv_scores(yv, pred[[name]][, j])
            rows[[length(rows)+1L]] <- data.frame(model = name, candidate_id = ids[j],
                auroc = score[["auroc"]], auprc = score[["auprc"]], log_loss = score[["log_loss"]],
                records = length(yv), events = sum(yv))
        }
    }
    list(validation_idx = validation_idx, predictions = pred, scores = do.call(rbind, rows),
         encoded_predictors = colnames(xt))
}

point_scores <- function(y, p) {
    s <- cv_scores(y, p)
    pp <- pmin(pmax(p, 1e-6), 1-1e-6); lp <- qlogis(pp)
    citl <- unname(coef(glm(y ~ 1, offset = lp, family = binomial()))[1L])
    slope <- if (length(unique(p)) > 1L) unname(coef(glm(y ~ lp, family = binomial()))[2L]) else NA_real_
    brier <- mean((y-p)^2); prevalence <- mean(y); null_brier <- prevalence*(1-prevalence)
    c(s, brier_score = brier, calibration_intercept = citl, calibration_slope = slope,
      event_prevalence = prevalence, descriptive_prevalence_brier = null_brier,
      scaled_brier_vs_descriptive_prevalence = 1-brier/null_brier)
}

run_window <- function(root, hour, contract, plan, run_hash) {
    tag <- paste0(hour, "h")
    completion <- file.path(root, "model_results", paste0("complete_", tag, ".json"))
    if (file.exists(completion)) {
        old <- jsonlite::read_json(completion, simplifyVector = TRUE)
        if (!identical(old$run_hash, run_hash)) stop("Completed model belongs to another configuration")
        log_progress(tag, " completed artifacts verified; not refitting")
        return(invisible(NULL))
    }
    mimic <- read_data(file.path(root, "private_data", "mimic_analysis_revision.csv"), contract, hour)
    split <- patient_partition(mimic, plan$seed+hour)
    write.csv(split$registry, file.path(root, "private_splits", paste0("patients_", tag, ".csv")), row.names = FALSE)
    x <- mimic$x[split$train, , drop = FALSE]; y <- mimic$y[split$train]
    log_progress(tag, " train: ", length(y), " records, ", sum(y), " deaths; holdout: ", length(split$test))
    folds <- vector("list", 5L)
    for (k in 1:5) {
        cache <- file.path(root, "private_cache", paste0("cv_", tag, "_fold", k, ".rds"))
        if (file.exists(cache)) {
            item <- readRDS(cache)
            if (!identical(item$run_hash, run_hash)) stop("CV cache configuration mismatch")
            folds[[k]] <- item$fold
            log_progress(tag, " using completed fold ", k)
        } else {
            validation_idx <- which(split$fold == k); training_idx <- which(split$fold != k)
            stopifnot(!length(intersect(mimic$metadata$group_id[split$train[training_idx]],
                                       mimic$metadata$group_id[split$train[validation_idx]])))
            log_progress(tag, " CV fold ", k, "/5")
            folds[[k]] <- run_fold(x, y, training_idx, validation_idx, plan, plan$seed+hour*100L+k*10000L)
            saveRDS(list(run_hash = run_hash, fold = folds[[k]]), cache)
        }
        folds[[k]]$scores$fold <- k
        gc(verbose = FALSE)
    }
    fold_scores <- do.call(rbind, lapply(folds, `[[`, "scores"))
    write.csv(fold_scores, file.path(root, "model_results", paste0("cv_fold_scores_", tag, ".csv")), row.names = FALSE)
    means <- aggregate(cbind(auroc, auprc, log_loss) ~ model+candidate_id, fold_scores, mean)
    write.csv(means, file.path(root, "model_results", paste0("cv_mean_scores_", tag, ".csv")), row.names = FALSE)
    selected <- list()
    for (name in c("xgboost", "random_forest", "elastic_net")) {
        candidates <- means[means$model == name, ]
        candidates <- candidates[order(-candidates$auroc, candidates$log_loss, candidates$candidate_id), ]
        selected[[name]] <- plan[[name]][match(candidates$candidate_id[1L], plan[[name]]$candidate_id), , drop = FALSE]
    }
    write_json(selected, file.path(root, "model_results", paste0("selected_parameters_", tag, ".json")))
    log_progress(tag, " training CV selection complete; freezing final training-partition fits")
    external <- read_data(file.path(root, "private_data", "eicu_analysis_revision.csv"), contract, hour)
    final <- prepare_matrix(x, list(mimic$x[split$test, , drop = FALSE], external$x))
    saveRDS(final$recipe, file.path(root, "private_models", paste0("recipe_", tag, ".rds")))
    write.csv(data.frame(encoded_predictor = colnames(final$train)),
        file.path(root, "model_results", paste0("encoded_predictors_", tag, ".csv")), row.names = FALSE)
    xt <- final$train; dtrain <- xgboost::xgb.DMatrix(xt, label = y, nthread = plan$threads)
    models <- list()
    models$xgboost <- xgboost::xgb.train(xgb_parameters(selected$xgboost, plan$seed+hour, plan$threads),
        dtrain, nrounds = as.integer(selected$xgboost$trees), verbose = 0)
    xgb_path <- file.path(root, "private_models", paste0("xgboost_", tag, ".ubj"))
    xgboost::xgb.save(models$xgboost, xgb_path)
    probe <- final$new[[1L]][seq_len(min(10L, nrow(final$new[[1L]]))), , drop = FALSE]
    before <- predict(models$xgboost, xgboost::xgb.DMatrix(probe))
    models$xgboost <- xgboost::xgb.load(xgb_path)
    stopifnot(max(abs(before-predict(models$xgboost, xgboost::xgb.DMatrix(probe)))) < 1e-7)
    models$random_forest <- fit_forest(selected$random_forest, xt, y, plan$seed+hour, plan$threads)
    models$elastic_net <- fit_enet_path(selected$elastic_net$alpha, xt, y, plan)
    models$clinical_reference <- fit_clinical(xt, y)
    models$constant_training_prevalence <- mean(y)
    for (name in setdiff(names(models), "xgboost"))
        saveRDS(models[[name]], file.path(root, "private_models", paste0(name, "_", tag, ".rds")))
    datasets <- list(MIMIC_internal_test = list(y = mimic$y[split$test],
        meta = mimic$metadata[split$test, ], x = final$new[[1L]]),
        eICU_external_validation = list(y = external$y, meta = external$metadata, x = final$new[[2L]]))
    results <- list(); predictions <- list()
    for (dataset in names(datasets)) {
        item <- datasets[[dataset]]
        probabilities <- list(
            xgboost = as.numeric(predict(models$xgboost, xgboost::xgb.DMatrix(item$x))),
            random_forest = predict(models$random_forest, as.data.frame(item$x), num.threads = plan$threads)$predictions[, "yes"],
            elastic_net = as.numeric(predict(models$elastic_net, newx = item$x, s = selected$elastic_net$lambda, type = "response")),
            clinical_reference = predict_clinical(models$clinical_reference, item$x),
            constant_training_prevalence = rep(models$constant_training_prevalence, length(item$y)))
        for (name in names(probabilities)) {
            p <- probabilities[[name]]
            score <- point_scores(item$y, p)
            results[[length(results)+1L]] <- data.frame(window_hours = hour, model = name,
                dataset = dataset, records = length(p), patients = length(unique(item$meta$group_id)),
                events = sum(item$y), as.list(score), check.names = FALSE)
            predictions[[length(predictions)+1L]] <- data.frame(model = name, dataset = dataset,
                window_hours = hour, item$meta, truth = ifelse(item$y == 1, "yes", "no"), probability = p)
        }
    }
    write.csv(do.call(rbind, results), file.path(root, "model_results", paste0("performance_points_", tag, ".csv")), row.names = FALSE)
    write.csv(do.call(rbind, predictions), file.path(root, "private_predictions", paste0("predictions_", tag, ".csv")), row.names = FALSE)
    selected_oof <- list()
    for (name in names(models)) {
        id <- if (name %in% names(selected)) selected[[name]]$candidate_id else name
        j <- if (name %in% names(selected)) match(id, plan[[name]]$candidate_id) else 1L
        p <- rep(NA_real_, length(y))
        for (k in 1:5) p[folds[[k]]$validation_idx] <- folds[[k]]$predictions[[name]][, j]
        stopifnot(!anyNA(p))
        selected_oof[[name]] <- data.frame(model = name, mimic$metadata[split$train, ],
            truth = ifelse(y == 1, "yes", "no"), fold = split$fold, probability = p)
    }
    write.csv(do.call(rbind, selected_oof), file.path(root, "private_predictions", paste0("selected_oof_", tag, ".csv")), row.names = FALSE)
    counts <- data.frame(window_hours = hour, partition = c("training", "internal_test", "external_validation"),
        records = c(length(y), length(split$test), length(external$y)),
        patients = c(length(unique(mimic$metadata$group_id[split$train])),
            length(unique(mimic$metadata$group_id[split$test])), length(unique(external$metadata$group_id))),
        events = c(sum(y), sum(mimic$y[split$test]), sum(external$y)))
    write.csv(counts, file.path(root, "model_results", paste0("partition_counts_", tag, ".csv")), row.names = FALSE)
    write_json(list(run_hash = run_hash, window_hours = hour,
        status = "training_and_point_evaluation_complete_not_submission_ready",
        checks = list(patient_partition_overlap = 0, folds_patient_overlap = 0,
            shared_recipe_internal_external = TRUE, shared_fitted_model_internal_external = TRUE,
            xgboost_save_reload_probability_check = TRUE),
        pending = c("confidence intervals and paired comparisons", "hospital clustering",
            "clinical score comparator clarification", "calibration/DCA/SHAP and sensitivity analyses")), completion)
    log_progress(tag, " model fitting and point evaluation complete")
}

run_revision <- function(root) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    for (dir in c("model_results", "private_cache", "private_models", "private_predictions", "private_splits"))
        dir.create(file.path(root, dir), showWarnings = FALSE)
    contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
    plan <- make_plan()
    script <- file.path(root, "R", "run_revision_models.R")
    fingerprints <- list(script_sha256 = sha_file(script), contract_sha256 = sha_file(file.path(root, "model_feature_contract.json")))
    for (source in c("mimic", "eicu")) {
        path <- file.path(root, "private_data", paste0(source, "_analysis_revision.csv"))
        saved <- jsonlite::read_json(file.path(root, "qc", paste0(source, "_extract_summary.json")), simplifyVector = TRUE)
        fingerprints[[paste0(source, "_sha256")]] <- sha_file(path)
        if (!identical(saved$output_sha256, fingerprints[[paste0(source, "_sha256")]])) stop("Extraction input hash mismatch")
    }
    run_hash <- digest::digest(list(plan, fingerprints), algo = "sha256")
    plan_path <- file.path(root, "model_results", "training_plan.json")
    if (file.exists(plan_path)) {
        previous <- jsonlite::read_json(plan_path, simplifyVector = FALSE)
        if (!identical(previous$run_hash, run_hash)) stop("Existing run has a different plan/input/source hash")
    } else write_json(list(run_hash = run_hash, fingerprints = fingerprints, plan = plan), plan_path)
    for (name in c("xgboost", "random_forest", "elastic_net"))
        write.csv(plan[[name]], file.path(root, "model_results", paste0("candidate_grid_", name, ".csv")), row.names = FALSE)
    packages <- unique(c(required, "rmarkdown", "knitr"))
    write.csv(data.frame(package = packages, version = vapply(packages, function(x) as.character(packageVersion(x)), character(1))),
        file.path(root, "model_results", "package_versions.csv"), row.names = FALSE)
    capture.output(sessionInfo(), file = file.path(root, "model_results", "sessionInfo.txt"))
    for (hour in plan$windows) run_window(root, hour, contract, plan, run_hash)
    points <- do.call(rbind, lapply(plan$windows, function(h) read.csv(file.path(root, "model_results", paste0("performance_points_", h, "h.csv")))))
    write.csv(points, file.path(root, "model_results", "performance_points_all.csv"), row.names = FALSE)
    log_progress("Both prediction windows completed; uncertainty analysis remains a separate step")
    invisible(points)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    root <- if (length(args)) args[1L] else file.path(getwd(), "revision_2026-11-04")
    run_revision(root)
}
