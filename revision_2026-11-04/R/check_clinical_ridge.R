# Bounded numerical verification of supplementary clinical ridge artifacts.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
root <- file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R/run_revision_uncertainty.R"), encoding = "UTF-8")
out <- file.path(root, "clinical_reference_results")
plan <- jsonlite::read_json(file.path(out, "clinical_ridge_plan.json"))
vp <- jsonlite::read_json(file.path(out, "validation_plan.json"))
stopifnot(identical(sha_file(file.path(root, "R/run_clinical_ridge.R")), plan$plan$script_sha256),
    identical(sha_file(file.path(root, "08_第五步临床对照规则.md")), plan$plan$rules_sha256),
    identical(sha_file(file.path(root, "R/validate_clinical_ridge.R")), vp$plan$script_sha256))
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
points <- read.csv(file.path(out, "performance_points.csv"))
maximum_prediction_difference <- 0
for (hour in c(24L, 6L)) {
    tag <- paste0(hour, "h")
    grid <- read.csv(file.path(out, paste0("training_oof_grid_", tag, ".csv")))
    selected <- jsonlite::read_json(file.path(out, paste0("selected_lambda_", tag, ".json")))
    stopifnot(nrow(grid) == 37L, abs(selected$lambda - grid$lambda[order(grid$log_loss, -grid$lambda)[1L]]) < 1e-14)
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", tag, ".csv")))
    stopifnot(!anyDuplicated(registry$group_id))
    oof <- read.csv(file.path(root, "private_predictions", paste0("clinical_ridge_oof_", tag, ".csv")))
    mapped <- match(oof$group_id, registry$group_id)
    stopifnot(!anyNA(mapped), all(registry$partition[mapped] == "train"), identical(oof$fold, registry$fold[mapped]),
        abs(cv_scores(as.integer(oof$truth == "yes"), oof$probability)["log_loss"] - selected$log_loss) < 1e-10)
    model <- readRDS(file.path(root, "private_models", paste0("clinical_ridge_", tag, ".rds")))
    stopifnot(identical(model$run_hash, plan$run_hash), abs(model$lambda-selected$lambda) < 1e-14,
        all(grepl("^(age_years|sex_|heart_rate_|respiratory_rate_|spo2_|sbp_|dbp_|map_)", model$columns)),
        all(is.finite(as.matrix(coef(model$fit, s = model$lambda)))))
    recpath <- file.path(root, "private_models", paste0("recipe_", tag, ".rds"))
    stopifnot(identical(sha_file(recpath), model$recipe_sha256))
    recipe <- readRDS(recpath)
    saved <- read.csv(file.path(root, "private_predictions", paste0("clinical_ridge_predictions_", tag, ".csv")))
    for (source_name in c("mimic", "eicu")) {
        data <- read_data(file.path(root, "private_data", paste0(source_name, "_analysis_revision.csv")), contract, hour)
        dataset <- if (source_name == "mimic") "MIMIC_internal_test" else "eICU_external_validation"
        if (source_name == "mimic") {
            indices <- which(registry$partition[match(data$metadata$group_id, registry$group_id)] == "test")
            stopifnot(!length(intersect(data$metadata$group_id[indices], oof$group_id)))
        } else indices <- seq_along(data$y)
        x <- as.matrix(recipes::bake(recipe, new_data = data$x[indices, ]))[, model$columns, drop = FALSE]
        probabilities <- as.numeric(predict(model$fit, newx = x, s = model$lambda, type = "response"))
        d <- saved[saved$dataset == dataset, ]
        d <- d[match(data$metadata$stay_id[indices], d$stay_id), ]
        stopifnot(!anyNA(d$probability), identical(d$truth, ifelse(data$y[indices] == 1, "yes", "no")))
        delta <- max(abs(probabilities - d$probability)); stopifnot(delta < 1e-12)
        maximum_prediction_difference <- max(maximum_prediction_difference, delta)
        s <- point_scores(data$y[indices], probabilities)
        r <- points[points$window_hours == hour & points$dataset == dataset, ]
        stopifnot(max(abs(unlist(r[, names(s)]) - s)) < 1e-8)
    }
}
ci <- function(v) as.numeric(quantile(v[is.finite(v)], c(.025, .975), type = 7))
metrics <- read.csv(file.path(out, "metric_intervals.csv"))
pairs <- read.csv(file.path(out, "paired_differences.csv"))
bins <- read.csv(file.path(out, "calibration_bins.csv"))
files <- list.files(file.path(root, "private_cache"), pattern = "^clinical_ridge_uncertainty_.*\\.rds$", full.names = TRUE)
stopifnot(length(files) == 8L, nrow(metrics) == 96L, nrow(pairs) == 32L, nrow(bins) == 80L)
for (f in files) {
    boot <- readRDS(f)
    stopifnot(identical(boot$run_hash, vp$run_hash), boot$completed == 1000L, boot$replay_checks == 10L)
    label <- sub("^clinical_ridge_", "", basename(f))
    old <- readRDS(file.path(root, "private_cache", label))
    select <- function(table) table[table$window_hours == old$hour & table$dataset == old$dataset & table$scheme == old$scheme, ]
    m <- select(metrics); p <- select(pairs); b <- select(bins)
    for (i in seq_len(nrow(m))) stopifnot(max(abs(ci(boot$metrics[, m$metric[i]]) - c(m$lower[i], m$upper[i]))) < 1e-10)
    for (i in seq_len(nrow(p))) {
        v <- old$metrics[, "xgboost", p$metric[i]] - boot$metrics[, p$metric[i]]
        stopifnot(max(abs(ci(v) - c(p$lower[i], p$upper[i]))) < 1e-10)
    }
    for (i in seq_len(nrow(b))) stopifnot(max(abs(ci(boot$bin_rates[, b$bin[i]]) - c(b$lower[i], b$upper[i]))) < 1e-10)
}
write_json(list(status = "passed_bounded_numerical_checks", models_reloaded = 2L,
    full_validation_prediction_sets_reproduced = 4L, maximum_prediction_difference = maximum_prediction_difference,
    candidate_count_each_window = 37L, schemes = 8L, paired_xgboost_weight_replay_checks = 80L,
    metric_intervals_verified = 96L, paired_intervals_verified = 32L, bin_intervals_verified = 80L,
    main_model_retraining = FALSE, validated_clinical_score_requirement_resolved = FALSE),
    file.path(root, "qc/clinical_ridge_checks.json"))
cat("CLINICAL_RIDGE_BOUNDED_CHECKS_PASSED\n")
