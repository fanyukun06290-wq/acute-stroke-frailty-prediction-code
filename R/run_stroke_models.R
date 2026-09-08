# Acute-stroke prognosis: R modeling pipeline
#
# Usage:
#   Rscript R/run_stroke_models.R data/mimic_analysis.csv data/eicu_analysis.csv results/model_results
#
# The input CSV files are the outputs exported from the SQL files in sql/.
# This script does not query or modify either
# database.

options(stringsAsFactors = FALSE)

# Keep the project-local library available to both RStudio and Rscript.
project_library <- Sys.getenv(
    "STROKE_R_LIBRARY",
    unset = ""
)
if (dir.exists(project_library)) {
    .libPaths(unique(c(project_library, .libPaths())))
}

# R 4.5 on Windows may install user packages outside the default library
# search path. Add the standard user library before checking dependencies.
user_library <- file.path(
    Sys.getenv("LOCALAPPDATA", unset = ""), "R", "win-library",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
)
if (nzchar(Sys.getenv("LOCALAPPDATA", unset = "")) && dir.exists(user_library)) {
    .libPaths(unique(c(user_library, .libPaths())))
}

args <- commandArgs(trailingOnly = TRUE)
env_mimic <- Sys.getenv("STROKE_MIMIC_CSV", unset = "")
env_eicu <- Sys.getenv("STROKE_EICU_CSV", unset = "")
env_out <- Sys.getenv("STROKE_OUTPUT_DIR", unset = "")
mimic_path <- if (length(args) >= 1) args[[1]] else
    if (nzchar(env_mimic)) env_mimic else "data/mimic_analysis.csv"
eicu_path <- if (length(args) >= 2) args[[2]] else
    if (nzchar(env_eicu)) env_eicu else "data/eicu_analysis.csv"
out_dir <- if (length(args) >= 3) args[[3]] else
    if (nzchar(env_out)) env_out else "results/model_results"

required_packages <- c("tidymodels", "xgboost", "glmnet", "ranger")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
    stop(
        "Missing R packages: ", paste(missing_packages, collapse = ", "),
        ". Install them before running this script; no packages are installed automatically."
    )
}

suppressPackageStartupMessages({
    library(tidymodels)
    library(xgboost)
    library(glmnet)
    library(ranger)
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260901)

tuning_grid_size <- suppressWarnings(as.integer(
    Sys.getenv("STROKE_TUNING_GRID_SIZE", unset = "5")
))
if (is.na(tuning_grid_size) || tuning_grid_size < 2) tuning_grid_size <- 20

labs <- c(
    "glucose", "sodium", "potassium", "chloride", "bicarbonate", "bun",
    "creatinine", "calcium", "wbc", "platelets", "hemoglobin", "hematocrit"
)
vitals <- c("heart_rate", "respiratory_rate", "spo2", "sbp", "dbp", "map")

to01 <- function(x) {
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(y))
    out[y %in% c("1", "true", "t", "yes")] <- 1L
    out[y %in% c("0", "false", "f", "no")] <- 0L
    out
}

collapse_stroke_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "AIS_STRICT", "AIS",
        ifelse(x == "ICH_STRICT", "ICH",
            ifelse(x == "SAH", "SAH", "MIXED_OR_AMBIGUOUS")))
}

prepare_dataset <- function(raw, window = c("landmark", "early")) {
    window <- match.arg(window)
    suffix <- if (window == "landmark") "mean_24h" else "first_6h"
    missing_suffix <- if (window == "landmark") "missing_24h" else "missing_6h"
    outcome_col <- if (window == "landmark") {
        "death_24h_to_30d_inhospital"
    } else {
        "death_30d_inhospital"
    }
    required <- c(
        "audit_phenotype", "age_years", "sex", "frailty_comorbidity_count",
        outcome_col, "outcome_missing_flag",
        paste0(labs, "_", suffix), paste0(labs, "_", missing_suffix),
        paste0(vitals, "_", suffix), paste0(vitals, "_", missing_suffix)
    )
    absent <- setdiff(required, names(raw))
    if (length(absent) > 0) {
        stop("Input is missing required columns: ", paste(absent, collapse = ", "))
    }
    outcome_value <- to01(raw[[outcome_col]])
    missing_value <- to01(raw$outcome_missing_flag)
    if (window == "landmark") {
        window_keep <- to01(raw$landmark_24h_eligible) == 1L
        window_keep[is.na(window_keep)] <- FALSE
        raw <- raw[window_keep, , drop = FALSE]
        outcome_value <- outcome_value[window_keep]
        missing_value <- missing_value[window_keep]
    }
    outcome_keep <- !is.na(outcome_value) & !is.na(missing_value) &
        missing_value == 0L
    raw <- raw[outcome_keep, , drop = FALSE]
    outcome_value <- outcome_value[outcome_keep]

    if ("subject_id" %in% names(raw)) {
        group_id <- as.character(raw$subject_id)
    } else if ("patienthealthsystemstayid" %in% names(raw)) {
        group_id <- as.character(raw$patienthealthsystemstayid)
    } else {
        stop("No patient-level grouping identifier was found in the input.")
    }

    feature_names <- c(
        "age_years", "sex", "stroke_group", "frailty_comorbidity_count",
        paste0(labs, "_", suffix), paste0(labs, "_", missing_suffix),
        paste0(vitals, "_", suffix), paste0(vitals, "_", missing_suffix)
    )
    include_missingness <- tolower(Sys.getenv(
        "STROKE_INCLUDE_MISSINGNESS", unset = "1"
    )) %in% c("1", "true", "t", "yes")
    if (!include_missingness) {
        feature_names <- setdiff(
            feature_names,
            c(paste0(labs, "_", missing_suffix),
              paste0(vitals, "_", missing_suffix))
        )
    }
    numeric_feature_names <- setdiff(
        feature_names,
        c("age_years", "sex", "stroke_group", "frailty_comorbidity_count")
    )
    numeric_features <- raw[, numeric_feature_names, drop = FALSE]
    numeric_features[] <- lapply(
        numeric_features,
        function(x) suppressWarnings(as.numeric(as.character(x)))
    )
    sex_value <- trimws(as.character(raw$sex))
    sex_value[is.na(sex_value) |
        tolower(sex_value) %in% c("", "unknown", "unk", "u", "na")] <- NA_character_
    dat <- data.frame(
        outcome = factor(outcome_value, levels = c(0, 1),
            labels = c("no", "yes")),
        split_group = group_id,
        age_years = suppressWarnings(as.numeric(as.character(raw$age_years))),
        sex = factor(sex_value),
        stroke_group = factor(collapse_stroke_group(raw$audit_phenotype),
            levels = c("AIS", "ICH", "SAH", "MIXED_OR_AMBIGUOUS")),
        frailty_comorbidity_count = suppressWarnings(as.numeric(
            as.character(raw$frailty_comorbidity_count)
        )),
        numeric_features,
        check.names = FALSE
    )
    dat
}

make_recipe <- function(train_data) {
    recipe(outcome ~ ., data = train_data) %>%
        update_role(split_group, new_role = "id") %>%
        update_role(split_stratum, new_role = "id") %>%
        step_unknown(all_nominal_predictors()) %>%
        step_novel(all_nominal_predictors()) %>%
        step_impute_median(all_numeric_predictors()) %>%
        step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
        step_zv(all_predictors()) %>%
        step_normalize(all_numeric_predictors())
}

calibration_summary <- function(truth, probability) {
    y <- as.integer(truth == "yes")
    p <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    lp <- qlogis(p)
    intercept_fit <- glm(y ~ 1, family = binomial(), offset = lp)
    slope_fit <- glm(y ~ lp, family = binomial())
    c(
        calibration_intercept = unname(coef(intercept_fit)[[1]]),
        calibration_slope = unname(coef(slope_fit)[[2]])
    )
}

score_predictions <- function(truth, probability) {
    p <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    cal <- calibration_summary(truth, p)
    tibble(
        auroc = roc_auc_vec(truth, p, event_level = "second"),
        auprc = pr_auc_vec(truth, p, event_level = "second"),
        brier_score = mean((as.integer(truth == "yes") - p)^2),
        calibration_intercept = cal[[1]],
        calibration_slope = cal[[2]]
    )
}

decision_curve <- function(truth, probability, thresholds = seq(0.05, 0.50, 0.05)) {
    y <- as.integer(truth == "yes")
    n <- length(y)
    tibble(
        threshold = thresholds,
        net_benefit = vapply(thresholds, function(t) {
            predicted_positive <- probability >= t
            tp <- sum(predicted_positive & y == 1)
            fp <- sum(predicted_positive & y == 0)
            tp / n - fp / n * t / (1 - t)
        }, numeric(1))
    )
}

tune_one <- function(model_name, recipe_obj, folds, grid_size = 20) {
    spec <- switch(model_name,
        xgboost = boost_tree(
            trees = tune(), tree_depth = tune(), learn_rate = tune(),
            min_n = tune(), loss_reduction = tune(), sample_size = tune()
        ) %>% set_engine("xgboost", objective = "binary:logistic",
            eval_metric = "logloss", verbosity = 0) %>% set_mode("classification"),
        elastic_net = logistic_reg(penalty = tune(), mixture = tune()) %>%
            set_engine("glmnet") %>% set_mode("classification"),
        random_forest = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
            set_engine("ranger", probability = TRUE, importance = "impurity") %>%
            set_mode("classification"),
        stop("Unknown model: ", model_name)
    )
    wf <- workflow() %>% add_recipe(recipe_obj) %>% add_model(spec)
    res <- tune_grid(
        wf, resamples = folds, grid = grid_size,
        metrics = metric_set(roc_auc, pr_auc),
        control = control_grid(save_pred = TRUE, verbose = FALSE)
    )
    best <- select_best(res, metric = "roc_auc")
    list(
        workflow = finalize_workflow(wf, best),
        tuning = res,
        best = best
    )
}

run_window <- function(mimic_raw, eicu_raw, window, out_dir) {
    message("Running window: ", window)
    mimic <- prepare_dataset(mimic_raw, window)
    eicu <- prepare_dataset(eicu_raw, window)
    # A patient may contribute more than one ICU stay, and outcomes can differ
    # between stays. Create a patient-level binary stratum before grouped
    # splitting so rsample does not reject non-constant strata within groups.
    mimic$split_stratum <- factor(
        ave(as.integer(mimic$outcome == "yes"), mimic$split_group, FUN = max),
        levels = c(0, 1), labels = c("no", "yes")
    )
    eicu$split_stratum <- factor(
        ave(as.integer(eicu$outcome == "yes"), eicu$split_group, FUN = max),
        levels = c(0, 1), labels = c("no", "yes")
    )
    set.seed(20260901)
    split <- group_initial_split(
        mimic, group = split_group, strata = split_stratum, prop = 0.80
    )
    train <- training(split)
    test <- testing(split)
    folds <- group_vfold_cv(
        train, group = split_group, v = 5, strata = split_stratum
    )
    rec <- make_recipe(train)

    model_set <- strsplit(Sys.getenv(
        "STROKE_MODEL_SET", unset = "xgboost,elastic_net,random_forest"
    ), ",", fixed = TRUE)[[1L]]
    model_set <- trimws(model_set)
    model_set <- model_set[model_set %in% c("xgboost", "elastic_net", "random_forest")]
    if (length(model_set) == 0L) model_set <- "xgboost"
    models <- lapply(model_set,
        tune_one, recipe_obj = rec, folds = folds,
        grid_size = tuning_grid_size)
    names(models) <- model_set

    metric_rows <- list()
    prediction_rows <- list()
    for (model_name in names(models)) {
        model <- models[[model_name]]
        fitted_train <- fit(model$workflow, data = train)
        fitted_full <- fit(model$workflow, data = mimic)
        p_test <- predict(fitted_train, test, type = "prob")$.pred_yes
        p_eicu <- predict(fitted_full, eicu, type = "prob")$.pred_yes

        test_score <- score_predictions(test$outcome, p_test) %>%
            mutate(model = model_name, dataset = "MIMIC_internal_test", window = window)
        eicu_score <- score_predictions(eicu$outcome, p_eicu) %>%
            mutate(model = model_name, dataset = "eICU_external_validation", window = window)
        metric_rows[[paste0(model_name, "_test")]] <- test_score
        metric_rows[[paste0(model_name, "_eicu")]] <- eicu_score

        prediction_rows[[paste0(model_name, "_test")]] <- tibble(
            model = model_name, dataset = "MIMIC_internal_test", window = window,
            group_id = test$split_group, truth = test$outcome, probability = p_test
        )
        prediction_rows[[paste0(model_name, "_eicu")]] <- tibble(
            model = model_name, dataset = "eICU_external_validation", window = window,
            group_id = eicu$split_group, truth = eicu$outcome, probability = p_eicu
        )
        dc_test <- decision_curve(test$outcome, p_test) %>%
            mutate(model = model_name, dataset = "MIMIC_internal_test", window = window)
        dc_eicu <- decision_curve(eicu$outcome, p_eicu) %>%
            mutate(model = model_name, dataset = "eICU_external_validation", window = window)
        write.csv(bind_rows(dc_test, dc_eicu),
            file.path(out_dir, paste0("decision_curve_", window, "_", model_name, ".csv")),
            row.names = FALSE)
        write.csv(collect_metrics(model$tuning),
            file.path(out_dir, paste0("cv_metrics_", window, "_", model_name, ".csv")),
            row.names = FALSE)
        saveRDS(fitted_full,
            file.path(out_dir, paste0("model_", window, "_", model_name, "_mimic_full.rds")))
    }
    write.csv(bind_rows(metric_rows),
        file.path(out_dir, paste0("model_metrics_", window, ".csv")), row.names = FALSE)
    write.csv(bind_rows(prediction_rows),
        file.path(out_dir, paste0("predictions_", window, ".csv")), row.names = FALSE)
    invisible(metric_rows)
}

if (!identical(Sys.getenv("STROKE_SKIP_RUN", unset = "0"), "1")) {
    mimic_raw <- read.csv(mimic_path, check.names = FALSE)
    eicu_raw <- read.csv(eicu_path, check.names = FALSE)
    run_window(mimic_raw, eicu_raw, "landmark", out_dir)
    run_window(mimic_raw, eicu_raw, "early", out_dir)
    message("Completed. Results written to: ", normalizePath(out_dir, mustWork = FALSE))
}
