# Publication-grade validation and explainability for the saved stroke models.
# This script consumes model predictions and saved workflows; it does not query
# or modify either database and does not refit the prediction models.

options(stringsAsFactors = FALSE)

safe_metric <- function(fun, truth, probability) {
    tryCatch(
        suppressWarnings(as.numeric(fun(truth, probability, event_level = "second"))),
        error = function(e) NA_real_
    )
}

safe_calibration <- function(truth, probability) {
    y <- as.integer(truth == "yes")
    p <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    if (length(unique(y)) < 2L) {
        return(c(calibration_intercept = NA_real_, calibration_slope = NA_real_))
    }
    lp <- qlogis(p)
    intercept <- tryCatch(
        suppressWarnings(coef(glm(y ~ 1, family = binomial(), offset = lp))[[1]]),
        error = function(e) NA_real_
    )
    slope <- tryCatch(
        suppressWarnings(coef(glm(y ~ lp, family = binomial()))[[2]]),
        error = function(e) NA_real_
    )
    c(calibration_intercept = unname(intercept), calibration_slope = unname(slope))
}

score_prediction_vectors <- function(truth, probability) {
    cal <- safe_calibration(truth, probability)
    tibble::tibble(
        auroc = safe_metric(yardstick::roc_auc_vec, truth, probability),
        auprc = safe_metric(yardstick::pr_auc_vec, truth, probability),
        brier_score = mean((as.integer(truth == "yes") - probability)^2),
        calibration_intercept = cal[[1]],
        calibration_slope = cal[[2]]
    )
}

bootstrap_group_metrics <- function(predictions, reps = 1000L, seed = 20260901L) {
    set.seed(seed)
    predictions$group_id <- as.character(predictions$group_id)
    groups <- unique(predictions$group_id)
    if (length(groups) < 2L) stop("At least two patient groups are required for bootstrap.")
    point <- score_prediction_vectors(predictions$truth, predictions$probability)
    boot <- vector("list", reps)
    for (i in seq_len(reps)) {
        drawn <- sample(groups, size = length(groups), replace = TRUE)
        row_index <- unlist(lapply(drawn, function(g) which(predictions$group_id == g)),
            use.names = FALSE)
        boot[[i]] <- score_prediction_vectors(
            predictions$truth[row_index], predictions$probability[row_index]
        )
    }
    boot_df <- dplyr::bind_rows(boot)
    metrics <- names(point)
    result <- dplyr::bind_rows(lapply(metrics, function(metric) {
        vals <- boot_df[[metric]]
        tibble::tibble(
            metric = metric,
            estimate = point[[metric]],
            ci_low = stats::quantile(vals, 0.025, na.rm = TRUE, names = FALSE),
            ci_high = stats::quantile(vals, 0.975, na.rm = TRUE, names = FALSE),
            bootstrap_reps = reps,
            patient_groups = length(groups),
            events = sum(predictions$truth == "yes")
        )
    }))
    result
}

make_calibration_bins <- function(predictions, bins = 10L) {
    predictions %>%
        dplyr::group_by(model, dataset, window) %>%
        dplyr::mutate(bin = dplyr::ntile(probability, bins)) %>%
        dplyr::group_by(model, dataset, window, bin) %>%
        dplyr::summarise(
            mean_predicted = mean(probability, na.rm = TRUE),
            observed_rate = mean(truth == "yes"),
            n = dplyr::n(),
            .groups = "drop"
        )
}

feature_domain <- function(feature) {
    dplyr::case_when(
        grepl("age_years|frailty_comorbidity_count|sex_|stroke_group_", feature) ~
            "Demographic/frailty/clinical",
        grepl("glucose|sodium|potassium|chloride|bicarbonate|bun|creatinine|calcium",
            feature) ~ "Metabolic/renal",
        grepl("wbc|platelets|hemoglobin|hematocrit", feature) ~
            "Inflammatory/hematologic",
        grepl("heart_rate|respiratory_rate|spo2|sbp|dbp|map", feature) ~
            "Vital signs",
        grepl("missing", feature) ~ "Missingness",
        TRUE ~ "Other"
    )
}

run_shap <- function(window, results_dir, validation_dir, max_rows = 10000L) {
    workflow_file <- file.path(
        results_dir, paste0("model_", window, "_xgboost_mimic_full.rds")
    )
    if (!file.exists(workflow_file)) stop("Missing workflow: ", workflow_file)
    workflow_fit <- readRDS(workflow_file)
    mold <- workflows::extract_mold(workflow_fit)
    engine <- workflows::extract_fit_engine(workflow_fit)
    x <- as.matrix(mold[["predictors"]])
    if (nrow(x) > max_rows) {
        set.seed(20260901)
        x <- x[sample(seq_len(nrow(x)), max_rows), , drop = FALSE]
    }
    shap <- predict(engine, newdata = x, predcontrib = TRUE)
    shap <- as.data.frame(shap, check.names = FALSE)
    bias_name <- tail(names(shap), 1L)
    shap_features <- shap[, setdiff(names(shap), bias_name), drop = FALSE]
    importance <- tibble::tibble(
        feature = names(shap_features),
        mean_abs_shap = vapply(shap_features, function(z) mean(abs(z)), numeric(1)),
        mean_shap = vapply(shap_features, mean, numeric(1)),
        domain = feature_domain(names(shap_features))
    ) %>%
        dplyr::arrange(dplyr::desc(mean_abs_shap))
    domain_summary <- importance %>%
        dplyr::group_by(domain) %>%
        dplyr::summarise(
            total_mean_abs_shap = sum(mean_abs_shap),
            mean_feature_abs_shap = mean(mean_abs_shap),
            n_features = dplyr::n(),
            .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(total_mean_abs_shap))
    shap_long <- shap_features %>%
        tibble::as_tibble() %>%
        dplyr::mutate(row_id = dplyr::row_number()) %>%
        tidyr::pivot_longer(
            cols = -row_id, names_to = "feature", values_to = "shap_value"
        )
    utils::write.csv(
        importance,
        file.path(validation_dir, paste0("shap_importance_", window, "_xgboost.csv")),
        row.names = FALSE
    )
    utils::write.csv(
        domain_summary,
        file.path(validation_dir, paste0("shap_domain_", window, "_xgboost.csv")),
        row.names = FALSE
    )
    utils::write.csv(
        shap_long,
        file.path(validation_dir, paste0("shap_values_", window, "_xgboost.csv")),
        row.names = FALSE
    )
    top <- importance %>% dplyr::slice_head(n = 20L)
    plot <- ggplot2::ggplot(
        top,
        ggplot2::aes(x = stats::reorder(feature, mean_abs_shap), y = mean_abs_shap,
            fill = domain)
    ) +
        ggplot2::geom_col(show.legend = TRUE) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = NULL, y = "Mean absolute SHAP value",
            title = paste0("Top SHAP features: ", window, " XGBoost")
        ) +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(
        file.path(validation_dir, paste0("shap_top20_", window, "_xgboost.png")),
        plot, width = 9, height = 7, dpi = 300
    )
    invisible(list(importance = importance, domain = domain_summary))
}

run_stroke_validation <- function(results_dir, validation_dir,
    bootstrap_reps = 1000L) {
    dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
    prediction_files <- list.files(
        results_dir, pattern = "^predictions_(landmark|early)\\.csv$",
        full.names = TRUE
    )
    if (length(prediction_files) != 2L) {
        stop("Expected two prediction files in ", results_dir)
    }
    predictions <- dplyr::bind_rows(lapply(prediction_files, function(f) {
        x <- utils::read.csv(f, check.names = FALSE)
        x$truth <- factor(x$truth, levels = c("no", "yes"))
        x$probability <- as.numeric(x$probability)
        x$group_id <- as.character(x$group_id)
        x
    }))

    groups <- split(
        predictions,
        interaction(predictions$model, predictions$dataset, predictions$window,
            drop = TRUE, sep = "||")
    )
    bootstrap <- dplyr::bind_rows(lapply(groups, function(x) {
        key <- paste(x$model[[1]], x$dataset[[1]], x$window[[1]], sep = "||")
        out <- bootstrap_group_metrics(x, reps = bootstrap_reps)
        out$model <- x$model[[1]]
        out$dataset <- x$dataset[[1]]
        out$window <- x$window[[1]]
        out
    })) %>%
        dplyr::select(model, dataset, window, dplyr::everything())
    utils::write.csv(
        bootstrap, file.path(validation_dir, "bootstrap_metric_ci.csv"), row.names = FALSE
    )

    calibration <- make_calibration_bins(predictions)
    utils::write.csv(
        calibration, file.path(validation_dir, "calibration_bins.csv"), row.names = FALSE
    )

    roc_data <- dplyr::bind_rows(lapply(groups, function(x) {
        curve <- yardstick::roc_curve(x, truth, probability, event_level = "second")
        curve$model <- x$model[[1]]
        curve$dataset <- x$dataset[[1]]
        curve$window <- x$window[[1]]
        curve
    }))
    pr_data <- dplyr::bind_rows(lapply(groups, function(x) {
        curve <- yardstick::pr_curve(x, truth, probability, event_level = "second")
        curve$model <- x$model[[1]]
        curve$dataset <- x$dataset[[1]]
        curve$window <- x$window[[1]]
        curve
    }))
    utils::write.csv(roc_data, file.path(validation_dir, "roc_curve_data.csv"), row.names = FALSE)
    utils::write.csv(pr_data, file.path(validation_dir, "pr_curve_data.csv"), row.names = FALSE)

    external_roc <- roc_data %>% dplyr::filter(dataset == "eICU_external_validation")
    roc_plot <- ggplot2::ggplot(
        external_roc,
        ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)
    ) +
        ggplot2::geom_path(linewidth = 0.9) +
        ggplot2::geom_abline(linetype = 2, color = "grey50") +
        ggplot2::facet_wrap(~window) +
        ggplot2::coord_equal() +
        ggplot2::labs(title = "eICU external validation: ROC curves",
            x = "1 - Specificity", y = "Sensitivity", color = "Model") +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(
        file.path(validation_dir, "roc_external.png"), roc_plot,
        width = 10, height = 5.5, dpi = 300
    )

    external_pr <- pr_data %>% dplyr::filter(dataset == "eICU_external_validation")
    pr_plot <- ggplot2::ggplot(
        external_pr, ggplot2::aes(x = recall, y = precision, color = model)
    ) +
        ggplot2::geom_path(linewidth = 0.9) +
        ggplot2::facet_wrap(~window) +
        ggplot2::labs(title = "eICU external validation: Precision-recall curves",
            x = "Recall", y = "Precision", color = "Model") +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(
        file.path(validation_dir, "pr_external.png"), pr_plot,
        width = 10, height = 5.5, dpi = 300
    )

    external_calibration <- calibration %>%
        dplyr::filter(dataset == "eICU_external_validation")
    calibration_plot <- ggplot2::ggplot(
        external_calibration,
        ggplot2::aes(x = mean_predicted, y = observed_rate, color = model)
    ) +
        ggplot2::geom_abline(linetype = 2, color = "grey50") +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::facet_wrap(~window) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "eICU external validation: Calibration curves",
            x = "Mean predicted risk", y = "Observed risk", color = "Model") +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(
        file.path(validation_dir, "calibration_external.png"), calibration_plot,
        width = 10, height = 5.5, dpi = 300
    )

    dca_files <- list.files(
        results_dir, pattern = "^decision_curve_(landmark|early)_.+\\.csv$",
        full.names = TRUE
    )
    dca <- dplyr::bind_rows(lapply(dca_files, utils::read.csv, check.names = FALSE))
    prevalences <- predictions %>%
        dplyr::group_by(dataset, window) %>%
        dplyr::summarise(prevalence = mean(truth == "yes"), .groups = "drop")
    dca_all <- dca %>%
        dplyr::filter(dataset == "eICU_external_validation") %>%
        dplyr::left_join(prevalences, by = c("dataset", "window"))
    treat_all <- dca_all %>%
        dplyr::distinct(dataset, window, threshold, prevalence) %>%
        dplyr::mutate(model = "treat_all",
            net_benefit = prevalence - (1 - prevalence) * threshold / (1 - threshold))
    treat_none <- dca_all %>%
        dplyr::distinct(dataset, window, threshold) %>%
        dplyr::mutate(model = "treat_none", net_benefit = 0)
    dca_plot_data <- dplyr::bind_rows(
        dca_all %>% dplyr::select(threshold, net_benefit, model, dataset, window),
        treat_all %>% dplyr::select(threshold, net_benefit, model, dataset, window),
        treat_none %>% dplyr::select(threshold, net_benefit, model, dataset, window)
    )
    utils::write.csv(
        dca_plot_data, file.path(validation_dir, "dca_external_data.csv"), row.names = FALSE
    )
    dca_plot <- ggplot2::ggplot(
        dca_plot_data, ggplot2::aes(x = threshold, y = net_benefit,
            color = model, linetype = model)
    ) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::facet_wrap(~window) +
        ggplot2::labs(title = "eICU external validation: Decision curve analysis",
            x = "Threshold probability", y = "Net benefit", color = "Model",
            linetype = "Model") +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(
        file.path(validation_dir, "dca_external.png"), dca_plot,
        width = 10, height = 5.5, dpi = 300
    )

    run_shap("landmark", results_dir, validation_dir)
    run_shap("early", results_dir, validation_dir)
    invisible(list(bootstrap = bootstrap, calibration = calibration,
        roc = roc_data, pr = pr_data, dca = dca_plot_data))
}
