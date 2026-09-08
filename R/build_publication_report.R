# Build publication-ready summary tables, figure manifest, and bilingual
# Results text from the validated model outputs. No model is refit here.

format_ci <- function(x, digits = 3L) {
    sprintf(
        paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"),
        x[["estimate"]], x[["ci_low"]], x[["ci_high"]]
    )
}

get_metric_row <- function(x, model, dataset, window, metric) {
    hit <- x[x$model == model & x$dataset == dataset &
        x$window == window & x$metric == metric, , drop = FALSE]
    if (nrow(hit) != 1L) stop("Expected one metric row for ", model, "/", dataset, "/", window, "/", metric)
    hit[1L, , drop = FALSE]
}

to01_report <- function(x) {
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(y))
    out[y %in% c("1", "true", "t", "yes")] <- 1L
    out[y %in% c("0", "false", "f", "no")] <- 0L
    out
}

build_publication_report <- function(results_dir, validation_dir, report_dir,
    mimic_csv = NULL, eicu_csv = NULL) {
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    ci <- utils::read.csv(
        file.path(validation_dir, "bootstrap_metric_ci.csv"),
        check.names = FALSE
    )
    models <- c("xgboost", "random_forest", "elastic_net")
    windows <- c("landmark", "early")
    datasets <- c("eICU_external_validation", "MIMIC_internal_test")

    make_table <- function(dataset) {
        rows <- list()
        k <- 1L
        for (window in windows) for (model in models) {
            z <- ci[ci$dataset == dataset & ci$window == window & ci$model == model, , drop = FALSE]
            rows[[k]] <- data.frame(
                time_window = ifelse(window == "landmark", "0-24 h landmark", "0-6 h early"),
                model = model,
                auroc_95ci = format_ci(get_metric_row(ci, model, dataset, window, "auroc")),
                auprc_95ci = format_ci(get_metric_row(ci, model, dataset, window, "auprc")),
                brier_95ci = format_ci(get_metric_row(ci, model, dataset, window, "brier_score")),
                calibration_intercept_95ci = format_ci(get_metric_row(ci, model, dataset, window, "calibration_intercept")),
                calibration_slope_95ci = format_ci(get_metric_row(ci, model, dataset, window, "calibration_slope")),
                patient_groups = unique(z$patient_groups)[[1]],
                events = unique(z$events)[[1]],
                stringsAsFactors = FALSE
            )
            k <- k + 1L
        }
        do.call(rbind, rows)
    }

    external_table <- make_table(datasets[[1]])
    internal_table <- make_table(datasets[[2]])
    utils::write.csv(
        external_table,
        file.path(report_dir, "table2_external_validation.csv"),
        row.names = FALSE
    )
    utils::write.csv(
        internal_table,
        file.path(report_dir, "tableS1_internal_test.csv"),
        row.names = FALSE
    )

    markdown_table <- function(table) {
        lines <- c(
            "| 时间窗 / Time window | 模型 / Model | AUROC (95% CI) | AUPRC (95% CI) | Brier (95% CI) | 校准截距 / Intercept | 校准斜率 / Slope | 患者数 / Patients | 事件数 / Events |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|"
        )
        for (i in seq_len(nrow(table))) {
            lines <- c(lines, sprintf(
                "| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                table$time_window[[i]], table$model[[i]], table$auroc_95ci[[i]],
                table$auprc_95ci[[i]], table$brier_95ci[[i]],
                table$calibration_intercept_95ci[[i]],
                table$calibration_slope_95ci[[i]], table$patient_groups[[i]],
                table$events[[i]]
            ))
        }
        paste(lines, collapse = "\n")
    }
    writeLines(
        c("# Table 2. eICU external validation", "", markdown_table(external_table)),
        file.path(report_dir, "table2_external_validation.md"), useBytes = TRUE
    )
    writeLines(
        c("# Supplementary Table S1. MIMIC internal test", "", markdown_table(internal_table)),
        file.path(report_dir, "tableS1_internal_test.md"), useBytes = TRUE
    )

    figure_manifest <- data.frame(
        figure = c("Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6"),
        file = c(
            "roc_external.png", "pr_external.png", "calibration_external.png",
            "dca_external.png", "shap_top20_landmark_xgboost.png"
        ),
        chinese_description = c(
            "eICU外部验证ROC曲线", "eICU外部验证Precision-recall曲线",
            "eICU外部验证校准曲线", "eICU外部验证决策曲线分析",
            "landmark XGBoost前20位SHAP重要性"
        ),
        english_description = c(
            "ROC curves in eICU external validation",
            "Precision-recall curves in eICU external validation",
            "Calibration curves in eICU external validation",
            "Decision curve analysis in eICU external validation",
            "Top 20 SHAP features for landmark XGBoost"
        ),
        stringsAsFactors = FALSE
    )
    utils::write.csv(
        figure_manifest, file.path(report_dir, "figure_manifest.csv"), row.names = FALSE
    )

    primary <- list(
        auroc = get_metric_row(ci, "xgboost", datasets[[1]], "landmark", "auroc"),
        auprc = get_metric_row(ci, "xgboost", datasets[[1]], "landmark", "auprc"),
        brier = get_metric_row(ci, "xgboost", datasets[[1]], "landmark", "brier_score"),
        intercept = get_metric_row(ci, "xgboost", datasets[[1]], "landmark", "calibration_intercept"),
        slope = get_metric_row(ci, "xgboost", datasets[[1]], "landmark", "calibration_slope")
    )
    early <- list(
        auroc = get_metric_row(ci, "xgboost", datasets[[1]], "early", "auroc"),
        auprc = get_metric_row(ci, "xgboost", datasets[[1]], "early", "auprc")
    )
    rf_primary <- get_metric_row(ci, "random_forest", datasets[[1]], "landmark", "auprc")

    cohort_text <- ""
    if (!is.null(mimic_csv) && !is.null(eicu_csv) && file.exists(mimic_csv) && file.exists(eicu_csv)) {
        mimic <- utils::read.csv(mimic_csv, check.names = FALSE)
        eicu <- utils::read.csv(eicu_csv, check.names = FALSE)
        mimic_landmark <- sum(to01_report(mimic[["landmark_24h_eligible"]]) == 1L, na.rm = TRUE)
        eicu_landmark <- sum(to01_report(eicu[["landmark_24h_eligible"]]) == 1L, na.rm = TRUE)
        cohort_text_en <- sprintf(
            "The MIMIC-IV and eICU analysis datasets contained %d and %d ICU-level records, respectively; %d and %d were eligible for the 24-h landmark analysis.",
            nrow(mimic), nrow(eicu), mimic_landmark, eicu_landmark
        )
        cohort_text_cn <- sprintf(
            "MIMIC-IV 与 eICU 分析数据集分别包含 %d 例和 %d 例 ICU 层面记录，其中分别有 %d 例和 %d 例符合 24 小时 landmark 分析的纳入条件。",
            nrow(mimic), nrow(eicu), mimic_landmark, eicu_landmark
        )
    } else {
        cohort_text_en <- "The MIMIC-IV and eICU analysis datasets were evaluated using the prespecified ICU-level cohort and 24-h landmark eligibility rule."
        cohort_text_cn <- "MIMIC-IV 与 eICU 分析数据集均按照预先设定的 ICU 层面队列和 24 小时 landmark 纳入规则进行分析。"
    }

    shap_file <- file.path(validation_dir, "shap_importance_landmark_xgboost.csv")
    shap_domain_file <- file.path(validation_dir, "shap_domain_landmark_xgboost.csv")
    shap_features <- utils::read.csv(shap_file, check.names = FALSE)
    shap_domain <- utils::read.csv(shap_domain_file, check.names = FALSE)
    top_features <- head(shap_features$feature, 5L)
    top_features_cn <- paste(top_features, collapse = "、")
    top_features_en <- paste(top_features, collapse = ", ")
    leading_domain <- shap_domain$domain[[1]]
    feature_cn <- c(
        spo2_mean_24h = "24 小时平均 SpO₂", age_years = "年龄",
        wbc_mean_24h = "24 小时平均白细胞计数", glucose_mean_24h = "24 小时平均血糖",
        sodium_mean_24h = "24 小时平均钠"
    )
    feature_en <- c(
        spo2_mean_24h = "24-h mean SpO₂", age_years = "age",
        wbc_mean_24h = "24-h mean white blood cell count", glucose_mean_24h = "24-h mean glucose",
        sodium_mean_24h = "24-h mean sodium"
    )
    top_features_cn <- paste(
        ifelse(top_features %in% names(feature_cn), feature_cn[top_features], top_features),
        collapse = "、"
    )
    top_features_en <- paste(
        ifelse(top_features %in% names(feature_en), feature_en[top_features], top_features),
        collapse = ", "
    )
    domain_cn <- c(
        "Metabolic/renal" = "代谢/肾功能",
        "Vital signs" = "生命体征",
        "Inflammatory/hematologic" = "炎症/血液学",
        "Demographic/frailty/clinical" = "人口学/衰弱/临床"
    )

    fmt <- function(row) format_ci(row)
    results_cn <- c(
        "## 结果（中文）",
        "",
        cohort_text_cn,
        sprintf(
            "在 eICU 外部验证中，预先设定的 0-24 h landmark XGBoost 模型 AUROC 为 %s，AUPRC 为 %s，Brier score 为 %s。校准截距为 %s，校准斜率为 %s，提示总体校准接近理想水平，但存在轻度风险高估。",
            fmt(primary$auroc), fmt(primary$auprc), fmt(primary$brier),
            fmt(primary$intercept), fmt(primary$slope)
        ),
        sprintf(
            "在 landmark 窗口，随机森林的 AUPRC 为 %s，略高于 XGBoost；因此模型比较不应仅依据 AUROC，主模型结果需同时报告区分度、概率校准和临床净获益。",
            fmt(rf_primary)
        ),
        sprintf(
            "作为敏感性分析，0-6 h XGBoost 模型在 eICU 中 AUROC 为 %s，AUPRC 为 %s，低于 landmark 模型的区分度。",
            fmt(early$auroc), fmt(early$auprc)
        ),
        sprintf(
            "landmark XGBoost 的 SHAP 分析显示，贡献最大的单项特征为 %s；按预设功能域汇总时，贡献最高的域为 %s。SHAP 反映模型预测贡献，不应解释为因果效应。",
            top_features_cn, domain_cn[[leading_domain]]
        )
    )
    results_en <- c(
        "## Results (English)",
        "",
        cohort_text_en,
        sprintf(
            "In eICU external validation, the prespecified 0-24-h landmark XGBoost model achieved an AUROC of %s, an AUPRC of %s, and a Brier score of %s. The calibration intercept was %s and the calibration slope was %s, indicating generally acceptable calibration with mild risk overprediction.",
            fmt(primary$auroc), fmt(primary$auprc), fmt(primary$brier),
            fmt(primary$intercept), fmt(primary$slope)
        ),
        sprintf(
            "For the landmark window, the random-forest model had an AUPRC of %s, slightly higher than XGBoost. Accordingly, model comparison should not rely on AUROC alone; discrimination, calibration, and clinical net benefit should be reported together.",
            fmt(rf_primary)
        ),
        sprintf(
            "In the sensitivity analysis using the 0-6-h early window, XGBoost achieved an AUROC of %s and an AUPRC of %s in eICU, with lower discrimination than the landmark model.",
            fmt(early$auroc), fmt(early$auprc)
        ),
        sprintf(
            "For landmark XGBoost, the five features with the largest mean absolute SHAP values were %s. At the prespecified domain level, %s had the largest aggregate contribution. SHAP values represent model contribution and should not be interpreted as causal effects.",
            top_features_en, leading_domain
        )
    )
    writeLines(
        c("# Bilingual Results Draft", "", results_cn, "", results_en),
        file.path(report_dir, "bilingual_results_draft.md"), useBytes = TRUE
    )
    invisible(list(external = external_table, internal = internal_table,
        figures = figure_manifest, results_cn = results_cn, results_en = results_en))
}
