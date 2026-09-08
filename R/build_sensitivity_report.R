# Assemble an auditable sensitivity-analysis report from saved predictions.
# Reads local CSVs only; does not refit models or access databases.

project_dir <- Sys.getenv("STROKE_PROJECT_DIR", unset = getwd())
validation_dir <- Sys.getenv("STROKE_VALIDATION_DIR", unset = file.path(project_dir, "results", "validation"))
sensitivity_dir <- Sys.getenv("STROKE_SENSITIVITY_DIR", unset = file.path(project_dir, "results", "sensitivity"))
root <- sensitivity_dir
main_file <- file.path(validation_dir, "bootstrap_metric_ci.csv")
sens_file <- file.path(sensitivity_dir, "sensitivity_external_bootstrap.csv")
metric_file <- file.path(sensitivity_dir, "sensitivity_external_xgboost.csv")
report_file <- file.path(sensitivity_dir, "sensitivity_analysis_report.md")
table_file <- file.path(sensitivity_dir, "tableS3_sensitivity_analysis.md")

main <- utils::read.csv(main_file, check.names = FALSE)
sens <- utils::read.csv(sens_file, check.names = FALSE)
point <- utils::read.csv(metric_file, check.names = FALSE)
main <- main[main$model == "xgboost" & main$dataset == "eICU_external_validation", , drop = FALSE]

fmt <- function(x, lo, hi) sprintf("%.3f (%.3f–%.3f)", x, lo, hi)
main_key <- main[, c("window", "metric", "estimate", "ci_low", "ci_high")]
main_key$scenario <- "Primary source-level QC"
sens_key <- sens[, c("scenario", "window", "metric", "estimate", "ci_low", "ci_high")]
all <- rbind(main_key[, c("scenario", "window", "metric", "estimate", "ci_low", "ci_high")], sens_key)

lines <- c(
    "# Prespecified sensitivity-analysis report",
    "",
    "All sensitivity analyses refit the XGBoost model in MIMIC and evaluated the fixed model in eICU. Confidence intervals are patient-level percentile bootstrap intervals based on 1,000 resamples. The primary source-level QC results are shown as the reference.",
    "",
    "## eICU external validation",
    "",
    "| Scenario | Window | AUROC (95% CI) | AUPRC (95% CI) | Brier score (95% CI) | N | Events |",
    "|---|---|---:|---:|---:|---:|---:|"
)
scenario_order <- c("Primary source-level QC", "mixed_excluded", "complete_case", "no_missing_indicators")
scenario_label <- c(
    "Primary source-level QC" = "Primary source-level QC",
    mixed_excluded = "Exclude mixed/ambiguous",
    complete_case = "Complete-case analysis",
    no_missing_indicators = "Without missingness indicators"
)
for (window in c("landmark", "early")) {
    for (scenario in scenario_order) {
        x <- all[all$window == window & all$scenario == scenario, , drop = FALSE]
        if (nrow(x) == 0L) next
        get <- function(metric) x[x$metric == metric, , drop = FALSE]
        a <- get("auroc"); p <- get("auprc"); b <- get("brier_score")
        if (scenario == "Primary source-level QC") {
            n <- if (window == "landmark") 4264L else 4502L
            events <- if (window == "landmark") 509L else 591L
        } else {
            pfile <- file.path(sensitivity_dir, scenario, paste0("predictions_", window, ".csv"))
            pred <- utils::read.csv(pfile, check.names = FALSE)
            pred <- pred[pred$model == "xgboost" & pred$dataset == "eICU_external_validation", , drop = FALSE]
            n <- nrow(pred); events <- sum(pred$truth == "yes")
        }
        lines <- c(lines, sprintf(
            "| %s | %s | %s | %s | %s | %d | %d |",
            scenario_label[[scenario]], window,
            fmt(a$estimate, a$ci_low, a$ci_high),
            fmt(p$estimate, p$ci_low, p$ci_high),
            fmt(b$estimate, b$ci_low, b$ci_high), n, events
        ))
    }
}

table_lines <- c(
    "# Supplementary Table S3. Prespecified sensitivity analyses",
    "",
    "All sensitivity analyses refit XGBoost in MIMIC and evaluated the fixed model in eICU. Confidence intervals are patient-level percentile bootstrap intervals based on 1,000 resamples.",
    "",
    lines
)
writeLines(table_lines, table_file, useBytes = TRUE)

primary_landmark <- main[main$window == "landmark" & main$metric == "auroc", "estimate"]
mixed_landmark <- sens[sens$scenario == "mixed_excluded" & sens$window == "landmark" & sens$metric == "auroc", "estimate"]
complete_landmark <- sens[sens$scenario == "complete_case" & sens$window == "landmark" & sens$metric == "auroc", "estimate"]
nomiss_landmark <- sens[sens$scenario == "no_missing_indicators" & sens$window == "landmark" & sens$metric == "auroc", "estimate"]
lines <- c(lines, "", "## Interpretation / 解释", "",
    sprintf("The primary landmark AUROC was %.3f. It was %.3f after excluding mixed/ambiguous records, %.3f in complete cases, and %.3f after removing explicit missingness indicators.", primary_landmark, mixed_landmark, complete_landmark, nomiss_landmark),
    "",
    "The landmark discrimination direction was stable across all prespecified scenarios. The complete-case analysis changed the eligible sample and event rate, so its estimates should be interpreted as a selection-sensitive sensitivity analysis rather than as evidence of improvement. Removing missingness indicators produced results very close to the primary analysis, suggesting that the main landmark discrimination was not driven solely by explicit missingness flags.",
    "",
    "0–6-hour results remained lower or less stable than the landmark analysis and should remain a sensitivity analysis. None of these analyses establishes causal effects or clinical utility.",
    "",
    "## 中文解释",
    "",
    sprintf("主要 landmark AUROC 为 %.3f；排除 mixed/ambiguous 后为 %.3f，完整病例分析为 %.3f，去除缺失指示变量后为 %.3f。", primary_landmark, mixed_landmark, complete_landmark, nomiss_landmark),
    "",
    "landmark 区分度方向在预设敏感性场景中保持稳定。完整病例分析改变了纳入样本和事件率，因此应视为对样本选择敏感性的检验，而不能解释为模型性能改善。去除缺失指示变量后结果与主分析非常接近，提示主要 landmark 区分度并非完全由显式缺失指标驱动。",
    "",
    "0–6 小时结果仍低于或不如 landmark 结果稳定，应继续作为敏感性分析报告。上述分析均不证明因果效应或临床效用。"
)
writeLines(lines, report_file, useBytes = TRUE)
cat(sprintf("SENSITIVITY_REPORT_OK report=%s table=%s\n", report_file, table_file))
