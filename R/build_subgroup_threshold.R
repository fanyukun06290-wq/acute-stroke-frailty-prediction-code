# Build subgroup performance and fixed-threshold operating characteristics.
# Reads saved predictions and analysis CSV files; does not query or modify databases.

project_dir <- Sys.getenv(
    "STROKE_PROJECT_DIR", unset = getwd()
)
results_dir <- Sys.getenv(
    "STROKE_RESULTS_DIR", unset = file.path(project_dir, "results", "model_results")
)
publication_dir <- Sys.getenv(
    "STROKE_PUBLICATION_DIR", unset = file.path(project_dir, "results", "publication")
)
dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)

mimic_csv <- Sys.getenv(
    "STROKE_MIMIC_CSV", unset = file.path(project_dir, "data", "mimic_analysis.csv")
)
eicu_csv <- Sys.getenv(
    "STROKE_EICU_CSV", unset = file.path(project_dir, "data", "eicu_analysis.csv")
)
if (!all(file.exists(c(mimic_csv, eicu_csv)))) stop("未找到 analysis CSV 文件。")

mimic_raw <- utils::read.csv(mimic_csv, check.names = FALSE, stringsAsFactors = FALSE)
eicu_raw <- utils::read.csv(eicu_csv, check.names = FALSE, stringsAsFactors = FALSE)

to01 <- function(x) {
    y <- tolower(trimws(as.character(x)))
    z <- rep(NA_integer_, length(y))
    z[y %in% c("1", "true", "t", "yes")] <- 1L
    z[y %in% c("0", "false", "f", "no")] <- 0L
    z
}

collapse_stroke_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "AIS_STRICT", "AIS",
        ifelse(x == "ICH_STRICT", "ICH",
            ifelse(x == "SAH", "SAH", "Mixed/ambiguous")))
}

make_metadata <- function(raw, dataset) {
    if (dataset == "MIMIC_internal_test") {
        group_id <- as.character(raw$subject_id)
    } else {
        group_id <- as.character(raw$patienthealthsystemstayid)
    }
    age <- suppressWarnings(as.numeric(as.character(raw$age_years)))
    age_group <- ifelse(is.na(age), "Unknown", ifelse(age < 65, "<65", ">=65"))
    sex_group <- trimws(as.character(raw$sex))
    sex_group[is.na(sex_group) | !sex_group %in% c("F", "M")] <- "Unknown"
    data.frame(
        group_id = group_id,
        age_group = age_group,
        sex_group = ifelse(sex_group == "F", "Female",
            ifelse(sex_group == "M", "Male", "Unknown")),
        stroke_group = collapse_stroke_group(raw$audit_phenotype),
        stringsAsFactors = FALSE
    )
}

fast_auc <- function(y, probability) {
    y <- as.integer(as.character(y) == "yes")
    probability <- as.numeric(probability)
    ok <- is.finite(probability) & !is.na(y)
    y <- y[ok]
    probability <- probability[ok]
    n_pos <- sum(y == 1L)
    n_neg <- sum(y == 0L)
    if (n_pos == 0L || n_neg == 0L) return(c(auroc = NA_real_, auprc = NA_real_))

    # Rank formula for AUROC, with average ranks for tied probabilities.
    ranks <- rank(probability, ties.method = "average")
    auroc <- (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) /
        (n_pos * n_neg)

    # Match yardstick's binary PR-AUC: trapezoidal integration over the
    # precision-recall curve after collapsing tied probability thresholds.
    ord <- order(probability, decreasing = TRUE)
    ys <- y[ord]
    ps <- probability[ord]
    loc_unique <- which(!duplicated(ps, fromLast = TRUE))
    tp <- cumsum(ys)[loc_unique]
    fp <- cumsum(1L - ys)[loc_unique]
    recall <- c(0, tp / n_pos)
    precision <- c(1, tp / (tp + fp))
    auprc <- sum((precision[-1] + precision[-length(precision)]) / 2 *
        diff(recall))
    c(auroc = auroc, auprc = auprc)
}

safe_auc <- function(df) {
    fast_auc(df$truth, df$probability)
}

bootstrap_rows <- function(df, metric_fun, reps = 1000L, seed = 20260901L) {
    df <- df[!is.na(df$truth) & is.finite(df$probability), , drop = FALSE]
    groups <- unique(as.character(df$group_id))
    group_indices <- split(seq_len(nrow(df)), as.character(df$group_id))
    point <- metric_fun(df)
    if (length(groups) == 0L) return(list(point = point, ci = matrix(NA_real_, nrow = 2, ncol = length(point))))
    set.seed(seed)
    boot <- matrix(NA_real_, nrow = reps, ncol = length(point))
    colnames(boot) <- names(point)
    for (b in seq_len(reps)) {
        sampled_groups <- sample(groups, length(groups), replace = TRUE)
        idx <- unlist(group_indices[sampled_groups], use.names = FALSE)
        boot[b, ] <- metric_fun(df[idx, , drop = FALSE])
    }
    ci <- apply(boot, 2, function(x) {
        if (all(is.na(x))) return(c(lower = NA_real_, upper = NA_real_))
        stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    })
    if (is.null(dim(ci))) {
        ci <- matrix(ci, nrow = 2, dimnames = list(c("lower", "upper"), names(point)))
    } else {
        rownames(ci) <- c("lower", "upper")
    }
    list(point = point, ci = ci, boot = boot)
}

threshold_metrics <- function(df, threshold) {
    y <- as.integer(df$truth == "yes")
    pred <- as.integer(df$probability >= threshold)
    tp <- sum(pred == 1L & y == 1L)
    tn <- sum(pred == 0L & y == 0L)
    fp <- sum(pred == 1L & y == 0L)
    fn <- sum(pred == 0L & y == 1L)
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    f1 <- if (is.finite(ppv) && is.finite(sensitivity) && (ppv + sensitivity) > 0) {
        2 * ppv * sensitivity / (ppv + sensitivity)
    } else NA_real_
    c(sensitivity = sensitivity, specificity = specificity, ppv = ppv,
        npv = npv, f1 = f1)
}

threshold_bootstrap <- function(df, threshold, reps = 1000L, seed = 20260901L) {
    metric_fun <- function(x) threshold_metrics(x, threshold)
    bootstrap_rows(df, metric_fun, reps = reps, seed = seed)
}

prediction_files <- list(
    landmark = file.path(results_dir, "predictions_landmark.csv"),
    early = file.path(results_dir, "predictions_early.csv")
)
if (!all(file.exists(unlist(prediction_files)))) stop("预测结果文件不完整。")

# Subgroup performance is reported for the primary XGBoost model in eICU,
# where the external validation grouping identifier is unique per record.
eicu_meta <- make_metadata(eicu_raw, "eICU_external_validation")
subgroup_specs <- list(
    list(name = "Age group", column = "age_group", levels = c("<65", ">=65")),
    list(name = "Sex", column = "sex_group", levels = c("Female", "Male")),
    list(name = "Stroke phenotype", column = "stroke_group",
        levels = c("AIS", "ICH", "SAH", "Mixed/ambiguous"))
)
subgroup_rows <- list()
for (window in names(prediction_files)) {
    p <- utils::read.csv(prediction_files[[window]], check.names = FALSE,
        stringsAsFactors = FALSE)
    p <- p[p$model == "xgboost" & p$dataset == "eICU_external_validation", , drop = FALSE]
    p$group_id <- as.character(p$group_id)
    p$truth <- as.character(p$truth)
    p$probability <- suppressWarnings(as.numeric(p$probability))
    p <- merge(p, eicu_meta, by = "group_id", all.x = FALSE, all.y = FALSE, sort = FALSE)
    for (spec in subgroup_specs) {
        for (level in spec$levels) {
            d <- p[p[[spec$column]] == level, , drop = FALSE]
            n <- nrow(d)
            events <- sum(d$truth == "yes", na.rm = TRUE)
            if (n < 20L || events < 5L || (n - events) < 5L) next
            b <- bootstrap_rows(d, safe_auc, reps = 1000L,
                seed = 20260901L + match(window, names(prediction_files)) + nchar(level))
            subgroup_rows[[length(subgroup_rows) + 1L]] <- data.frame(
                window = window, dataset = "eICU_external_validation", model = "xgboost",
                subgroup = spec$name, level = level, n = n, events = events,
                auroc = unname(b$point[["auroc"]]),
                auroc_lower = unname(b$ci["lower", "auroc"]),
                auroc_upper = unname(b$ci["upper", "auroc"]),
                auprc = unname(b$point[["auprc"]]),
                auprc_lower = unname(b$ci["lower", "auprc"]),
                auprc_upper = unname(b$ci["upper", "auprc"]),
                bootstrap_reps = 1000L, stringsAsFactors = FALSE
            )
        }
    }
}
subgroup_results <- do.call(rbind, subgroup_rows)
utils::write.csv(subgroup_results,
    file.path(publication_dir, "tableS2_subgroup_performance.csv"),
    row.names = FALSE)

subgroup_md <- c(
    "# Supplementary Table S2. eICU external validation performance by prespecified subgroup",
    "",
    "Results are for XGBoost. Confidence intervals are patient-level percentile bootstrap intervals based on 1,000 resamples. Subgroups with fewer than 20 records, fewer than 5 events, or fewer than 5 non-events were not estimated.",
    "",
    "| Window | Subgroup | Level | N | Events | AUROC (95% CI) | AUPRC (95% CI) |",
    "|---|---|---|---:|---:|---:|---:|",
    apply(subgroup_results, 1, function(r) sprintf(
        "| %s | %s | %s | %s | %s | %.3f (%.3f-%.3f) | %.3f (%.3f-%.3f) |",
        r[["window"]], r[["subgroup"]], r[["level"]], r[["n"]], r[["events"]],
        as.numeric(r[["auroc"]]), as.numeric(r[["auroc_lower"]]), as.numeric(r[["auroc_upper"]]),
        as.numeric(r[["auprc"]]), as.numeric(r[["auprc_lower"]]), as.numeric(r[["auprc_upper"]])
    ))
)
writeLines(subgroup_md, file.path(publication_dir, "tableS2_subgroup_performance.md"),
    useBytes = TRUE)

# Fixed thresholds are prespecified operating points, not optimized on eICU.
# The 0.10 threshold is the primary operating point; other thresholds are a
# sensitivity grid aligned with the decision-curve range.
thresholds <- c(0.05, 0.10, 0.20, 0.30, 0.50)
threshold_rows <- list()
for (window in names(prediction_files)) {
    p <- utils::read.csv(prediction_files[[window]], check.names = FALSE,
        stringsAsFactors = FALSE)
    p <- p[p$model == "xgboost" & p$dataset %in% c("MIMIC_internal_test", "eICU_external_validation"), , drop = FALSE]
    p$group_id <- as.character(p$group_id)
    p$truth <- as.character(p$truth)
    p$probability <- suppressWarnings(as.numeric(p$probability))
    for (dataset in unique(p$dataset)) {
        d <- p[p$dataset == dataset, , drop = FALSE]
        for (threshold in thresholds) {
            b <- threshold_bootstrap(d, threshold, reps = 1000L,
                seed = 20260901L + round(threshold * 100) + match(window, names(prediction_files)))
            point <- b$point
            threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
                window = window, dataset = dataset, model = "xgboost",
                threshold = threshold, primary_threshold = threshold == 0.10,
                n_records = nrow(d), patient_groups = length(unique(d$group_id)),
                events = sum(d$truth == "yes", na.rm = TRUE),
                sensitivity = point[["sensitivity"]], sensitivity_lower = b$ci["lower", "sensitivity"], sensitivity_upper = b$ci["upper", "sensitivity"],
                specificity = point[["specificity"]], specificity_lower = b$ci["lower", "specificity"], specificity_upper = b$ci["upper", "specificity"],
                ppv = point[["ppv"]], ppv_lower = b$ci["lower", "ppv"], ppv_upper = b$ci["upper", "ppv"],
                npv = point[["npv"]], npv_lower = b$ci["lower", "npv"], npv_upper = b$ci["upper", "npv"],
                f1 = point[["f1"]], f1_lower = b$ci["lower", "f1"], f1_upper = b$ci["upper", "f1"],
                bootstrap_reps = 1000L, stringsAsFactors = FALSE
            )
        }
    }
}
threshold_results <- do.call(rbind, threshold_rows)
utils::write.csv(threshold_results,
    file.path(publication_dir, "table3_threshold_metrics.csv"), row.names = FALSE)

fmt_ci <- function(x, lo, hi) sprintf("%.3f (%.3f-%.3f)", x, lo, hi)
threshold_md <- c(
    "# Table 3. Fixed-threshold operating characteristics of XGBoost",
    "",
    "Thresholds were prespecified and were not optimized using eICU external-validation results. The 0.10 threshold is the primary operating point; the remaining thresholds form a sensitivity grid. Confidence intervals are patient-level percentile bootstrap intervals based on 1,000 resamples.",
    "",
    "| Window | Dataset | Threshold | N | Events | Sensitivity | Specificity | PPV | NPV | F1 |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    apply(threshold_results, 1, function(r) sprintf(
        "| %s | %s | %.2f%s | %s | %s | %s | %s | %s | %s | %s |",
        r[["window"]], r[["dataset"]], as.numeric(r[["threshold"]]),
        ifelse(r[["primary_threshold"]] == "TRUE", "*", ""), r[["n_records"]], r[["events"]],
        fmt_ci(as.numeric(r[["sensitivity"]]), as.numeric(r[["sensitivity_lower"]]), as.numeric(r[["sensitivity_upper"]])),
        fmt_ci(as.numeric(r[["specificity"]]), as.numeric(r[["specificity_lower"]]), as.numeric(r[["specificity_upper"]])),
        fmt_ci(as.numeric(r[["ppv"]]), as.numeric(r[["ppv_lower"]]), as.numeric(r[["ppv_upper"]])),
        fmt_ci(as.numeric(r[["npv"]]), as.numeric(r[["npv_lower"]]), as.numeric(r[["npv_upper"]])),
        fmt_ci(as.numeric(r[["f1"]]), as.numeric(r[["f1_lower"]]), as.numeric(r[["f1_upper"]]))
    )),
    "",
    "* Primary prespecified operating threshold: 0.10."
)
writeLines(threshold_md, file.path(publication_dir, "table3_threshold_metrics.md"),
    useBytes = TRUE)

readme <- c(
    "# Subgroup and threshold analysis notes",
    "",
    "- Subgroup analysis: primary XGBoost model in eICU external validation, stratified by age (<65 vs >=65), sex, and collapsed stroke phenotype.",
    "- Subgroup confidence intervals: patient-level percentile bootstrap, 1,000 resamples; unstable groups (<20 records, <5 events, or <5 non-events) are omitted.",
    "- Threshold analysis: XGBoost predictions from the MIMIC internal test set and eICU external validation, at prespecified thresholds 0.05, 0.10, 0.20, 0.30, and 0.50.",
    "- The 0.10 threshold is the primary operating point. No threshold was selected or optimized using eICU results.",
    "- Threshold confidence intervals also use patient-level percentile bootstrap with 1,000 resamples.",
    "- These analyses estimate predictive performance and operating characteristics; they do not establish causal effects or clinical benefit."
)
writeLines(readme, file.path(publication_dir, "subgroup_threshold_analysis_readme.md"), useBytes = TRUE)

cat(sprintf("SUBGROUP_THRESHOLD_OK subgroup_rows=%d; threshold_rows=%d; primary_threshold=0.10; bootstrap_reps=1000\n",
    nrow(subgroup_results), nrow(threshold_results)))
