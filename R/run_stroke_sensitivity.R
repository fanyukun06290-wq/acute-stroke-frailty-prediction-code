# Prespecified sensitivity analyses for the source-level QC stroke dataset.
# Each scenario refits XGBoost in MIMIC and evaluates the fixed model in eICU.
# No database is queried or modified.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- Sys.getenv("STROKE_PROJECT_DIR", unset = getwd())
source_mimic_path <- if (length(args) >= 1L) args[[1L]] else
    Sys.getenv("STROKE_MIMIC_CSV", unset = file.path(project_dir, "data", "mimic_analysis.csv"))
source_eicu_path <- if (length(args) >= 2L) args[[2L]] else
    Sys.getenv("STROKE_EICU_CSV", unset = file.path(project_dir, "data", "eicu_analysis.csv"))
out_root <- if (length(args) >= 3L) args[[3L]] else
    Sys.getenv("STROKE_SENSITIVITY_DIR", unset = file.path(project_dir, "results", "sensitivity"))
model_script <- Sys.getenv("STROKE_MODEL_SCRIPT", unset = file.path(project_dir, "R", "run_stroke_models.R"))

if (!all(file.exists(c(source_mimic_path, source_eicu_path, model_script)))) {
    stop("Sensitivity-analysis input or model script is missing.")
}

Sys.setenv(
    STROKE_SKIP_RUN = "1",
    STROKE_MODEL_SET = "xgboost",
    STROKE_TUNING_GRID_SIZE = Sys.getenv("STROKE_TUNING_GRID_SIZE", unset = "5"),
    STROKE_R_LIBRARY = Sys.getenv(
        "STROKE_R_LIBRARY", unset = ""
    )
)
source(model_script, local = TRUE)

mimic_raw <- utils::read.csv(source_mimic_path, check.names = FALSE)
eicu_raw <- utils::read.csv(source_eicu_path, check.names = FALSE)

labs <- c(
    "glucose", "sodium", "potassium", "chloride", "bicarbonate", "bun",
    "creatinine", "calcium", "wbc", "platelets", "hemoglobin", "hematocrit"
)
vitals <- c("heart_rate", "respiratory_rate", "spo2", "sbp", "dbp", "map")

complete_case_window <- function(raw, window) {
    suffix <- if (window == "landmark") "mean_24h" else "first_6h"
    outcome <- if (window == "landmark") {
        "death_24h_to_30d_inhospital"
    } else "death_30d_inhospital"
    cols <- c(
        "age_years", "sex", "frailty_comorbidity_count", "outcome_missing_flag",
        outcome, paste0(labs, "_", suffix), paste0(vitals, "_", suffix)
    )
    keep <- stats::complete.cases(raw[, cols, drop = FALSE])
    if (window == "landmark") {
        eligible <- to01(raw$landmark_24h_eligible) == 1L
        keep <- keep & eligible
    }
    raw[keep, , drop = FALSE]
}

strict_only <- function(raw) {
    raw[raw$audit_phenotype %in% c("AIS_STRICT", "ICH_STRICT", "SAH"), , drop = FALSE]
}

scenarios <- list(
    mixed_excluded = list(
        description = "Exclude mixed, ambiguous, and CVST records; retain strict AIS/ICH/SAH.",
        include_missingness = TRUE,
        filter = function(raw, window) strict_only(raw)
    ),
    complete_case = list(
        description = "Restrict each window to complete cases for all core predictors and outcome fields.",
        include_missingness = TRUE,
        filter = complete_case_window
    ),
    no_missing_indicators = list(
        description = "Retain the primary cohort but remove explicit laboratory/vital-sign missingness indicators.",
        include_missingness = FALSE,
        filter = function(raw, window) raw
    )
)

all_metrics <- list()
for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    scenario_dir <- file.path(out_root, scenario_name)
    dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(
        c(paste0("# ", scenario_name), "", scenario$description),
        file.path(scenario_dir, "README.md"), useBytes = TRUE
    )
    Sys.setenv(
        STROKE_INCLUDE_MISSINGNESS = ifelse(scenario$include_missingness, "1", "0")
    )
    for (window in c("landmark", "early")) {
        m <- scenario$filter(mimic_raw, window)
        e <- scenario$filter(eicu_raw, window)
        message("Running sensitivity scenario=", scenario_name, " window=", window,
            " MIMIC=", nrow(m), " eICU=", nrow(e))
        run_window(m, e, window, scenario_dir)
        metric_file <- file.path(scenario_dir, paste0("model_metrics_", window, ".csv"))
        x <- utils::read.csv(metric_file, check.names = FALSE)
        x$scenario <- scenario_name
        all_metrics[[paste(scenario_name, window, sep = "||")]] <- x
    }
}

summary <- do.call(rbind, all_metrics)
summary <- summary[, c("scenario", "window", "model", "dataset",
    "auroc", "auprc", "brier_score", "calibration_intercept", "calibration_slope")]
utils::write.csv(summary, file.path(out_root, "sensitivity_model_metrics.csv"), row.names = FALSE)
utils::write.csv(
    subset(summary, model == "xgboost" & dataset == "eICU_external_validation"),
    file.path(out_root, "sensitivity_external_xgboost.csv"), row.names = FALSE
)

Sys.unsetenv(c("STROKE_SKIP_RUN", "STROKE_INCLUDE_MISSINGNESS", "STROKE_MODEL_SET"))
cat(sprintf("SENSITIVITY_OK scenarios=%d windows=2 models=xgboost output=%s\n",
    length(scenarios), normalizePath(out_root, mustWork = FALSE)))
