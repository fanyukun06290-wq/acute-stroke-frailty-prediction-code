# Diagnose the existing reference without changing or re-estimating any model.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R", "run_revision_models.R"), encoding = "UTF-8")
dir.create(file.path(root, "validation_results"), showWarnings = FALSE)
contract <- jsonlite::read_json(file.path(root, "model_feature_contract.json"), simplifyVector = TRUE)
summaries <- list()
for (hour in c(24L, 6L)) {
    tag <- paste0(hour, "h")
    mimic <- read_data(file.path(root, "private_data", "mimic_analysis_revision.csv"), contract, hour)
    external <- read_data(file.path(root, "private_data", "eicu_analysis_revision.csv"), contract, hour)
    registry <- read.csv(file.path(root, "private_splits", paste0("patients_", tag, ".csv")))
    keep <- mimic$metadata$group_id %in% registry$group_id[registry$partition == "train"]
    model <- readRDS(file.path(root, "private_models", paste0("clinical_reference_", tag, ".rds")))
    rec <- readRDS(file.path(root, "private_models", paste0("recipe_", tag, ".rds")))
    matrix <- as.matrix(bake(rec, new_data = external$x)); storage.mode(matrix) <- "double"
    probability <- predict_clinical(model, matrix)
    low <- probability < 1e-6
    flags <- grep("^(heart_rate|respiratory_rate|spo2|sbp|dbp|map)_missing", names(mimic$x), value = TRUE)
    detailed <- list()
    for (flag in flags) {
        observed <- mimic$x[[flag]][keep]
        y <- mimic$y[keep]
        positive <- observed == 1
        detailed[[flag]] <- data.frame(window_hours = hour, flag = flag,
            training_flag_positive_records = sum(positive),
            training_flag_positive_events = sum(y[positive]),
            training_flag_negative_records = sum(!positive),
            training_flag_negative_events = sum(y[!positive]),
            external_flag_positive_records = sum(external$x[[flag]] == 1),
            external_low_predictions_with_flag = sum(low & external$x[[flag]] == 1),
            external_low_predictions_total = sum(low),
            no_training_events_in_flag_positive_stratum = sum(positive) > 0L && sum(y[positive]) == 0L)
    }
    details <- do.call(rbind, detailed)
    write.csv(details, file.path(root, "private_cache", paste0("clinical_missingness_diagnostics_", tag, ".csv")), row.names = FALSE)
    write.csv(data.frame(encoded_term = names(model$coefficients), coefficient = unname(model$coefficients)),
        file.path(root, "private_cache", paste0("clinical_coefficients_", tag, ".csv")), row.names = FALSE)
    flagged <- details$flag[details$no_training_events_in_flag_positive_stratum]
    summaries[[tag]] <- list(model_changed = FALSE, model_retrained = FALSE,
        external_records = length(probability), external_low_probability_records = sum(low),
        external_events_in_low_probability_records = sum(external$y[low]),
        zero_training_event_missingness_flags = flagged,
        low_prediction_overlap = as.list(setNames(details$external_low_predictions_with_flag, details$flag)),
        interpretation = if (length(flagged))
            "A missingness-positive training stratum has no events, supporting quasi-separation and unstable unpenalized coefficients; assess workflow shift rather than infer protection from missingness"
        else "No zero-event positive missingness stratum identified by this bounded diagnostic",
        limitation = "This routine-vital-sign comparator is not a validated stroke/ICU clinical score; do not use weak reference performance to claim clinical superiority")
}
write_json(summaries, file.path(root, "validation_results", "clinical_reference_diagnostic_summary.json"))
cat(jsonlite::toJSON(summaries, pretty = TRUE, auto_unbox = TRUE), "\n")
