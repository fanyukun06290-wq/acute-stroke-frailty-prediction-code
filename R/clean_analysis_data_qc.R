# Predefined, outcome-blind numerical quality control for the exported analysis CSVs.
# The original CSV files are never overwritten. Invalid values are set to NA and
# the corresponding missingness flags are updated. This is a conservative interim
# QC step; source-level SQL/item/unit auditing remains required.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- Sys.getenv("STROKE_PROJECT_DIR", unset = getwd())
input_dir <- if (length(args) >= 1L) args[[1L]] else
    Sys.getenv("STROKE_INPUT_DIR", unset = file.path(project_dir, "data"))
output_dir <- if (length(args) >= 2L) args[[2L]] else
    Sys.getenv("STROKE_QC_DIR", unset = file.path(project_dir, "results", "qc"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rules <- data.frame(
    variable = c(
        "glucose", "sodium", "potassium", "chloride", "bicarbonate", "bun",
        "creatinine", "calcium", "wbc", "platelets", "hemoglobin", "hematocrit",
        "heart_rate", "respiratory_rate", "spo2", "sbp", "dbp", "map"
    ),
    lower = c(20, 100, 1.5, 50, 5, 1, 0.1, 4, 0.1, 10, 2, 5, 20, 4, 50, 50, 20, 30),
    upper = c(1000, 200, 10, 150, 60, 300, 30, 15, 300, 2000, 25, 75, 250, 100, 100, 300, 200, 250),
    rationale = rep("Outcome-blind broad physiologic screening bound; confirm against source units/item dictionary before final analysis.", 18),
    stringsAsFactors = FALSE
)

clean_one <- function(path, dataset) {
    dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    log_rows <- list()
    k <- 1L
    for (i in seq_len(nrow(rules))) {
        v <- rules$variable[[i]]
        for (window in c("first_6h", "mean_24h")) {
            value_col <- paste0(v, "_", window)
            flag_col <- paste0(v, "_missing_", sub("^.*_", "", window))
            if (!value_col %in% names(dat)) next
            value <- suppressWarnings(as.numeric(as.character(dat[[value_col]])))
            bad <- !is.na(value) & (value < rules$lower[[i]] | value > rules$upper[[i]])
            invalid_n <- sum(bad)
            if (invalid_n > 0L) {
                value[bad] <- NA_real_
                dat[[value_col]] <- value
                if (flag_col %in% names(dat)) {
                    flag <- suppressWarnings(as.numeric(as.character(dat[[flag_col]])))
                    # A 24-hour missing flag is only defined for landmark-eligible rows.
                    flag[bad & (window == "first_6h" | dat$landmark_24h_eligible == 1)] <- 1
                    dat[[flag_col]] <- flag
                }
            }
            log_rows[[k]] <- data.frame(
                dataset = dataset,
                variable = v,
                window = window,
                lower_bound = rules$lower[[i]],
                upper_bound = rules$upper[[i]],
                invalid_n = invalid_n,
                invalid_pct_all_rows = round(100 * invalid_n / nrow(dat), 4),
                stringsAsFactors = FALSE
            )
            k <- k + 1L
        }
    }
    clean_path <- file.path(output_dir, paste0(tolower(dataset), "_analysis_qc.csv"))
    log_path <- file.path(output_dir, paste0(tolower(dataset), "_qc_log.csv"))
    utils::write.csv(dat, clean_path, row.names = FALSE, na = "")
    utils::write.csv(do.call(rbind, log_rows), log_path, row.names = FALSE)
    list(data = dat, log = do.call(rbind, log_rows), path = clean_path, log_path = log_path)
}

mimic <- clean_one(file.path(input_dir, "mimic_analysis.csv"), "MIMIC")
eicu <- clean_one(file.path(input_dir, "eicu_analysis.csv"), "eICU")
qc_log <- rbind(mimic$log, eicu$log)
utils::write.csv(qc_log, file.path(output_dir, "qc_log_all.csv"), row.names = FALSE)

summary_lines <- c(
    "# Outcome-blind numerical QC summary",
    "",
    "This interim QC does not overwrite the original analysis CSVs. Values outside the predefined broad screening bounds were set to missing, and matching missingness flags were updated. The bounds are not a substitute for source-level unit/item verification.",
    "",
    sprintf("MIMIC invalid cells set to missing: %d", sum(mimic$log$invalid_n)),
    sprintf("eICU invalid cells set to missing: %d", sum(eicu$log$invalid_n)),
    "",
    "| Dataset | Variable | Window | Bound | Invalid values |",
    "|---|---|---|---|---:|"
)
hits <- qc_log[qc_log$invalid_n > 0, , drop = FALSE]
if (nrow(hits) > 0L) {
    summary_lines <- c(summary_lines, vapply(seq_len(nrow(hits)), function(i) {
        sprintf("| %s | %s | %s | %s–%s | %d |", hits$dataset[[i]], hits$variable[[i]],
            hits$window[[i]], hits$lower_bound[[i]], hits$upper_bound[[i]], hits$invalid_n[[i]])
    }, character(1)))
} else {
    summary_lines <- c(summary_lines, "| None | — | — | — | 0 |")
}
summary_lines <- c(
    summary_lines,
    "",
    "## Required before final reporting",
    "",
    "1. Verify every bound against the source variable dictionary and units.",
    "2. Prefer raw-observation filtering before aggregation; this interim CSV-level QC cannot reconstruct the valid raw values that contributed to an already contaminated mean.",
    "3. Rerun the complete modeling and validation pipeline on the QC data, and treat the pre-QC model results as provisional."
)
writeLines(summary_lines, file.path(output_dir, "qc_summary.md"), useBytes = TRUE)
cat(sprintf("QC_OK mimic_invalid_cells=%d eicu_invalid_cells=%d output_dir=%s\n",
    sum(mimic$log$invalid_n), sum(eicu$log$invalid_n), output_dir))
