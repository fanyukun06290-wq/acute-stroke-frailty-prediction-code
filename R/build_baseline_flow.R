# Build publication-oriented Table 1 and cohort flow diagram.
# Reads the already exported analysis CSV files; does not query or modify databases.

project_dir <- Sys.getenv(
    "STROKE_PROJECT_DIR", unset = getwd()
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
if (!file.exists(mimic_csv) || !file.exists(eicu_csv)) {
    stop("未找到 MIMIC/eICU analysis CSV 文件。")
}

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

prepare_landmark <- function(raw) {
    eligible <- to01(raw$landmark_24h_eligible) == 1L
    outcome_available <- to01(raw$outcome_missing_flag) == 0L
    keep <- eligible & outcome_available
    raw <- raw[keep, , drop = FALSE]
    raw$death <- to01(raw$death_24h_to_30d_inhospital)
    raw$stroke_group <- collapse_stroke_group(raw$audit_phenotype)
    raw$sex_label <- trimws(as.character(raw$sex))
    raw$sex_label[is.na(raw$sex_label) |
        tolower(raw$sex_label) %in% c("", "unknown", "unk", "u", "na")] <- "Unknown"
    raw
}

mimic <- prepare_landmark(mimic_raw)
eicu <- prepare_landmark(eicu_raw)

fmt_cont <- function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    ok <- is.finite(x)
    if (!any(ok)) return("NA")
    q <- stats::quantile(x[ok], probs = c(0.25, 0.50, 0.75), names = FALSE)
    sprintf("%.1f [%.1f-%.1f]; missing %.1f%%", q[[2]], q[[1]], q[[3]],
        100 * mean(!ok))
}

fmt_cat <- function(x, level) {
    x <- as.character(x)
    x[is.na(x) | trimws(x) == ""] <- "Unknown"
    n <- sum(x == level)
    sprintf("%d (%.1f%%)", n, 100 * n / length(x))
}

groups <- list(
    MIMIC_overall = mimic,
    MIMIC_survivors = mimic[mimic$death == 0L, , drop = FALSE],
    MIMIC_deaths = mimic[mimic$death == 1L, , drop = FALSE],
    eICU_overall = eicu,
    eICU_survivors = eicu[eicu$death == 0L, , drop = FALSE],
    eICU_deaths = eicu[eicu$death == 1L, , drop = FALSE]
)

continuous_vars <- c(
    age_years = "Age, years",
    frailty_comorbidity_count = "Frailty-related comorbidity count",
    glucose_mean_24h = "Glucose, 24-h mean",
    sodium_mean_24h = "Sodium, 24-h mean",
    potassium_mean_24h = "Potassium, 24-h mean",
    chloride_mean_24h = "Chloride, 24-h mean",
    bicarbonate_mean_24h = "Bicarbonate, 24-h mean",
    bun_mean_24h = "BUN, 24-h mean",
    creatinine_mean_24h = "Creatinine, 24-h mean",
    calcium_mean_24h = "Total calcium, 24-h mean",
    wbc_mean_24h = "White blood cell count, 24-h mean",
    platelets_mean_24h = "Platelets, 24-h mean",
    hemoglobin_mean_24h = "Hemoglobin, 24-h mean",
    hematocrit_mean_24h = "Hematocrit, 24-h mean",
    heart_rate_mean_24h = "Heart rate, 24-h mean",
    respiratory_rate_mean_24h = "Respiratory rate, 24-h mean",
    spo2_mean_24h = "SpO2, 24-h mean",
    sbp_mean_24h = "Systolic blood pressure, 24-h mean",
    dbp_mean_24h = "Diastolic blood pressure, 24-h mean",
    map_mean_24h = "Mean arterial pressure, 24-h mean"
)

categorical_vars <- list(
    list(name = "Sex", column = "sex_label", levels = c("Female", "Male", "Unknown"),
        recode = function(x) ifelse(x == "F", "Female", ifelse(x == "M", "Male", "Unknown"))),
    list(name = "Stroke phenotype", column = "stroke_group",
        levels = c("AIS", "ICH", "SAH", "Mixed/ambiguous"), recode = identity)
)

for (g in names(groups)) {
    groups[[g]]$sex_label <- categorical_vars[[1]]$recode(groups[[g]]$sex_label)
}

rows <- list()
add_row <- function(variable, category, values) {
    rows[[length(rows) + 1L]] <<- data.frame(
        variable = variable, category = category,
        MIMIC_overall = values[[1]], MIMIC_survivors = values[[2]],
        MIMIC_deaths = values[[3]], eICU_overall = values[[4]],
        eICU_survivors = values[[5]], eICU_deaths = values[[6]],
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

for (v in names(continuous_vars)) {
    add_row(continuous_vars[[v]], "",
        lapply(groups, function(d) fmt_cont(d[[v]])))
}
for (v in categorical_vars) {
    for (level in v$levels) {
        add_row(v$name, level,
            lapply(groups, function(d) fmt_cat(d[[v$column]], level)))
    }
}

table1 <- do.call(rbind, rows)
utils::write.csv(table1, file.path(publication_dir, "table1_landmark.csv"),
    row.names = FALSE, fileEncoding = "UTF-8")

group_labels <- c(
    "MIMIC_overall" = sprintf("MIMIC overall (n=%d)", nrow(groups$MIMIC_overall)),
    "MIMIC_survivors" = sprintf("MIMIC survivors (n=%d)", nrow(groups$MIMIC_survivors)),
    "MIMIC_deaths" = sprintf("MIMIC deaths (n=%d)", nrow(groups$MIMIC_deaths)),
    "eICU_overall" = sprintf("eICU overall (n=%d)", nrow(groups$eICU_overall)),
    "eICU_survivors" = sprintf("eICU survivors (n=%d)", nrow(groups$eICU_survivors)),
    "eICU_deaths" = sprintf("eICU deaths (n=%d)", nrow(groups$eICU_deaths))
)
md_lines <- c(
    "# Table 1. Baseline characteristics of the 24-hour landmark modeling cohorts",
    "",
    "Continuous variables are shown as median [IQR] with the percentage of missing values; categorical variables are shown as n (%). The table uses records with an eligible 24-hour landmark and an available primary outcome.",
    "",
    paste0("| ", paste(c("Variable", "Category", unname(group_labels)), collapse = " | "), " |"),
    paste0("| ", paste(c("---", "---", rep("---:", length(group_labels))), collapse = " | "), " |"),
    apply(table1, 1, function(r) paste0("| ", paste(c(
        r[["variable"]], r[["category"]], unname(r[names(group_labels)])
    ), collapse = " | "), " |"))
)
writeLines(md_lines, file.path(publication_dir, "table1_landmark.md"), useBytes = TRUE)

flow_rows <- list()
add_flow <- function(dataset, stage, label, n) {
    flow_rows[[length(flow_rows) + 1L]] <<- data.frame(
        dataset = dataset, stage = stage, label = label, n = as.integer(n),
        stringsAsFactors = FALSE
    )
}
for (entry in list(list(name = "MIMIC-IV", raw = mimic_raw),
                   list(name = "eICU-CRD", raw = eicu_raw))) {
    raw <- entry$raw
    eligible <- to01(raw$landmark_24h_eligible) == 1L
    outcome_available <- to01(raw$outcome_missing_flag) == 0L
    model_ready <- eligible & outcome_available
    death_primary <- to01(raw$death_24h_to_30d_inhospital)
    death_early <- to01(raw$death_30d_inhospital)
    add_flow(entry$name, "total", "Analysis records", nrow(raw))
    add_flow(entry$name, "landmark", "24-h landmark eligible", sum(eligible, na.rm = TRUE))
    add_flow(entry$name, "outcome", "Primary outcome available", sum(eligible & outcome_available, na.rm = TRUE))
    add_flow(entry$name, "primary", "Primary model-ready cohort", sum(model_ready, na.rm = TRUE))
    add_flow(entry$name, "primary_event", "Primary landmark events", sum(death_primary[model_ready] == 1L, na.rm = TRUE))
    add_flow(entry$name, "primary_nonevent", "Primary landmark non-events", sum(death_primary[model_ready] == 0L, na.rm = TRUE))
    add_flow(entry$name, "early", "Early model-ready cohort", sum(outcome_available, na.rm = TRUE))
    add_flow(entry$name, "early_event", "Early-window events", sum(death_early[outcome_available] == 1L, na.rm = TRUE))
}
flow <- do.call(rbind, flow_rows)
utils::write.csv(flow, file.path(publication_dir, "flow_diagram_counts.csv"),
    row.names = FALSE, fileEncoding = "UTF-8")

box <- function(x, y, label, fill = "#EAF2F8") {
    rect(x - 0.32, y - 0.055, x + 0.32, y + 0.055,
        col = fill, border = "#2C3E50", lwd = 1.5)
    text(x, y, label, cex = 0.78)
}
arrow <- function(y1, y2) arrows(0.5, y1 - 0.07, 0.5, y2 + 0.07,
    length = 0.08, lwd = 1.2, col = "#566573")

png(file.path(publication_dir, "flow_diagram.png"), width = 2200, height = 1500, res = 220)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 1.8, 0.5))
for (entry in c("MIMIC-IV", "eICU-CRD")) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    dat <- flow[flow$dataset == entry, , drop = FALSE]
    getn <- function(stage) dat$n[match(stage, dat$stage)]
    title(main = entry, cex.main = 1.15)
    box(0.5, 0.86, sprintf("Analysis records\nn=%d", getn("total")))
    box(0.5, 0.62, sprintf("24-h landmark eligible\nn=%d", getn("landmark")))
    box(0.5, 0.38, sprintf("Primary model-ready cohort\nn=%d", getn("primary")), fill = "#D5F5E3")
    arrow(0.86, 0.62); arrow(0.62, 0.38)
    text(0.82, 0.50, sprintf("Excluded: missing outcome n=%d",
        getn("landmark") - getn("outcome")), cex = 0.62, col = "#566573")
    text(0.5, 0.22, sprintf("Landmark events: %d\nNon-events: %d",
        getn("primary_event"), getn("primary_nonevent")), cex = 0.82)
    text(0.5, 0.08, sprintf("Early model-ready: %d; events: %d",
        getn("early"), getn("early_event")), cex = 0.72, col = "#566573")
}
par(oldpar)
dev.off()

cat(sprintf("TABLE1_FLOW_OK table1_nrows=%d; MIMIC_primary=%d/%d; eICU_primary=%d/%d\n",
    nrow(table1), nrow(groups$MIMIC_overall), sum(mimic$death == 1L),
    nrow(groups$eICU_overall), sum(eicu$death == 1L)))
