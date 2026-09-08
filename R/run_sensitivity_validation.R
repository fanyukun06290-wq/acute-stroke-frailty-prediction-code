# Fast patient-level bootstrap for prespecified sensitivity predictions.
# Focuses on AUROC, AUPRC, and Brier score; no model refitting or database access.

options(stringsAsFactors = FALSE)
project_library <- Sys.getenv("STROKE_R_LIBRARY", unset = "")
if (nzchar(project_library) && dir.exists(project_library)) {
    .libPaths(unique(c(project_library, .libPaths())))
}
suppressPackageStartupMessages(library(yardstick))

args <- commandArgs(trailingOnly = TRUE)
project_dir <- Sys.getenv("STROKE_PROJECT_DIR", unset = getwd())
root <- if (length(args) >= 1L) args[[1L]] else
    Sys.getenv("STROKE_SENSITIVITY_DIR", unset = file.path(project_dir, "results", "sensitivity"))
bootstrap_reps <- suppressWarnings(as.integer(
    Sys.getenv("STROKE_BOOTSTRAP_REPS", unset = "1000")
))
if (is.na(bootstrap_reps) || bootstrap_reps < 100L) bootstrap_reps <- 1000L
scenario_names <- c("mixed_excluded", "complete_case", "no_missing_indicators")
windows <- c("landmark", "early")

score_fast <- function(truth, probability) {
    truth <- factor(as.character(truth), levels = c("no", "yes"))
    probability <- pmin(pmax(as.numeric(probability), 1e-6), 1 - 1e-6)
    c(
        auroc = as.numeric(roc_auc_vec(truth, probability, event_level = "second")),
        auprc = as.numeric(pr_auc_vec(truth, probability, event_level = "second")),
        brier_score = mean((as.integer(truth == "yes") - probability)^2)
    )
}

bootstrap_one <- function(predictions, reps = 1000L, seed = 20260901L) {
    predictions$group_id <- as.character(predictions$group_id)
    predictions$truth <- factor(predictions$truth, levels = c("no", "yes"))
    predictions$probability <- as.numeric(predictions$probability)
    groups <- unique(predictions$group_id)
    index_by_group <- split(seq_len(nrow(predictions)), predictions$group_id)
    point <- score_fast(predictions$truth, predictions$probability)
    set.seed(seed)
    boot <- matrix(NA_real_, nrow = reps, ncol = length(point))
    colnames(boot) <- names(point)
    for (b in seq_len(reps)) {
        sampled <- sample(groups, length(groups), replace = TRUE)
        idx <- unlist(index_by_group[sampled], use.names = FALSE)
        boot[b, ] <- score_fast(predictions$truth[idx], predictions$probability[idx])
    }
    out <- do.call(rbind, lapply(names(point), function(metric) {
        data.frame(
            metric = metric,
            estimate = unname(point[[metric]]),
            ci_low = unname(stats::quantile(boot[, metric], 0.025, na.rm = TRUE)),
            ci_high = unname(stats::quantile(boot[, metric], 0.975, na.rm = TRUE)),
            bootstrap_reps = reps,
            n_records = nrow(predictions),
            patient_groups = length(groups),
            events = sum(predictions$truth == "yes"),
            stringsAsFactors = FALSE
        )
    }))
    out
}

rows <- list()
for (scenario in scenario_names) {
    for (window in windows) {
        f <- file.path(root, scenario, paste0("predictions_", window, ".csv"))
        if (!file.exists(f)) stop("Missing sensitivity predictions: ", f)
        p <- utils::read.csv(f, check.names = FALSE)
        p <- p[p$model == "xgboost" & p$dataset == "eICU_external_validation", , drop = FALSE]
        x <- bootstrap_one(p, reps = bootstrap_reps,
            seed = 20260901L + match(scenario, scenario_names) + match(window, windows))
        x$scenario <- scenario
        x$window <- window
        rows[[length(rows) + 1L]] <- x
    }
}
result <- do.call(rbind, rows)
result <- result[, c("scenario", "window", "metric", "estimate", "ci_low", "ci_high",
    "bootstrap_reps", "n_records", "patient_groups", "events")]
utils::write.csv(result, file.path(root, "sensitivity_external_bootstrap.csv"), row.names = FALSE)
cat(sprintf("SENSITIVITY_BOOTSTRAP_OK rows=%d reps=%d\n", nrow(result), bootstrap_reps))
