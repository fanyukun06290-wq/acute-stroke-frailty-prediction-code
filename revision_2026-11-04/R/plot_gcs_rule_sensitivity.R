Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "gcs_rule_sensitivity_results")
labels <- c(allow_reverse_recording = "Allow reversed recording order", first_valid = "First valid score", last_valid = "Last valid score",
    motor_requires_complete = "Require concurrent full GCS", airway_60 = "Airway lookback 60 min", airway_240 = "Airway lookback 240 min")
changes <- read.csv(file.path(out, "score_changes.csv"))
x <- changes[changes$score_type == "gcs_motor" & changes$variant %in% names(labels)[1:4], ]
pieces <- do.call(rbind, lapply(c("changed_among_both_observed", "newly_missing", "newly_available"), function(name) {
    a <- x; a$kind <- name; a$percent <- 100*a[[name]]/a$records; a
}))
pieces$rule <- factor(unname(labels[pieces$variant]), levels = rev(unname(labels[1:4])))
pieces$database <- factor(pieces$database, levels = c("mimic", "eicu"), labels = c("MIMIC", "eICU"))
pieces$window <- factor(paste0(pieces$window_hours, "-hour landmark"), levels = c("24-hour landmark", "6-hour landmark"))
theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.spacing = grid::unit(1, "lines"),
    plot.caption = element_text(hjust = 0, size = 9), legend.position = "bottom")
fig <- ggplot(pieces, aes(percent, rule, fill = kind)) + geom_col(width = .65) + facet_grid(database ~ window) +
    scale_fill_manual(values = c(changed_among_both_observed = "#2166AC", newly_missing = "#B2182B", newly_available = "#4D9221"),
        labels = c(changed_among_both_observed = "Changed score", newly_missing = "Newly missing", newly_available = "Newly available")) +
    labs(title = "Documented motor-score sensitivity to construction rules", x = "Affected records (% of unchanged eligible cohort)", y = NULL, fill = NULL,
        caption = "Score changes and availability switches are disjoint counts. These are measurement-rule sensitivities, not clinical correctness rankings.") + theme
dir.create(file.path(out, "figures"), showWarnings = FALSE)
ggsave(file.path(out, "figures/GCS_rule_score_changes.png"), fig, width = 12, height = 6.5, dpi = 300)
ggsave(file.path(out, "figures/GCS_rule_score_changes.svg"), fig, width = 12, height = 6.5, device = grDevices::svg)
p <- read.csv(file.path(out, "paired_input_rule_differences.csv"))
p <- p[p$dataset == "eICU_external_validation" & p$model == "omit_all_flags_BMI" & p$metric == "auroc", ]
p$rule <- factor(unname(labels[p$variant]), levels = rev(unname(labels)))
p$window <- factor(paste0(p$window_hours, "-hour landmark"), levels = c("24-hour landmark", "6-hour landmark"))
fig <- ggplot(p, aes(difference, rule)) + geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
    geom_segment(aes(x = lower, xend = upper, yend = rule), color = "#2166AC", linewidth = .7) +
    geom_point(color = "#2166AC", size = 2.2) + facet_wrap(~window, nrow = 1) +
    labs(title = "Frozen combination model: external AUROC changes under score substitution", x = "Paired AUROC difference from strict-baseline inputs", y = NULL,
        caption = paste("Fixed omit-all-flags B+M+I coefficients and preprocessing; no model was retrained. 95% conditional paired bootstrap intervals.",
            "\nChanges do not estimate performance after training with the new definition, or rank the clinical validity of the definitions.",
            "\nUnchanged predictions are retained as structural zero differences.")) + theme
ggsave(file.path(out, "figures/GCS_rule_fixed_model_AUROC_stress.png"), fig, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out, "figures/GCS_rule_fixed_model_AUROC_stress.svg"), fig, width = 12, height = 5.5, device = grDevices::svg)
cat("GCS_RULE_SENSITIVITY_PLOTS_WRITTEN\n")
