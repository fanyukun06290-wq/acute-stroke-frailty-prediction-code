# Draw calibration summaries from fixed predictions and cluster-bootstrap bin intervals.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else file.path(getwd(), "revision_2026-11-04")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""), .libPaths())))
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(root, "validation_results")
dir.create(file.path(out, "figures"), showWarnings = FALSE)
bins <- read.csv(file.path(out, "calibration_bins_intervals.csv"))
metrics <- read.csv(file.path(out, "metric_confidence_intervals.csv"))
pairs <- read.csv(file.path(out, "paired_model_differences.csv"))
primary <- function(x) x[(x$dataset == "MIMIC_internal_test" & x$scheme == "patient") |
    (x$dataset == "eICU_external_validation" & x$scheme == "hospital_patient_component"), , drop = FALSE]
write.csv(primary(metrics), file.path(out, "primary_metric_intervals.csv"), row.names = FALSE)
write.csv(primary(pairs), file.path(out, "primary_paired_differences.csv"), row.names = FALSE)
write.csv(primary(bins), file.path(out, "primary_calibration_bins.csv"), row.names = FALSE)
external <- primary(bins)
external <- external[external$dataset == "eICU_external_validation" & external$model != "constant_training_prevalence", ]
external$window_label <- factor(ifelse(external$window_hours == 24, "24-hour landmark", "6-hour landmark"),
    levels = c("24-hour landmark", "6-hour landmark"))
model_names <- c(xgboost = "XGBoost", random_forest = "Random forest", elastic_net = "Elastic net",
    clinical_reference = "Clinical reference\n(unpenalized)")
external$model_label <- factor(unname(model_names[external$model]), levels = unname(model_names))
theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2"),
    plot.caption = element_text(hjust = 0, size = 9), plot.title = element_text(size = 14),
    panel.spacing.x = grid::unit(1.1, "lines"), panel.spacing.y = grid::unit(1, "lines"))
external$label_x <- external$mean_predicted
external$label_y <- external$observed_rate + .025
low <- external$bin <= 5
external$label_y[low] <- .045 + .034 * external$bin[low]
main <- ggplot(external[external$model == "xgboost", ], aes(mean_predicted, observed_rate)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#777777") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .012, color = "#2166AC") +
    geom_line(color = "#2166AC", linewidth = .55) +
    geom_point(shape = 21, size = 2.4, fill = "white", color = "#2166AC") +
    geom_segment(aes(xend = label_x, yend = label_y - .009), color = "#9AAABA", linewidth = .25) +
    geom_text(aes(x = label_x, y = label_y, label = bin), size = 3, color = "#173D66") +
    facet_wrap(~ window_label, nrow = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, .2)) + scale_y_continuous(breaks = seq(0, 1, .2)) +
    labs(title = "External calibration of the fixed XGBoost models",
        x = "Mean predicted probability", y = "Observed event proportion",
        caption = paste("Numbers identify prediction-risk groups (up to 10 deciles; ties remain together).",
            "\nBars: 95% percentile intervals from 1,000 hospital/shared-patient block resamples.",
            "\nThe diagonal denotes ideal calibration. The two landmarks refer to different risk populations.")) + theme
ggsave(file.path(out, "figures", "Figure_4_calibration_REVISED.png"), main, width = 10, height = 5.7, dpi = 300)
ggsave(file.path(out, "figures", "Figure_4_calibration_REVISED.svg"), main, width = 10, height = 5.7,
    device = grDevices::svg)
supplement <- ggplot(external, aes(mean_predicted, observed_rate)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#777777") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .012, color = "#2166AC") +
    geom_line(color = "#2166AC", linewidth = .5) +
    geom_point(shape = 21, size = 1.8, fill = "white", color = "#2166AC") +
    facet_grid(model_label ~ window_label) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = "External calibration by fixed model and landmark",
        x = "Mean predicted probability", y = "Observed event proportion",
        caption = paste("Prediction-defined groups; 95% intervals from hospital/shared-patient block bootstrap.",
                        "\nThe provisional unpenalized clinical reference has a diagnosed missingness/separation problem.")) + theme
ggsave(file.path(out, "figures", "Supplement_calibration_all_models.png"), supplement,
    width = 10.5, height = 14, dpi = 300)
ggsave(file.path(out, "figures", "Supplement_calibration_all_models.svg"), supplement,
    width = 10.5, height = 14, device = grDevices::svg)
cat("PRIMARY_TABLES_AND_CALIBRATION_PLOTS_WRITTEN\n")
