# Publication-style exploratory increment plot, separate from the submitted manuscript.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "neuro_reference_results")
p <- read.csv(file.path(out, "paired_differences.csv"))
p <- p[p$dataset == "eICU_external_validation" & p$metric == "auroc" & p$neuro == "motor", ]
labels <- c(BM_minus_B = "Add metabolic/renal to reference", BI_minus_B = "Add hematologic to reference",
    BMI_minus_B = "Add both domains to reference", BMI_minus_BM = "Add hematologic after metabolic/renal",
    BMI_minus_BI = "Add metabolic/renal after hematologic")
p$comparison <- factor(unname(labels[p$contrast]), levels = rev(unname(labels)))
p$window <- factor(paste0(p$window_hours, "-hour landmark"), levels = c("24-hour landmark", "6-hour landmark"))
p$population_label <- factor(ifelse(p$population == "all_eligible", "All eligible records", "Motor score observed"),
    levels = c("All eligible records", "Motor score observed"))
fig <- ggplot(p, aes(difference, comparison)) + geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
    geom_segment(aes(x = lower, xend = upper, yend = comparison), color = "#2166AC", linewidth = .7) +
    geom_point(color = "#2166AC", size = 2) + facet_grid(population_label ~ window) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
        panel.spacing = grid::unit(1, "lines")) +
    labs(x = "Paired AUROC difference (augmented minus reference)", y = NULL,
        title = "External increments over the clinical/motor reference",
        caption = paste("Points and 95% percentile intervals from 1,000 hospital/shared-patient block resamples.",
            "\nAll models use ridge logistic regression; exploratory comparisons without multiplicity adjustment.",
            "\nObserved-score panels concern a selected population; 24-hour all-record calibration is compromised by missingness shift."))
dir.create(file.path(out, "figures"), showWarnings = FALSE)
ggsave(file.path(out, "figures/Neurological_domain_AUROC_increments.png"), fig, width = 12, height = 6.5, dpi = 300)
ggsave(file.path(out, "figures/Neurological_domain_AUROC_increments.svg"), fig, width = 12, height = 6.5, device = grDevices::svg)
cat("NEURO_INCREMENT_PLOT_WRITTEN\n")
