Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
suppressPackageStartupMessages(library(ggplot2))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "missingness_sensitivity_results")
old <- read.csv(file.path(root, "neuro_reference_results/missingness_diagnostics.csv"))
old <- old[old$window_hours == 24 & old$dataset == "eICU_external_validation" & old$model == "motor_B", ]
new <- read.csv(file.path(out, "missingness_diagnostics.csv"))
new <- new[new$window_hours == 24 & new$dataset == "eICU_external_validation" & new$domain == "B", ]
plotdata <- rbind(data.frame(strategy = "original_flags", missing = old$neuro_missing,
    mean_predicted = old$mean_predicted, observed_rate = old$event_prevalence),
    data.frame(strategy = new$strategy, missing = new$neuro_missing,
    mean_predicted = new$mean_predicted, observed_rate = new$observed_rate))
strategy_names <- c(original_flags = "Original flags", omit_neuro_flag = "Omit neurological flag", omit_all_flags = "Omit all flags")
plotdata$strategy_label <- factor(unname(strategy_names[plotdata$strategy]), levels = unname(strategy_names))
plotdata$group <- ifelse(plotdata$missing, "Motor score missing (n = 420; deaths = 87)", "Motor score observed (n = 3,310; deaths = 401)")
theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9), panel.spacing = grid::unit(1, "lines"))
risk <- ggplot(plotdata, aes(strategy_label, mean_predicted, group = 1)) +
    geom_hline(aes(yintercept = observed_rate), color = "#B2182B", linetype = "dashed", linewidth = .6) +
    geom_line(color = "#2166AC", linewidth = .6) + geom_point(color = "#2166AC", size = 2.4) +
    facet_wrap(~group, nrow = 1) + scale_y_continuous(limits = c(0, 1)) +
    labs(title = "24-hour clinical/motor reference: external risk estimates", x = NULL,
        y = "Mean predicted risk or observed event rate",
        caption = paste("Blue: mean predicted risk; dashed red: observed in-hospital event rate. Descriptive means; no intervals shown.",
            "\nPost-hoc sensitivity after external diagnosis: omission reduces extreme overprediction but leaves underprediction in the missing-score group.")) + theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
dir.create(file.path(out, "figures"), showWarnings = FALSE)
ggsave(file.path(out, "figures/Missingness_risk_diagnostic.png"), risk, width = 11, height = 5.3, dpi = 300)
ggsave(file.path(out, "figures/Missingness_risk_diagnostic.svg"), risk, width = 11, height = 5.3, device = grDevices::svg)
old <- read.csv(file.path(root, "neuro_reference_results/paired_differences.csv"))
old <- old[old$dataset == "eICU_external_validation" & old$neuro == "motor" & old$contrast == "BMI_minus_B" & old$metric == "auroc", ]
old <- old[, c("window_hours", "population", "difference", "lower", "upper")]; old$strategy <- "original_flags"
new <- read.csv(file.path(out, "paired_domain_differences.csv"))
new <- new[new$dataset == "eICU_external_validation" & new$contrast == "BMI_minus_B" & new$metric == "auroc", ]
d <- rbind(old, new[, names(old)])
d$strategy_label <- factor(unname(strategy_names[d$strategy]), levels = rev(unname(strategy_names)))
d$window <- factor(paste0(d$window_hours, "-hour landmark"), levels = c("24-hour landmark", "6-hour landmark"))
d$population_label <- factor(ifelse(d$population == "all_eligible", "All eligible records", "Motor score observed"),
    levels = c("All eligible records", "Motor score observed"))
increments <- ggplot(d, aes(difference, strategy_label)) +
    geom_vline(xintercept = 0, color = "#777777", linetype = "dashed") +
    geom_segment(aes(x = lower, xend = upper, yend = strategy_label), color = "#2166AC", linewidth = .7) +
    geom_point(color = "#2166AC", size = 2.2) + facet_grid(population_label ~ window) +
    labs(title = "Increment from adding both laboratory domains: missingness sensitivity", y = NULL,
        x = "Paired AUROC difference (B+M+I minus B)",
        caption = paste("Points and 95% conditional percentile intervals from paired hospital/shared-patient block bootstrap (1,000 replicates).",
            "\nNew strategies are post-hoc sensitivities, not newly independent external validations; original analyses are retained.",
            "\nIntervals describe within-strategy increments; interval overlap does not test differences between strategies.")) + theme
ggsave(file.path(out, "figures/Missingness_domain_AUROC_sensitivity.png"), increments, width = 11, height = 6.2, dpi = 300)
ggsave(file.path(out, "figures/Missingness_domain_AUROC_sensitivity.svg"), increments, width = 11, height = 6.2, device = grDevices::svg)
cat("MISSINGNESS_SENSITIVITY_FIGURES_WRITTEN\n")
