Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "gcs_rule_sensitivity_results")
plan <- jsonlite::read_json(file.path(out, "validation_plan.json"))
sha <- function(p) digest::digest(file = p, algo = "sha256")
stopifnot(identical(sha(file.path(root, "R/validate_gcs_rule_stress.R")), plan$plan$script_sha256))
for (h in c("24", "6")) stopifnot(identical(sha(file.path(root, "private_predictions", paste0("gcs_rule_stress_", h, "h.csv"))), plan$plan$prediction_hashes[[h]]))
p <- read.csv(file.path(out, "paired_input_rule_differences.csv"))
d <- read.csv(file.path(out, "domain_increment_stress.csv"))
pop <- read.csv(file.path(out, "evaluation_populations.csv"))
stopifnot(nrow(p) == 384L, nrow(d) == 112L, nrow(pop) == 4L)
ci <- function(v) as.numeric(quantile(v[is.finite(v)], c(.025, .975), type = 7))
count_p <- count_d <- count_change <- replay <- 0L
for (i in seq_len(nrow(pop))) {
    r <- pop[i, ]; label <- paste(r$window_hours, r$dataset, sep = "_")
    boot <- readRDS(file.path(root, "private_cache", paste0("gcs_rules_validation_", label, ".rds")))
    stopifnot(identical(boot$run_hash, plan$run_hash), boot$completed == 1000L)
    replay <- replay+boot$baseline_vectors_checked
    a <- p[p$window_hours == r$window_hours & p$dataset == r$dataset, ]
    b <- d[d$window_hours == r$window_hours & d$dataset == r$dataset, ]
    for (j in seq_len(nrow(a))) {
        v <- boot$metrics[, paste(a$variant[j], a$model[j], sep = "/"), a$metric[j]]-
            boot$metrics[, paste("strict_baseline", a$model[j], sep = "/"), a$metric[j]]
        stopifnot(max(abs(ci(v)-c(a$lower[j], a$upper[j]))) < 1e-10, sum(is.finite(v)) == a$valid_reps[j])
        if (a$unchanged_predictions[j]) stopifnot(max(abs(v)) == 0, a$difference[j] == 0)
        count_p <- count_p+1L
    }
    strict <- boot$metrics[, "strict_baseline/omit_all_flags_BMI", ]-boot$metrics[, "strict_baseline/omit_all_flags_B", ]
    for (j in seq_len(nrow(b))) {
        v <- boot$metrics[, paste(b$variant[j], "omit_all_flags_BMI", sep = "/"), b$metric[j]]-
            boot$metrics[, paste(b$variant[j], "omit_all_flags_B", sep = "/"), b$metric[j]]
        change <- v-strict[, b$metric[j]]
        stopifnot(max(abs(ci(v)-c(b$lower[j], b$upper[j]))) < 1e-10,
            max(abs(ci(change)-c(b$change_lower[j], b$change_upper[j]))) < 1e-10,
            sum(is.finite(v)) == b$valid_reps[j])
        count_d <- count_d+1L; count_change <- count_change+1L
    }
}
stopifnot(count_p == 384L, count_d == 112L, count_change == 112L, replay == 160L)
jsonlite::write_json(list(status = "passed_all_fixed_model_stress_interval_checks",
    input_rule_difference_intervals = count_p, domain_increment_intervals = count_d,
    change_in_increment_intervals = count_change, original_bootstrap_vectors_reproduced = replay,
    configurations = 4L, replicates_each = 1000L, models_refitted = FALSE, posthoc = TRUE),
    file.path(root, "qc/gcs_rules_interval_checks.json"), auto_unbox = TRUE, pretty = TRUE)
cat("GCS_RULE_STRESS_INTERVAL_CHECKS_PASSED\n")
