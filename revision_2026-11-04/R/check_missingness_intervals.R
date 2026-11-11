Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "missingness_sensitivity_results")
plan <- jsonlite::read_json(file.path(out, "validation_plan.json"))
sha <- function(p) digest::digest(file = p, algo = "sha256")
stopifnot(identical(sha(file.path(root, "R/validate_missingness_sensitivity.R")), plan$plan$script_sha256))
for (hour in c("24", "6")) {
    stopifnot(identical(sha(file.path(root, "private_predictions", paste0("missing_sensitivity_predictions_", hour, "h.csv"))), plan$plan$prediction_hashes[[hour]]),
        identical(sha(file.path(root, "private_predictions", paste0("neuro_predictions_", hour, "h.csv"))), plan$plan$prior_prediction_hashes[[hour]]))
}
m <- read.csv(file.path(out, "metric_intervals.csv"))
p <- read.csv(file.path(out, "paired_domain_differences.csv"))
o <- read.csv(file.path(out, "paired_original_differences.csv"))
pop <- read.csv(file.path(out, "evaluation_populations.csv"))
stopifnot(nrow(m) == 576L, nrow(p) == 240L, nrow(o) == 96L, nrow(pop) == 6L)
ci <- function(v) as.numeric(quantile(v[is.finite(v)], c(.025, .975), type = 7))
count_m <- count_p <- count_o <- replay <- 0L
for (i in seq_len(nrow(pop))) {
    r <- pop[i, ]; label <- paste(r$window_hours, r$dataset, r$population, sep = "_")
    boot <- readRDS(file.path(root, "private_cache", paste0("missing_validation_", label, ".rds")))
    old_label <- paste(label, if (r$population == "all_eligible") "motor_total" else "motor", sep = "_")
    old <- readRDS(file.path(root, "private_cache", paste0("neuro_validation_", old_label, ".rds")))
    stopifnot(identical(boot$run_hash, plan$run_hash), boot$completed == 1000L,
        identical(old$run_hash, plan$plan$prior_validation_hash))
    replay <- replay+boot$replay_checks
    subset_rows <- function(x) x[x$window_hours == r$window_hours & x$dataset == r$dataset & x$population == r$population, ]
    a <- subset_rows(m); b <- subset_rows(p); c <- subset_rows(o)
    for (j in seq_len(nrow(a))) {
        v <- boot$metrics[, a$model[j], a$metric[j]]
        stopifnot(max(abs(ci(v)-c(a$lower[j], a$upper[j]))) < 1e-10, sum(is.finite(v)) == a$valid_reps[j]); count_m <- count_m+1L
    }
    for (j in seq_len(nrow(b))) {
        v <- boot$metrics[, b$augmented[j], b$metric[j]]-boot$metrics[, b$reference[j], b$metric[j]]
        stopifnot(max(abs(ci(v)-c(b$lower[j], b$upper[j]))) < 1e-10, sum(is.finite(v)) == b$valid_reps[j]); count_p <- count_p+1L
    }
    for (j in seq_len(nrow(c))) {
        v <- boot$metrics[, c$model[j], c$metric[j]]-old$metrics[, c$original_model[j], c$metric[j]]
        stopifnot(max(abs(ci(v)-c(c$lower[j], c$upper[j]))) < 1e-10, sum(is.finite(v)) == c$valid_reps[j]); count_o <- count_o+1L
    }
}
stopifnot(count_m == 576L, count_p == 240L, count_o == 96L, replay == 60L)
jsonlite::write_json(list(status = "passed_all_cached_interval_checks", metric_intervals_verified = count_m,
    domain_intervals_verified = count_p, original_comparison_intervals_verified = count_o,
    weight_replay_checks = replay, configurations = 6L, replicates_each = 1000L, posthoc = TRUE),
    file.path(root, "qc/missing_sensitivity_interval_checks.json"), auto_unbox = TRUE, pretty = TRUE)
cat("MISSINGNESS_SENSITIVITY_INTERVAL_CHECKS_PASSED\n")
