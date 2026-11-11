# Recompute every reported interval from frozen validation arrays.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "neuro_reference_results")
plan <- jsonlite::read_json(file.path(out, "validation_plan.json"))
sha <- function(p) digest::digest(file = p, algo = "sha256")
stopifnot(identical(sha(file.path(root, "R/validate_neuro_nested.R")), plan$plan$script_sha256))
for (h in c("24", "6")) stopifnot(identical(sha(file.path(root, "private_predictions", paste0("neuro_predictions_", h, "h.csv"))), plan$plan$prediction_hashes[[h]]))
m <- read.csv(file.path(out, "metric_intervals.csv"))
p <- read.csv(file.path(out, "paired_differences.csv"))
pop <- read.csv(file.path(out, "evaluation_populations.csv"))
stopifnot(nrow(m) == 576L, nrow(p) == 240L, nrow(pop) == 8L,
    all(m$valid_reps >= 950L), all(p$valid_reps >= 950L))
ci <- function(v) as.numeric(quantile(v[is.finite(v)], c(.025, .975), type = 7))
count_m <- count_p <- replay <- 0L
for (i in seq_len(nrow(pop))) {
    r <- pop[i, ]
    label <- paste(r$window_hours, r$dataset, r$population, gsub(",", "_", r$neuro_types), sep = "_")
    boot <- readRDS(file.path(root, "private_cache", paste0("neuro_validation_", label, ".rds")))
    stopifnot(identical(boot$run_hash, plan$run_hash), boot$completed == 1000L)
    replay <- replay+boot$weight_replay_checks
    models <- dimnames(boot$metrics)[[2L]]
    a <- m[m$window_hours == r$window_hours & m$dataset == r$dataset & m$population == r$population & m$model %in% models, ]
    b <- p[p$window_hours == r$window_hours & p$dataset == r$dataset & p$population == r$population & p$augmented %in% models, ]
    for (j in seq_len(nrow(a))) {
        v <- boot$metrics[, a$model[j], a$metric[j]]
        stopifnot(max(abs(ci(v)-c(a$lower[j], a$upper[j]))) < 1e-10, sum(is.finite(v)) == a$valid_reps[j])
        count_m <- count_m+1L
    }
    for (j in seq_len(nrow(b))) {
        v <- boot$metrics[, b$augmented[j], b$metric[j]]-boot$metrics[, b$reference[j], b$metric[j]]
        stopifnot(max(abs(ci(v)-c(b$lower[j], b$upper[j]))) < 1e-10, sum(is.finite(v)) == b$valid_reps[j])
        count_p <- count_p+1L
    }
}
stopifnot(count_m == 576L, count_p == 240L, replay == 40L)
jsonlite::write_json(list(status = "passed_all_cached_interval_checks", metric_intervals_verified = count_m,
    paired_intervals_verified = count_p, original_xgboost_weight_replay_checks = replay,
    resampling_configurations = 8L, replicates_each = 1000L, retraining = FALSE),
    file.path(root, "qc/neuro_interval_checks.json"), auto_unbox = TRUE, pretty = TRUE)
cat("NEURO_INTERVAL_CHECKS_PASSED\n")
