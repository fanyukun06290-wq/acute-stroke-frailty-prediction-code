# Verify saved paired bootstrap quantiles without retraining or resampling.
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
.libPaths(unique(c(Sys.getenv("STROKE_R_LIBRARY", ""),.libPaths())))
root <- file.path(getwd(), "revision_2026-11-04")
out <- file.path(root, "validation_results")
plan <- jsonlite::read_json(file.path(out, "uncertainty_plan.json"))
pairs <- read.csv(file.path(out, "paired_model_differences.csv"))
metrics <- read.csv(file.path(out, "metric_confidence_intervals.csv"))
bins <- read.csv(file.path(out, "calibration_bins_intervals.csv"))
files <- list.files(file.path(root, "private_cache"), pattern = "^uncertainty_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_checkpoint", files)]
stopifnot(length(files) == 8L)
checked_pairs <- checked_metrics <- checked_bins <- 0L
same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10))
ci <- function(v) as.numeric(quantile(v[is.finite(v)], c(.025, .975), type = 7))
for (f in files) {
    x <- readRDS(f)
    stopifnot(x$completed == 1000L, dim(x$metrics)[1L] == 1000L, identical(x$run_hash, plan$run_hash))
    p <- pairs[pairs$window_hours == x$hour & pairs$dataset == x$dataset & pairs$scheme == x$scheme, ]
    for (i in seq_len(nrow(p))) {
        r <- p[i, ]; v <- x$metrics[, "xgboost", r$metric] - x$metrics[, r$comparator, r$metric]
        stopifnot(same(ci(v), c(r$lower, r$upper)), sum(is.finite(v)) == r$valid_reps)
        checked_pairs <- checked_pairs + 1L
    }
    m <- metrics[metrics$window_hours == x$hour & metrics$dataset == x$dataset & metrics$scheme == x$scheme, ]
    for (i in seq_len(nrow(m))) {
        r <- m[i, ]; v <- x$metrics[, r$model, r$metric]
        if (r$status == "computed") {
            stopifnot(same(ci(v), c(r$lower, r$upper)), sum(is.finite(v)) == r$valid_reps)
            checked_metrics <- checked_metrics + 1L
        }
    }
    b <- bins[bins$window_hours == x$hour & bins$dataset == x$dataset & bins$scheme == x$scheme, ]
    for (i in seq_len(nrow(b))) {
        r <- b[i, ]; v <- x$bin_rates[[r$model]][, r$bin]
        stopifnot(same(ci(v), c(r$lower, r$upper)), sum(is.finite(v)) == r$valid_reps)
        checked_bins <- checked_bins + 1L
    }
}
result <- list(status = "passed", schemes = length(files), paired_intervals_checked = checked_pairs,
    metric_intervals_checked = checked_metrics, bin_intervals_checked = checked_bins,
    model_retraining = FALSE, new_resampling = FALSE)
jsonlite::write_json(result, file.path(root, "qc", "uncertainty_cache_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(result)
