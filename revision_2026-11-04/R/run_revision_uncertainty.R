# Conditional validation uncertainty for frozen models; no retraining/calibration updating.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else file.path(getwd(), "revision_2026-11-04")
source(file.path(root, "R", "run_revision_models.R"), encoding = "UTF-8")

metric_names <- c("auroc", "auprc", "brier_score", "log_loss", "calibration_intercept",
    "calibration_slope", "mean_predicted", "observed_rate", "risk_bias",
    "observed_expected_ratio", "descriptive_null_brier", "scaled_brier")

prepare_probability <- function(y, p, constant = FALSE) {
    order <- order(-p, seq_along(p))
    sorted <- p[order]
    groups <- cumsum(c(TRUE, diff(sorted) != 0))
    cutoffs <- unique(as.numeric(quantile(p, probs = seq(0.1, 0.9, 0.1), names = FALSE)))
    cutoffs <- cutoffs[cutoffs > min(p) & cutoffs < max(p)]
    bins <- cut(p, breaks = c(-Inf, cutoffs, Inf), labels = FALSE, right = TRUE)
    list(y = y, p = p, lp = qlogis(pmin(pmax(p, 1e-6), 1-1e-6)),
        order = order, score_group = groups, bins = bins, cutoffs = cutoffs, constant = constant)
}

weighted_metrics <- function(object, weights) {
    y <- object$y; p <- object$p; n <- sum(weights)
    if (n == 0) return(setNames(rep(NA_real_, length(metric_names)), metric_names))
    events <- sum(weights*y); nonevents <- n-events
    observed <- events/n; predicted <- sum(weights*p)/n
    brier <- sum(weights*(y-p)^2)/n
    pp <- pmin(pmax(p, 1e-6), 1-1e-6)
    logloss <- -sum(weights*(y*log(pp)+(1-y)*log(1-pp)))/n
    auc <- pr <- citl <- slope <- NA_real_
    if (events > 0 && nonevents > 0) {
        idx <- object$order
        by_score <- rowsum(cbind(weights[idx]*y[idx], weights[idx]*(1-y[idx])),
                          object$score_group, reorder = FALSE)
        by_score <- by_score[rowSums(by_score) > 0, , drop = FALSE]
        tp <- cumsum(by_score[, 1L]); fp <- cumsum(by_score[, 2L])
        recall <- tp/events; fpr <- fp/nonevents; precision <- tp/(tp+fp)
        auc <- sum(diff(c(0, fpr))*(c(0, head(recall, -1L))+recall)/2)
        if (!object$constant)
            pr <- sum(diff(c(0, recall))*(c(1, head(precision, -1L))+precision)/2)
        citl <- tryCatch(uniroot(function(a) sum(weights*(y-plogis(object$lp+a))),
                                interval = c(-50, 50), tol = 1e-10)$root,
                         error = function(e) NA_real_)
        if (!object$constant) {
            keep <- weights > 0
            fit <- tryCatch(suppressWarnings(glm.fit(cbind(1, object$lp[keep]), y[keep],
                weights = weights[keep], family = binomial(), control = glm.control(maxit = 100))),
                error = function(e) NULL)
            if (!is.null(fit) && isTRUE(fit$converged) && length(fit$coefficients) == 2L)
                slope <- unname(fit$coefficients[2L])
        }
    }
    null_brier <- observed*(1-observed)
    setNames(c(auc, pr, brier, logloss, citl, slope, predicted, observed, predicted-observed,
        if (predicted > 0) observed/predicted else NA_real_, null_brier,
        if (null_brier > 0) 1-brier/null_brier else NA_real_), metric_names)
}

hospital_components <- function(group_id, hospital_id) {
    hospitals <- sort(unique(as.character(hospital_id)))
    parent <- seq_along(hospitals)
    find <- function(x) {
        while (parent[x] != x) x <- parent[x]
        x
    }
    patient_hospitals <- split(as.character(hospital_id), group_id)
    for (h in patient_hospitals) {
        idx <- match(unique(h), hospitals)
        if (length(idx) > 1L) for (j in idx[-1L]) {
            a <- find(idx[1L]); b <- find(j)
            if (a != b) parent[b] <- a
        }
    }
    roots <- vapply(seq_along(hospitals), find, integer(1))
    component <- match(roots, unique(roots))
    mapped <- component[match(as.character(hospital_id), hospitals)]
    stopifnot(all(vapply(split(mapped, group_id), function(x) length(unique(x)), integer(1)) == 1L),
              all(vapply(split(mapped, hospital_id), function(x) length(unique(x)), integer(1)) == 1L))
    list(row_component = mapped, map = data.frame(hospitalid = hospitals, component = component),
         cross_hospital_patients = sum(vapply(patient_hospitals, function(x) length(unique(x)), integer(1)) > 1L))
}

align_predictions <- function(path, dataset) {
    raw <- read.csv(path)
    raw <- raw[raw$dataset == dataset, , drop = FALSE]
    models <- c("xgboost", "random_forest", "elastic_net", "clinical_reference", "constant_training_prevalence")
    reference <- raw[raw$model == "xgboost", ]
    stopifnot(!anyDuplicated(reference$stay_id))
    reference <- reference[order(reference$stay_id), ]
    objects <- list()
    for (name in models) {
        data <- raw[raw$model == name, ]
        data <- data[match(reference$stay_id, data$stay_id), ]
        stopifnot(nrow(data) == nrow(reference), !anyNA(data$probability),
            identical(data$truth, reference$truth), identical(data$group_id, reference$group_id),
            identical(as.character(data$hospitalid), as.character(reference$hospitalid)))
        objects[[name]] <- prepare_probability(as.integer(data$truth == "yes"), data$probability,
            constant = name == "constant_training_prevalence")
    }
    list(reference = reference, objects = objects)
}

validate_weighted_engine <- function(aligned) {
    n <- nrow(aligned$reference)
    for (seed in c(20261104L, 20261105L)) {
        set.seed(seed)
        counts <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
        duplicated_rows <- rep(seq_len(n), counts)
        for (name in names(aligned$objects)) {
            object <- aligned$objects[[name]]
            fast <- weighted_metrics(object, counts)
            library <- point_scores(object$y[duplicated_rows], object$p[duplicated_rows])
            metrics <- c("auroc", "brier_score", "log_loss", "calibration_intercept")
            if (!object$constant) metrics <- c(metrics, "auprc", "calibration_slope")
            if (max(abs(fast[metrics]-library[metrics])) > 2e-6)
                stop("Weighted metric engine differs from explicit duplicated-record calculation")
        }
    }
    TRUE
}

quantile_ci <- function(x, structural = FALSE) {
    x <- x[is.finite(x)]
    if (structural) return(c(lower = NA_real_, upper = NA_real_, valid_reps = 0))
    if (!length(x)) return(c(lower = NA_real_, upper = NA_real_, valid_reps = 0))
    q <- quantile(x, probs = c(0.025, 0.975), names = FALSE, type = 7)
    c(lower = q[1L], upper = q[2L], valid_reps = length(x))
}

bootstrap_scheme <- function(aligned, hour, dataset, scheme, cluster, run_hash, root, seed, reps = 1000L) {
    models <- names(aligned$objects); n <- nrow(aligned$reference)
    group_index <- match(cluster, unique(cluster)); groups <- length(unique(cluster))
    stopifnot(!anyNA(group_index), groups >= 2L)
    label <- paste(hour, dataset, scheme, sep = "_")
    cache <- file.path(root, "private_cache", paste0("uncertainty_", label, ".rds"))
    checkpoint <- file.path(root, "private_cache", paste0("uncertainty_", label, "_checkpoint.rds"))
    if (file.exists(cache)) {
        result <- readRDS(cache)
        if (!identical(result$run_hash, run_hash)) stop("Completed bootstrap hash mismatch")
        log_progress(label, " completed cache verified")
        return(result)
    }
    if (file.exists(checkpoint)) {
        result <- readRDS(checkpoint)
        if (!identical(result$run_hash, run_hash)) stop("Bootstrap checkpoint hash mismatch")
        assign(".Random.seed", result$rng, envir = .GlobalEnv)
        start <- result$completed + 1L
    } else {
        set.seed(seed)
        result <- list(run_hash = run_hash, completed = 0L,
            metrics = array(NA_real_, c(reps, length(models), length(metric_names)),
                            dimnames = list(NULL, models, metric_names)),
            bin_rates = lapply(aligned$objects, function(x)
                matrix(NA_real_, reps, max(x$bins))),
            records = numeric(reps), events = numeric(reps), clusters = groups,
            hour = hour, dataset = dataset, scheme = scheme)
        start <- 1L
    }
    points <- lapply(aligned$objects, weighted_metrics, weights = rep(1, n))
    result$point_metrics <- points
    if (start <= reps) for (b in start:reps) {
        multiplicity <- tabulate(sample.int(groups, groups, replace = TRUE), nbins = groups)
        weights <- multiplicity[group_index]
        result$records[b] <- sum(weights)
        result$events[b] <- sum(weights*aligned$objects[[1L]]$y)
        for (name in models) {
            object <- aligned$objects[[name]]
            result$metrics[b, name, ] <- weighted_metrics(object, weights)
            by_bin <- rowsum(cbind(weights*object$y, weights), object$bins, reorder = TRUE)
            rates <- by_bin[, 1L]/by_bin[, 2L]
            result$bin_rates[[name]][b, ] <- rates
        }
        if (b %% 100L == 0L || b == reps) {
            result$completed <- b
            result$rng <- .Random.seed
            saveRDS(result, checkpoint)
            log_progress(label, " bootstrap ", b, "/", reps)
        }
    }
    saveRDS(result, cache)
    result
}

summarize_scheme <- function(result, aligned) {
    models <- names(aligned$objects); metrics <- list(); pairs <- list(); bins <- list()
    for (name in models) {
        object <- aligned$objects[[name]]
        for (metric in metric_names) {
            structural <- object$constant && metric %in% c("auprc", "calibration_slope")
            ci <- quantile_ci(result$metrics[, name, metric], structural)
            metrics[[length(metrics)+1L]] <- data.frame(window_hours = result$hour,
                dataset = result$dataset, scheme = result$scheme, clusters = result$clusters,
                model = name, metric = metric, estimate = result$point_metrics[[name]][metric],
                lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]],
                status = if (structural) "not_applicable_constant_score" else
                    if (ci[["valid_reps"]] < 950L) "unstable_insufficient_valid_replicates" else "computed")
        }
        for (bin in seq_len(max(object$bins))) {
            keep <- object$bins == bin
            ci <- quantile_ci(result$bin_rates[[name]][, bin])
            bins[[length(bins)+1L]] <- data.frame(window_hours = result$hour,
                dataset = result$dataset, scheme = result$scheme, model = name, bin = bin,
                records = sum(keep), patients = length(unique(aligned$reference$group_id[keep])),
                hospitals = if (all(is.na(aligned$reference$hospitalid))) NA_integer_
                    else length(unique(aligned$reference$hospitalid[keep])),
                mean_predicted = mean(object$p[keep]), observed_rate = mean(object$y[keep]),
                events = sum(object$y[keep]), minimum_predicted = min(object$p[keep]),
                maximum_predicted = max(object$p[keep]), lower = ci[["lower"]], upper = ci[["upper"]],
                valid_reps = ci[["valid_reps"]], zero_event_bin = sum(object$y[keep]) == 0L)
        }
    }
    for (comparator in setdiff(models, "xgboost")) for (metric in c("auroc", "auprc", "brier_score", "log_loss")) {
        if (comparator == "constant_training_prevalence" && metric == "auprc") next
        ci <- quantile_ci(result$metrics[, "xgboost", metric]-result$metrics[, comparator, metric])
        pairs[[length(pairs)+1L]] <- data.frame(window_hours = result$hour,
            dataset = result$dataset, scheme = result$scheme, comparator = comparator, metric = metric,
            difference_xgboost_minus_comparator = result$point_metrics$xgboost[metric]-result$point_metrics[[comparator]][metric],
            lower = ci[["lower"]], upper = ci[["upper"]], valid_reps = ci[["valid_reps"]],
            preferred_difference_direction = if (metric %in% c("auroc", "auprc")) "positive" else "negative")
    }
    list(metrics = do.call(rbind, metrics), pairs = do.call(rbind, pairs), bins = do.call(rbind, bins))
}

run_uncertainty <- function(root) {
    out <- file.path(root, "validation_results")
    dir.create(out, showWarnings = FALSE)
    script <- file.path(root, "R", "run_revision_uncertainty.R")
    plan <- list(version = "conditional_validation_1.0", replicates = 1000L,
        seed_base = 20261104L, interval_method = "2.5 and 97.5 percentile; type 7 quantiles",
        target = "Pooled validation performance conditional on frozen fitted models and preprocessing; not development uncertainty or a prediction interval for a new hospital",
        resampling = "MIMIC: patient. eICU: patient, hospital, and connected hospital blocks linked by shared patients",
        primary_reporting = "MIMIC patient bootstrap; eICU hospital/shared-patient connected-component bootstrap; other two eICU schemes sensitivity",
        pairing = "Same cluster multiplicities for all models and all paired differences within each replicate",
        calibration_bins = "Fixed prediction-only empirical decile cutpoints; tied predictions never split; intervals for observed rates conditional on fixed bins",
        multiplicity = "Exploratory multiple comparisons; no multiplicity-adjusted confirmatory inference",
        script_sha256 = sha_file(script),
        predictions_sha256 = setNames(lapply(c(24L, 6L), function(h)
            sha_file(file.path(root, "private_predictions", paste0("predictions_", h, "h.csv")))), c("24", "6")))
    run_hash <- digest::digest(plan, algo = "sha256")
    plan_path <- file.path(out, "uncertainty_plan.json")
    if (file.exists(plan_path)) {
        previous <- jsonlite::read_json(plan_path)
        if (!identical(previous$run_hash, run_hash)) stop("Existing uncertainty run has another configuration")
    } else write_json(list(run_hash = run_hash, plan = plan), plan_path)
    metric_tables <- pair_tables <- bin_tables <- list(); structure <- list()
    engine_validated <- FALSE
    counter <- 0L
    for (hour in c(24L, 6L)) for (dataset in c("MIMIC_internal_test", "eICU_external_validation")) {
        aligned <- align_predictions(file.path(root, "private_predictions", paste0("predictions_", hour, "h.csv")), dataset)
        if (!engine_validated) {
            engine_validated <- validate_weighted_engine(aligned)
            log_progress("Weighted ROC/PR/calibration checked against explicit record replication")
        }
        reference <- aligned$reference
        schemes <- list(patient = reference$group_id)
        if (dataset == "eICU_external_validation") {
            stopifnot(!anyNA(reference$hospitalid))
            components <- hospital_components(reference$group_id, reference$hospitalid)
            schemes$hospital <- reference$hospitalid
            schemes$hospital_patient_component <- components$row_component
            write.csv(components$map, file.path(root, "private_splits", paste0("hospital_components_", hour, "h.csv")), row.names = FALSE)
            structure[[paste0(hour, "h")]] <- list(hospitals = length(unique(reference$hospitalid)),
                patients = length(unique(reference$group_id)),
                patients_crossing_hospitals = components$cross_hospital_patients,
                independent_resampling_blocks = length(unique(components$row_component)),
                largest_block_hospitals = max(table(components$map$component)))
            log_progress(hour, "h external: ", length(unique(reference$hospitalid)), " hospitals; ",
                components$cross_hospital_patients, " cross-hospital patients; ",
                length(unique(components$row_component)), " connected blocks")
        }
        for (scheme in names(schemes)) {
            counter <- counter+1L
            result <- bootstrap_scheme(aligned, hour, dataset, scheme, schemes[[scheme]], run_hash, root,
                seed = plan$seed_base+counter*1000L, reps = plan$replicates)
            summary <- summarize_scheme(result, aligned)
            metric_tables[[length(metric_tables)+1L]] <- summary$metrics
            pair_tables[[length(pair_tables)+1L]] <- summary$pairs
            bin_tables[[length(bin_tables)+1L]] <- summary$bins
        }
    }
    write.csv(do.call(rbind, metric_tables), file.path(out, "metric_confidence_intervals.csv"), row.names = FALSE)
    write.csv(do.call(rbind, pair_tables), file.path(out, "paired_model_differences.csv"), row.names = FALSE)
    write.csv(do.call(rbind, bin_tables), file.path(out, "calibration_bins_intervals.csv"), row.names = FALSE)
    write_json(structure, file.path(out, "cluster_structure.json"))
    write_json(list(run_hash = run_hash, weighted_engine_validation = engine_validated,
        status = "conditional_validation_uncertainty_computed_not_final_paper_ready",
        model_retraining = FALSE, calibration_updating = FALSE,
        pending = c("clinical score reference clarification", "DCA, SHAP and remaining sensitivity analyses", "manuscript and response revision")),
        file.path(out, "uncertainty_complete.json"))
    log_progress("Conditional uncertainty, paired comparisons and calibration-bin intervals complete")
}

if (sys.nframe() == 0L) run_uncertainty(root)
