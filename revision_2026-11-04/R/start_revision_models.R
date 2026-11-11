# Windows launcher: ensure UTF-8 support before reading/writing Chinese paths.
if (.Platform$OS.type == "windows") {
    active <- suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))
    if (!nzchar(active) || !isTRUE(l10n_info()[["UTF-8"]]))
        stop("An active UTF-8 Windows locale is required for the project path")
}
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else {
    if (file.exists(file.path(getwd(), "model_feature_contract.json"))) getwd()
    else file.path(getwd(), "revision_2026-11-04")
}
source(file.path(root, "R", "run_revision_models.R"), encoding = "UTF-8")
run_revision(root)
