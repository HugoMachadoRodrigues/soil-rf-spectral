# =============================================================================
# main_rf_ossl.R — Soil Spectral Modeling with Random Forest (OSSL / ASD)
# =============================================================================
# Pipeline:
#   0. Setup & configuration
#   1. Data download / load from OSSL (ASD VisNIR sensor, random fraction)
#   2. Spectral preprocessing (SG smooth, SG 1st deriv., SNV)
#   3. Hyperparameter tuning  (grid search, OOB-based, parallelised)
#   4. k-fold cross-validation (parallelised)
#   5. Final model fit on full data
#   6. Results export (metrics CSV, model .rds, figures)
#
# Author  : Hugo Machado Rodrigues
#           University of Florida — BMP Project
#           ORCID    : https://orcid.org/0000-0002-8070-8126
#           GitHub   : https://github.com/HugoMachadoRodrigues
#           ResearchGate: https://www.researchgate.net/profile/Hugo-Rodrigues-12
#           X/Twitter: https://x.com/Hugo_MRodrigues
#           LinkedIn : https://www.linkedin.com/in/hugo-rodrigues-52b535119/
#
# Inspired by: Clingensmith, C.M. & Grunwald, S. (2022). Predicting soil
#              properties and interpreting Vis-NIR models from across
#              continental United States. Sensors, 22(9), 3187.
#              https://doi.org/10.3390/s22093187
# Date    : 2026
# =============================================================================

# ── Auto working directory ─────────────────────────────────────────────────────
# Automatically sets the working directory to the folder containing this script.
# Works whether you source() the file or run it from the RStudio editor.
local({
  src <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
    error = function(e) ""
  )
  if (!nzchar(src) && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    src <- rstudioapi::getActiveDocumentContext()$path
  }
  if (nzchar(src)) {
    d <- dirname(src)
    if (!identical(normalizePath(d), normalizePath(getwd()))) {
      setwd(d)
      message("Working directory set to: ", d)
    }
  }
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════\n")
cat("  Soil RF Spectral Modeling — OSSL / ASD Pipeline\n")
cat("══════════════════════════════════════════════════\n\n")

# Install any missing packages from CRAN
pkg_required <- c(
  "ranger", "prospectr", "foreach", "doSNOW",
  "ggplot2", "ggpubr", "viridis", "moments", "httr",
  "dplyr", "tidyr", "readr"
)
pkg_missing <- pkg_required[!sapply(pkg_required, requireNamespace, quietly = TRUE)]
if (length(pkg_missing) > 0) {
  message("Installing missing packages: ", paste(pkg_missing, collapse = ", "))
  install.packages(pkg_missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(httr)
})

# Source utility modules (paths relative to this script's directory)
source("R/config.R")
source("R/utils_preprocessing.R")
source("R/utils_metrics.R")
source("R/utils_rf.R")
source("R/utils_visualization.R")

# Set parallel workers
if (is.null(N_CORES)) N_CORES <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Using %d parallel cores.\n\n", N_CORES))

# Output directories
dir.create(file.path(OUTPUT_DIR, "models"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# ── 1. OSSL Data —————————————————————————————————————————————————————————————
cat("── Step 1: Loading OSSL data ──────────────────────\n")

# OSSL public data — version 1.2
# Files can also be downloaded manually from https://soilspectroscopy.org
# CSV.gz format: no special serialisation package required; readr handles
# gzip decompression transparently. VisNIR is at L0; lab data at L1.
ossl_visnir_url  <- "https://storage.googleapis.com/soilspec4gg-public/ossl_visnir_L0_v1.2.csv.gz"
ossl_soillab_url <- "https://storage.googleapis.com/soilspec4gg-public/ossl_soillab_L1_v1.2.csv.gz"

ossl_visnir_local  <- file.path(DATA_DIR, "ossl_visnir_L0_v1.2.csv.gz")
ossl_soillab_local <- file.path(DATA_DIR, "ossl_soillab_L1_v1.2.csv.gz")

download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading: ", basename(dest))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    resp <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE),
                      httr::progress())
    httr::stop_for_status(resp)
  } else {
    message("Using cached: ", basename(dest))
  }
}

download_if_missing(ossl_visnir_url,  ossl_visnir_local)
download_if_missing(ossl_soillab_url, ossl_soillab_local)

cat("Reading spectral data...\n")
visnir_raw <- readr::read_csv(ossl_visnir_local,  show_col_types = FALSE)
soillab    <- readr::read_csv(ossl_soillab_local, show_col_types = FALSE)

# ── Filter to ASD sensor ──────────────────────────────────────────────────────
# ASD FieldSpec columns: scan_visnir.350_ref … scan_visnir.2500_ref (1-nm steps)
# Identify wavelength columns in the range [WAVE_MIN, WAVE_MAX]
wave_cols <- grep("^scan_visnir\\.(\\d+)_ref$", names(visnir_raw), value = TRUE)
wave_nm   <- as.numeric(gsub("scan_visnir\\.(\\d+)_ref", "\\1", wave_cols))

# Retain only wavelengths within the configured range
keep_cols <- wave_cols[wave_nm >= WAVE_MIN & wave_nm <= WAVE_MAX]
wave_nm   <- wave_nm[wave_nm  >= WAVE_MIN & wave_nm <= WAVE_MAX]

cat(sprintf("  VisNIR columns found: %d  |  After range filter [%d–%d nm]: %d\n",
            length(wave_cols), WAVE_MIN, WAVE_MAX, length(keep_cols)))

# Keep only rows that have all wavelength columns populated (ASD-like completeness)
spec_df <- visnir_raw %>%
  dplyr::select(id.layer_uuid_txt, all_of(keep_cols)) %>%
  tidyr::drop_na()

n_rows <- nrow(spec_df)
n_uuid <- dplyr::n_distinct(spec_df$id.layer_uuid_txt)
if (n_uuid < n_rows) {
  cat(sprintf("  Complete spectral observations: %d (%d unique layers; %.1f scans/layer avg)\n",
              n_rows, n_uuid, n_rows / n_uuid))
  cat("  Replicate scans retained — CV folds are assigned by layer UUID to prevent leakage.\n")
} else {
  cat(sprintf("  Complete spectral observations: %d (1 scan per layer)\n", n_rows))
}

# ── Join with lab data ────────────────────────────────────────────────────────
if (!TARGET_PROPERTY %in% names(soillab)) {
  stop(sprintf(
    "Target property '%s' not found in soillab. Available columns:\n%s",
    TARGET_PROPERTY,
    paste(grep("_usda|_iso", names(soillab), value = TRUE), collapse = "\n")
  ))
}

lab_sub <- soillab %>%
  dplyr::select(id.layer_uuid_txt, all_of(TARGET_PROPERTY)) %>%
  dplyr::filter(!is.na(.data[[TARGET_PROPERTY]]))

joined <- dplyr::inner_join(spec_df, lab_sub, by = "id.layer_uuid_txt")
cat(sprintf("  Observations after joining with lab data: %d\n", nrow(joined)))

# ── Random fraction + hard cap ────────────────────────────────────────────────
set.seed(SEED)
sampled <- dplyr::slice_sample(joined, prop = OSSL_FRACTION)
if (!is.null(N_MAX_SAMPLES) && is.finite(N_MAX_SAMPLES) &&
    nrow(sampled) > N_MAX_SAMPLES) {
  sampled <- dplyr::slice_sample(sampled, n = N_MAX_SAMPLES)
  cat(sprintf("  Sampled fraction (%.0f%%) capped at %d observations\n",
              OSSL_FRACTION * 100, nrow(sampled)))
} else {
  cat(sprintf("  Sampled fraction (%.0f%%): %d observations\n",
              OSSL_FRACTION * 100, nrow(sampled)))
}

# ── Matrices ──────────────────────────────────────────────────────────────────
X_raw   <- as.matrix(sampled[, keep_cols])
colnames(X_raw) <- paste0("w", wave_nm)

y_orig  <- sampled[[TARGET_PROPERTY]]
uuid_id <- sampled$id.layer_uuid_txt   # group key for leakage-free CV folds

# Summary statistics of target
cat("\nTarget property summary:\n")
print(summary_stats(y_orig))

# ── 2. Preprocessing ──────────────────────────────────────────────────────────
# Preprocessing of X is independent of y transformation — built once.
# Pipeline: raw reflectance → SG smooth → {sg_smooth, sg_deriv1, snv}
cat("\n── Step 2: Spectral preprocessing ────────────────\n")
cat(sprintf("  SG parameters: window = %d | poly order = %d\n", SG_WINDOW, SG_POLY))

spec_list <- build_preprocessing_stack(
  X_raw     = X_raw,
  sg_window = SG_WINDOW,
  sg_poly   = SG_POLY
)

for (nm in names(spec_list)) {
  cat(sprintf("  %-15s → %d × %d\n", nm, nrow(spec_list[[nm]]), ncol(spec_list[[nm]])))
}

# ── 3-6. Per-log-mode × per-preprocessing RF pipeline ────────────────────────
all_log_results <- list()

for (log_mode in LOG_MODES) {

  use_log  <- identical(log_mode, "log")
  mode_tag <- if (use_log) "log" else "raw"

  cat(sprintf(
    "\n╔══════════════════════════════════════════════════╗\n  Log mode : %s\n╚══════════════════════════════════════════════════╝\n\n",
    if (use_log) "log1p(y)  — back-transformed for all metrics"
    else         "original scale (no transformation)"
  ))

  if (use_log) {
    y_model       <- log1p(y_orig)
    back_trans_fn <- expm1
  } else {
    y_model       <- y_orig
    back_trans_fn <- NULL
  }

  # Folds are assigned by layer UUID so replicate scans of the same soil layer
  # always fall in the same fold — prevents leakage when n_uuid < n_rows.
  folds <- make_folds(y_model, k = CV_FOLDS, repeats = CV_REPEATS,
                      stratified = CV_STRATIFIED, seed = SEED,
                      groups = uuid_id)

  all_results    <- list()
  all_cv_metrics <- data.frame()

  for (prep_name in names(spec_list)) {

    cat(sprintf("\n══ [%s] Preprocessing: %s ══\n", toupper(mode_tag), toupper(prep_name)))
    X_prep <- spec_list[[prep_name]]

    # Remove any columns with zero variance (can occur with derivatives)
    col_var <- apply(X_prep, 2, var, na.rm = TRUE)
    X_prep  <- X_prep[, col_var > 1e-12, drop = FALSE]
    cat(sprintf("  Predictors after zero-var filter: %d\n", ncol(X_prep)))

    # ── Step 3: Tuning ────────────────────────────────────────────────────────
    cat("\n[3] Hyperparameter Tuning\n")
    tune_res <- tune_rf(
      X                = X_prep,
      y                = y_model,
      ntree_candidates = NTREE_CANDIDATES,
      mtry_fracs       = MTRY_FRACS,
      nodesize_vals    = NODESIZE_CANDIDATES,
      n_cores          = N_CORES,
      seed             = SEED
    )

    cat(sprintf("\nFull 3D grid: %d combinations evaluated\n", nrow(tune_res$grid)))
    cat("Top 10 combinations (ordered by OOB RMSE):\n")
    print(head(tune_res$grid[, c("ntree","mtry","min_nodesize","oob_rmse","oob_r2")], 10),
          row.names = FALSE, digits = 2)
    cat(sprintf("\nBest: ntree = %d | mtry = %d | min.node.size = %d | OOB RMSE = %.2f | OOB R² = %.2f\n\n",
                tune_res$best$ntree, tune_res$best$mtry, tune_res$best$min_nodesize,
                tune_res$best$oob_rmse, tune_res$best$oob_r2))

    p_tune <- plot_tuning_grid(
      tune_res$grid,
      title = sprintf("Tuning Grid — %s [%s | %s]", TARGET_LABEL, prep_name, mode_tag)
    )
    save_plot(p_tune,
              filename = sprintf("tuning_%s_%s_%s", TARGET_PROPERTY, mode_tag, prep_name),
              dir      = file.path(OUTPUT_DIR, "figures"),
              devices  = PLOT_DEVICES, dpi = PLOT_DPI)

    best_params <- list(
      ntree        = tune_res$best$ntree,
      mtry         = tune_res$best$mtry,
      min_nodesize = tune_res$best$min_nodesize
    )

    # ── Step 4: Cross-validation ──────────────────────────────────────────────
    cat("[4] Cross-Validation\n")
    cv_res <- cv_rf(
      X             = X_prep,
      y             = y_model,
      best_params   = best_params,
      folds         = folds,
      n_cores       = N_CORES,
      seed          = SEED,
      back_trans_fn = back_trans_fn
    )

    resamp_summ <- metrics_summary(cv_res$resample_metrics)
    cat(sprintf("\nCV Summary (%d resamples) [%s | %s]:\n",
                length(cv_res$resample_metrics), prep_name, mode_tag))
    print_metrics(resamp_summ)

    row_summ <- data.frame(
      preprocessing = prep_name,
      t(setNames(resamp_summ$mean, resamp_summ$metric)),
      stringsAsFactors = FALSE
    )
    for (m in resamp_summ$metric) {
      row_summ[[paste0(m, "_sd")]] <- resamp_summ$sd[resamp_summ$metric == m]
    }
    all_cv_metrics <- dplyr::bind_rows(all_cv_metrics, row_summ)

    p_op <- plot_obs_pred(
      cv_res$obs_pred,
      title = sprintf("Obs vs Pred — %s [%s | %s]", TARGET_LABEL, prep_name, mode_tag),
      label = TARGET_LABEL
    )
    save_plot(p_op,
              filename = sprintf("obs_pred_%s_%s_%s", TARGET_PROPERTY, mode_tag, prep_name),
              dir      = file.path(OUTPUT_DIR, "figures"),
              devices  = PLOT_DEVICES, dpi = PLOT_DPI)

    imp_wl <- as.numeric(gsub("w", "", names(cv_res$avg_importance)))
    p_imp  <- plot_importance(
      importance  = cv_res$avg_importance,
      wavelengths = imp_wl,
      title = sprintf("Variable Importance — %s [%s | %s]", TARGET_LABEL, prep_name, mode_tag)
    )
    save_plot(p_imp,
              filename = sprintf("importance_%s_%s_%s", TARGET_PROPERTY, mode_tag, prep_name),
              dir      = file.path(OUTPUT_DIR, "figures"),
              devices  = PLOT_DEVICES, dpi = PLOT_DPI)

    # ── Step 5: Final model on full data ──────────────────────────────────────
    cat("[5] Final Model (full dataset)\n")
    final_model <- fit_final_rf(
      X           = X_prep,
      y           = y_model,
      best_params = best_params,
      n_cores     = N_CORES,
      seed        = SEED
    )

    cat(sprintf("  OOB R²: %.2f | OOB RMSE: %.2f\n\n",
                final_model$r.squared, sqrt(final_model$prediction.error)))

    all_results[[prep_name]] <- list(
      tune        = tune_res,
      cv          = cv_res,
      resamp_summ = resamp_summ,
      final_model = final_model
    )

    if (SAVE_MODELS) {
      saveRDS(final_model,
              file.path(OUTPUT_DIR, "models",
                        sprintf("rf_final_%s_%s_%s.rds",
                                TARGET_PROPERTY, mode_tag, prep_name)))
    }
  }

  # ── Step 6: Cross-preprocessing comparison (within this log mode) ──────────
  cat(sprintf("\n── Step 6: Preprocessing comparison [%s] ─────\n\n", mode_tag))
  tmp <- all_cv_metrics[, c("preprocessing", "r2", "rmse", "mae", "me", "ccc", "rpd", "rpiq")]
  tmp[, -1] <- lapply(tmp[, -1], round, 2)
  print(tmp, row.names = FALSE)

  metrics_file <- file.path(OUTPUT_DIR, "metrics",
                             sprintf("cv_metrics_%s_%s.csv", TARGET_PROPERTY, mode_tag))
  readr::write_csv(all_cv_metrics, metrics_file)
  cat(sprintf("\nMetrics saved: %s\n", metrics_file))

  if ("rpd_sd" %in% names(all_cv_metrics)) {
    p_comp_rpd <- plot_cv_comparison(
      all_cv_metrics, metric = "rpd",
      title = sprintf("CV — RPD by Preprocessing [%s | %s]", TARGET_LABEL, mode_tag)
    )
    save_plot(p_comp_rpd,
              filename = sprintf("comparison_rpd_%s_%s", TARGET_PROPERTY, mode_tag),
              dir      = file.path(OUTPUT_DIR, "figures"),
              devices  = PLOT_DEVICES, dpi = PLOT_DPI)

    p_comp_r2 <- plot_cv_comparison(
      all_cv_metrics, metric = "r2",
      title = sprintf("CV — R² by Preprocessing [%s | %s]", TARGET_LABEL, mode_tag)
    )
    save_plot(p_comp_r2,
              filename = sprintf("comparison_r2_%s_%s", TARGET_PROPERTY, mode_tag),
              dir      = file.path(OUTPUT_DIR, "figures"),
              devices  = PLOT_DEVICES, dpi = PLOT_DPI)
  }

  saveRDS(all_results,
          file.path(OUTPUT_DIR, "models",
                    sprintf("all_results_%s_%s.rds", TARGET_PROPERTY, mode_tag)))

  all_log_results[[log_mode]] <- list(
    cv_metrics = all_cv_metrics,
    results    = all_results
  )
}

# ── Step 7: Log-mode comparison ───────────────────────────────────────────────
if (length(LOG_MODES) > 1) {
  cat("\n── Step 7: Log vs. No-Log comparison ─────────────\n\n")

  comp_rows <- lapply(names(all_log_results), function(lm) {
    cm      <- all_log_results[[lm]]$cv_metrics
    best_i  <- which.max(cm$r2)
    best    <- cm[best_i, c("preprocessing", "r2", "rmse", "mae", "rpd", "rpiq", "ccc")]
    best[, -1] <- lapply(best[, -1], round, 2)
    cbind(log_mode = lm, best)
  })
  comp_df <- do.call(rbind, comp_rows)
  print(comp_df, row.names = FALSE)

  readr::write_csv(comp_df,
    file.path(OUTPUT_DIR, "metrics",
              sprintf("log_comparison_%s.csv", TARGET_PROPERTY)))
}

cat("\n══════════════════════════════════════════════════\n")
cat("  Pipeline complete. Outputs in:", OUTPUT_DIR, "\n")
cat("══════════════════════════════════════════════════\n\n")
