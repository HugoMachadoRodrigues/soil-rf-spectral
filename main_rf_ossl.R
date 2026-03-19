# =============================================================================
# main_rf_ossl.R — Soil Spectral Modeling with Random Forest (OSSL / ASD)
# =============================================================================
# Pipeline:
#   0. Setup & configuration
#   1. Data download / load from OSSL (ASD VisNIR sensor, random fraction)
#   2. Spectral preprocessing (SG smooth, SG 1st deriv., SNV, absorbance)
#   3. Hyperparameter tuning  (grid search, OOB-based, parallelised)
#   4. k-fold cross-validation (parallelised)
#   5. Final model fit on full data
#   6. Results export (metrics CSV, model .rds, figures)
#
# Authors : [Your Name]
# Inspired by: Clingensmith, C. & Grunwald, S. — NRCS Soil Spectral Modeling
# Date    : 2024
# =============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════\n")
cat("  Soil RF Spectral Modeling — OSSL / ASD Pipeline\n")
cat("══════════════════════════════════════════════════\n\n")

# Install any missing packages from CRAN
pkg_required <- c(
  "ranger", "prospectr", "foreach", "doParallel",
  "ggplot2", "ggpubr", "viridis", "moments", "httr",
  "dplyr", "tidyr", "readr", "qs2"
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
  library(qs2)
  library(httr)
})

# Source utility modules — paths are relative to the working directory.
# Run setwd("/path/to/soil-rf-spectral") before sourcing this script.
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
ossl_visnir_url <- "https://storage.googleapis.com/soilspec4gg-public/ossl_visnir_L1_v1.2.qs"
ossl_soillab_url <- "https://storage.googleapis.com/soilspec4gg-public/ossl_soillab_L1_v1.2.qs"

ossl_visnir_local  <- file.path(DATA_DIR, "ossl_visnir_L1_v1.2.qs")
ossl_soillab_local <- file.path(DATA_DIR, "ossl_soillab_L1_v1.2.qs")

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
# qread_qs1() reads files written by the original qs package (OSSL v1.2 format)
visnir_raw <- qs2::qread_qs1(ossl_visnir_local)
soillab    <- qs2::qread_qs1(ossl_soillab_local)

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

cat(sprintf("  Complete spectral observations: %d\n", nrow(spec_df)))

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

# ── Random fraction ───────────────────────────────────────────────────────────
set.seed(SEED)
sampled <- dplyr::slice_sample(joined, prop = OSSL_FRACTION)
cat(sprintf("  Sampled fraction (%.0f%%): %d observations\n",
            OSSL_FRACTION * 100, nrow(sampled)))

# ── Matrices ──────────────────────────────────────────────────────────────────
X_raw <- as.matrix(sampled[, keep_cols])
colnames(X_raw) <- paste0("w", wave_nm)

y_orig <- sampled[[TARGET_PROPERTY]]

# Summary statistics of target
cat("\nTarget property summary:\n")
print(summary_stats(y_orig))

# ── Target transformation ─────────────────────────────────────────────────────
if (USE_LOG_TARGET) {
  # Shift to ensure positivity before log: log1p handles zeros
  y_shift <- min(y_orig, na.rm = TRUE)
  if (y_shift <= 0) y_shift <- 0  # log1p(0) = 0 is fine
  y_model       <- log1p(y_orig)
  back_trans_fn <- function(x) expm1(x)
  cat(sprintf("\nTarget log-transformed: log1p(%s)\n", TARGET_PROPERTY))
} else {
  y_model       <- y_orig
  back_trans_fn <- NULL
  cat(sprintf("\nNo log transformation applied to target.\n"))
}

# ── 2. Preprocessing ──────────────────────────────────────────────────────────
# Pipeline: raw reflectance → SG smooth → {sg_smooth, sg_deriv1, snv}
# SG smooth is applied once as the shared base; sg_deriv1 and snv are derived
# from the smoothed spectra, not from raw reflectance.
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

# ── CV folds (shared across preprocessing methods for fair comparison) ────────
folds <- make_folds(y_model, k = CV_FOLDS, repeats = CV_REPEATS,
                    stratified = CV_STRATIFIED, seed = SEED)

# ── 3-5. Per-preprocessing RF pipeline ───────────────────────────────────────
all_results   <- list()
all_cv_metrics <- data.frame()

for (prep_name in names(spec_list)) {

  cat(sprintf("\n══ Preprocessing: %s ══\n", toupper(prep_name)))
  X_prep <- spec_list[[prep_name]]

  # Remove any columns with zero variance (can occur with derivatives)
  col_var <- apply(X_prep, 2, var, na.rm = TRUE)
  X_prep  <- X_prep[, col_var > 1e-12, drop = FALSE]
  cat(sprintf("  Predictors after zero-var filter: %d\n", ncol(X_prep)))

  # ── Step 3: Tuning ──────────────────────────────────────────────────────────
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

  cat(sprintf("\nntree convergence — selected: %d\n", tune_res$ntree_sel))
  cat("\nStage 2 grid (ordered by OOB RMSE):\n")
  print(tune_res$stage2[, c("ntree","mtry","min_nodesize","oob_rmse","oob_r2")],
        row.names = FALSE, digits = 5)
  cat(sprintf("\nBest: ntree = %d | mtry = %d | min.node.size = %d | OOB RMSE = %.5f | OOB R² = %.4f\n\n",
              tune_res$best$ntree, tune_res$best$mtry, tune_res$best$min_nodesize,
              tune_res$best$oob_rmse, tune_res$best$oob_r2))

  # Save tuning grid plot
  p_tune <- plot_tuning_grid(tune_res$stage2,
                              title = sprintf("Tuning Grid — %s [%s]",
                                             TARGET_LABEL, prep_name))
  save_plot(p_tune,
            filename = sprintf("tuning_%s_%s", TARGET_PROPERTY, prep_name),
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

  # Per-resample metrics (30 resamples = 10 folds × 3 repeats)
  resamp_summ <- metrics_summary(cv_res$resample_metrics)
  cat(sprintf("\nCV Summary (%d resamples) [%s]:\n",
              length(cv_res$resample_metrics), prep_name))
  print_metrics(resamp_summ)

  # Aggregate for comparison across methods
  row_summ <- data.frame(
    preprocessing = prep_name,
    t(setNames(resamp_summ$mean, resamp_summ$metric)),
    stringsAsFactors = FALSE
  )
  # Add SD columns for error bars in plots
  for (m in resamp_summ$metric) {
    row_summ[[paste0(m, "_sd")]] <-
      resamp_summ$sd[resamp_summ$metric == m]
  }
  all_cv_metrics <- dplyr::bind_rows(all_cv_metrics, row_summ)

  # Obs vs Pred scatter
  p_op <- plot_obs_pred(
    cv_res$obs_pred,
    title = sprintf("Obs vs Pred — %s [%s]", TARGET_LABEL, prep_name),
    label = TARGET_LABEL
  )
  save_plot(p_op,
            filename = sprintf("obs_pred_%s_%s", TARGET_PROPERTY, prep_name),
            dir      = file.path(OUTPUT_DIR, "figures"),
            devices  = PLOT_DEVICES, dpi = PLOT_DPI)

  # Variable importance
  imp_wl  <- as.numeric(gsub("w", "", names(cv_res$avg_importance)))
  p_imp <- plot_importance(
    importance  = cv_res$avg_importance,
    wavelengths = imp_wl,
    title = sprintf("Variable Importance — %s [%s]", TARGET_LABEL, prep_name)
  )
  save_plot(p_imp,
            filename = sprintf("importance_%s_%s", TARGET_PROPERTY, prep_name),
            dir      = file.path(OUTPUT_DIR, "figures"),
            devices  = PLOT_DEVICES, dpi = PLOT_DPI)

  # ── Step 5: Final model on full data ─────────────────────────────────────
  cat("[5] Final Model (full dataset)\n")
  final_model <- fit_final_rf(
    X           = X_prep,
    y           = y_model,
    best_params = best_params,
    n_cores     = N_CORES,
    seed        = SEED
  )

  cat(sprintf("  OOB R²: %.4f | OOB RMSE: %.5f\n\n",
              final_model$r.squared, sqrt(final_model$prediction.error)))

  # Save results
  all_results[[prep_name]] <- list(
    tune        = tune_res,
    cv          = cv_res,
    resamp_summ = resamp_summ,
    final_model = final_model
  )

  if (SAVE_MODELS) {
    saveRDS(final_model,
            file.path(OUTPUT_DIR, "models",
                      sprintf("rf_final_%s_%s.rds",
                              TARGET_PROPERTY, prep_name)))
  }
}

# ── 6. Cross-method comparison ───────────────────────────────────────────────
cat("\n── Step 6: Cross-method performance comparison ─────\n\n")
print(all_cv_metrics[, c("preprocessing", "r2", "rmse", "mae", "me", "ccc", "rpd", "rpiq")],
      row.names = FALSE, digits = 4)

# Write full metrics to CSV
metrics_file <- file.path(OUTPUT_DIR, "metrics",
                           sprintf("cv_metrics_%s.csv", TARGET_PROPERTY))
readr::write_csv(all_cv_metrics, metrics_file)
cat(sprintf("\nMetrics saved: %s\n", metrics_file))

# Comparison bar charts (RPD and R²)
if ("rpd_sd" %in% names(all_cv_metrics)) {
  p_comp_rpd <- plot_cv_comparison(all_cv_metrics, metric = "rpd",
                                    title = sprintf("CV — RPD by Preprocessing [%s]",
                                                   TARGET_LABEL))
  save_plot(p_comp_rpd,
            filename = sprintf("comparison_rpd_%s", TARGET_PROPERTY),
            dir      = file.path(OUTPUT_DIR, "figures"),
            devices  = PLOT_DEVICES, dpi = PLOT_DPI)

  p_comp_r2 <- plot_cv_comparison(all_cv_metrics, metric = "r2",
                                   title = sprintf("CV — R² by Preprocessing [%s]",
                                                  TARGET_LABEL))
  save_plot(p_comp_r2,
            filename = sprintf("comparison_r2_%s", TARGET_PROPERTY),
            dir      = file.path(OUTPUT_DIR, "figures"),
            devices  = PLOT_DEVICES, dpi = PLOT_DPI)
}

# Save full results object
saveRDS(all_results,
        file.path(OUTPUT_DIR, "models",
                  sprintf("all_results_%s.rds", TARGET_PROPERTY)))

cat("\n══════════════════════════════════════════════════\n")
cat("  Pipeline complete. Outputs in:", OUTPUT_DIR, "\n")
cat("══════════════════════════════════════════════════\n\n")
