# =============================================================================
# config.R — Global Configuration for Soil Spectral RF Modeling
# =============================================================================
# Edit this file to customize the modeling pipeline.
# All other scripts source this file at startup.
#
# Project: BMP Soil Spectral Modeling (ASD FieldSpec 4, 350–2500 nm)
# Targets: pH, TP, M3-P, WEP, TC, TN, Fe, Al, Ca, Mg
# =============================================================================

# ── Target soil property ──────────────────────────────────────────────────────
# Column name in the lab dataset.
# For OSSL prototype — change to project column names when your data are ready.
#
#   OSSL examples:
#     "oc_usda.c729_w.pct"          → Total organic carbon (%)
#     "ph.h2o_usda.a268_index"      → pH in water
#     "n.tot_usda.a623_w.pct"       → Total nitrogen (%)
#
#   BMP project targets (use exact column names from your lab database):
#     pH, TP, M3P, WEP, TC, TN, Fe, Al, Ca, Mg
TARGET_PROPERTY <- "oc_usda.c729_w.pct"
TARGET_LABEL    <- "Soil Organic Carbon (%)"   # used in figures and outputs

# To loop over all BMP targets, set this vector and iterate in main_rf_ossl.R:
# BMP_TARGETS <- c("pH", "TP", "M3P", "WEP", "TC", "TN", "Fe", "Al", "Ca", "Mg")

# ── Target transformation ─────────────────────────────────────────────────────
# Controls whether to model y on log1p scale, original scale, or both.
# Allowed values (single or combined):
#   "log"    — model log1p(y), back-transform predictions for all metrics
#   "no_log" — model y on original scale
#   c("log", "no_log") — run both and print a side-by-side comparison
# Recommended: c("log", "no_log") for right-skewed properties (OC, TN, TP…)
LOG_MODES <- c("log", "no_log")

# ── Spectral preprocessing methods ───────────────────────────────────────────
# Three treatments, all derived from a common SG-smooth base (see pipeline):
#
#   Treatment 1 — sg_smooth : SG smoothed reflectance (base)
#   Treatment 2 — sg_deriv1 : SG 1st derivative of the smoothed spectra
#   Treatment 3 — snv       : SNV applied to the smoothed spectra
#
# Do not edit this vector — the three treatments are fixed by design.
# Adjust SG parameters (window / poly) below if needed.
PREPROCESSING_METHODS <- c("sg_smooth", "sg_deriv1", "snv")

# Savitzky–Golay parameters (reported in model documentation)
SG_WINDOW <- 11   # window size (must be odd integer; e.g. 11 = ±5 bands)
SG_POLY   <- 3    # polynomial order (< SG_WINDOW)

# ── OSSL data (prototype) ─────────────────────────────────────────────────────
# Fraction of OSSL observations to use for prototyping (0 < FRACTION ≤ 1).
# Set to 1.0 when running with your own complete dataset.
OSSL_FRACTION  <- 0.20
# Hard cap on sample size after fraction sampling (set to Inf to disable).
# Prototype default: 1000 samples are sufficient to evaluate the pipeline.
N_MAX_SAMPLES  <- 1000
# ASD FieldSpec 4 usable range (nm) — edges trimmed to reduce noise
WAVE_MIN      <- 400
WAVE_MAX      <- 2450
# Local cache directory (OSSL files downloaded once and reused)
DATA_DIR      <- "data"

# ── Cross-validation ──────────────────────────────────────────────────────────
# Per proposal: 10-fold CV with 3 repetitions.
# Performance reported as mean ± SD across all (folds × repeats) resamples.
CV_FOLDS      <- 10
CV_REPEATS    <- 3      # 3 repeats × 10 folds = 30 resamples total
CV_STRATIFIED <- TRUE   # stratify folds by quantile bins of target
SEED          <- 2024L

# ── Random Forest hyperparameter tuning ───────────────────────────────────────
# Tuning is performed via OOB RMSE on the training set (grid search).
# All grid combinations are evaluated in parallel.
#
# ntree candidates: 500–2000 in steps of 100 (per proposal)
NTREE_CANDIDATES    <- seq(500, 2000, by = 100)   # evaluated during tuning
NUM_TREES_FINAL     <- 1000   # trees for final model after tuning confirms convergence

# mtry candidates: fractions 1/9, 1/6, 1/3, 1 of total predictors (per proposal)
MTRY_FRACS          <- c(1/9, 1/6, 1/3, 1)

# min.node.size candidates (not in original proposal; kept for completeness)
NODESIZE_CANDIDATES <- c(3, 5, 10)

# ── Parallelisation ───────────────────────────────────────────────────────────
# NULL = auto-detect (all available cores minus 1)
N_CORES <- NULL

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR   <- "output"
SAVE_MODELS  <- TRUE          # save fitted ranger objects (.qs)
PLOT_DEVICES <- c("png")      # "png", "pdf", or c("png", "pdf")
PLOT_DPI     <- 300
