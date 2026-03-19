# =============================================================================
# config.R — Global Configuration for Soil Spectral RF Modeling
# =============================================================================
# Edit this file to customize the modeling pipeline.
# All other scripts source this file at startup.
# =============================================================================

# ── Target soil property ──────────────────────────────────────────────────────
# Column name in the OSSL lab dataset (ossl_soillab)
# Examples: "oc_usda.c729_w.pct"  (organic carbon, %)
#           "clay.tot_usda.a334_w.pct" (total clay, %)
#           "ph.h2o_usda.a268_index" (pH in water)
TARGET_PROPERTY <- "oc_usda.c729_w.pct"
TARGET_LABEL    <- "Soil Organic Carbon (%)"   # used in figures and outputs

# ── Target transformation ─────────────────────────────────────────────────────
# Apply log1p() to the target before modeling and exp()-1 back-transform
# afterward. Useful for skewed properties (e.g., OC, TN).
USE_LOG_TARGET <- TRUE

# ── Spectral preprocessing methods ───────────────────────────────────────────
# Any combination of: "raw", "sg_smooth", "sg_deriv1", "snv", "absorbance"
# "absorbance" = log(1/R); can be applied alone or as a base for SG/SNV.
PREPROCESSING_METHODS <- c("sg_smooth", "sg_deriv1", "snv", "absorbance")

# Savitzky–Golay parameters (applied to all SG methods)
SG_WINDOW   <- 11   # window size (odd integer)
SG_POLY     <- 3    # polynomial order

# ── OSSL data ─────────────────────────────────────────────────────────────────
# Fraction of OSSL observations to use (0 < FRACTION ≤ 1)
OSSL_FRACTION  <- 0.20
# Minimum wavelength (nm) and maximum wavelength (nm) — ASD FieldSpec range
WAVE_MIN       <- 400
WAVE_MAX       <- 2450
# Local cache directory (data files are downloaded once and reused)
DATA_DIR       <- here::here("data")

# ── Cross-validation ──────────────────────────────────────────────────────────
CV_FOLDS       <- 10
CV_REPEATS     <- 1      # increase for more stable estimates (costs time)
CV_STRATIFIED  <- TRUE   # stratify folds by quantile of target
SEED           <- 2024L

# ── Random Forest hyperparameter tuning ───────────────────────────────────────
# Tuning is performed via OOB error on the training set using a grid search.
# Each combination is evaluated in parallel.
NUM_TREES_TUNE  <- 500            # trees used during grid search (fast)
NUM_TREES_FINAL <- 1000           # trees for the final model
MTRY_CANDIDATES <- c(0.05, 0.10, 0.20, 0.33)  # fractions of total predictors
NODESIZE_CANDIDATES <- c(3, 5, 10, 20)         # min samples per terminal node

# ── Parallelisation ───────────────────────────────────────────────────────────
# Number of CPU cores. Set to NULL to use all available minus 1.
N_CORES <- NULL

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR     <- here::here("output")
SAVE_MODELS    <- TRUE   # save fitted ranger objects (can be large)
PLOT_DEVICES   <- c("png")  # "png", "pdf", or both
PLOT_DPI       <- 300
