# =============================================================================
# prepare_report_data.R
# =============================================================================
# Extracts the 1,000-sample modelling dataset from cached OSSL files,
# builds the preprocessing stack, and saves a compact report_data.rds
# for fast loading by the Rmd report.
#
# Run this ONCE from the project root before knitting the report:
#   source("reports/report_01_soc_ossl_prototype/prepare_report_data.R")
# =============================================================================

message("── Preparing report data ───────────────────────────────────────────")

# ── 0. Resolve project root (works regardless of caller's working directory) ──
.script_dir  <- normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE)
# Fallback when sourced interactively without a frame (e.g. source() with no env)
if (!nzchar(.script_dir) || .script_dir == ".") {
  .script_dir <- normalizePath(
    dirname(rstudioapi::getSourceEditorContext()$path), mustWork = FALSE
  )
}
proj_root    <- normalizePath(file.path(.script_dir, "../.."), mustWork = TRUE)
message("  Project root: ", proj_root)

# ── 1. Source project config and utilities ────────────────────────────────────
source(file.path(proj_root, "R/config.R"))
source(file.path(proj_root, "R/utils_preprocessing.R"))

# ── 2. Load OSSL from local cache ─────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(readr)

abs_data_dir <- file.path(proj_root, DATA_DIR)
message("  Loading spectral data (", abs_data_dir, ") …")
visnir_raw <- readr::read_csv(
  file.path(abs_data_dir, "ossl_visnir_L0_v1.2.csv.gz"),
  show_col_types = FALSE
)
soillab <- readr::read_csv(
  file.path(abs_data_dir, "ossl_soillab_L1_v1.2.csv.gz"),
  show_col_types = FALSE
)

# ── 3. Wavelength filter ───────────────────────────────────────────────────────
wave_cols <- grep("^scan_visnir\\.(\\d+)_ref$", names(visnir_raw), value = TRUE)
wave_nm   <- as.numeric(gsub("scan_visnir\\.(\\d+)_ref", "\\1", wave_cols))
keep_cols <- wave_cols[wave_nm >= WAVE_MIN & wave_nm <= WAVE_MAX]
wave_nm   <- wave_nm[wave_nm  >= WAVE_MIN & wave_nm <= WAVE_MAX]

# ── 4. Completeness filter + UUID counts ──────────────────────────────────────
spec_df <- visnir_raw %>%
  dplyr::select(id.layer_uuid_txt, all_of(keep_cols)) %>%
  tidyr::drop_na()

n_rows_total <- nrow(spec_df)
n_uuid_total <- dplyr::n_distinct(spec_df$id.layer_uuid_txt)
message(sprintf("  Complete spectral rows: %d (%d unique layers)", n_rows_total, n_uuid_total))

# ── 5. Join with lab data ─────────────────────────────────────────────────────
lab_sub <- soillab %>%
  dplyr::select(id.layer_uuid_txt, all_of(TARGET_PROPERTY)) %>%
  dplyr::filter(!is.na(.data[[TARGET_PROPERTY]]))

joined       <- dplyr::inner_join(spec_df, lab_sub, by = "id.layer_uuid_txt")
n_joined     <- nrow(joined)
n_uuid_joined <- dplyr::n_distinct(joined$id.layer_uuid_txt)
message(sprintf("  After SOC join: %d rows (%d unique layers)", n_joined, n_uuid_joined))

# ── 6. Reproducible subsample (same seed as main script) ──────────────────────
set.seed(SEED)
sampled <- dplyr::slice_sample(joined, prop = OSSL_FRACTION)
if (!is.null(N_MAX_SAMPLES) && is.finite(N_MAX_SAMPLES) &&
    nrow(sampled) > N_MAX_SAMPLES) {
  sampled <- dplyr::slice_sample(sampled, n = N_MAX_SAMPLES)
}
n_final      <- nrow(sampled)
n_uuid_final <- dplyr::n_distinct(sampled$id.layer_uuid_txt)
message(sprintf("  Modelling sample: %d rows (%d unique layers)", n_final, n_uuid_final))

# ── 7. Build matrices and preprocessing stack ─────────────────────────────────
X_raw   <- as.matrix(sampled[, keep_cols])
colnames(X_raw) <- paste0("w", wave_nm)
y_orig  <- sampled[[TARGET_PROPERTY]]
uuid_id <- sampled$id.layer_uuid_txt

message("  Building preprocessing stack …")
spec_list <- build_preprocessing_stack(
  X_raw     = X_raw,
  sg_window = SG_WINDOW,
  sg_poly   = SG_POLY
)

# ── 8. Save compact report data ───────────────────────────────────────────────
out_path <- file.path(.script_dir, "report_data.rds")
saveRDS(
  list(
    X_raw         = X_raw,
    spec_list     = spec_list,
    y_orig        = y_orig,
    uuid_id       = uuid_id,
    wave_nm       = wave_nm,
    n_wave        = length(wave_nm),
    n_rows_total  = n_rows_total,
    n_uuid_total  = n_uuid_total,
    avg_scans     = n_rows_total / n_uuid_total,
    n_joined      = n_joined,
    n_uuid_joined = n_uuid_joined,
    n_final       = n_final,
    n_uuid_final  = n_uuid_final,
    WAVE_MIN      = WAVE_MIN,
    WAVE_MAX      = WAVE_MAX,
    OSSL_FRACTION = OSSL_FRACTION,
    N_MAX_SAMPLES = N_MAX_SAMPLES,
    SG_WINDOW     = SG_WINDOW,
    SG_POLY       = SG_POLY,
    SEED          = SEED,
    TARGET_PROPERTY = TARGET_PROPERTY,
    TARGET_LABEL    = TARGET_LABEL
  ),
  file = out_path
)

sz <- round(file.size(out_path) / 1024^2, 1)
message(sprintf("  Saved: %s  (%.1f MB)", out_path, sz))
message("── Done. Now knit the .Rmd report. ────────────────────────────────")
