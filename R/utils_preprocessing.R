# =============================================================================
# utils_preprocessing.R — Spectral Preprocessing Functions
# =============================================================================
# Three preprocessing treatments, all derived from a common SG-smoothed base:
#
#   Raw reflectance (input)
#         │
#         ▼
#   SG Smooth  ──────────────────► Treatment 1: sg_smooth
#         │
#         ├──► SG 1st Derivative ► Treatment 2: sg_deriv1
#         │
#         └──► SNV               ► Treatment 3: snv
#
# The SG smoothing step is always applied first to reduce high-frequency noise
# before any further transformation. sg_deriv1 and snv are computed from the
# smoothed spectra, not from raw reflectance.
#
# All functions accept and return a numeric matrix (n × p):
#   rows    = samples
#   columns = wavelengths
# SG derivatives reduce the number of columns by (w − 1) due to edge removal.
# =============================================================================

library(prospectr)

# ── Primitive transforms ──────────────────────────────────────────────────────

#' Savitzky–Golay smoothing (0th derivative).
#' @param X Numeric matrix (n × p).
#' @param w Window size (odd integer).
#' @param p Polynomial order (p < w).
sg_smooth <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  prospectr::savitzkyGolay(X, w = w, p = p, m = 0)
}

#' Savitzky–Golay 1st derivative.
#' @param X Numeric matrix already smoothed.
#' @param w Window size (odd integer).
#' @param p Polynomial order (p < w).
sg_deriv1 <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  prospectr::savitzkyGolay(X, w = w, p = p, m = 1)
}

#' Standard Normal Variate.
#' @param X Numeric matrix already smoothed.
snv <- function(X) {
  stopifnot(is.matrix(X))
  prospectr::standardNormalVariate(X)
}

# ── Pipeline: builds all three treatments from raw reflectance ────────────────

#' Build the three spectral preprocessing treatments from raw reflectance.
#'
#' Applies SG smoothing once as the shared base, then derives each treatment:
#'   - sg_smooth : smoothed spectra (base)
#'   - sg_deriv1 : SG 1st derivative applied to smoothed spectra
#'   - snv       : Standard Normal Variate applied to smoothed spectra
#'
#' @param X_raw    Numeric matrix of raw reflectance (n × p).
#' @param sg_window Savitzky–Golay window size (odd integer; default 11).
#' @param sg_poly   Savitzky–Golay polynomial order (default 3).
#' @return Named list with three elements: sg_smooth, sg_deriv1, snv.
build_preprocessing_stack <- function(X_raw, sg_window = 11, sg_poly = 3) {
  stopifnot(is.matrix(X_raw))

  # Step 1 — shared SG smooth base
  X_smooth <- sg_smooth(X_raw, w = sg_window, p = sg_poly)

  # Step 2 — three treatments derived from the smooth base
  list(
    sg_smooth = X_smooth,
    sg_deriv1 = sg_deriv1(X_smooth, w = sg_window, p = sg_poly),
    snv       = snv(X_smooth)
  )
}
