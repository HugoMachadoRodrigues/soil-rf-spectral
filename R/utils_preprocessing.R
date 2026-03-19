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
#
# References:
#   Savitzky & Golay (1964)  — original SG filter:
#     Savitzky, A. & Golay, M.J.E. Smoothing and differentiation of data by
#     simplified least squares procedures. Analytical Chemistry, 36(8),
#     1627–1639. https://doi.org/10.1021/ac60214a047
#
#   Barnes et al. (1989)  — Standard Normal Variate (SNV):
#     Barnes, R.J., Dhanoa, M.S., & Lister, S.J. Standard normal variate
#     transformation and de-trending of near-infrared diffuse reflectance
#     spectra. Applied Spectroscopy, 43(5), 772–777.
#     https://doi.org/10.1366/0003702894202201
#
#   Rinnan et al. (2009)  — comprehensive review of NIR preprocessing:
#     Rinnan, Å., van den Berg, F., & Engelsen, S.B. Review of the most common
#     pre-processing techniques for near-infrared spectra. TrAC Trends in
#     Analytical Chemistry, 28(10), 1201–1222.
#     https://doi.org/10.1016/j.trac.2009.07.007
#
#   Viscarra Rossel & Behrens (2010)  — RF + derivative preprocessing for soil:
#     Viscarra Rossel, R.A. & Behrens, T. Using data mining to model and
#     interpret soil diffuse reflectance spectra. Geoderma, 158(1–2), 46–54.
#     https://doi.org/10.1016/j.geoderma.2009.12.025
# =============================================================================

library(prospectr)

# ── Primitive transforms ──────────────────────────────────────────────────────

#' Savitzky–Golay smoothing (0th derivative).
#'
#' Fits a polynomial of degree p through w consecutive points by least squares,
#' using the central fitted value as the smoothed output. This removes
#' high-frequency instrument noise while preserving band shape better than
#' a moving average (Savitzky & Golay 1964).
#'
#' @param X Numeric matrix (n × p).
#' @param w Window size (odd integer).
#' @param p Polynomial order (p < w).
sg_smooth <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  prospectr::savitzkyGolay(X, w = w, p = p, m = 0)
}

#' Savitzky–Golay 1st derivative.
#'
#' Computes the first derivative of the smoothed polynomial fit (m = 1),
#' removing additive baseline offsets and slope effects that arise from
#' variable path length or particle size differences (Rinnan et al. 2009).
#' Applied to already-smoothed spectra to avoid amplifying noise.
#'
#' @param X Numeric matrix already SG-smoothed.
#' @param w Window size (odd integer).
#' @param p Polynomial order (p < w).
sg_deriv1 <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  prospectr::savitzkyGolay(X, w = w, p = p, m = 1)
}

#' Standard Normal Variate (SNV).
#'
#' Centres and scales each spectrum independently to zero mean and unit
#' variance: xSNV = (x − mean(x)) / sd(x). This corrects multiplicative
#' scatter effects and path-length variation without requiring a reference
#' spectrum (Barnes et al. 1989). Preferred over MSC when a stable mean
#' spectrum cannot be defined (Rinnan et al. 2009).
#'
#' @param X Numeric matrix already SG-smoothed.
snv <- function(X) {
  stopifnot(is.matrix(X))
  prospectr::standardNormalVariate(X)
}

# ── Pipeline: builds all three treatments from raw reflectance ────────────────

#' Build the three spectral preprocessing treatments from raw reflectance.
#'
#' Applies SG smoothing once as the shared base, then derives each treatment.
#' The shared-base design (Rinnan et al. 2009) ensures that differences in
#' model performance are attributable to the transformation step alone, not to
#' differences in the underlying smoothing.
#'
#'   Treatment 1 — sg_smooth : SG-smoothed reflectance (base; Savitzky & Golay 1964)
#'   Treatment 2 — sg_deriv1 : SG 1st derivative of smoothed spectra (Rinnan et al. 2009)
#'   Treatment 3 — snv       : SNV of smoothed spectra (Barnes et al. 1989)
#'
#' @param X_raw    Numeric matrix of raw reflectance (n × p).
#' @param sg_window Savitzky–Golay window size (odd integer; default 11).
#' @param sg_poly   Savitzky–Golay polynomial order (default 3).
#' @return Named list with three elements: sg_smooth, sg_deriv1, snv.
build_preprocessing_stack <- function(X_raw, sg_window = 11, sg_poly = 3) {
  stopifnot(is.matrix(X_raw))

  # Step 1 — shared SG smooth base (Savitzky & Golay 1964)
  X_smooth <- sg_smooth(X_raw, w = sg_window, p = sg_poly)

  # Step 2 — three treatments derived from the smooth base
  list(
    sg_smooth = X_smooth,                                       # Treatment 1
    sg_deriv1 = sg_deriv1(X_smooth, w = sg_window, p = sg_poly), # Treatment 2
    snv       = snv(X_smooth)                                  # Treatment 3
  )
}
