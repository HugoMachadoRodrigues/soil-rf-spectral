# =============================================================================
# utils_preprocessing.R — Spectral Preprocessing Functions
# =============================================================================
# Implements five preprocessing strategies on a numeric matrix X (n × p)
# where rows = samples and columns = wavelengths.
#
#   1. raw        – raw reflectance (identity, no transform)
#   2. sg_smooth  – Savitzky–Golay smoothing (0th derivative)
#   3. sg_deriv1  – Savitzky–Golay 1st derivative
#   4. snv        – Standard Normal Variate
#   5. absorbance – log(1 / R) — pseudo-absorbance from reflectance
#
# All functions accept and return a numeric matrix with the same dimnames,
# except SG derivatives which reduce the number of columns by (w - 1).
# =============================================================================

library(prospectr)

# ── 1. Raw reflectance ────────────────────────────────────────────────────────
preprocess_raw <- function(X) {
  stopifnot(is.matrix(X))
  X
}

# ── 2. Savitzky–Golay smoothing (m = 0) ──────────────────────────────────────
preprocess_sg_smooth <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  out <- prospectr::savitzkyGolay(X, w = w, p = p, m = 0)
  out
}

# ── 3. Savitzky–Golay 1st derivative (m = 1) ─────────────────────────────────
preprocess_sg_deriv1 <- function(X, w = 11, p = 3) {
  stopifnot(is.matrix(X), w %% 2 == 1, p < w)
  out <- prospectr::savitzkyGolay(X, w = w, p = p, m = 1)
  out
}

# ── 4. Standard Normal Variate ────────────────────────────────────────────────
preprocess_snv <- function(X) {
  stopifnot(is.matrix(X))
  out <- prospectr::standardNormalVariate(X)
  out
}

# ── 5. Absorbance: log(1 / R) ─────────────────────────────────────────────────
preprocess_absorbance <- function(X) {
  stopifnot(is.matrix(X))
  if (any(X <= 0, na.rm = TRUE)) {
    warning("Non-positive reflectance values found; replacing with small epsilon before log transform.")
    X[X <= 0] <- .Machine$double.eps
  }
  out <- log(1 / X)
  dimnames(out) <- dimnames(X)
  out
}

# ── Dispatcher ────────────────────────────────────────────────────────────────
#' Apply a named preprocessing method to a reflectance matrix.
#'
#' @param X         Numeric matrix (n × p), rows = samples, cols = wavelengths.
#' @param method    Character: one of "raw", "sg_smooth", "sg_deriv1",
#'                  "snv", "absorbance".
#' @param sg_window Savitzky–Golay window size (odd integer).
#' @param sg_poly   Savitzky–Golay polynomial order.
#' @return          Preprocessed numeric matrix.
apply_preprocessing <- function(X, method,
                                sg_window = 11, sg_poly = 3) {
  method <- match.arg(method,
                      c("raw", "sg_smooth", "sg_deriv1", "snv", "absorbance"))
  switch(method,
    raw        = preprocess_raw(X),
    sg_smooth  = preprocess_sg_smooth(X, w = sg_window, p = sg_poly),
    sg_deriv1  = preprocess_sg_deriv1(X, w = sg_window, p = sg_poly),
    snv        = preprocess_snv(X),
    absorbance = preprocess_absorbance(X)
  )
}

#' Apply all requested preprocessing methods and return a named list.
#'
#' @param X       Numeric matrix of raw reflectance.
#' @param methods Character vector of method names.
#' @param ...     Forwarded to apply_preprocessing().
#' @return Named list, one element per method.
preprocess_all <- function(X, methods, ...) {
  out <- lapply(methods, function(m) apply_preprocessing(X, m, ...))
  names(out) <- methods
  out
}
