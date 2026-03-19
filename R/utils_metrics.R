# =============================================================================
# utils_metrics.R — Performance Metrics for Soil Spectral Models
# =============================================================================
# Implements the full set of metrics required by the BMP project proposal:
#   RMSE, MAE, R², ME (bias), RPD, RPIQ, CCC
# plus the MSD decomposition framework (Gauch et al. 2003) inspired by
# Clingensmith & Grunwald (NRCS Soil Spectral Modeling Project):
#   MSE = SB + NU + LC
#
# Formulas
# ─────────────────────────────────────────────────────────────────────────────
#
#  Let:  n   = number of observations
#        yᵢ  = observed value for sample i
#        ŷᵢ  = predicted value for sample i
#        ȳ   = mean of observed values
#        ȳ̂   = mean of predicted values
#        s_y = SD of observed values   [n-1 denominator]
#        s_ŷ = SD of predicted values  [n-1 denominator]
#        r   = Pearson correlation coefficient
#
#  ── Primary metrics ──────────────────────────────────────────────────────────
#
#  ME (Mean Error / Bias):
#    ME = (1/n) Σ (ŷᵢ − yᵢ)
#    Signed; positive = over-prediction, negative = under-prediction.
#
#  RMSE (Root Mean Squared Error):
#    RMSE = √[ (1/n) Σ (ŷᵢ − yᵢ)² ]
#
#  MAE (Mean Absolute Error):
#    MAE = (1/n) Σ |ŷᵢ − yᵢ|
#
#  R² (Coefficient of Determination):
#    R² = 1 − [ Σ (yᵢ − ŷᵢ)² / Σ (yᵢ − ȳ)² ]
#    (computed from regression of observed on predicted)
#
#  RPD (Residual Prediction Deviation):
#    RPD = s_y / RMSE
#    Interpretation (Chang et al. 2001):
#      RPD < 1.5  → poor model
#      1.5–2.0    → fair
#      2.0–2.5    → good
#      > 2.5      → excellent
#
#  RPIQ (Ratio of Performance to InterQuartile range; Bellon-Maurel et al. 2010):
#    RPIQ = (Q3 − Q1) / RMSE
#    Preferred over RPD for skewed distributions.
#
#  CCC (Lin's Concordance Correlation Coefficient; Lin 1989):
#    CCC = 2 · s_yŷ / (s²_y + s²_ŷ + (ȳ − ȳ̂)²)
#    where s_yŷ = covariance(y, ŷ)  [n denominator for consistency]
#    CCC ∈ [−1, 1]; combines precision (r) and accuracy (distance from 1:1).
#
#  ── MSD decomposition (Gauch et al. 2003) ────────────────────────────────────
#
#  MSE = SB + NU + LC
#
#  SB  (Squared Bias):         SB  = (ȳ̂ − ȳ)²
#  NU  (Non-Unity slope):      NU  = (1 − b)² · σ²_ŷ      [n denominator]
#  LC  (Lack of Correlation):  LC  = (1 − r²) · σ²_y       [n denominator]
#
#  where b is the slope of the OLS regression  ŷ = a + b · y.
#  A well-calibrated model has SB ≈ 0, NU ≈ 0, LC ≈ 0.
#
# =============================================================================

library(moments)   # skewness / kurtosis

# ── CCC helper ────────────────────────────────────────────────────────────────
#' Lin's Concordance Correlation Coefficient.
#' @param obs  Numeric vector of observed values.
#' @param pred Numeric vector of predicted values.
#' @return Scalar CCC in [−1, 1].
lin_ccc <- function(obs, pred) {
  n     <- length(obs)
  mu_y  <- mean(obs);   mu_p  <- mean(pred)
  s2_y  <- var(obs)  * (n - 1) / n   # population variance (n denominator)
  s2_p  <- var(pred) * (n - 1) / n
  s_yp  <- cov(obs, pred) * (n - 1) / n
  2 * s_yp / (s2_y + s2_p + (mu_y - mu_p)^2)
}

# ── Main metrics function ─────────────────────────────────────────────────────
#' Compute the full metric suite for a single observed/predicted pair.
#'
#' @param obs      Numeric vector of observed values (original scale).
#' @param pred     Numeric vector of predicted values (same scale as obs).
#' @param r2_oob   OOB R² from ranger (scalar; NA if not available).
#' @param rmse_oob OOB RMSE from ranger (scalar; NA if not available).
#' @return Named numeric vector:
#'   me, rmse, mae, r2, ccc, rpd, rpiq,
#'   r, intercept, slope, r2_oob, rmse_oob,
#'   sb, nu, lc, rmse_c, n
calc_metrics <- function(obs, pred,
                         r2_oob   = NA_real_,
                         rmse_oob = NA_real_) {

  stopifnot(length(obs) == length(pred))
  idx  <- complete.cases(obs, pred)
  obs  <- obs[idx];  pred <- pred[idx]
  n    <- length(obs)

  if (n < 3) {
    warning("Too few observations to compute metrics (n < 3).")
    nms <- c("me","rmse","mae","r2","ccc","rpd","rpiq",
             "r","intercept","slope","r2_oob","rmse_oob",
             "sb","nu","lc","rmse_c","n")
    return(setNames(rep(NA_real_, length(nms)), nms))
  }

  # ── Regression of obs on pred ──────────────────────────────────────────────
  lm_fit    <- lm(obs ~ pred)
  r         <- cor(obs, pred, method = "pearson")
  intercept <- unname(coef(lm_fit)[1])
  slope     <- unname(coef(lm_fit)[2])
  r2        <- summary(lm_fit)$r.squared

  # ── Primary metrics ────────────────────────────────────────────────────────
  me    <- mean(pred - obs)                           # Mean Error (bias)
  mse   <- mean((pred - obs)^2)
  rmse  <- sqrt(mse)
  mae   <- mean(abs(pred - obs))
  ccc   <- lin_ccc(obs, pred)

  # ── Relative performance ───────────────────────────────────────────────────
  rpd  <- sd(obs) / rmse
  iqr  <- as.numeric(diff(quantile(obs, c(0.25, 0.75), na.rm = TRUE)))
  rpiq <- iqr / rmse

  # ── MSD decomposition (Gauch et al. 2003) ──────────────────────────────────
  # Population variances (n denominator) for additive consistency
  sigma2_y <- var(obs)  * (n - 1) / n
  sigma2_p <- var(pred) * (n - 1) / n

  sb    <- (mean(pred) - mean(obs))^2           # Squared Bias
  nu    <- (1 - slope)^2 * sigma2_p             # Non-Unity slope
  lc    <- (1 - r^2)    * sigma2_y              # Lack of Correlation
  rmse_c <- sqrt(mse - sb)                      # bias-corrected RMSE

  round(c(
    me        = me,
    rmse      = rmse,
    mae       = mae,
    r2        = r2,
    ccc       = ccc,
    rpd       = rpd,
    rpiq      = rpiq,
    r         = r,
    intercept = intercept,
    slope     = slope,
    r2_oob    = r2_oob,
    rmse_oob  = rmse_oob,
    sb        = sb,
    nu        = nu,
    lc        = lc,
    rmse_c    = rmse_c,
    n         = n
  ), 2)
}

# ── Aggregate fold metrics ─────────────────────────────────────────────────────
#' Aggregate per-resample metric vectors into mean ± SD summary.
#'
#' @param resample_metrics List of named numeric vectors from calc_metrics().
#' @return Data frame with columns: metric, mean, sd.
metrics_summary <- function(resample_metrics) {
  mat  <- do.call(rbind, resample_metrics)
  mn   <- colMeans(mat, na.rm = TRUE)
  sds  <- apply(mat, 2, sd, na.rm = TRUE)
  data.frame(
    metric = names(mn),
    mean   = round(mn,  2),
    sd     = round(sds, 2),
    row.names = NULL
  )
}

#' Print a formatted metric summary to the console.
#' @param df Data frame from metrics_summary() or a named numeric vector.
print_metrics <- function(df) {
  if (is.numeric(df)) df <- data.frame(metric = names(df), value = df)
  cat("\n", strrep("─", 56), "\n", sep = "")
  print(df, row.names = FALSE, digits = 2)
  cat(strrep("─", 56), "\n\n", sep = "")
}

# ── Descriptive statistics ────────────────────────────────────────────────────
#' Compute summary statistics for a soil property vector.
#'
#' @param y Numeric vector (original scale).
#' @return Named numeric vector.
summary_stats <- function(y) {
  y <- na.omit(y)
  q <- quantile(y, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  c(
    min    = unname(q[1]),
    q1     = unname(q[2]),
    median = unname(q[3]),
    mean   = mean(y),
    q3     = unname(q[4]),
    max    = unname(q[5]),
    sd     = sd(y),
    cv     = sd(y) / mean(y),
    skew   = moments::skewness(y),
    kurt   = moments::kurtosis(y),
    n      = length(y),
    n_na   = sum(is.na(y))
  )
}
