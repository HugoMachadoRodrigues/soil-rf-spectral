# =============================================================================
# utils_metrics.R — Performance Metrics for Soil Spectral Models
# =============================================================================
# Inspired by the MSD decomposition approach of Clingensmith & Grunwald (NRCS).
#
# calc_metrics()  — full metric suite for a single observed/predicted pair
# metrics_table() — aggregate metrics across CV folds
# =============================================================================

library(moments)  # skewness / kurtosis

#' Compute a comprehensive suite of regression metrics.
#'
#' @param obs  Numeric vector of observed values (original scale after
#'             back-transformation if log was used).
#' @param pred Numeric vector of predicted values (same scale as obs).
#' @param r2_oob  OOB R² from ranger (optional; NA if not available).
#' @param rmse_oob OOB RMSE from ranger (optional; NA if not available).
#' @return Named numeric vector with the following metrics:
#'   r, int, slope, r2, r2_oob, rmse_oob, bias, rmse, mae,
#'   sb, nu, lc, rmse_c, rpd, rpiq, n
calc_metrics <- function(obs, pred, r2_oob = NA_real_, rmse_oob = NA_real_) {
  stopifnot(length(obs) == length(pred))

  idx <- complete.cases(obs, pred)
  obs  <- obs[idx]
  pred <- pred[idx]
  n    <- length(obs)

  if (n < 3) {
    warning("Too few observations to compute metrics.")
    return(rep(NA_real_, 16))
  }

  # ── Core regression ────────────────────────────────────────────────────────
  lm_fit <- lm(obs ~ pred)
  r       <- cor(obs, pred, method = "pearson")
  int     <- coef(lm_fit)[1]
  slope   <- coef(lm_fit)[2]
  r2      <- summary(lm_fit)$r.squared

  # ── Error statistics ───────────────────────────────────────────────────────
  bias   <- mean(pred) - mean(obs)
  mse    <- mean((pred - obs)^2)
  rmse   <- sqrt(mse)
  mae    <- mean(abs(pred - obs))
  mse_c  <- mean((pred - bias - obs)^2)   # bias-corrected MSE
  rmse_c <- sqrt(mse_c)

  # ── MSD decomposition (Gauch et al. 2003) ──────────────────────────────────
  sb  <- bias^2                                        # squared bias
  nu  <- ((1 - slope)^2) * (var(pred) * (n - 1) / n)  # non-unity
  lc  <- (1 - r^2) * (var(obs) * (n - 1) / n)         # lack of correlation

  # ── Relative performance indices ───────────────────────────────────────────
  rpd  <- sd(obs) / rmse
  iqr  <- as.numeric(diff(quantile(obs, c(0.25, 0.75))))
  rpiq <- iqr / rmse

  metrics <- round(c(
    r        = r,
    int      = int,
    slope    = slope,
    r2       = r2,
    r2_oob   = r2_oob,
    rmse_oob = rmse_oob,
    bias     = bias,
    rmse     = rmse,
    mae      = mae,
    sb       = sb,
    nu       = nu,
    lc       = lc,
    rmse_c   = rmse_c,
    rpd      = rpd,
    rpiq     = rpiq,
    n        = n
  ), 6)

  metrics
}

#' Aggregate per-fold metric lists into a summary data frame.
#'
#' @param fold_metrics List of named numeric vectors from calc_metrics().
#' @return Data frame with mean ± SD across folds.
metrics_summary <- function(fold_metrics) {
  mat  <- do.call(rbind, fold_metrics)
  mn   <- colMeans(mat, na.rm = TRUE)
  sds  <- apply(mat, 2, sd, na.rm = TRUE)
  data.frame(
    metric = names(mn),
    mean   = round(mn, 4),
    sd     = round(sds, 4),
    row.names = NULL
  )
}

#' Print a formatted metric table to the console.
#'
#' @param df  Data frame returned by metrics_summary() or a named numeric vector.
print_metrics <- function(df) {
  if (is.numeric(df)) df <- data.frame(metric = names(df), value = df)
  cat("\n", strrep("─", 50), "\n", sep = "")
  print(df, row.names = FALSE, digits = 4)
  cat(strrep("─", 50), "\n\n", sep = "")
}

#' Compute summary statistics for a soil property vector.
#'
#' @param y Numeric vector.
#' @return Named numeric vector (min, Q1, median, mean, Q3, max, sd, cv, skew, kurt).
summary_stats <- function(y) {
  y <- na.omit(y)
  q <- quantile(y, c(0, 0.25, 0.5, 0.75, 1))
  c(
    min    = q[1], q1 = q[2], median = q[3], mean  = mean(y),
    q3     = q[4], max = q[5],
    sd     = sd(y),
    cv     = sd(y) / mean(y),
    skew   = skewness(y),
    kurt   = kurtosis(y),
    n      = length(y),
    n_na   = sum(is.na(y))
  )
}
