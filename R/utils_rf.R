# =============================================================================
# utils_rf.R — Random Forest Tuning and Cross-Validation Utilities
# =============================================================================
# Implements the two-stage workflow specified in the BMP project proposal:
#
#   1. tune_rf()   — grid search over ntree × mtry × min.node.size using OOB
#                    error. Fully parallelised with doParallel / foreach.
#                    ntree range: 500–2000 (step 100)
#                    mtry fracs:  1/9, 1/6, 1/3, 1  of total predictors
#
#   2. cv_rf()     — 10-fold × 3-repeat stratified cross-validation with the
#                    best hyperparameters. Resamples run in parallel.
#                    Performance reported as mean ± SD across 30 resamples.
#
# =============================================================================

library(ranger)
library(foreach)
library(doParallel)

# ── CV fold construction ──────────────────────────────────────────────────────

#' Create (optionally stratified, optionally grouped) k-fold × repeat CV indices.
#'
#' When \code{groups} is supplied, fold assignment is performed at the group
#' level (e.g. layer UUID), so every row that belongs to the same group is
#' always assigned to the same fold.  This prevents leakage when multiple
#' spectral replicates share a single lab measurement.
#'
#' @param y          Numeric response vector.
#' @param k          Number of folds (default 10).
#' @param repeats    Number of repetitions (default 3).
#' @param stratified Logical; stratify by quantile bins of y (or group-mean y).
#' @param seed       Integer seed for reproducibility.
#' @param groups     Optional character/factor vector of group IDs (same length
#'                   as y).  When provided, folds are assigned by group.
#' @return List of length k×repeats. Each element is an integer vector of
#'         *validation* row indices for that resample.
make_folds <- function(y, k = 10, repeats = 3,
                       stratified = TRUE, seed = 2024L,
                       groups = NULL) {
  set.seed(seed)
  n     <- length(y)
  folds <- vector("list", k * repeats)
  idx   <- 0L

  for (r in seq_len(repeats)) {

    if (!is.null(groups)) {
      # ── Grouped fold assignment ──────────────────────────────────────────
      ugrp <- unique(as.character(groups))
      ng   <- length(ugrp)

      if (stratified) {
        grp_y  <- as.numeric(
          tapply(y, as.character(groups), mean, na.rm = TRUE)[ugrp]
        )
        bins   <- cut(grp_y,
                      breaks = quantile(grp_y, seq(0, 1, length.out = k + 1),
                                        na.rm = TRUE),
                      include.lowest = TRUE, labels = FALSE)
        bins[is.na(bins)] <- 1L
        grp_fold <- integer(ng)
        for (b in unique(bins)) {
          ii <- which(bins == b)
          grp_fold[ii] <- sample(rep_len(seq_len(k), length(ii)))
        }
      } else {
        grp_fold <- sample(rep_len(seq_len(k), ng))
      }
      names(grp_fold) <- ugrp
      fold_id <- unname(grp_fold[as.character(groups)])

    } else {
      # ── Row-level fold assignment ────────────────────────────────────────
      if (stratified) {
        bins <- cut(y,
                    breaks = quantile(y, seq(0, 1, length.out = k + 1),
                                      na.rm = TRUE),
                    include.lowest = TRUE, labels = FALSE)
        bins[is.na(bins)] <- 1L
        fold_id <- integer(n)
        for (b in unique(bins)) {
          ii <- which(bins == b)
          fold_id[ii] <- sample(rep_len(seq_len(k), length(ii)))
        }
      } else {
        fold_id <- sample(rep_len(seq_len(k), n))
      }
    }

    for (f in seq_len(k)) {
      idx          <- idx + 1L
      folds[[idx]] <- which(fold_id == f)
    }
  }
  folds
}

# ── Hyperparameter tuning ─────────────────────────────────────────────────────

#' Full 3D grid search over ntree × mtry × min.node.size using OOB RMSE.
#'
#' All combinations are evaluated in a single parallel sweep — no sequential
#' staging, no assumption of independence between parameters.
#'
#'   ntree candidates : seq(500, 2000, by = 100)  → 16 values  (configurable)
#'   mtry fracs       : 1/9, 1/6, 1/3, 1          → 4 values   (configurable)
#'   min.node.size    : 3, 5, 10                   → 3 values   (configurable)
#'   ──────────────────────────────────────────────────────────
#'   Total combinations (default): 16 × 4 × 3 = 192
#'
#' Each combination trains one ranger model (OOB only, importance = "none")
#' and records OOB RMSE and R². All combinations run in parallel via foreach.
#' The combination with the lowest OOB RMSE is selected.
#'
#' @param X               Numeric matrix of predictors (n × p).
#' @param y               Numeric response vector (length n).
#' @param ntree_candidates Integer vector of tree counts to evaluate.
#' @param mtry_fracs      Fractions of p to try as mtry (1/9, 1/6, 1/3, 1).
#' @param nodesize_vals   Integer vector of min.node.size candidates.
#' @param n_cores         Number of parallel workers (NULL = auto).
#' @param seed            Random seed.
#' @return List:
#'   $grid      — full grid data frame ordered by OOB RMSE (ascending)
#'   $ntree_sel — ntree of the best combination
#'   $best      — single-row data frame with the best hyperparameters
tune_rf <- function(X, y,
                    ntree_candidates = seq(500, 2000, by = 100),
                    mtry_fracs       = c(1/9, 1/6, 1/3, 1),
                    nodesize_vals    = c(3, 5, 10),
                    n_cores          = NULL,
                    seed             = 2026L) {

  p <- ncol(X)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)

  # Translate mtry fractions → actual integer values (at least 1)
  mtry_vals <- pmax(1L, unique(round(mtry_fracs * p)))

  # Full 3D grid
  grid <- expand.grid(
    ntree        = ntree_candidates,
    mtry         = mtry_vals,
    min_nodesize = nodesize_vals,
    stringsAsFactors = FALSE
  )
  n_comb <- nrow(grid)

  message(sprintf(
    "[Tuning] Full 3D grid — %d ntree × %d mtry × %d nodesize = %d combinations | Cores: %d",
    length(ntree_candidates), length(mtry_vals), length(nodesize_vals),
    n_comb, n_cores
  ))

  df_all <- as.data.frame(X)
  df_all[[".y"]] <- y

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ parallel::stopCluster(cl); doParallel::stopImplicitCluster() },
          add = TRUE)

  results <- foreach::foreach(
    i         = seq_len(n_comb),
    .packages = "ranger",
    .combine  = rbind
  ) %dopar% {
    set.seed(seed + i)
    fit <- ranger::ranger(
      formula       = .y ~ .,
      data          = df_all,
      num.trees     = grid$ntree[i],
      mtry          = grid$mtry[i],
      min.node.size = grid$min_nodesize[i],
      importance    = "none",
      verbose       = FALSE
    )
    c(oob_rmse = sqrt(fit$prediction.error),
      oob_r2   = fit$r.squared)
  }

  grid$oob_rmse <- results[, "oob_rmse"]
  grid$oob_r2   <- results[, "oob_r2"]
  grid          <- grid[order(grid$oob_rmse), ]
  best          <- grid[1L, ]

  message(sprintf(
    "  → Best: ntree = %d | mtry = %d | min.node.size = %d | OOB RMSE = %.4f | OOB R² = %.4f",
    best$ntree, best$mtry, best$min_nodesize, best$oob_rmse, best$oob_r2
  ))

  list(
    grid      = grid,
    ntree_sel = best$ntree,
    best      = best
  )
}

# ── Cross-validation ──────────────────────────────────────────────────────────

#' k-fold × repeated cross-validation of a ranger Random Forest.
#'
#' Each resample trains on (k−1)/k of the data and predicts the held-out fold.
#' All resamples run in parallel. Performance is reported per resample so that
#' mean ± SD can be computed across the full resampling distribution.
#'
#' @param X             Numeric matrix of predictors (n × p).
#' @param y             Numeric response vector (model scale; length n).
#' @param best_params   Named list: ntree, mtry, min_nodesize (from tune_rf).
#' @param folds         List of validation-set indices from make_folds().
#' @param n_cores       Number of parallel workers (NULL = auto).
#' @param seed          Random seed.
#' @param back_trans_fn Function to convert model-scale predictions back to
#'                      the original measurement scale (or NULL).
#' @return List:
#'   $obs_pred       — data frame (row_id, fold_id, repeat_id, obs, pred,
#'                     obs_mod, pred_mod)
#'   $resample_metrics — list of per-resample metric vectors (calc_metrics)
#'   $avg_importance — named numeric vector: mean variable importance
cv_rf <- function(X, y,
                  best_params,
                  folds,
                  n_cores       = NULL,
                  seed          = 2024L,
                  back_trans_fn = NULL) {

  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  n_resamples <- length(folds)

  message(sprintf(
    "[CV] %d resamples | ntree = %d | mtry = %d | min.node.size = %d | Cores = %d",
    n_resamples,
    best_params$ntree,
    best_params$mtry,
    best_params$min_nodesize,
    n_cores
  ))

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ parallel::stopCluster(cl); doParallel::stopImplicitCluster() },
          add = TRUE)

  resample_results <- foreach::foreach(
    f         = seq_len(n_resamples),
    .packages = "ranger",
    .combine  = "list",
    .multicombine = TRUE
  ) %dopar% {
    val_idx   <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(X)), val_idx)

    df_train <- as.data.frame(X[train_idx, , drop = FALSE])
    df_train[[".y"]] <- y[train_idx]

    set.seed(seed + f)
    fit <- ranger::ranger(
      formula       = .y ~ .,
      data          = df_train,
      num.trees     = best_params$ntree,
      mtry          = best_params$mtry,
      min.node.size = best_params$min_nodesize,
      importance    = "impurity",
      verbose       = FALSE
    )

    pred_val <- predict(fit,
                        data = as.data.frame(X[val_idx, , drop = FALSE]))$predictions

    list(
      obs       = y[val_idx],
      pred      = pred_val,
      fold_id   = f,
      row_idx   = val_idx,
      r2_oob    = fit$r.squared,
      rmse_oob  = sqrt(fit$prediction.error),
      importance = fit$variable.importance
    )
  }

  # ── Assemble full prediction table ────────────────────────────────────────
  obs_pred <- data.frame(
    row_id   = unlist(lapply(resample_results, `[[`, "row_idx")),
    fold_id  = unlist(lapply(resample_results, `[[`, "fold_id")),
    obs_mod  = unlist(lapply(resample_results, `[[`, "obs")),
    pred_mod = unlist(lapply(resample_results, `[[`, "pred"))
  )
  obs_pred <- obs_pred[order(obs_pred$row_id), ]

  # Back-transform to original scale
  if (!is.null(back_trans_fn)) {
    obs_pred$obs  <- back_trans_fn(obs_pred$obs_mod)
    obs_pred$pred <- back_trans_fn(obs_pred$pred_mod)
  } else {
    obs_pred$obs  <- obs_pred$obs_mod
    obs_pred$pred <- obs_pred$pred_mod
  }

  # ── Per-resample metrics (original scale) ─────────────────────────────────
  resample_metrics <- lapply(resample_results, function(rr) {
    obs_bt  <- if (!is.null(back_trans_fn)) back_trans_fn(rr$obs)  else rr$obs
    pred_bt <- if (!is.null(back_trans_fn)) back_trans_fn(rr$pred) else rr$pred
    calc_metrics(obs_bt, pred_bt,
                 r2_oob   = rr$r2_oob,
                 rmse_oob = rr$rmse_oob)
  })

  # ── Average variable importance across resamples ──────────────────────────
  imp_mat <- do.call(rbind, lapply(resample_results, `[[`, "importance"))
  avg_imp <- colMeans(imp_mat, na.rm = TRUE)

  list(
    obs_pred          = obs_pred,
    resample_metrics  = resample_metrics,
    avg_importance    = avg_imp
  )
}

# ── Final model ───────────────────────────────────────────────────────────────

#' Fit a final ranger model on the full dataset with the tuned hyperparameters.
#'
#' @param X           Numeric matrix of predictors (n × p).
#' @param y           Numeric response vector (model scale).
#' @param best_params Named list: ntree, mtry, min_nodesize.
#' @param n_cores     Number of ranger threads (NULL = auto).
#' @param seed        Random seed.
#' @return A fitted ranger object.
fit_final_rf <- function(X, y, best_params,
                         n_cores = NULL,
                         seed    = 2024L) {
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  df <- as.data.frame(X)
  df[[".y"]] <- y
  set.seed(seed)
  ranger::ranger(
    formula        = .y ~ .,
    data           = df,
    num.trees      = best_params$ntree,
    mtry           = best_params$mtry,
    min.node.size  = best_params$min_nodesize,
    importance     = "impurity",
    num.threads    = n_cores,
    verbose        = TRUE
  )
}
