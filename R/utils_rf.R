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

#' Create (optionally stratified) k-fold × repeat CV indices.
#'
#' @param y          Numeric response vector.
#' @param k          Number of folds (default 10).
#' @param repeats    Number of repetitions (default 3).
#' @param stratified Logical; stratify by quantile bins of y.
#' @param seed       Integer seed for reproducibility.
#' @return List of length k×repeats. Each element is an integer vector of
#'         *validation* row indices for that resample.
make_folds <- function(y, k = 10, repeats = 3,
                       stratified = TRUE, seed = 2024L) {
  set.seed(seed)
  n     <- length(y)
  folds <- vector("list", k * repeats)
  idx   <- 0L

  for (r in seq_len(repeats)) {
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
    for (f in seq_len(k)) {
      idx          <- idx + 1L
      folds[[idx]] <- which(fold_id == f)
    }
  }
  folds
}

# ── Hyperparameter tuning ─────────────────────────────────────────────────────

#' Two-stage grid search: (1) ntree convergence, (2) mtry × min.node.size.
#'
#' Stage 1 evaluates ntree candidates (500–2000, step 100 per proposal) at
#' default mtry. The smallest ntree whose OOB RMSE is within 0.1% of the
#' minimum is selected as the operating ntree for Stage 2.
#'
#' Stage 2 evaluates all mtry × min.node.size combinations at the chosen
#' ntree. OOB RMSE is the selection criterion throughout.
#'
#' @param X               Numeric matrix of predictors (n × p).
#' @param y               Numeric response vector (length n).
#' @param ntree_candidates Integer vector of tree counts to evaluate.
#' @param mtry_fracs      Fractions of p to try as mtry (1/9, 1/6, 1/3, 1).
#' @param nodesize_vals   Integer vector of min.node.size candidates.
#' @param n_cores         Number of parallel workers (NULL = auto).
#' @param seed            Random seed.
#' @return List:
#'   $stage1      — data frame of ntree vs. OOB RMSE
#'   $ntree_sel   — selected ntree
#'   $stage2      — data frame of full grid ordered by OOB RMSE
#'   $best        — single-row data frame with best hyperparameters
tune_rf <- function(X, y,
                    ntree_candidates = seq(500, 2000, by = 100),
                    mtry_fracs       = c(1/9, 1/6, 1/3, 1),
                    nodesize_vals    = c(3, 5, 10),
                    n_cores          = NULL,
                    seed             = 2024L) {

  p <- ncol(X)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ parallel::stopCluster(cl); doParallel::stopImplicitCluster() },
          add = TRUE)

  df_all <- as.data.frame(X)
  df_all[[".y"]] <- y

  # ── Stage 1: ntree convergence ─────────────────────────────────────────────
  message(sprintf(
    "[Tuning | Stage 1] ntree convergence — %d candidates | Cores: %d",
    length(ntree_candidates), n_cores
  ))

  s1 <- foreach::foreach(
    nt        = ntree_candidates,
    .packages = "ranger",
    .combine  = rbind
  ) %dopar% {
    set.seed(seed)
    fit <- ranger::ranger(
      formula       = .y ~ .,
      data          = df_all,
      num.trees     = nt,
      importance    = "none",
      verbose       = FALSE
    )
    c(ntree    = nt,
      oob_rmse = sqrt(fit$prediction.error),
      oob_r2   = fit$r.squared)
  }
  s1 <- as.data.frame(s1)

  # Select smallest ntree within 0.1% of minimum OOB RMSE
  min_rmse  <- min(s1$oob_rmse, na.rm = TRUE)
  ntree_sel <- min(s1$ntree[s1$oob_rmse <= min_rmse * 1.001])
  message(sprintf("  → Selected ntree = %d (OOB RMSE = %.6f)", ntree_sel, min_rmse))

  # ── Stage 2: mtry × min.node.size grid ────────────────────────────────────
  mtry_vals <- pmax(1L, unique(round(mtry_fracs * p)))
  grid <- expand.grid(
    ntree        = ntree_sel,
    mtry         = mtry_vals,
    min_nodesize = nodesize_vals,
    oob_rmse     = NA_real_,
    oob_r2       = NA_real_,
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "[Tuning | Stage 2] mtry × min.node.size — %d combinations | ntree = %d",
    nrow(grid), ntree_sel
  ))

  s2 <- foreach::foreach(
    i         = seq_len(nrow(grid)),
    .packages = "ranger",
    .combine  = rbind
  ) %dopar% {
    set.seed(seed + i)
    fit <- ranger::ranger(
      formula       = .y ~ .,
      data          = df_all,
      num.trees     = ntree_sel,
      mtry          = grid$mtry[i],
      min.node.size = grid$min_nodesize[i],
      importance    = "none",
      verbose       = FALSE
    )
    c(oob_rmse = sqrt(fit$prediction.error),
      oob_r2   = fit$r.squared)
  }

  grid$oob_rmse <- s2[, "oob_rmse"]
  grid$oob_r2   <- s2[, "oob_r2"]
  grid          <- grid[order(grid$oob_rmse), ]

  list(
    stage1    = s1,
    ntree_sel = ntree_sel,
    stage2    = grid,
    best      = grid[1, ]
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
