# =============================================================================
# utils_rf.R — Random Forest Tuning and Cross-Validation Utilities
# =============================================================================
# Two-stage workflow:
#   1. tune_rf()      — grid search over mtry × min.node.size using OOB error.
#                       Parallelised with doParallel / foreach.
#   2. cv_rf()        — k-fold (optionally repeated) cross-validation with
#                       the best hyperparameters. Folds run in parallel.
# =============================================================================

library(ranger)
library(foreach)
library(doParallel)

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Create stratified CV fold indices.
#'
#' @param y          Numeric response vector.
#' @param k          Number of folds.
#' @param repeats    Number of repeats.
#' @param stratified Logical; stratify by quantile bins.
#' @param seed       Random seed.
#' @return List of length (k × repeats), each element a vector of row indices
#'         for the *validation* set of that fold.
make_folds <- function(y, k = 10, repeats = 1,
                       stratified = TRUE, seed = 2024L) {
  set.seed(seed)
  n    <- length(y)
  folds <- vector("list", k * repeats)
  idx   <- 0L

  for (r in seq_len(repeats)) {
    if (stratified) {
      # Assign quantile bins then sample within bins
      bins <- cut(y, breaks = quantile(y, seq(0, 1, length.out = k + 1),
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
      idx <- idx + 1L
      folds[[idx]] <- which(fold_id == f)
    }
  }
  folds
}

# ── Hyperparameter tuning ─────────────────────────────────────────────────────

#' Grid search over mtry and min.node.size using OOB RMSE.
#'
#' @param X              Numeric matrix of predictors (n × p).
#' @param y              Numeric response vector (length n).
#' @param mtry_fracs     Numeric vector of mtry as fractions of p.
#' @param nodesize_vals  Integer vector of min.node.size candidates.
#' @param num_trees      Number of trees used during tuning.
#' @param n_cores        Number of parallel workers.
#' @param seed           Random seed.
#' @return List with elements:
#'   $grid  — full grid with OOB_RMSE and OOB_R2 for each combination,
#'   $best  — row of grid with lowest OOB_RMSE.
tune_rf <- function(X, y,
                    mtry_fracs    = c(0.05, 0.10, 0.20, 0.33),
                    nodesize_vals = c(3, 5, 10, 20),
                    num_trees     = 500L,
                    n_cores       = NULL,
                    seed          = 2024L) {

  p    <- ncol(X)
  grid <- expand.grid(
    mtry         = pmax(1L, round(mtry_fracs * p)),
    min_nodesize = nodesize_vals,
    oob_rmse     = NA_real_,
    oob_r2       = NA_real_,
    stringsAsFactors = FALSE
  )

  # Remove duplicate mtry values that round to the same integer
  grid <- grid[!duplicated(grid[, c("mtry", "min_nodesize")]), ]

  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ parallel::stopCluster(cl); doParallel::stopImplicitCluster() },
          add = TRUE)

  message(sprintf(
    "[Tuning] Grid: %d combinations | Cores: %d | Trees: %d",
    nrow(grid), n_cores, num_trees
  ))

  results <- foreach::foreach(
    i      = seq_len(nrow(grid)),
    .packages = "ranger",
    .combine  = rbind
  ) %dopar% {
    set.seed(seed + i)
    df <- as.data.frame(X)
    df[[".y"]] <- y
    fit <- ranger::ranger(
      formula      = .y ~ .,
      data         = df,
      num.trees    = num_trees,
      mtry         = grid$mtry[i],
      min.node.size = grid$min_nodesize[i],
      importance   = "none",
      verbose      = FALSE
    )
    c(oob_rmse = sqrt(fit$prediction.error),
      oob_r2   = fit$r.squared)
  }

  grid$oob_rmse <- results[, "oob_rmse"]
  grid$oob_r2   <- results[, "oob_r2"]
  grid          <- grid[order(grid$oob_rmse), ]

  list(grid = grid, best = grid[1, ])
}

# ── Cross-validation ──────────────────────────────────────────────────────────

#' k-fold cross-validation of a ranger Random Forest.
#'
#' @param X             Numeric matrix of predictors (n × p).
#' @param y             Numeric response vector (length n).
#' @param best_params   Named list with mtry and min_nodesize (from tune_rf).
#' @param num_trees     Number of trees in the final per-fold model.
#' @param folds         List of validation-set indices (from make_folds).
#' @param n_cores       Number of parallel workers.
#' @param seed          Random seed.
#' @param back_trans_fn Function to back-transform predictions (or NULL).
#' @return List with elements:
#'   $obs_pred — data frame with columns obs, pred, fold_id (original scale),
#'   $fold_metrics — list of per-fold metric vectors (calc_metrics output).
cv_rf <- function(X, y,
                  best_params,
                  num_trees     = 1000L,
                  folds,
                  n_cores       = NULL,
                  seed          = 2024L,
                  back_trans_fn = NULL) {

  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ parallel::stopCluster(cl); doParallel::stopImplicitCluster() },
          add = TRUE)

  k <- length(folds)
  message(sprintf(
    "[CV] %d folds | mtry = %d | min.node.size = %d | Trees = %d | Cores = %d",
    k, best_params$mtry, best_params$min_nodesize, num_trees, n_cores
  ))

  fold_results <- foreach::foreach(
    f         = seq_len(k),
    .packages = c("ranger"),
    .combine  = "list",
    .multicombine = TRUE
  ) %dopar% {
    val_idx   <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(X)), val_idx)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx,   , drop = FALSE]
    y_val   <- y[val_idx]

    df_train <- as.data.frame(X_train)
    df_train[[".y"]] <- y_train

    set.seed(seed + f)
    fit <- ranger::ranger(
      formula       = .y ~ .,
      data          = df_train,
      num.trees     = num_trees,
      mtry          = best_params$mtry,
      min.node.size = best_params$min_nodesize,
      importance    = "impurity",
      verbose       = FALSE
    )

    pred_val <- predict(fit, data = as.data.frame(X_val))$predictions

    list(
      obs       = y_val,
      pred      = pred_val,
      fold_id   = rep(f, length(val_idx)),
      row_idx   = val_idx,
      r2_oob    = fit$r.squared,
      rmse_oob  = sqrt(fit$prediction.error),
      importance = fit$variable.importance
    )
  }

  # Assemble predictions (model scale)
  obs_pred <- data.frame(
    row_id  = unlist(lapply(fold_results, `[[`, "row_idx")),
    fold_id = unlist(lapply(fold_results, `[[`, "fold_id")),
    obs_mod = unlist(lapply(fold_results, `[[`, "obs")),
    pred_mod = unlist(lapply(fold_results, `[[`, "pred"))
  )
  obs_pred <- obs_pred[order(obs_pred$row_id), ]

  # Back-transform if needed
  if (!is.null(back_trans_fn)) {
    obs_pred$obs  <- back_trans_fn(obs_pred$obs_mod)
    obs_pred$pred <- back_trans_fn(obs_pred$pred_mod)
  } else {
    obs_pred$obs  <- obs_pred$obs_mod
    obs_pred$pred <- obs_pred$pred_mod
  }

  # Per-fold metrics (original scale)
  fold_metrics <- lapply(fold_results, function(fr) {
    obs_bt  <- if (!is.null(back_trans_fn)) back_trans_fn(fr$obs)  else fr$obs
    pred_bt <- if (!is.null(back_trans_fn)) back_trans_fn(fr$pred) else fr$pred
    calc_metrics(obs_bt, pred_bt,
                 r2_oob   = fr$r2_oob,
                 rmse_oob = fr$rmse_oob)
  })

  # Average variable importance across folds
  imp_mat  <- do.call(rbind, lapply(fold_results, `[[`, "importance"))
  avg_imp  <- colMeans(imp_mat, na.rm = TRUE)

  list(
    obs_pred     = obs_pred,
    fold_metrics = fold_metrics,
    avg_importance = avg_imp
  )
}

# ── Final model ───────────────────────────────────────────────────────────────

#' Fit a final ranger model on the full dataset with the best hyperparameters.
#'
#' @param X           Numeric matrix of predictors.
#' @param y           Numeric response vector.
#' @param best_params Named list with mtry and min_nodesize.
#' @param num_trees   Number of trees.
#' @param n_cores     Number of parallel threads for ranger.
#' @param seed        Random seed.
#' @return A fitted ranger object.
fit_final_rf <- function(X, y, best_params,
                         num_trees = 1000L,
                         n_cores   = NULL,
                         seed      = 2024L) {
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  df <- as.data.frame(X)
  df[[".y"]] <- y
  set.seed(seed)
  ranger::ranger(
    formula        = .y ~ .,
    data           = df,
    num.trees      = num_trees,
    mtry           = best_params$mtry,
    min.node.size  = best_params$min_nodesize,
    importance     = "impurity",
    num.threads    = n_cores,
    verbose        = TRUE
  )
}
