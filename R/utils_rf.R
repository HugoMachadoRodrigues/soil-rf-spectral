# =============================================================================
# utils_rf.R — Random Forest Tuning and Cross-Validation Utilities
# =============================================================================
# Implements the full modeling workflow for the BMP spectral project:
#
#   1. make_folds() — grouped stratified k-fold × repeat CV index builder.
#                     Fold assignment at the layer-UUID level prevents leakage
#                     when multiple spectral replicates share one lab value
#                     (Roberts et al. 2017).
#
#   2. tune_rf()    — full 3D grid search over ntree × mtry × min.node.size
#                     using OOB RMSE (Probst et al. 2019; Liaw & Wiener 2002).
#                     All 192 combinations run in parallel via doSNOW/foreach.
#                     No parameter-independence assumption (unlike two-stage).
#
#   3. cv_rf()      — 10-fold × 3-repeat stratified cross-validation with the
#                     tuned hyperparameters (Viscarra Rossel & Behrens 2010).
#                     Performance reported as mean ± SD across 30 resamples.
#
#   4. fit_final_rf() — final ranger model fitted on the full dataset with
#                       the best hyperparameters (Breiman 2001).
#
# Key references:
#   Breiman (2001)           — original Random Forest algorithm:
#     Breiman, L. Random forests. Machine Learning, 45(1), 5–32.
#     https://doi.org/10.1023/A:1010933404324
#
#   Liaw & Wiener (2002)     — OOB-based hyperparameter selection:
#     Liaw, A. & Wiener, M. Classification and regression by randomForest.
#     R News, 2(3), 18–22. https://cran.r-project.org/doc/Rnews/
#
#   Probst et al. (2019)     — joint 3D tuning rationale:
#     Probst, P., Wright, M.N., & Boulesteix, A.-L. Hyperparameters and
#     tuning strategies for random forest. WIREs Data Mining and Knowledge
#     Discovery, 9(3), e1301. https://doi.org/10.1002/widm.1301
#
#   Viscarra Rossel & Behrens (2010) — RF for soil VisNIR spectroscopy:
#     Viscarra Rossel, R.A. & Behrens, T. Using data mining to model and
#     interpret soil diffuse reflectance spectra. Geoderma, 158(1–2), 46–54.
#     https://doi.org/10.1016/j.geoderma.2009.12.025
#
#   Roberts et al. (2017)    — grouped / spatially-aware CV:
#     Roberts, D.R. et al. Cross-validation strategies for data with temporal,
#     spatial, hierarchical, or phylogenetic structure. Ecography, 40(8),
#     913–929. https://doi.org/10.1111/ecog.02881
#
# =============================================================================

library(ranger)
library(foreach)
library(doSNOW)     # progress-bar-aware parallel backend

# ── Internal helpers ──────────────────────────────────────────────────────────

# Format elapsed/ETA seconds as "Xm Ys" or "Xs".
fmt_time <- function(secs) {
  secs <- max(0, round(as.numeric(secs)))
  if (secs < 60L) return(sprintf("%ds", secs))
  sprintf("%dm %02ds", secs %/% 60L, secs %% 60L)
}

# Build a doSNOW progress-callback that writes an ASCII progress bar + ETA
# to stderr(), overwriting the same console line with \r.
#
# @param n_total  Total number of parallel tasks.
# @param start_t  proc.time()["elapsed"] captured just before foreach().
# @return A function(n_done) suitable for .options.snow = list(progress = ...).
.make_progress_fn <- function(n_total, start_t) {
  bar_w <- 38L
  function(n_done) {
    elapsed <- proc.time()["elapsed"] - start_t
    rate    <- if (n_done > 0L) elapsed / n_done else 0
    eta     <- if (rate > 0) rate * (n_total - n_done) else NA_real_
    pct     <- n_done / n_total
    filled  <- round(bar_w * pct)
    bar <- paste0(
      strrep("=", max(0L, filled - 1L)),
      if (filled > 0L && n_done < n_total) ">" else if (filled > 0L) "=" else "",
      strrep(" ", bar_w - filled)
    )
    cat(sprintf("\r  [%s] %d/%d (%.0f%%) | %s elapsed | ETA %s   ",
                bar, n_done, n_total, pct * 100,
                fmt_time(elapsed),
                if (is.na(eta)) "..." else fmt_time(eta)),
        file = stderr())
    flush(stderr())
    if (n_done == n_total) cat("\n", file = stderr())
  }
}

# ── CV fold construction ──────────────────────────────────────────────────────

#' Create (optionally stratified, optionally grouped) k-fold × repeat CV indices.
#'
#' When \code{groups} is supplied, fold assignment is performed at the group
#' level (e.g. layer UUID), so every row that belongs to the same group is
#' always assigned to the same fold.  This prevents leakage when multiple
#' spectral replicates share a single lab measurement — the grouped CV design
#' recommended by Roberts et al. (2017) for datasets with hierarchical or
#' replicated structure.
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
#' staging, no assumption of independence between parameters.  Probst et al.
#' (2019) show that mtry and min.node.size interact: lower mtry (more
#' diversity per tree) typically requires more trees to converge, so fixing
#' ntree before tuning mtry (two-stage design) can miss the true optimum.
#' The full 3D grid eliminates this assumption.
#'
#' OOB error is used as the tuning criterion following Liaw & Wiener (2002):
#' each ranger model computes OOB predictions as a by-product of bootstrap
#' sampling, providing an unbiased estimate without a separate validation set.
#' importance = "none" is set for all tuning models to reduce computation
#' (Viscarra Rossel & Behrens 2010 apply the same strategy).
#'
#'   ntree candidates : seq(500, 2000, by = 100)  → 16 values  (configurable)
#'   mtry fracs       : 1/9, 1/6, 1/3, 1          → 4 values   (configurable)
#'   min.node.size    : 3, 5, 10                   → 3 values   (configurable)
#'   ──────────────────────────────────────────────────────────
#'   Total combinations (default): 16 × 4 × 3 = 192
#'
#' Each combination trains one ranger model (OOB only, importance = "none")
#' and records OOB RMSE and R². All combinations run in parallel via foreach.
#' The combination with the lowest OOB RMSE is selected (Breiman 2001).
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

  # Full 3D grid: all combinations evaluated simultaneously
  grid <- expand.grid(
    ntree        = ntree_candidates,
    mtry         = mtry_vals,
    min_nodesize = nodesize_vals,
    stringsAsFactors = FALSE
  )
  n_comb <- nrow(grid)

  # ── Grid summary ───────────────────────────────────────────────────────────
  cat(sprintf(
    "\n  [Tuning] Full 3D grid search\n  %s\n  ntree  : %d values — %d to %d (step %d)\n  mtry   : %d values — %s\n  node   : %d values — %s\n  %s\n  Total  : %d × %d × %d = %d combinations | %d cores\n\n",
    strrep("─", 55),
    length(ntree_candidates),
    min(ntree_candidates), max(ntree_candidates),
    if (length(ntree_candidates) > 1L)
      as.integer(diff(range(ntree_candidates)) / (length(ntree_candidates) - 1L))
    else 0L,
    length(mtry_vals),
    paste(mtry_vals, collapse = ", "),
    length(nodesize_vals),
    paste(nodesize_vals, collapse = ", "),
    strrep("─", 55),
    length(ntree_candidates), length(mtry_vals), length(nodesize_vals),
    n_comb, n_cores
  ))

  df_all <- as.data.frame(X)
  df_all[[".y"]] <- y

  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  on.exit({ parallel::stopCluster(cl) }, add = TRUE)

  start_t  <- proc.time()["elapsed"]
  prog_fn  <- .make_progress_fn(n_comb, start_t)
  opts     <- list(progress = prog_fn)

  results <- foreach::foreach(
    i             = seq_len(n_comb),
    .packages     = "ranger",
    .combine      = rbind,
    .options.snow = opts
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

  total_t       <- proc.time()["elapsed"] - start_t
  grid$oob_rmse <- results[, "oob_rmse"]
  grid$oob_r2   <- results[, "oob_r2"]
  grid          <- grid[order(grid$oob_rmse), ]
  best          <- grid[1L, ]

  # ── Results table (top 10) ─────────────────────────────────────────────────
  cat(sprintf("  Completed in %s\n\n", fmt_time(total_t)))
  cat("  Top 10 combinations (ordered by OOB RMSE):\n")
  top10 <- head(grid, 10L)
  top10$oob_rmse <- round(top10$oob_rmse, 4)
  top10$oob_r2   <- round(top10$oob_r2,   4)
  rownames(top10) <- NULL
  print(top10[, c("ntree", "mtry", "min_nodesize", "oob_rmse", "oob_r2")],
        row.names = FALSE)
  cat(sprintf(
    "\n  \u2192 Best: ntree = %d | mtry = %d | node = %d | OOB RMSE = %.4f | OOB R\u00b2 = %.4f\n\n",
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

  cat(sprintf(
    "  [CV] %d resamples (%d-fold × %d-repeat) | ntree = %d | mtry = %d | node = %d | %d cores\n\n",
    n_resamples,
    length(unique(seq_len(n_resamples) %% (n_resamples %/% 3L))),  # approx folds
    3L,
    best_params$ntree,
    best_params$mtry,
    best_params$min_nodesize,
    n_cores
  ))

  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)
  on.exit({ parallel::stopCluster(cl) }, add = TRUE)

  start_t <- proc.time()["elapsed"]
  prog_fn <- .make_progress_fn(n_resamples, start_t)
  opts    <- list(progress = prog_fn)

  resample_results <- foreach::foreach(
    f             = seq_len(n_resamples),
    .packages     = "ranger",
    .combine      = "list",
    .multicombine = TRUE,
    .options.snow = opts
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

  cat(sprintf("  CV completed in %s\n\n", fmt_time(proc.time()["elapsed"] - start_t)))

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
