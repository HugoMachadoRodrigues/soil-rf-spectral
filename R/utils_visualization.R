# =============================================================================
# utils_visualization.R — Plotting Functions
# =============================================================================

library(ggplot2)
library(ggpubr)    # stat_cor / ggarrange
library(viridis)

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_spectral <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#2C3E50", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold"),
      legend.position   = "bottom",
      plot.title        = element_text(face = "bold", hjust = 0),
      plot.subtitle     = element_text(colour = "grey40", hjust = 0)
    )
}

# ── Observed vs Predicted scatter ─────────────────────────────────────────────
#' @param obs_pred  Data frame with columns obs, pred, (optionally) fold_id.
#' @param title     Plot title.
#' @param label     Axis label for the soil property.
plot_obs_pred <- function(obs_pred, title = "Observed vs Predicted",
                          label = "Property") {
  lims <- range(c(obs_pred$obs, obs_pred$pred), na.rm = TRUE)

  p <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = fold_id), alpha = 0.55, size = 1.8) +
    scale_colour_viridis_c(name = "Fold") +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title    = title,
         subtitle = "Dashed line = 1:1; colour = CV fold",
         x        = paste("Observed", label),
         y        = paste("Predicted", label)) +
    ggpubr::stat_cor(aes(label = paste(after_stat(rr.label),
                                       after_stat(..p.label..), sep = "~`,`~")),
                     size = 3.5, colour = "#C0392B") +
    theme_spectral()
  p
}

# ── Variable importance ───────────────────────────────────────────────────────
#' @param importance Named numeric vector (from ranger).
#' @param wavelengths Numeric vector of wavelengths.
#' @param top_n      How many top variables to label.
plot_importance <- function(importance, wavelengths = NULL, top_n = 20L,
                            title = "Variable Importance") {
  if (is.null(wavelengths)) {
    wavelengths <- seq_along(importance)
  }
  df <- data.frame(
    wavelength = wavelengths,
    importance = importance / max(importance, na.rm = TRUE)
  )

  p <- ggplot(df, aes(x = wavelength, y = importance)) +
    geom_area(fill = "#2980B9", alpha = 0.4) +
    geom_line(colour = "#2980B9", linewidth = 0.6) +
    geom_point(
      data = df[order(df$importance, decreasing = TRUE)[seq_len(top_n)], ],
      colour = "#C0392B", size = 1.8
    ) +
    labs(title = title,
         x = "Wavelength (nm)",
         y = "Normalised Importance") +
    theme_spectral()
  p
}

# ── Tuning grid heatmap ───────────────────────────────────────────────────────
#' @param tune_grid  Data frame from tune_rf()$grid.
plot_tuning_grid <- function(tune_grid, title = "Tuning Grid — OOB RMSE") {
  ggplot(tune_grid,
         aes(x = factor(mtry), y = factor(min_nodesize), fill = oob_rmse)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = round(oob_rmse, 4)),
              size = 3, colour = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "plasma", name = "OOB RMSE",
                         direction = -1) +
    labs(title = title,
         x = "mtry", y = "min.node.size") +
    theme_spectral()
}

# ── CV metrics bar chart ──────────────────────────────────────────────────────
#' @param metrics_df  Data frame from metrics_summary(), one row per preprocessing.
#' @param metric      Column name to plot (e.g., "rmse", "rpd").
plot_cv_comparison <- function(metrics_df, metric = "rpd",
                               title = "CV Performance by Preprocessing") {
  ggplot(metrics_df, aes(x = reorder(preprocessing, .data[[metric]]),
                          y = .data[[metric]],
                          fill = preprocessing)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = .data[[metric]] - .data[[paste0(metric, "_sd")]],
                       ymax = .data[[metric]] + .data[[paste0(metric, "_sd")]]),
                  width = 0.25) +
    scale_fill_viridis_d(option = "turbo") +
    coord_flip() +
    labs(title = title, x = NULL, y = toupper(metric)) +
    theme_spectral() +
    theme(legend.position = "none")
}

# ── Save helper ───────────────────────────────────────────────────────────────
#' Save a ggplot to disk in one or more formats.
save_plot <- function(plot, filename, dir = "output/figures",
                      devices = "png", width = 8, height = 6, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (dev in devices) {
    fp <- file.path(dir, paste0(filename, ".", dev))
    ggsave(fp, plot = plot, device = dev,
           width = width, height = height, dpi = dpi, bg = "white")
    message("Saved: ", fp)
  }
  invisible(NULL)
}
