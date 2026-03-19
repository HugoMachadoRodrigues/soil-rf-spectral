# Soil Spectral Modeling with Random Forest

> **A reproducible, parallelised pipeline for predicting soil properties from VisNIR spectroscopy using tuned Random Forests and multiple spectral preprocessing strategies.**

---

## Overview

This repository implements an end-to-end soil spectral modeling workflow using the [Open Soil Spectral Library (OSSL)](https://soilspectroscopy.org) and ASD FieldSpec VisNIR data (350–2500 nm). The pipeline is inspired by the work of Clingensmith & Grunwald (NRCS Soil Spectral Modeling Project) and extends it with:

- **Modular spectral preprocessing** — four strategies evaluated in parallel
- **Intelligent hyperparameter tuning** — OOB-based grid search (no data leakage)
- **Stratified k-fold cross-validation** — robust performance estimation
- **Parallelisation** throughout via `doParallel` / `foreach`
- **Comprehensive diagnostic output** — metrics, figures, and serialised models

---

## Preprocessing Strategies

| Method | Description |
|---|---|
| `sg_smooth` | Savitzky–Golay smoothing (0th derivative) — noise reduction |
| `sg_deriv1` | Savitzky–Golay 1st derivative — removes additive baseline effects |
| `snv` | Standard Normal Variate — corrects for scatter and path-length variation |
| `absorbance` | log(1/R) — pseudo-absorbance, linearises Beer–Lambert relationships |

All SG methods use configurable window size and polynomial order (see `R/config.R`).

---

## Repository Structure

```
soil-rf-spectral/
├── main_rf_ossl.R          # Main pipeline — entry point
├── R/
│   ├── config.R            # All user-configurable parameters
│   ├── utils_preprocessing.R  # Spectral preprocessing functions
│   ├── utils_metrics.R        # MSD-based performance metrics
│   ├── utils_rf.R             # Tuning, CV, and final model fitting
│   └── utils_visualization.R  # ggplot2-based figures
├── data/                   # OSSL files downloaded automatically (gitignored)
├── output/
│   ├── models/             # Serialised ranger objects (.qs)
│   ├── metrics/            # CSV metric tables
│   └── figures/            # Plots (PNG / PDF)
└── README.md
```

---

## Modeling Workflow

```
OSSL VisNIR (ASD)
       │
       ▼
  Random fraction (configurable %)
       │
       ▼
  Preprocessing ──► sg_smooth ──┐
                ──► sg_deriv1 ──┤
                ──► snv ────────┤
                ──► absorbance ─┘
                                │ (for each method)
                                ▼
                   Hyperparameter Tuning
                   (OOB grid: mtry × min.node.size)
                                │
                                ▼
                   k-fold Cross-Validation
                   (stratified, parallelised)
                                │
                                ▼
                   Final Model (full dataset)
                                │
                                ▼
              Metrics · Figures · Saved Models
```

---

## Performance Metrics

Inspired by the MSD (Mean Squared Deviation) decomposition framework:

| Metric | Description |
|---|---|
| R² | Coefficient of determination |
| RMSE | Root mean squared error |
| Bias | Systematic prediction offset |
| RPD | Residual Prediction Deviation (SD / RMSE) |
| RPIQ | Ratio of Performance to InterQuartile range |
| SB | Squared Bias (MSD component) |
| NU | Non-Unity slope (MSD component) |
| LC | Lack of Correlation (MSD component) |

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/<your-org>/soil-rf-spectral.git
cd soil-rf-spectral
```

### 2. Configure

Edit `R/config.R` to set:

```r
TARGET_PROPERTY       <- "oc_usda.c729_w.pct"   # soil property
USE_LOG_TARGET        <- TRUE                    # log-transform target?
PREPROCESSING_METHODS <- c("sg_smooth", "sg_deriv1", "snv", "absorbance")
OSSL_FRACTION         <- 0.20                    # fraction of OSSL to use
CV_FOLDS              <- 10
N_CORES               <- NULL                    # NULL = all cores minus 1
```

### 3. Run

```r
source("main_rf_ossl.R")
```

OSSL data files are downloaded automatically on first run and cached in `data/`.

---

## Dependencies

| Package | Purpose |
|---|---|
| `ranger` | Fast Random Forest implementation |
| `prospectr` | Savitzky–Golay, SNV spectral preprocessing |
| `doParallel` / `foreach` | Parallel tuning and CV |
| `qs` | Fast serialisation of R objects and OSSL files |
| `httr` | OSSL data download |
| `ggplot2` / `ggpubr` / `viridis` | Figures |
| `dplyr` / `tidyr` / `readr` | Data wrangling |
| `moments` | Skewness / kurtosis |
| `here` | Reproducible file paths |

Install all dependencies at once:

```r
install.packages("pak")
pak::pkg_install(c(
  "ranger", "prospectr", "doParallel", "foreach",
  "qs", "httr", "ggplot2", "ggpubr", "viridis",
  "dplyr", "tidyr", "readr", "moments", "here"
))
```

---

## Data Source

This pipeline uses the **Open Soil Spectral Library (OSSL)** — a global, harmonised collection of soil spectral and laboratory data.

> Safanelli, J.L., Hengl, T., Parente, L., et al. (2023). *Open Soil Spectral Library (OSSL): Building reproducible soil calibration models through open development and community engagement.* PLOS ONE. https://doi.org/10.1371/journal.pone.0296545

OSSL data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{soil_rf_spectral_2024,
  author  = {[Your Name]},
  title   = {Soil Spectral Modeling with Random Forest},
  year    = {2024},
  url     = {https://github.com/<your-org>/soil-rf-spectral}
}
```

---

## Acknowledgements

Preprocessing and metric design inspired by:
- Clingensmith, C.M. & Grunwald, S. — NRCS Soil Spectral Modeling Project
- Gauch, H.G., Hwang, J.T.G., & Fick, G.W. (2003). Model evaluation by comparison of model-based predictions and measured values. *Agronomy Journal*, 95(6), 1442–1446.

---

*Built with R · ranger · OSSL*
