# Soil Spectral Modeling with Random Forest

> **A reproducible, parallelised pipeline for predicting soil properties from VisNIR spectroscopy (ASD FieldSpec 4, 350–2500 nm) using tuned Random Forests and multiple spectral preprocessing strategies.**

---

## Project Context

This repository implements the soil spectral modeling workflow described in the **BMP Mini Proposal: Soil Spectral Modeling Assignments (Random Forest Prototype)**. The goal is to develop prototype Random Forest (RF) models linking ASD VisNIR spectra to soil laboratory measurements for Years 1–5 of the BMP project, generating report-ready figures and performance tables for the annual report.

The pipeline is inspired by and adapted from the work of **Clingensmith, C. & Grunwald, S.** (NRCS Soil Spectral Modeling Project), extended with intelligent hyperparameter tuning, repeated cross-validation, and a complete metric suite including Lin's CCC.

### Target soil properties

| Variable | Description |
|---|---|
| pH | Soil pH |
| TP | Total phosphorus |
| M3-P | Mehlich-3 extractable phosphorus |
| WEP | Water-extractable phosphorus |
| TC | Total carbon |
| TN | Total nitrogen |
| Fe | Iron |
| Al | Aluminum |
| Ca | Calcium |
| Mg | Magnesium |

> **Prototype stage:** The pipeline uses the Open Soil Spectral Library (OSSL) with ASD sensor data. Replace the data loading section with your project database when the Years 1–5 dataset is ready.

---

## Repository Structure

```
soil-rf-spectral/
├── main_rf_ossl.R               # Entry point — run this script
├── R/
│   ├── config.R                 # All tunable parameters (edit before running)
│   ├── utils_preprocessing.R    # Spectral preprocessing functions
│   ├── utils_metrics.R          # Performance metrics with formulas
│   ├── utils_rf.R               # Hyperparameter tuning and cross-validation
│   └── utils_visualization.R    # ggplot2 figures
├── data/                        # OSSL files (auto-downloaded; gitignored)
├── output/
│   ├── models/                  # Serialised ranger objects (.qs)
│   ├── metrics/                 # CSV performance tables
│   └── figures/                 # Plots (PNG / PDF)
└── README.md
```

---

## Modeling Workflow

```
OSSL VisNIR / ASD (350–2500 nm)
          │
          ▼
   Wavelength trimming + fraction sampling
          │
          ▼
   Raw reflectance
          │
          ▼
   SG Smooth (shared base) ──────────────────► Treatment 1: sg_smooth
          │
          ├──► SG 1st Derivative ────────────► Treatment 2: sg_deriv1
          │
          └──► SNV ──────────────────────────► Treatment 3: snv
                                                           │ (each treatment evaluated)
                                                           ▼
                                          Stage 1: ntree convergence
                                          (500–2000, step 100, OOB RMSE)
                                                           │
                                                           ▼
                                          Stage 2: mtry × min.node.size grid
                                          (fracs: 1/9, 1/6, 1/3, 1 of p)
                                                           │
                                                           ▼
                                   10-fold × 3-repeat stratified CV
                                   (30 resamples; parallel execution)
                                                           │
                                                           ▼
                                   Final model (full dataset)
                                                           │
                                                           ▼
                          Metrics · Figures · Saved models · CSV tables
```

---

## Spectral Preprocessing

SG smoothing is applied to raw reflectance first and serves as the common base for all three treatments. The other two treatments are derived from the smoothed spectra, not from raw reflectance.

| Treatment | Input | Transform | Removes / Linearises |
|---|---|---|---|
| `sg_smooth` | Raw reflectance | SG smooth (m = 0) | High-frequency instrument noise |
| `sg_deriv1` | SG smooth | SG 1st derivative (m = 1) | Additive baseline offsets and slope effects |
| `snv` | SG smooth | Standard Normal Variate | Multiplicative scatter and path-length variation |

SG parameters (`SG_WINDOW`, `SG_POLY`) apply uniformly across all treatments and are reported in model outputs. All preprocessing is applied identically within each cross-validation resample to prevent data leakage.

---

## Performance Metrics

All metrics are computed on the **original measurement scale** after back-transformation (when `USE_LOG_TARGET = TRUE`). For repeated CV, values are reported as **mean ± SD across all 30 resamples** (10 folds × 3 repeats).

---

### Notation

| Symbol | Meaning |
|---|---|
| n | Number of observations |
| yᵢ | Observed value for sample i |
| ŷᵢ | Predicted value for sample i |
| ȳ | Mean of observed values |
| ȳ̂ | Mean of predicted values |
| sᵧ | Sample SD of observed values |
| sŷ | Sample SD of predicted values |
| r | Pearson correlation coefficient |
| b | Slope of OLS regression ŷ = a + b·y |

---

### Primary Metrics

**ME — Mean Error (signed bias)**

$$\text{ME} = \frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)$$

Positive values indicate systematic over-prediction; negative values indicate under-prediction.

---

**RMSE — Root Mean Squared Error**

$$\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)^2}$$

Penalises large errors more heavily than MAE.

---

**MAE — Mean Absolute Error**

$$\text{MAE} = \frac{1}{n} \sum_{i=1}^{n} |\hat{y}_i - y_i|$$

Robust to outliers; same units as the response variable.

---

**R² — Coefficient of Determination**

$$R^2 = 1 - \frac{\sum_{i=1}^{n}(y_i - \hat{y}_i)^2}{\sum_{i=1}^{n}(y_i - \bar{y})^2}$$

Proportion of variance in y explained by ŷ. Computed from the regression of observed on predicted.

---

**RPD — Residual Prediction Deviation** *(Chang et al. 2001)*

$$\text{RPD} = \frac{s_y}{\text{RMSE}}$$

| RPD range | Model quality |
|---|---|
| < 1.5 | Poor — not recommended |
| 1.5 – 2.0 | Fair — rough quantitative estimates |
| 2.0 – 2.5 | Good — quantitative predictions |
| > 2.5 | Excellent |

---

**RPIQ — Ratio of Performance to InterQuartile range** *(Bellon-Maurel et al. 2010)*

$$\text{RPIQ} = \frac{Q_3 - Q_1}{\text{RMSE}}$$

Preferred over RPD for skewed or non-normal distributions (e.g., P fractions, Fe, Al).

---

**CCC — Lin's Concordance Correlation Coefficient** *(Lin 1989)*

$$\text{CCC} = \frac{2 \, s_{y\hat{y}}}{s_y^2 + s_{\hat{y}}^2 + (\bar{y} - \bar{\hat{y}})^2}$$

where $s_{y\hat{y}}$ is the covariance of observed and predicted values (n denominator).

CCC ∈ [−1, 1]. It combines **precision** (captured by r) with **accuracy** (distance from the 1:1 line). A model with high r but systematic bias will have CCC < r.

---

### MSD Decomposition *(Gauch et al. 2003)*

The Mean Squared Deviation decomposes into three additive components:

$$\text{MSE} = \underbrace{(\bar{\hat{y}} - \bar{y})^2}_{\text{SB}} \;+\; \underbrace{(1 - b)^2 \, \sigma_{\hat{y}}^2}_{\text{NU}} \;+\; \underbrace{(1 - r^2) \, \sigma_y^2}_{\text{LC}}$$

| Component | Name | Meaning |
|---|---|---|
| SB | Squared Bias | Systematic offset — model consistently over- or under-predicts |
| NU | Non-Unity slope | Scale error — predictions spread too wide or too narrow |
| LC | Lack of Correlation | Random error — irreducible noise not captured by any trend |

A well-calibrated model has SB ≈ 0, NU ≈ 0, with most error in LC.

> *σ² denotes population variance (n denominator) for additive consistency of the decomposition.*

---

## Hyperparameter Tuning

Tuning is performed entirely on the **training set using OOB error**, avoiding any information leakage from the validation folds.

### Stage 1 — ntree convergence

ntree is evaluated from 500 to 2000 in steps of 100. The smallest ntree whose OOB RMSE is within 0.1% of the global minimum is selected as the operating value.

### Stage 2 — mtry × min.node.size grid

| Parameter | Candidates |
|---|---|
| `mtry` | ⌊p/9⌋, ⌊p/6⌋, ⌊p/3⌋, p (fractions of total predictors p) |
| `min.node.size` | 3, 5, 10 |

All combinations are evaluated in parallel. The combination with the lowest OOB RMSE is used for cross-validation and the final model.

---

## Cross-Validation

| Setting | Value |
|---|---|
| Strategy | k-fold repeated CV |
| Folds | 10 |
| Repeats | 3 |
| Total resamples | 30 |
| Fold assignment | Stratified by quantile bins of the target |
| Execution | Parallel (all resamples simultaneously) |

Performance is summarised as **mean ± SD** across all 30 resamples, providing uncertainty estimates for each metric.

---

## Quick Start

### 1. Clone

```bash
git clone https://github.com/HugoMachadoRodrigues/soil-rf-spectral.git
cd soil-rf-spectral
```

### 2. Configure

Edit `R/config.R`:

```r
TARGET_PROPERTY       <- "oc_usda.c729_w.pct"              # soil property column
USE_LOG_TARGET        <- TRUE                               # log1p transform?
PREPROCESSING_METHODS <- c("raw", "sg_deriv1", "snv")      # up to 3 methods
OSSL_FRACTION         <- 0.20                               # prototype fraction
CV_FOLDS              <- 10
CV_REPEATS            <- 3
NTREE_CANDIDATES      <- seq(500, 2000, by = 100)
MTRY_FRACS            <- c(1/9, 1/6, 1/3, 1)
N_CORES               <- NULL                               # NULL = auto
```

### 3. Run

```r
source("main_rf_ossl.R")
```

OSSL data files are downloaded automatically on first run and cached in `data/`.

---

## Dependencies

```r
install.packages(c(
  "ranger", "prospectr", "doParallel", "foreach",
  "qs2", "httr", "ggplot2", "ggpubr", "viridis",
  "dplyr", "tidyr", "readr", "moments"
), repos = "https://cloud.r-project.org")
```

| Package | Purpose |
|---|---|
| `ranger` | Fast Random Forest (parallelised C++) |
| `prospectr` | Savitzky–Golay smoothing/derivatives, SNV |
| `doParallel` / `foreach` | Parallel tuning and CV loops |
| `qs2` | Read OSSL data files (`.qs` legacy format via `qread_qs1()`); models saved as `.rds` |
| `httr` | OSSL data download |
| `ggplot2` / `ggpubr` / `viridis` | Publication-quality figures |
| `dplyr` / `tidyr` / `readr` | Data wrangling |
| `moments` | Skewness / kurtosis for summary statistics |

---

## Deliverables (per proposal)

- [ ] Analysis-ready Years 1–5 modeling dataset with documented QA/QC
- [x] Reproducible R scripts for preprocessing, tuning, and CV
- [x] Performance summary table (RMSE, MAE, R², ME, RPD, RPIQ, CCC ± SD)
- [x] Figures: raw vs. processed spectra, variable importance, observed vs. predicted

---

## Data Source (prototype)

**Open Soil Spectral Library (OSSL)**

> Safanelli, J.L., Hengl, T., Parente, L., et al. (2023). Open Soil Spectral Library (OSSL): Building reproducible soil calibration models through open development and community engagement. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0296545

OSSL data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## References

- Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J.-M., & McBratney, A. (2010). Critical review of chemometric indicators commonly used for assessing the quality of the prediction of soil attributes by NIR spectroscopy. *TrAC Trends in Analytical Chemistry*, 29(9), 1073–1081.
- Chang, C.-W., Laird, D.A., Mausbach, M.J., & Hurburgh, C.R. (2001). Near-infrared reflectance spectroscopy–principal components regression analyses of soil properties. *Soil Science Society of America Journal*, 65(2), 480–490.
- Gauch, H.G., Hwang, J.T.G., & Fick, G.W. (2003). Model evaluation by comparison of model-based predictions and measured values. *Agronomy Journal*, 95(6), 1442–1446.
- Lin, L.I.-K. (1989). A concordance correlation coefficient to evaluate reproducibility. *Biometrics*, 45(1), 255–268.
- Clingensmith, C.M. & Grunwald, S. — NRCS Soil Spectral Modeling Project (internal reference).

---

## Citation

```bibtex
@software{soil_rf_spectral_2024,
  author  = {Rodrigues, Hugo M.},
  title   = {Soil Spectral Modeling with Random Forest — BMP Project},
  year    = {2024},
  url     = {https://github.com/HugoMachadoRodrigues/soil-rf-spectral}
}
```

---

*Built with R · ranger · OSSL — BMP Project, University of Florida*
