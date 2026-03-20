# Soil Spectral Modeling with Random Forest

> **A reproducible, parallelised pipeline for predicting soil properties from VisNIR spectroscopy (ASD FieldSpec 4, 350–2500 nm) using tuned Random Forests and multiple spectral preprocessing strategies.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![ranger](https://img.shields.io/badge/ranger-Random%20Forest-2E7D32)](https://github.com/imbs-hl/ranger)
[![OSSL](https://img.shields.io/badge/Data-OSSL%20v1.2-8B4513)](https://soilspectroscopy.org)

---

## Author

**Hugo Rodrigues** — University of Florida, BMP Project

[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--8070--8126-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0000-0002-8070-8126)
[![GitHub](https://img.shields.io/badge/GitHub-HugoMachadoRodrigues-181717?logo=github&logoColor=white)](https://github.com/HugoMachadoRodrigues)
[![ResearchGate](https://img.shields.io/badge/ResearchGate-Hugo--Rodrigues-00CCBB?logo=researchgate&logoColor=white)](https://www.researchgate.net/profile/Hugo-Rodrigues-12)
[![X](https://img.shields.io/badge/X-Hugo__MRodrigues-000000?logo=x&logoColor=white)](https://x.com/Hugo_MRodrigues)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Hugo%20Rodrigues-0A66C2?logo=linkedin&logoColor=white)](https://www.linkedin.com/in/hugo-rodrigues-52b535119/)

---

## Project Context

This repository implements the soil spectral modeling workflow described in the **BMP Mini Proposal: Soil Spectral Modeling Assignments (Random Forest Prototype)**. The goal is to develop prototype Random Forest (RF) models linking ASD VisNIR spectra to soil laboratory measurements for Years 1–5 of the BMP project, generating report-ready figures and performance tables for the annual report.

The pipeline is inspired by and adapted from **Clingensmith & Grunwald (2022)**, who developed RF models linking VisNIR spectra to soil properties across the continental United States, extended here with full 3D hyperparameter tuning, grouped repeated cross-validation, and a complete metric suite including Lin's CCC and MSD decomposition.

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
│   ├── models/                  # Serialised ranger objects (.rds)
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
                                   Full 3D grid search (OOB RMSE, parallel)
                                   ntree × mtry × min.node.size
                                   16 × 4 × 3 = 192 combinations
                                                           │
                                                           ▼
                                   10-fold × 3-repeat stratified CV
                                   (30 resamples; grouped by layer UUID)
                                                           │
                                                           ▼
                                   Final model (full dataset)
                                                           │
                                                           ▼
                          Metrics · Figures · Saved models · CSV tables
```

---

## Spectral Preprocessing

Raw VisNIR reflectance spectra (350–2500 nm, ASD FieldSpec 4) were first subjected to wavelength range trimming (400–2450 nm) to remove edge-of-range noise, then processed through three candidate preprocessing treatments. All three treatments share a common Savitzky–Golay (SG) smoothed base, following the shared-base design of Rinnan et al. (2009): SG smoothing is applied once to raw reflectance, and the two additional treatments are derived from that smoothed output rather than from raw reflectance. This design ensures that observed differences in model performance are attributable to the transformation step alone, not to differences in the underlying noise reduction.

| Treatment | Input | Transform | Physical justification | Key reference |
|---|---|---|---|---|
| `sg_smooth` | Raw reflectance | SG polynomial smooth (m = 0, w = 11, p = 3) | Reduces high-frequency instrument noise while preserving spectral band shape | Savitzky & Golay (1964) |
| `sg_deriv1` | SG smooth | SG 1st derivative (m = 1) | Removes additive baseline offsets and multiplicative slope effects arising from variable path length and particle size | Rinnan et al. (2009) |
| `snv` | SG smooth | Standard Normal Variate: xSNV = (x − μ) / σ per spectrum | Corrects multiplicative scatter and path-length variation without requiring a reference spectrum; preferred over MSC for heterogeneous sample sets | Barnes et al. (1989) |

SG window size and polynomial order (`SG_WINDOW`, `SG_POLY`) are reported in all model outputs for reproducibility. All preprocessing transformations are computed exclusively within each cross-validation training partition and applied to the corresponding held-out fold, ensuring no information from validation samples influences the preprocessing step. Viscarra Rossel & Behrens (2010) demonstrated that SG derivatives and SNV consistently improve RF predictive performance for soil VisNIR spectra relative to raw reflectance across a wide range of soil properties.

---

## Performance Metrics

All metrics are computed on the **original measurement scale** after back-transformation (when `LOG_MODES = "log"`). For repeated CV, values are reported as **mean ± SD across all 30 resamples** (10 folds × 3 repeats).

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

### Full 3D grid search

All three hyperparameters are optimised jointly in a single parallel sweep — no sequential staging, no assumption that any parameter is independent of the others.

| Parameter | Candidates | Values |
|---|---|---|
| `num.trees` | 500 to 2000 in steps of 100 | 16 |
| `mtry` | ⌊p/9⌋, ⌊p/6⌋, ⌊p/3⌋, p (fractions of total predictors p) | 4 |
| `min.node.size` | 3, 5, 10 | 3 |
| **Total combinations** | **16 × 4 × 3** | **192** |

Every combination trains one `ranger` model using only OOB predictions (importance disabled for speed). All 192 models run simultaneously across all available CPU cores via `foreach %dopar%`. The combination with the **lowest OOB RMSE** is selected for cross-validation and the final model.

**Why OOB and not a separate validation set?** Each `ranger` model computes OOB predictions as a free by-product of bootstrap sampling (Breiman 2001; Liaw & Wiener 2002): trees are not evaluated on the ~37% of samples excluded from their bootstrap draw, providing an unbiased error estimate without requiring additional data splits. This is the standard tuning criterion for Random Forests (Probst et al. 2019).

**Why a full 3D grid and not two stages?** A sequential design that fixes `ntree` first and then searches `mtry × min.node.size` assumes that the optimal `ntree` is independent of `mtry` and `min.node.size`. Probst et al. (2019) show this assumption fails in practice: lower `mtry` (more per-tree diversity) typically requires more trees to converge, so the optimal `ntree` shifts with `mtry`. The full 3D grid eliminates this assumption entirely at a modest computational cost (~192 OOB models vs 28 for a two-stage approach; each model is fast with `importance = "none"`).

---

## Cross-Validation

| Setting | Value |
|---|---|
| Strategy | k-fold repeated CV |
| Folds | 10 |
| Repeats | 3 |
| Total resamples | 30 |
| Fold assignment | **Grouped by layer UUID** — all spectral replicates of the same soil layer always fall in the same fold, preventing leakage when multiple scans share a single lab measurement. Stratified by quantile bins of group-mean y. Grouped assignment follows the leakage-prevention framework of Kaufman & Rosset (2012). |
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
LOG_MODES             <- c("log", "no_log")                 # "log", "no_log", or both
PREPROCESSING_METHODS <- c("sg_smooth", "sg_deriv1", "snv") # three treatments
OSSL_FRACTION         <- 0.20                               # prototype fraction
N_MAX_SAMPLES         <- 1000                               # cap total rows
CV_FOLDS              <- 10
CV_REPEATS            <- 3
NTREE_CANDIDATES      <- seq(500, 2000, by = 100)
MTRY_FRACS            <- c(1/9, 1/6, 1/3, 1)
NODESIZE_CANDIDATES   <- c(3, 5, 10)
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
  "ranger", "prospectr", "doSNOW", "foreach",
  "httr", "ggplot2", "ggpubr", "viridis",
  "dplyr", "tidyr", "readr", "moments"
), repos = "https://cloud.r-project.org")
```

| Package | Purpose |
|---|---|
| `ranger` | Fast Random Forest (parallelised C++) |
| `prospectr` | Savitzky–Golay smoothing/derivatives, SNV |
| `doSNOW` / `foreach` | Parallel tuning and CV loops with ASCII progress bar |
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

> Safanelli et al. (2023) — see References below.

OSSL data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## References

### Spectral preprocessing
- Barnes, R.J., Dhanoa, M.S., & Lister, S.J. (1989). Standard normal variate transformation and de-trending of near-infrared diffuse reflectance spectra. *Applied Spectroscopy*, 43(5), 772–777. https://doi.org/10.1366/0003702894202201
- Rinnan, Å., van den Berg, F., & Engelsen, S.B. (2009). Review of the most common pre-processing techniques for near-infrared spectra. *TrAC Trends in Analytical Chemistry*, 28(10), 1201–1222. https://doi.org/10.1016/j.trac.2009.07.007
- Savitzky, A. & Golay, M.J.E. (1964). Smoothing and differentiation of data by simplified least squares procedures. *Analytical Chemistry*, 36(8), 1627–1639. https://doi.org/10.1021/ac60214a047

### Random Forest and hyperparameter tuning
- Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
- Liaw, A. & Wiener, M. (2002). Classification and regression by randomForest. *R News*, 2(3), 18–22. https://cran.r-project.org/doc/Rnews/
- Probst, P., Wright, M.N., & Boulesteix, A.-L. (2019). Hyperparameters and tuning strategies for random forest. *WIREs Data Mining and Knowledge Discovery*, 9(3), e1301. https://doi.org/10.1002/widm.1301

### Soil spectroscopy
- Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J.-M., & McBratney, A. (2010). Critical review of chemometric indicators commonly used for assessing the quality of the prediction of soil attributes by NIR spectroscopy. *TrAC Trends in Analytical Chemistry*, 29(9), 1073–1081. https://doi.org/10.1016/j.trac.2010.05.006
- Chang, C.-W., Laird, D.A., Mausbach, M.J., & Hurburgh, C.R. (2001). Near-infrared reflectance spectroscopy–principal components regression analyses of soil properties. *Soil Science Society of America Journal*, 65(2), 480–490. https://doi.org/10.2136/sssaj2001.652480x
- Viscarra Rossel, R.A. & Behrens, T. (2010). Using data mining to model and interpret soil diffuse reflectance spectra. *Geoderma*, 158(1–2), 46–54. https://doi.org/10.1016/j.geoderma.2009.12.025

### Performance metrics and cross-validation
- Gauch, H.G., Hwang, J.T.G., & Fick, G.W. (2003). Model evaluation by comparison of model-based predictions and measured values. *Agronomy Journal*, 95(6), 1442–1446. https://doi.org/10.2134/agronj2003.1442
- Kaufman, S. & Rosset, S. (2012). Leakage in data mining: Formulation, detection, and avoidance. *ACM Transactions on Knowledge Discovery from Data*, 6(4), 15. https://doi.org/10.1145/2382577.2382579
- Lin, L.I.-K. (1989). A concordance correlation coefficient to evaluate reproducibility. *Biometrics*, 45(1), 255–268. https://doi.org/10.2307/2532051

### OSSL data
- Safanelli, J.L., Hengl, T., Parente, L., et al. (2023). Open Soil Spectral Library (OSSL): Building reproducible soil calibration models through open development and community engagement. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0296545

### Pipeline inspiration
- Clingensmith, C.M. & Grunwald, S. (2022). Predicting soil properties and interpreting Vis-NIR models from across continental United States. *Sensors*, 22(9), 3187. https://doi.org/10.3390/s22093187

---

## Citation

```bibtex
@software{soil_rf_spectral_2026,
  author  = {Rodrigues, Hugo},
  title   = {Soil Spectral Modeling with Random Forest — BMP Project},
  year    = {2026},
  orcid   = {0000-0002-8070-8126},
  url     = {https://github.com/HugoMachadoRodrigues/soil-rf-spectral}
}
```

---

*Built with R · ranger · OSSL — BMP Project, University of Florida*
*© Hugo Rodrigues · [ORCID](https://orcid.org/0000-0002-8070-8126) · [GitHub](https://github.com/HugoMachadoRodrigues)*
