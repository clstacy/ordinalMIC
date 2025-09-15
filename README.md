# ordinalMIC

Minimum Inhibitory Concentration (MIC*) Estimation and groupwise comparisons from proportional-odds ordinal regression models.

> Package website: <https://clstacy.github.io/ordinalMIC/>

> Companion Website: <https://LewisLabUARK.github.io/MICalculator/>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17123352.svg)](https://doi.org/10.5281/zenodo.17123352)

## Purpose

Microbiological assays often record growth on an ordered score scale (e.g., 0–4) across concentrations. Treating these scores as binary or continuous can lead to inaccurate conclusions. `ordinalMIC` turns cumulative-link models (`ordinal::clm`) into **interpretable concentration-scale summaries**:

- **MIC\*** per group (the concentration where predicted probability of “no growth” is 50%)
- **Δ-MIC\*** (additive differences)
- **MIC\* ratio** (fold-change)
- **Difference-of-differences (DoD)** on **additive** and **ratio** scales to test interactions

It also provides compact output tables, a `ggplot2::autoplot()` for quick figures generation.


---

## Installation

### Stable (from r-universe)

```r
install.packages("ordinalMIC", repos = c("https://clstacy.r-universe.dev"))
```

### From source tarball
```r
install.packages("ordinalMIC_1.0.0.tar.gz", repos = NULL, type = "source")
```
---

## Dependencies
 - Required: `ggplot2`, `stats`, `utils`, `cli`, `tibble`, R (>= 4.1)
 - Suggested: `ordinal` for modeling

---

## Core Idea

Fit a proportional-odds ordinal regression to your scored outcomes across concentration. Then use `mic_solve()` to:

 1. Compute group-specific MIC\* on the original concentration scale (using your chosen transform/inverse).
 2. Compare groups by ΔMIC\* and MIC ratios with CIs and p-values.
 3. Test interactions via difference-of-differences (additive and ratio scales).

Default transform is `log1p` for concentration; override if your model uses a different link-side transform.

---

## Quick Start

```r
library(ordinalMIC)

# Example data included with the package
data(yeast_df)

# Fit a proportional-odds model (requires 'ordinal' package)
fit <- ordinal::clm(score ~ strain * treatment + log1p(conc), data = yeast_df)

# Estimate MIC* and compare groups
res <- mic_solve(
  clm_fit    = fit,
  conc_name  = "conc",
  transform_fun     = log1p,
  inv_transform_fun = expm1,
  alpha      = 0.05,
  compare_pairs = "all"
)

# Tables
head(tidy.mic_solve(res, table = "mic"))
head(tidy.mic_solve(res, table = "ratio"))

# Plots
ggplot2::autoplot(res, type = "mic")
```

---

## Main functions
`mic_solve()` 

 - Inputs: fitted `ordinal::clm`, newdata, conc_name, transform/inverse functions
 - Outputs: MICs, Δ-MIC, ratios, DoD summaries, all with CIs and p-values

`ggplot2::autoplot.mic_solve()`

 - Quick visualization of MICs, Δ-MICs, ratios, DoDs
 
`tidy.mic_solve()`

 - Extracts summary tables for MICs, Δ-MICs, ratios, DoDs as tibble objects
 
 ---

## Data
Example dataset `yeast_df` included with the package:


| Column     | Type    | Description                                   |
|:-----------|:--------|:----------------------------------------------|
| `strain`   | factor  | Yeast strain (WT, Mut)                        |
| `treatment`| factor  | Treatment type (None, Salt)                   |
| `conc`     | numeric | Concentration (e.g., drug concentration)      |
| `score`    | ordered | Growth score (0 = no growth, 4 = full growth) |
| `rep`      | factor  | Replicate identifier                          |


---

## Modeling guidance
 - Outcome must be an ordered factor with levels increasing in severity
 - Transform concentration (default `log1p`); supply both transform and inverse
 - `newdata` gives groups to estimate MIC*, automatically generated if not given.
 - Use DoD (Difference of Difference) summaries to interpret interaction effects on MIC*

---

## MIC definition
MIC is the concentration where the predicted probability of “no growth” equals 0.5, computed from the fitted clm.

 - MIC ratios: differences in predicted log-MIC

 - DoD ratios: ratio-of-ratios (interaction on log scale)

---

## Troubleshooting
 - **Outcome not ordered**: convert with `ordered()`
 - **Transform mismatch**: supply correct transform/inverse
 - **No Estimate for MIC**: reconsider dose range (censored results/no change observed)
 - **Comparisons**: restrict with `compare_pairs = "share_any"`
 
---

## Citation
If you use `ordinalMIC` in published research, please cite:
Stacy C (2025). ordinalMIC: Minimum Inhibitory Concentration Estimation from Ordinal Regression Models. R package version 1.0.0, https://clstacy.github.io/ordinalMIC/.

---

## Contributing
 - Issues/PRs: https://github.com/clstacy/ordinalMIC/
 - Include `sessionInfo()` and minimal example for bug reports

---

## License
MIT License. See `LICENSE` file for details.
