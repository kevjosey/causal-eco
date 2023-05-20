PM2.5 Exposure Response Curves Stratified by Age, Race, and Medicaid Eligibility
================================================================================

[TODO] short summary about the research

Below is a brief summary of the files contained within this repository: 

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights to be used in [`erf_models.R`](https://github.com/kevjosey/erc-strata/tree/main/R/erf_models.R).
- [`gam_models.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_models.R): Wrapper function for estimating generalized additive outcome model with a Poisson distribution. Output used by [`erf_models.R`](https://github.com/kevjosey/erc-strata/tree/main/R/erf_models.R).
- [`erf_models.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/erf_models.R): Functions for estimating doubly-robust exposure response function with stratified, aggregated binary (count) observations. This includes functions for performing kernel-weighted least squares, a wrapper function for estimating the nuisance parameters appearing in the pseudo-outcome, and a function that pools the pseudo-outcomes in preparation for least squares regression.

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

### Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript.

### Models

- [`dr_models.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/dr_models.R): Script that produces required strata-specific nuisance parameters and elements for [`fit_erf.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/fit_erf.R). The generalized propensity score is fit using 1) linear regression, 2) xgboost (via SuperLearner), and 3) covariate balance (calibration) weighting methods. Care should be taken in commenting out unwanted GPS implementations both in this script and in [`erf.R`](https://github.com/kevjosey/erc-strata/tree/main/R/erf.R). A generalized additive model of the outcome is fit for the outcome and is controlled by the weights, family, and df (degrees of freedom) arguments.
- [`fit_erf.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/fit_erf.R): Script for pooling pseudo-outcomes and fitting doubly-robust estimates of the ERF for a targeted strata. This function requires arguments from objects constructed by [`dr_models.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/dr_models.R).

### Plots

- [`Plots/plot_main.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_main.R): Script for generating the main plots that appear in the manuscript.
- [`Plots/plot_age.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`Plots/plot_sens.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment and the covariate balance statistics associated with the different weighting methods highlighted in [`erf.R`](https://github.com/kevjosey/erc-strata/tree/main/R/erf.R).
