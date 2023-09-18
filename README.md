PM2.5 Exposure Response Functions (ERFs) Stratified by Age, Race, and Medicaid Eligibility
================================================================================

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights with an entropy loss function.
- [`kwls.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/kwls.R): Functions for estimating exposure response function with kernel weighted least squares. Includes a sandwich variance estimator for the point-wise confidence intervals of the ERF, which incorporates/propagates uncertainty from the fitted calibration weights.
- [`gam.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam.R):  A GAM implementation of the ERF estimate similar to the KWLS implementation. Instead of pointwise estimates/confidence intervals, a GAM implementation yields estimates and a sandwich variance for the simultaneous confidence bands as opposed to the pointwise confidence bands found in [`kwls.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/kwls.R).

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

### Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript. Additional code for generating the covariate balance plots is also provided here.

### Models

- [`models_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_ipw.R): Script that fits the strata-specific exposure response functions using estimated inverse probability weights. We use entropy balancing to estimate the IPWs - think of entropy balancing as a type of method of moments estimator whereas the more traditional way of estimating IPWs is with a plug-in estimator. Uses code from [`kwls.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/kwls.R).
- [`models_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_dr.R): A doubly-robust estimator of the exposure response curves for strata specidic covariates. These models are used for our final analysis as they better accommodate the ecological design of our study.
- [`excess_deaths.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/excess_deaths.R) Predicts excess deaths attributable to heightened levels of PM2.5 from the fitted ERFs obtained in [`models_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_dr.R)

### Plots

- [`Plots/plot_main.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_main.R): Script for generating the main plots that appear in the manuscript.
- [`Plots/plot_age.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`Plots/plot_sens.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment.
