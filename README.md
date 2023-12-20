PM2.5 Exposure Response Functions (ERFs) Stratified by Region, Race, and Medicaid Eligibility
================================================================================

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights with an entropy loss function.
- [`gam_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_ipw.R): Functions for estimating exposure response function with generalized additive models and inverse probability weights. Includes a sandwich variance estimator for the point-wise confidence intervals of the ERF, which incorporates/propagates uncertainty from the fitted calibration weights.
- [`gam_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_dr.R):  A doubly-robust implementation of the GAM ERF estimate similar to the estimator in [`gam_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_ipw.R) which incorporates an additional outcome nuisance model.

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

### Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript. Additional code for generating the covariate balance plots is also provided here.

### Models

- [`models_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_ipw.R): Script that fits the strata-specific exposure response functions using estimated inverse probability weights. We use entropy balancing to estimate the IPWs - think of entropy balancing as a type of method of moments estimator whereas the more traditional way of estimating IPWs is with a plug-in estimator. Uses code from [`gam_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_ipw.R).
- [`models_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_dr.R): A doubly-robust implementation to estimate the strata-specific exposure response curves. Uses code from [`gam_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_dr.R).
- [`excess_deaths.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/excess_deaths.R) Predicts excess deaths attributable to heightened levels of PM2.5 from the fitted ERFs obtained in [`models_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_dr.R)
- [`excess_deaths.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/simulations.R) Simulation study comparing the implementations featured in [`models_ipw.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_ipw.R) and [`models_dr.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models_dr.R).

### Plots

- [`Plots/plot_main.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_main.R): Script for generating the main plots that appear in the manuscript.
- [`Plots/plot_age.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`Plots/plot_sens.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment.
