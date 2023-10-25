PM2.5 Exposure Response Functions (ERFs) and Shift Response Functions (SRFs) Stratified by Age, Race, and Medicaid Eligibility
==========================================================================================

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights with an entropy loss function.
- [`gam_std.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_ipw.R): A doubly-robust implementation of a GAM ERF which incorporates an additional outcome nuisance model.

## Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript. Additional code for generating the covariate balance plots is also provided here.

## Models

- [`model_erc.R`](https://github.com/kevjosey/erc-strata/tree/main/model_erc.R): A doubly-robust implementation to estimate the strata-specific exposure response curves. Uses code from [`gam_std.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/gam_std.R). We use entropy balancing to estimate the inverse probability weights (IPWs) - think of entropy balancing as a type of method of moments estimator whereas the more traditional way of estimating IPWs is with a plug-in estimator. Uses code from [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R).
- [`model_si.R`](https://github.com/kevjosey/erc-strata/tree/main/model_si.R): Evaluation of the effects of a stochastic shift response function where we cap the total exposure to PM2.5 to 12,11,10,9,8 mug/m^3. Also included in this code is an estimation method for estimating the excess deaths avoided.

## [`Plots`](https://github.com/kevjosey/erc-strata/tree/main/Plots)

- [`Plots/plot_erc.R`](https://github.com/kevjosey/erc-strata/tree/main/Plots/plot_erc.R): Script for generating exposure response curve plots that appear in the manuscript.
- [`Plots/plot_si.R`](https://github.com/kevjosey/erc-strata/tree/main/Plots/plot_si.R): Script for generating  stochastic intervention plots and tables (excess deaths).
- [`Plots/plot_age.R`](https://github.com/kevjosey/erc-strata/tree/main/Plots/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`Plots/plot_sens.R`](https://github.com/kevjosey/erc-strata/tree/main/Plots/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment.
