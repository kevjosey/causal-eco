PM2.5 Exposure Response Curves Stratified by Age, Race, and Medicaid Eligibility
================================================================================

[TODO] short summary about the research

Below is a brief summary of the files contained within this repository: 

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights with an entropy loss function.
- [`kwls.R`](https://github.com/kevjosey/erc-strata/tree/main/Functions/kwls.R): Functions for estimating exposure response function with kernel weighted least squares. Included is a function that incorporates the ecological study design weights wherease the other ``standard" function ignores the ecological design.

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

### Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript.

### Models

- [`models.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models.R): Script that produces inverse probability weights and then subsequently fitting the strata-specific exposure response functions using these weights. We use entropy balancing to estimate the IPWs - think of entropy balancing as a type of method of moments estimator whereas the more traditional way of estimating IPWs is a plug-in estimator.
- [`excess_deaths.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/excess_deaths.R) Predicts excess deaths attributable to heigtened levels of PM2.5 from the fitted ERFs obtained in [`models.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/models.R)

### Plots

- [`Plots/plot_main.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_main.R): Script for generating the main plots that appear in the manuscript.
- [`Plots/plot_age.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`Plots/plot_sens.R`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment and the covariate balance statistics associated with the different weighting methods highlighted in [`erf.R`](https://github.com/kevjosey/erc-strata/tree/main/R/erf.R).
