PM2.5 Exposure Response Curves Stratified by Age, Race, and Medicaid Eligibility
================================================================================

[TODO] short summary about the research

Below is a brief summary of the files contained within this repository: 

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/R)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/R/calibrate.R): Calibration function for estimating covariate balance weights to be used in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R).
- [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R): Functions for estimating doubly-robust exposure response function with stratified, aggregated binary (count) observations. This includes functions for performing kernel-weighted least squares, a wrapper function for estimating the nuisance parameters appearing in the pseudo-outcome, and a function that pools the pseudo-outcomes in preparation for least squares regression.
- [`match.R`](https://github.com/kevjosey//erc-strata/tree/main/R/match.R): ![#f03c15](https://via.placeholder.com/15/f03c15/f03c15.png)`DEPRECATED` - code not used in manuscript. Function for estimating exposure response function via the [`CausalGPS`](https://github.com/NSAPH-Software/CausalGPS) package with aggregated observations.
- [`tmle.R`](https://github.com/kevjosey//erc-strata/tree/main/R/tmle.R): ![#f03c15](https://via.placeholder.com/15/f03c15/f03c15.png)`DEPRECATED` - code not used in manuscript. A reworked function that adapts the methods coded in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R) using targeted maximum likelihood estimation techniques. 

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

### Data Cleaning and Descriptives

- [`desriptives.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript.
- [`data_process.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.

### Models

- [`dr_models.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/dr_models.R): Script that produces required strata-specific nuisance parameters and elements for [`fit_erf.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/fit_erf.R). The generalized propensity score is fit using 1) linear regression, 2) xgboost (via SuperLearner), and 3) covariate balance (calibration) weighting methods. Care should be taken in commenting out unwanted GPS implementations both in this script and in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R). A generalized additive model of the outcome is fit for the outcome and is controlled by the weights, family, and df (degrees of freedom) arguments.
- [`fit_erf.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/fit_erf.R): Script for pooling pseudo-outcomes and fitting doubly-robust estimates of the ERF for a targetted strata. This function requires arguments from objects constructed by [`dr_models.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/dr_models.R).

### Plots

- [`plot_main.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_dr.R): Script for generating the main plots that appear in the manuscript.
- [`plot_age.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves with the associated strata-specific sample sizes.
- [`plot_sens.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_sens.R): Script for plotting sensitivity analysis results. Plots include a comparison of the overall ERF using an alternative exposure assessment and the covariate balance statistics associated with the different weighting methods highlighted in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R).
