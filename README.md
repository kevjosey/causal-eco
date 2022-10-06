PM2.5 Exposure Response Curves Stratified by Race, Age, and Income
==================================================================

[TODO] short summary about the research

Below is a brief summary of the files contained within this repository: 

## [`Functions`](https://github.com/kevjosey/erc-strata/tree/main/R)

- [`calibrate.R`](https://github.com/kevjosey/erc-strata/tree/main/R/calibrate.R): Calibration function for estimating covariate balance weights to be used in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R).
- [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R): Functions for estimating doubly-robust exposure response function with stratified, aggregated binary (count) observations.
- [`match.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R): REDUNDANT - code not used in manuscript. Function for estimating exposure response function via the [`CausalGPS`](https://github.com/NSAPH-Software/CausalGPS) package with aggregated observations.
- [`tmle.R`](https://github.com/kevjosey//erc-strata/tree/main/R/tmle.R): REDUNDANT - code not used in manuscript. A reworked function that adapts the methods coded in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R) using targeted maximum likelihood estimation techniques. 

## [`Analysis`](https://github.com/kevjosey/erc-strata/tree/main/Analysis)

- [`desriptives.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript.
- [`data_process.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.

### [`Modelling Scripts`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Models)

- [`boot_gam.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/boot_gam.R): REDUNDANT - code not used in manuscript. Script for fitting G-computation estimates of the ERF with m-out-of-n (m = n/log(n)) bootstrap confidence intervals.
- [`boot_match.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/boot_match.R): REDUNDANT - code not used in manuscript. Script for fitting estimates of the ERF with a generalized propensity score matched pseudo-population and with m-out-of-n (m = 2*sqrt(n)) bootstrap confidence intervals.
- [`boot_tmle.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/boot_tmle.R): REDUNDANT - code not used in manuscript. Script for fitting targeted maximum likelihood estimates of the ERF with m-out-of-n (m = n/log(n)) bootstrap confidence intervals.
- [`fit_dr.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Model/fit_dr.R): Script for fitting doubly-robust estimates of the ERF while respecting the stratified structures in the outcome data. The generalized propensity score is fit using 1) generalized propensity scores, 2) xgboost (via SuperLearner), and 3) covariate balancing (calibration) weighting methods. Care should be taken in commenting out unwanted GPS implementations both in this script and in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R).

### [`Plotting Scripts`](https://github.com/kevjosey/erc-strata/tree/main/Analysis/Plots)

- [`plot_dr.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_dr.R): Script for generating the main plots that appear in the manuscript.
- [`plot_age.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_age.R): Script for generating the full age x race x dual Medicare/Medicaid eligibility curves.
- [`plot_sens.R`](https://github.com/kevjosey//erc-strata/tree/main/Analysis/Plot/plot_sens.R): Script for plotting sensitivity analyses, including a comparison of the overall ERF using an alternative exposure assessment and the covariate balance plots associated with the different weighting methods highlighted in [`erf.R`](https://github.com/kevjosey//erc-strata/tree/main/R/erf.R).
