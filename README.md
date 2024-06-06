Stratified Fine Particulate Matter Exposure Response Functions (ERFs)
=====================================================================

## [`Functions`](https://github.com/kevjosey/causal-eco/tree/main/Functions)

- [`calibrate.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/calibrate.R): Calibration function for estimating covariate balance weights with an entropy loss function.
- [`gam_ipw.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/gam_ipw.R): Functions for estimating exposure response function with generalized additive models and inverse probability weights. Includes a sandwich variance estimator for the point-wise confidence intervals of the ERF, which incorporates/propagates uncertainty from the fitted calibration weights.
- [`gam_dr.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/gam_dr.R):  A doubly-robust implementation of the GAM ERF estimate similar to the estimator in [`gam_ipw.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/gam_ipw.R) that requires an initial outcome (GAM) outcome model fit. Should be used in tandem with [`erf_fun.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/erf_fun.R).
-  [`erf_fun.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/erf_fun.R): Implementation of the doubly-robust exposure response function specific to the ecological PM2.5 regression.
-  [`srf_fun.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/srf_fun.R): Implementation of the stochastic intervention analysis (useful for obtaining excess events associated with different NAAQS policies).
-  [`simulate.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/srf_fun.R): Functions for running simulations in the [`Experiments`](https://github.com/kevjosey/causal-eco/tree/main/Experiments/) directory.

## [`Analysis`](https://github.com/kevjosey/causal-eco/tree/main/)

### Data Cleaning and Descriptives

- [`data_process.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/data_process.R): Data processing script which aggregates binary responses into strata within ZIP-code years.
- [`desriptives.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/descriptives.R): Descriptive statistics appearing in Tables 1 and 2 of the manuscript. Additional code for generating the covariate balance plots is also provided here.

### Models

- [`region_erc.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/region_erc.R) and [`region_srf`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/region_srf.R): Scripts that fit the region-specific (and whole US) exposure response functions and stochastic intervention curves using doubly-robust methods. We use entropy balancing to estimate the IPWs - think of entropy balancing as a type of method of moments estimator whereas the more traditional way of estimating IPWs is with a plug-in estimator. Uses code from [`gam_dr.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/gam_dr.R).
- [`state_erc.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/state_erc.R) and [`state_srf`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/state_srf.R): A doubly-robust implementation to estimate the state-specific exposure response curves and shift responses (stochastic interventions). Uses code from [`gam_dr.R`](https://github.com/kevjosey/causal-eco/tree/main/Functions/gam_dr.R).

## [`Experiments`](https://github.com/kevjosey/causal-eco/tree/main/Experiments/)

-  [`erf_sim.R`](https://github.com/kevjosey/causal-eco/tree/main/Experiments/erf_sim.R): Simulation study to demonstrate the effectiveness of the exposure response estimator.
-  [`srf_sim.R`](https://github.com/kevjosey/causal-eco/tree/main/Experiments/srf_sim.R): Simulation study examining accuracy of the stochastic interventions.
- [`test.R`](https://github.com/kevjosey/causal-eco/tree/main/Experiments/test.R): Test whether the functions necessary for simulation work.


## Plots [`Plots`](https://github.com/kevjosey/causal-eco/tree/main/Plots/)

- [`Plots/plot_region.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/Plots/plot_region.R): Script for generating the main plots that appear in the manuscript.
- [`Plots/plot_state.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/Plots/plot_state.R): Script for generating state-specific outcomes.
- [`Plots/plot_eco.R`](https://github.com/kevjosey/causal-eco/tree/main/Analysis/Plots/plot_eco.R): Script for plotting sensitivity analysis results of the exposure response curve. Plots include a comparison of the overall ERF using an alternative exposure assessment.
