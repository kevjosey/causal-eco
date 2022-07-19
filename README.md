# Subpopulation-Specific Exposure Responses between Fine Particulate Matter and Mortality in the Medicare Population

Code for fitting exposure response curves stratified by medicaid status and race.

## R

- match_models.R: Function for fitting exposure response curves with GPS matching.
- tmle_glm.R: Function for fitting targeted maximum likelihood model of the exposure response curves.

## Analyses

- data_process.R: Processes individual health data into more computationally friendly counts.
- descriptives.R: Descriptive analyses of the data.
- boot_tmle.R: Code for fitting/bootstrapping TMLEs of the exposure response curve.
- boot_gam.R: Code for fitting/bootstrapping GAM outcome model.
- plot_tmle.R: Plot TMLE output.
- plot_tmle.R: Plot GAM output.


## Experimental

- boot_gam_fix.R: Experimenting with GAM fits.
- example.R: Example file.