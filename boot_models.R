library(data.table)
library(tidyr)
library(dplyr)
library(zipcode)
library(parallel)
library(doParallel)
library(CausalGPS)
library(splines)
library(ranger)
library(xgboost)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/dr_fun.R')

# scenarios
scenarios <- expand.grid(sex = c("male", "female"), race = c("white", "black", "hispanic", "asian"))
a.vals <- seq(3, 18, length.out = 76)

# Load Poisson model
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd_strata/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd_strata/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

load(paste0(dir_data_qd, "both_all_qd.RData"))
boot_data <- wrapper(data = new_data, a.vals = a.vals, n.boot = 1000)
save(boot_data, file = paste0(dir_out_qd, "both_all_qd.RData"))

lapply(1:nrow(scenarios), function(i, ...) {}) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  boot_data <- wrapper(data = new_data, a.vals = a.vals, n.boot = 1000)
  save(boot_data, file = paste0(dir_out_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

load(paste0(dir_data_rm, "both_all_rm.RData"))
boot_data <- wrapper(data = new_data, a.vals = a.vals, n.boot = 1000)
save(boot_data, file =  paste0(dir_out_rm, "both_all_rm.RData"))

lapply(1:nrow(scenarios), function(i, ...) {}) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  boot_data <- wrapper(data = new_data, a.vals = a.vals, n.boot = 1000)
  save(boot_data, file = paste0(dir_out_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  
}