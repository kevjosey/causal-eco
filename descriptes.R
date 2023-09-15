library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(fst)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"),
                         race = c("all", "white", "black", "hispanic", "asian"))
# scenarios <- expand.grid(dual = c("high", "low"),
#                          race = c("white", "black", "asian", "hispanic", "other"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(2, 31, length.out = 146)

### Fit Outcome Models

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'

options(digits = 5)

for (i in nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  # load ZIP-code data
  load(paste0(dir_data, scenario$dual, "_", scenario$race, "_both_all.RData"))
  zip_data <- new_data$wx
  wts <- zip_data$n 
  
  weighted.mean(zip_data$mean_bmi, w = wts, na.rm=T); sd(zip_data$mean_bmi,na.rm=T)
  weighted.mean(zip_data$smoke_rate, w = wts, na.rm=T)*100;  sd(zip_data$smoke_rate*100,na.rm=T)
  weighted.mean(zip_data$poverty, w = wts, na.rm=T)*100; sd(zip_data$poverty*100,na.rm=T)
  weighted.mean(zip_data$education, w = wts, na.rm=T)*100; sd(zip_data$education*100,na.rm=T)
  weighted.mean(zip_data$pct_owner_occ, w = wts, na.rm=T)*100; sd(zip_data$pct_owner_occ*100,na.rm=T)
  weighted.mean(zip_data$pct_blk, w = wts, na.rm=T)*100; sd(zip_data$pct_blk*100,na.rm=T)
  weighted.mean(zip_data$hispanic, w = wts, na.rm=T)*100; sd(zip_data$hispanic*100,na.rm=T)
  weighted.mean(zip_data$medhouseholdincome, w = wts, na.rm=T)/1000; sd(zip_data$medhouseholdincome/1000,na.rm=T)
  weighted.mean(zip_data$medianhousevalue, w = wts, na.rm=T)/1000; sd(zip_data$medianhousevalue/1000,na.rm=T)
  weighted.mean(zip_data$popdensity, w = wts, na.rm=T); sd(zip_data$popdensity,na.rm=T)
  weighted.mean(zip_data$summer_tmmx, w = wts, na.rm=T)-273.15; sd(zip_data$summer_tmmx,na.rm=T)
  weighted.mean(zip_data$winter_tmmx, w = wts, na.rm=T)-273.15; sd(zip_data$winter_tmmx,na.rm=T)
  weighted.mean(zip_data$summer_rmax, w = wts, na.rm=T); sd(zip_data$summer_rmax,na.rm=T)
  weighted.mean(zip_data$winter_rmax, w = wts, na.rm=T); sd(zip_data$winter_rmax,na.rm=T)
  weighted.mean(zip_data$pm25, w = wts, na.rm=T); sd(zip_data$pm25,na.rm=T)
  weighted.mean(zip_data$regionNORTHEAST,  w = wts)*100
  weighted.mean(zip_data$regionSOUTH, w = wts)*100
  (1-weighted.mean(zip_data$regionNORTHEAST, w = wts)-weighted.mean(zip_data$regionSOUTH, w = wts)-weighted.mean(zip_data$regionWEST, w = wts))*100
  weighted.mean(zip_data$regionWEST, w = wts)*100
  
}
