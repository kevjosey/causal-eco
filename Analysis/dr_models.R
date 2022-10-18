
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(SuperLearner)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/erf.R')
source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1), race = c("white","black"),
                         age_break = c("[65,75)","[75,85)","[85,95)","[95,125)"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(4, 16, length.out = 121)
n.boot <- 1000

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_mod_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_mod/'

for (i in c(13:16)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(subset(new_data$w, age_break == scenario$age_break))
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  y <- wx.tmp$dead
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count, age_break))
  
  model_data <- gam_est(a_w = a_w, y = y, w = w, log.pop = log.pop,
                        a_x = a_x, x = x, w.id = w.id, x.id = x.id,
                        a.vals = a.vals, sl.lib = c("SL.mean", "SL.glm", "SL.xgboost"))
  
  individual_data <- data.frame(wx.tmp,
                                resid.lm = model_data$resid.lm,
                                weights.lm = model_data$weights.lm_w,
                                resid.sl = model_data$resid.sl,
                                weights.sl = model_data$weights.sl_w,
                                resid.cal = model_data$resid.cal,
                                weights.cal = model_data$weights.cal_w)
  
  zip_data <- data.frame(x.tmp, weights.lm = model_data$weights.lm_x,
                         weights.sl = model_data$weights.sl_x,
                         weights.cal = model_data$weights.cal_x)
  
  save(model_data, individual_data, zip_data, 
       file = paste0(dir_mod_qd, scenario$dual, "_",
                     scenario$race, "_", scenario$age_break, "_qd.RData"))
  
}
