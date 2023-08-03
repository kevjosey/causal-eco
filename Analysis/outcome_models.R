
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(ggplot2)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c("high", "low"), race = c("white","black","hispanic","asian","other"),
                         sex = c("female","male"), age_break = c("[65,75)","[75,85)","[85,95)","[95,125)"), )
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$sex <- as.numeric(scenarios$sex)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(2, 31, length.out = 146)

# Load/Save models
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data/'

for (i in 1:nrow(scenarios)) {
  
  print(i)
  
  scenario <- scenarios[i,]
  load(paste0(dir_data, scenario$dual, "_", scenario$race, ".RData"))
  
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  phat.vals <- new_data$phat.vals
  
  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")
  
  # data
  y <- wx.tmp$dead
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, age_break, sex,
                                  id, dead, time_count, ipw, cal, cal_trunc))
  
  # fitted weights
  weights.lm <- wx.tmp$ipw
  weights.cal <- wx.tmp$cal
  weights.cal_trunc <- wx.tmp$cal_trunc
  
  model_data <- count_erf(a = a_w, y = y, w = w, weights = weights.cal_trunc,
                          id = w.id, a.vals = a.vals, log.pop = log.pop)
  
  individual_data <- data.frame(wx.tmp, resid = model_data$resid)
  zip_data <- data.frame(x.tmp)
  
  save(model_data, individual_data, zip_data, phat.vals, 
       file = paste0(dir_mod, scenario$dual, "_", scenario$race, "_", 
                     scenario$sex, "_", scenario$age_break, ".RData"))
  
}
