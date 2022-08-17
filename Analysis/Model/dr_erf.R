
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(SuperLearner)
library(ranger)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/erf.R')
source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(3, 17, length.out = 106)

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(w.tmp$zip, w.tmp$year, sep = "-")
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  log.pop <- log(wx.tmp$time_count)

  x <- subset(x.tmp, select = -c(zip, pm25))
  
  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count))
  } else if (scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  }
  
  target <- count_erf(a_w = a_w, y = y, w = w, log.pop = log.pop,
                      a_x = a_x, x = x, w.id = w.id, x.id = x.id,
                      a.vals = a.vals, bw = 1, se.fit = TRUE,
                      sl.lib = c("SL.mean", "SL.glm", "SL.ranger"))
  
  individual_data <- data.frame(wx.tmp, weights.lm = target$weights.lm_w,
                                weights.sl = target$weights.sl_w,
                                weights.cal = target$weights.cal_w)
  zip_data <- data.frame(x.tmp, weights.lm = target$weights.lm_x,
                         weights.sl = target$weights.sl_x, 
                         weights.cal = target$weights.cal_x)
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.sl = target$estimate.sl, se.sl = sqrt(target$variance.sl),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal))
  
  print(paste0("Fit Complete: Scenario ", i, " QD"))
  save(individual_data, zip_data, est_data, n.zip,
       file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for (i in 1:nrow(scenarios)) {

  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))

  x.tmp <- setDF(new_data$x)
  w.tmp <- setDF(new_data$w)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(w.tmp$zip, w.tmp$year, sep = "-")

  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)

  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))

  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count))
  } else if (scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  }

  target <- count_erf(a_w = a_w, y = y, w = w, log.pop = log.pop,
                      a_x = a_x, x = x, w.id = w.id, x.id = x.id,
                      a.vals = a.vals, bw = 1, se.fit = TRUE,
                      sl.lib = c("SL.mean", "SL.glm", "SL.ranger"))
  
  individual_data <- data.frame(wx.tmp, weights.lm = target$weights.lm_w,
                                weights.sl = target$weights.sl_w,
                                weights.cal = target$weights.cal_w)
  zip_data <- data.frame(x.tmp, weights.lm = target$weights.lm_x,
                         weights.sl = target$weights.sl_x, 
                         weights.cal = target$weights.cal_x)
  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.sl = target$estimate.sl, se.sl = sqrt(target$variance.sl),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal))

  print(paste0("Fit Complete: Scenario ", i, " RM"))
  save(individual_data, zip_data, est_data, n.zip,
       file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))

}
