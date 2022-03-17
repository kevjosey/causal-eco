library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(ranger)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/tmle_glm.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(dual = 2, race = "all"), scenarios)
a.vals <- seq(2, 18, length.out = 81)
n.boot <- 1000
df <- 4

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  w.tmp <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  offset <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  if (i == 1){
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  }
  
  nsa <- ns(a_x, knots = c(8, 12), intercept = TRUE)
  
  target <- tmle_glm(a_w = a_w, a_x = a_x, w = w, x = x,
                     y = y, offset = offset, df = df,
                     nsa = nsa, family = poisson(link = "log"), 
                     a.vals = a.vals, trunc = 0.01)
  
  print(paste0("Initial Fit Complete: Scenario ", i, " QD"))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, FUN = function(j, ...){
    
    # print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE)
    aa <- u.zip[idx]
    bb <- table(aa)
    wx.boot.tmp <- NULL
    x.boot.tmp <- NULL
    
    for (k in 1:max(bb)) {
      cc <- wx.tmp[wx.tmp$zip %in% names(bb[which(bb == k)]),]
      dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
      for (l in 1:k) {
        wx.boot.tmp <- rbind(wx.boot.tmp, cc)
        x.boot.tmp <- rbind(x.boot.tmp, dd)
      }
    }
    
    a_x.boot <- x.boot.tmp$pm25
    a_w.boot <- wx.boot.tmp$pm25
    y.boot <- wx.boot.tmp$dead
    offset.boot <- log(wx.boot.tmp$time_count)
    x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
    
    if (i == 1){
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
    } else {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    }
    
    boot_target <- tmle_glm(a_w = a_w.boot, a_x = a_x.boot, 
                            w = w.boot, x = x.boot, df = df,
                            y = y.boot, offset = offset.boot, 
                            nsa = nsa, family = poisson(link = "log"), 
                            a.vals = a.vals, trunc = 0.01)
    return(boot_target$estimate)
    
  })
  
  individual_data <- data.frame(wx.tmp, weights = target$weights_w)
  zip_data <- data.frame(x.tmp, weights = target$weights_x)
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i, " QD"))
  save(individual_data, zip_data, boot_data, n.zip, 
       file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
  w.tmp <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  offset <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  if (i == 1){
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  }
  
  nsa <- ns(a_x, knots = c(8, 13), intercept = TRUE)
  
  target <- tmle_glm(a_w = a_w, a_x = a_x, w = w, x = x,
                     y = y, offset = offset, df = df,
                     family = poisson(link = "log"), 
                     a.vals = a.vals, trunc = 0.01)
  
  print(paste0("Initial Fit Complete: Scenario ", i, " RM"))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, FUN = function(j, ...){
    
    # print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE)
    aa <- u.zip[idx]
    bb <- table(aa)
    wx.boot.tmp <- NULL
    x.boot.tmp <- NULL
    
    for (k in 1:max(bb)) {
      cc <- wx.tmp[wx.tmp$zip %in% names(bb[which(bb == k)]),]
      dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
      for (l in 1:k) {
        wx.boot.tmp <- rbind(wx.boot.tmp, cc)
        x.boot.tmp <- rbind(x.boot.tmp, dd)
      }
    }
    
    a_x.boot <- x.boot.tmp$pm25
    a_w.boot <- wx.boot.tmp$pm25
    y.boot <- wx.boot.tmp$dead
    offset.boot <- log(wx.boot.tmp$time_count)
    x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
    
    if (i == 1){
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
    } else {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    }
    
    boot_target <- tmle_glm(a_w = a_w.boot, a_x = a_x.boot, 
                            w = w.boot, x = x.boot, df = df,
                            y = y.boot, offset = offset.boot, 
                            nsa = nsa, family = poisson(link = "log"), 
                            a.vals = a.vals, trunc = 0.01)
    return(boot_target$estimate)
    
  })
  
  individual_data <- data.frame(wx.tmp, weights = target$weights_w)
  zip_data <- data.frame(x.tmp, weights = target$weights_x)
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i, " RM"))
  save(individual_data, zip_data, boot_data, n.zip, 
       file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
}