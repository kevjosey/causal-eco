library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/tmle_glm.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1), race = c("white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(dual = 2, race = "all"), scenarios)
a.vals <- seq(3, 18, length.out = 76)
n.boot <- 1000

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  w.tmp <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  w.list <- split(wx.tmp, list(wx.tmp$zip))
  x.list <- split(x.tmp, list(x.tmp$zip))
  n.zip <- length(unique(x.tmp$zip))
  w.ord <- order(names(w.list))
  x.ord <- order(names(x.list))
  
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  offset <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  
  target <- tmle_glm(a_w = a_w, a_x = a_x, w = w, x = x,
                     y = y, offset = offset, df = 4,
                     family = poisson(link = "log"), 
                     a.vals = a.vals, trunc = 0.01)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- lapply(1:n.boot, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE) 
    wx.tmp <- data.frame(Reduce(rbind, w.list[w.ord[idx]]))
    x.tmp <- data.frame(Reduce(rbind, x.list[x.ord[idx]]))
    
    a_x.boot <- x.tmp$pm25
    a_w.boot <- wx.tmp$pm25
    y.boot <- wx.tmp$dead
    offset.boot <- log(wx.tmp$time_count)
    x.boot <- subset(x.tmp, select = -c(zip, pm25))
    w.boot <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
    
    boot_target <- tmle_glm(a_w = a_w.boot, a_x = a_x.boot, 
                            w = w.boot, x = x.boot, df = 4,
                            y = y.boot, offset = offset.boot, 
                            family = poisson(link = "log"), 
                            a.vals = a.vals, trunc = 0.01)
    return(boot_target$estimate)
    
  })
  
  individual_data <- data.frame(wx.tmp, weights = target$weights[-(1:nrow(x))])
  zip_data <- data.frame(x.tmp, weights = target$weights[(1:nrow(x))])
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(individual_data, zip_data, boot_data, n.zip, file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
  w.tmp <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
  w.list <- split(wx.tmp, list(wx.tmp$zip))
  x.list <- split(x.tmp, list(x.tmp$zip))
  n.zip <- length(unique(x.tmp$zip))
  w.ord <- order(names(w.list))
  x.ord <- order(names(x.list))
  
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  offset <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))
  w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  
  target <- tmle_glm(a_w = a_w, a_x = a_x, w = w, x = x,
                     y = y, offset = offset, df = 4,
                     family = poisson(link = "log"), 
                     a.vals = a.vals, trunc = 0.01)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- lapply(1:n.boot, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE) 
    wx.tmp <- data.frame(Reduce(rbind, w.list[w.ord[idx]]))
    x.tmp <- data.frame(Reduce(rbind, x.list[x.ord[idx]]))
    
    a_x.boot <- x.tmp$pm25
    a_w.boot <- wx.tmp$pm25
    y.boot <- wx.tmp$dead
    offset.boot <- log(wx.tmp$time_count)
    x.boot <- subset(x.tmp, select = -c(zip, pm25))
    w.boot <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
    
    boot_target <- tmle_glm(a_w = a_w.boot, a_x = a_x.boot, 
                            w = w.boot, x = x.boot, df = 4,
                            y = y.boot, offset = offset.boot, 
                            family = poisson(link = "log"), 
                            a.vals = a.vals, trunc = 0.01)
    return(boot_target$estimate)
    
  })
  
  individual_data <- data.frame(wx.tmp, weights = target$weights[-(1:nrow(x))])
  zip_data <- data.frame(x.tmp, weights = target$weights[(1:nrow(x))])
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(individual_data, zip_data, boot_data, n.zip, file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
}