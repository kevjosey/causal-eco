library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(CausalGPS)
library(splines)
library(ranger)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/match_estimate.R')
set_logger(logger_file_path = "CausalGPS.log", logger_level = "DEBUG")
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1), race = c("white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(dual = 2, race = "all"), scenarios)
a.vals <- seq(3, 18, length.out = 76)
knot.list <- list(pm25 = c(0,10))
n.boot <- 1000

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_rm/'

## Run Models QD

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  w <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  zip <- x.tmp$zip
  a <- x.tmp$pm25
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  if (i == 1) {
    fmla <- formula(dead ~ s(pm25, bs = 'cr') + factor(sex) + factor(race) + factor(dual) + factor(age_break))
  } else {
    fmla <- formula(dead ~ s(pm25, bs = 'tp') + factor(sex) + factor(age_break))
  }
  
  target <- match_estimate(a = a, w = w, x = x, zip = zip, knot.list,
                           fmla = fmla, a.vals = a.vals, trim = 0.05)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, 2*sqrt(n.zip), replace = TRUE)
    aa <- u.zip[idx]
    bb <- table(aa)
    w.boot <- NULL
    x.boot.tmp <- NULL
    
    for (k in 1:max(bb)) {
      cc <- w[w$zip %in% names(bb[which(bb == k)]),]
      dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
      for (l in 1:k) {
        w.boot <- rbind(w.boot, cc)
        x.boot.tmp <- rbind(x.boot.tmp, dd)
      }
    }
    
    zip.boot <- x.boot.tmp$zip
    a.boot <- x.boot.tmp$pm25
    x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
    
    boot_target <- match_estimate(a = a.boot, w = w.boot, x = x.boot, zip = zip.boot,
                                  attempts = 1, fmla = fmla, a.vals = a.vals, trim = 0.05)
    return(boot_target$estimate)
    
  })
  
  match_data <- target$match_data
  corr_data <- data.frame(original = target$original_corr_results, 
                          adjusted = target$adjusted_corr_results)
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(boot_data, match_data, corr_data, n.zip, file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
  w <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  
  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)
  
  zip <- x.tmp$zip
  a <- x.tmp$pm25
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  if (i == 1) {
    fmla <- formula(dead ~ s(pm25, bs = 'cr', k = 4) + factor(sex) + factor(race) + factor(dual) + factor(age_break))
  } else {
    fmla <- formula(dead ~ s(pm25, bs = 'cr', k = 4) + factor(sex) + factor(age_break))
  }
  
  target <- match_estimate(a = a, w = w, x = x, zip = zip,
                           fmla = fmla, a.vals = a.vals, trim = 0.05)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, 2*sqrt(n.zip), replace = TRUE)
    aa <- u.zip[idx]
    bb <- table(aa)
    w.boot <- NULL
    x.boot.tmp <- NULL
    
    for (k in 1:max(bb)) {
      cc <- w[w$zip %in% names(bb[which(bb == k)]),]
      dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
      for (l in 1:k) {
        w.boot <- rbind(w.boot, cc)
        x.boot.tmp <- rbind(x.boot.tmp, dd)
      }
    }
    
    zip.boot <- x.boot.tmp$zip
    a.boot <- x.boot.tmp$pm25
    x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
    
    boot_target <- match_estimate(a = a.boot, w = w.boot, x = x.boot, zip = zip.boot,
                                  attempts = 1, a.vals = a.vals, fmla = fmla, trim = 0.05)
    return(boot_target$estimate)
    
  })
  
  match_data <- target$match_data
  corr_data <- data.frame(original = target$original_corr_results$absolute_corr, 
                          adjusted = target$adjusted_corr_results$absolute_corr)
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(boot_data, match_data, corr_data, n.zip, file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
}