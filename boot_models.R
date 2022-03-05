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

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/match_fun.R')
set_logger(logger_file_path = "CausalGPS.log", logger_level = "DEBUG")
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
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_rm/'

## Run Models QD

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  w <- setDF(new_data$w)
  x.tmp <- setDF(new_data$x)
  
  w.list <- split(w, list(w$zip))
  x.list <- split(x.tmp, list(x.tmp$zip))
  n.zip <- length(unique(x.tmp$zip))
  w.ord <- order(names(w.list))
  x.ord <- order(names(x.list))
  
  zip <- x.tmp$zip
  a <- x.tmp$pm25
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  rm(x.tmp); gc()
  
  if (i == 1) {
    fmla <- formula(dead ~ ns(pm25, df = 4) + factor(sex) + factor(race) + factor(dual) + factor(age_break))
  } else {
    fmla <- formula(dead ~ ns(pm25, df = 4) + factor(sex) + factor(age_break))
  }
  
  target <- match_models(a = a, w = w, x = x, zip = zip, 
                         fmla = fmla, a.vals = a.vals, trim = 0.05)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, 2*sqrt(n.zip), replace = TRUE) 
    w.boot <- data.frame(Reduce(rbind, w.list[w.ord[idx]]))
    x.tmp <- data.frame(Reduce(rbind, x.list[x.ord[idx]]))
    
    zip.boot <- x.tmp$zip
    a.boot <- x.tmp$pm25
    x.boot <- subset(x.tmp, select = -c(zip, pm25))
    
    boot_target <- match_models(a = a.boot, w = w.boot, x = x.boot, zip = zip.boot,
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
  
  w.list <- split(w, list(w$zip))
  x.list <- split(x.tmp, list(x.tmp$zip))
  n.zip <- length(unique(x.tmp$zip))
  w.ord <- order(names(w.list))
  x.ord <- order(names(x.list))
  
  zip <- x.tmp$zip
  a <- x.tmp$pm25
  x <- subset(x.tmp, select = -c(zip, pm25))
  
  rm(x.tmp); gc()
  
  if (i == 1) {
    fmla <- formula(dead ~ ns(pm25, df = 4) + factor(sex) + factor(race) + factor(dual) + factor(age_break))
  } else {
    fmla <- formula(dead ~ ns(pm25, df = 4) + factor(sex) + factor(age_break))
  }
  
  target <- match_models(a = a, w = w, x = x, zip = zip,
                         fmla = fmla, a.vals = a.vals, trim = 0.05)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, 2*sqrt(n.zip), replace = TRUE) 
    w.boot <- data.frame(Reduce(rbind, w.list[w.ord[idx]]))
    x.tmp <- data.frame(Reduce(rbind, x.list[x.ord[idx]]))
    
    zip.boot <- x.tmp$zip
    a.boot <- x.tmp$pm25
    x.boot <- subset(x.tmp, select = -c(zip, pm25))
    
    boot_target <- match_models(a = a.boot, w = w.boot, x = x.boot, zip = zip.boot,
                                attempts = 1, a.vals = a.vals, fmla = fmla, trim = 0.05)
    return(boot_target$estimate)
    
  })
  
  match_data <- target$match_data
  corr_data <- data.frame(original = target$original_corr_results, 
                          adjusted = target$adjusted_corr_results)
  boot_data <- data.frame(a.vals = a.vals, estimate = target$estimate, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(boot_data, match_data, corr_data, n.zip, file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
}
