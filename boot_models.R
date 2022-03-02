.libPaths("~/lib/R")

library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(CausalGPS)
library(splines)
library(pscl)
library(ranger)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/dr_fun.R')

# scenarios
scenarios <- expand.grid(sex = c("male", "female"), race = c("white", "black", "hispanic", "asian"))
scenarios$sex <- as.character(scenarios$sex)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(sex = "both", race = "all"), scenarios)
a.vals <- seq(3, 18, length.out = 151)
n.boot <- 1000

# Load Poisson model
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD
for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  new_data.list <- split(new_data, list(new_data$zip))
  n.zip <- length(unique(new_data$zip))
  
  w <- setDF(subset(new_data, select = -c(zip, dead, time_count, pm25)))
  x <- setDF(subset(new_data, select = c(2,6:22)))
  a <- new_data$pm25
  y <- new_data$dead
  offset <- log(new_data$time_count)
  
  target <- tmle_glm(a = a, w = w, x = x, y = y, offset = offset, a.vals = a.vals)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE) 
    boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    
    w.boot <- setDF(subset(boot_data, select = -c(zip, dead, time_count, pm25)))
    x.boot <- setDF(boot_data[,c(2, 6:22)])
    a.boot <- boot_data$pm25
    y.boot <- boot_data$dead
    offset.boot <- log(boot_data$time_count)
    
    boot_target <- tmle_glm(a = a.boot, w = w.boot, x = x.boot,
                            y = y.boot, offset = offset.boot, a.vals = a.vals)
    return(boot_target$estimate)
    
  })
  
  estimate <- target$estimate
  out_data <- data.frame(y = y, a = a, x, offset = offset, weights = target$weights)
  boot_out <- data.frame(a.vals = a.vals, estimate = estimate, Reduce(cbind, boot_list))
  colnames(boot_out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(boot_out, out_data, n.zip, file = paste0(dir_out_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  new_data.list <- split(new_data, list(new_data$zip))
  n.zip <- length(unique(new_data$zip))
  
  w <- setDF(subset(new_data, select = -c(zip, dead, time_count, pm25)))
  x <- setDF(subset(new_data, select = c(2, 6:22)))
  a <- new_data$pm25
  y <- new_data$dead
  offset <- log(new_data$time_count)
  
  target <- tmle_glm(a = a, w = w, x = x, y = y, offset = offset, a.vals = a.vals)
  
  print(paste0("Initial Fit Complete: Scenario ", i))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 1, function(j, ...){
    
    print(j)
    
    idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE) 
    boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    
    w.boot <- setDF(subset(boot_data, select = -c(zip, dead, time_count, pm25)))
    x.boot <- setDF(boot_data[,c(2, 6:22)])
    a.boot <- boot_data$pm25
    y.boot <- boot_data$dead
    offset.boot <- log(boot_data$time_count)
    
    boot_target <- tmle_glm(a = a.boot, w = w.boot, x = x.boot,
                            y = y.boot, offset = offset.boot, a.vals = a.vals)
    return(boot_target$estimate)
    
  })
  
  estimate <- target$estimate
  out_data <- data.frame(y = y, a = a, x, offset = offset, weights = target$weights)
  boot_out <- data.frame(a.vals = a.vals, estimate = estimate, Reduce(cbind, boot_list))
  colnames(boot_out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  print(paste0("Bootstrap Complete: Scenario ", i))
  save(boot_out, out_data, n.zip, file = paste0(dir_out_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  
}
