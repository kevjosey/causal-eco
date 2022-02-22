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
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/dr_fun.R')

# scenarios
scenarios <- expand.grid(sex = c("male", "female"), race = c("white", "black", "hispanic", "asian"))
scenarios$sex <- as.character(scenarios$sex)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(sex = "both", race = "all"), scenarios)
a.vals <- seq(3, 18, length.out = 76)
n.boot <- 500

# Load Poisson model
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd_strata/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm_strata/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  new_data.list <- split(new_data, list(new_data$zip))
  n_zip <- length(unique(new_data$zip))
  
  x <- setDF(subset(new_data, select = -c(zip, dead, time_count, pm25)))
  a <- new_data$pm25
  y <- new_data$dead
  offset <- log(new_data$time_count)
  
  target <- tmle_glm(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
  # target_ranger <- tmle_ranger(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
  # plot(a.vals, target$estimate, ylab = "Risk", xlab = "PM2.5", main = "All QD", 
  #      type = "l", ylim = c(0.045, 0.05), col = "blue")
  # lines(a.vals, target_ranger$estimate, col = "red")
  # legend(x = 3, y = 0.05, lty = 1, legend = c("GLM GPS", "RF GPS"), col = c("blue", "red"))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 24, function(j, new_data, new_data.list, a.vals, ...) {

    idx <- sample(1:n_zip, n_zip, replace = TRUE) 
    boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    
    x <- setDF(subset(boot_data, select = -c(zip, dead, time_count, pm25)))
    a <- boot_data$pm25
    y <- boot_data$dead
    offset <- log(boot_data$time_count)
    
    boot_target <- tmle_glm(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
    return(boot_target$estimate)
    
  }, new_data = new_data, new_data.list = new_data.list, a.vals = a.vals)
  
  estimate <- target$estimate
  out_data <- data.frame(y = y, a = a, x, offset = offset, weights = target$weights)
  boot_out <- data.frame(a.vals = a.vals, estimate = estimate, Reduce(cbind, boot_list))
  colnames(boot_out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  save(boot_out, out_data, file = paste0(dir_out_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  
}

## Run Models RM

for(i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  new_data.list <- split(new_data, list(new_data$zip))
  n_zip <- length(unique(new_data$zip))
  
  x <- setDF(subset(new_data, select = -c(zip, dead, time_count, pm25)))
  a <- new_data$pm25
  y <- new_data$dead
  offset <- log(new_data$time_count)
  
  target <- tmle_glm(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
  # target_ranger <- tmle_ranger(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
  # plot(a.vals, target$estimate, ylab = "Risk", xlab = "PM2.5", main = "All RM", 
  #      type = "l", ylim = c(0.045, 0.05), col = "blue")
  # lines(a.vals, target_ranger$estimate, col = "red")
  # legend(x = 3, y = 0.05, lty = 1, legend = c("GLM GPS", "RF GPS"), col = c("blue", "red"))
  
  boot_list <- mclapply(1:n.boot, mc.cores = 24, function(j, new_data, new_data.list, a.vals, ...) {
    
    idx <- sample(1:n_zip, n_zip, replace = TRUE) 
    boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    
    x <- setDF(subset(boot_data, select = -c(zip, dead, time_count, pm25)))
    a <- boot_data$pm25
    y <- boot_data$dead
    offset <- log(boot_data$time_count)
    
    boot_target <- tmle_glm(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
    return(boot_target$estimate)
    
  }, new_data = new_data, new_data.list = new_data.list, a.vals = a.vals)
  
  estimate <- target$estimate
  out_data <- data.frame(y = y, a = a, x, offset = offset, weights = target$weights)
  boot_out <- data.frame(a.vals = a.vals, estimate = estimate, Reduce(cbind, boot_list))
  colnames(boot_out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  save(boot_out, out_data, file = paste0(dir_out_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  
}
