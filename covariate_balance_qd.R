library("wCorr")
library("parallel")
library(fst)
library(data.table)
library("xgboost")
require(polycor)
require(dplyr)
require(tidyr)

process <- c(0:49)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
load(paste0(dir_data,"aggregate_data_qd.RData"))
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])
rm(aggregate_data_qd)

# Matching function based on Xgboost GPS
matching.fun.dose.l1.caliper_xgb <- function(simulated.data,
                                             GPS_mod,
                                             a,
                                             delta_n = 1,
                                             scale) {
  ## cosmetic changes only
  simulated.data[["treat"]] <- simulated.data[["pm25_ensemble"]]
  p.a <- dnorm(a,
               mean = predict(GPS_mod2, data.matrix(simulated.data[, c(4:19)])),
               sd = sd(simulated.data$pm25_ensemble - predict(GPS_mod2, data.matrix(simulated.data[, c(4:19)]))))
  
  ## calculate min and max once, cache result
  treat.min <- min(simulated.data[["treat"]], na.rm = TRUE)
  treat.max <- max(simulated.data[["treat"]], na.rm = TRUE)
  GPS.min <- min(simulated.data[["GPS2"]], na.rm = TRUE)
  GPS.max <- max(simulated.data[["GPS2"]], na.rm = TRUE)
  ## using transform instead of $ is mostly cosmetic
  simulated.data <- transform(simulated.data,
                              std.treat = (treat - treat.min) / (treat.max - treat.min),
                              std.GPS = (GPS2 - GPS.min) / (GPS.max - GPS.min))
  std.a <- (a - treat.min) / (treat.max - treat.min)
  std.p.a <- (p.a - GPS.min) / (GPS.max - GPS.min)
  ## this subsetting doesn't depend on i, and therefore doesn't need to be done on each iteration
  simulated.data.subset <- simulated.data[abs(simulated.data[["treat"]] - a) <= (delta_n / 2), ]
  ## doing the subtraction with `outer` is faster than looping over with sapply or parSapply
  if (nrow(simulated.data.subset) >= 1) {
    wm <- sapply(1:length(std.p.a), function(iter) {
      return(apply(abs(outer(simulated.data.subset[["std.GPS"]], std.p.a[iter], `-`)) * scale,
                   2,
                   function(x) which.min(abs(simulated.data.subset[["std.treat"]] - std.a) * (1 - scale) + x))
      )
    })
    dp <- simulated.data.subset[wm, ]
  } else {dp <- data.frame(matrix(rep(NA, ncol(simulated.data)), nrow = 1))}
  dp$pm25 <- a
  return(dp)
  gc()
}


#White female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_female_qd.RData")
covariates<-covariates_white_female_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_female/CB_", process, ".rds"))
rm(match_data_xgb, covariates, covariates_white_female_qd, GPS_mod2)
gc()


#White male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_male_qd.RData")
covariates<-covariates_white_male_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_male/CB_", process, ".rds"))
rm(match_data_xgb, covariates, covariates_white_male_qd, GPS_mod2)
gc()


#Black female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_female_qd.RData")
covariates<-covariates_black_female_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_female/CB_", process, ".rds"))
rm(match_data_xgb, covariates, covariates_black_female_qd, GPS_mod2)
gc()

#Black male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_male_qd.RData")
covariates<-covariates_black_male_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_male/CB_", process, ".rds"))

rm(covariates, match_data_xgb, GPS_mod2, covariates_black_male_qd)
gc()


#Hispanic female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_female_qd.RData")
covariates<-covariates_hispanic_female_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/hispanic_female/CB_", process, ".rds"))

rm(covariates, match_data_xgb, GPS_mod2, covariates_hispanic_female_qd)
gc()


#Hispanic male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_male_qd.RData")
covariates<-covariates_hispanic_male_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/hispanic_male/CB_", process, ".rds"))

rm(covariates, match_data_xgb, GPS_mod2, covariates_hispanic_male_qd)
gc()



#Asian female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_female_qd.RData")
covariates<-covariates_asian_female_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_female/CB_", process, ".rds"))

rm(covariates, match_data_xgb, GPS_mod2, covariates_asian_female_qd)
gc()


#Asian male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_male_qd.RData")
covariates<-covariates_asian_male_qd

# Fit GPS model using Xgboost machine
GPS_mod2 <- xgboost(data = data.matrix(covariates[, c(4:19)]),
                    label = covariates$pm25_ensemble,
                    nrounds = 50)
covariates$GPS2 <- dnorm(covariates$pm25_ensemble,
                         mean = predict(GPS_mod2, data.matrix(covariates[, c(4:19)])),
                         sd = sd(covariates$pm25_ensemble - predict(GPS_mod2, data.matrix(covariates[, c(4:19)]))))


match_data_xgb <- matching.fun.dose.l1.caliper_xgb(a.vals[process+1] + delta_n/2,
                                                   simulated.data = covariates,
                                                   GPS_mod = GPS_mod2,
                                                   delta_n = delta_n,
                                                   scale = 1)
# Saved  matched results
saveRDS(match_data_xgb, 
        file = paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_male/CB_", process, ".rds"))

rm(covariates, match_data_xgb, GPS_mod2, covariates_asian_male_qd)
gc()
