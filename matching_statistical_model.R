library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
library(fst)
require(xgboost)
require(parallel)


dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
load(paste0(dir_data,"national_merged2016_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-aggregate_data
rm(aggregate_data)
load(paste0(dir_data,"national_merged2016_qd.RData"))
#national_merged2016_qd<-national_merged2016
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

# Matching by GPS
#QD
# We utilize parallel computing to accelarate matching
process <- c(0:49)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]
mc_cores=10
aggregate_data.list<-split(aggregate_data_qd, f=aggregate_data_qd$year)
mod_sd<-mod_sd_qd
feature_names<-feature_names_qd
# Matching on single exposure level a, a.vals selects the caliper
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

# Matching on single exposure level a
matching.fun.dose.l1.caliper2 <- function(simulated.data,
                                          GPS_mod,
                                          a,
                                          delta_n=1,
                                          scale=1) {
  ## cosmetic changes only
  simulated.data[["treat"]] <- simulated.data[["pm25_ensemble"]]
  simulated.data[["year_fac"]] <- as.factor(simulated.data[["year"]])
  simulated.data[["region"]] <- as.factor(simulated.data[["region"]])
  p.a <- dnorm(a,
               mean = predict(GPS_mod,data.matrix(simulated.data[feature_names])),
               sd = mod_sd)
  
  ## calculate min and max once, cache result
  treat.min <- min(simulated.data[["treat"]], na.rm = TRUE)
  treat.max <- max(simulated.data[["treat"]], na.rm = TRUE)
  GPS.min <- min(simulated.data[["GPS"]], na.rm = TRUE)
  GPS.max <- max(simulated.data[["GPS"]], na.rm = TRUE)
  if (nrow(simulated.data) > 1) {
    simulated.data <- transform(simulated.data,
                                std.treat = (treat - treat.min) / (treat.max - treat.min),
                                std.GPS = (GPS - GPS.min) / (GPS.max - GPS.min))
    std.a <- (a - treat.min) / (treat.max - treat.min)
    std.p.a <- (p.a - GPS.min) / (GPS.max - GPS.min)
  } else {
    simulated.data <- transform(simulated.data,
                                std.treat = treat,
                                std.GPS = GPS )
    std.a <- a 
    std.p.a <- p.a 
  }
  simulated.data.subset <- simulated.data[abs(simulated.data[["treat"]] - a) <= (delta_n / 2), ]
  ## doing the subtraction with `outer` is faster than looping over with sapply or parSapply
  if (nrow(simulated.data.subset) >= 1) {
    wm <- sapply(1:length(std.p.a), function(iter) {
      return(apply(abs(outer(simulated.data.subset[["std.GPS"]], std.p.a[iter], `-`)) * scale,
                   2,
                   function(x) which.min(abs(simulated.data.subset[["std.treat"]] - std.a) * (1 - scale) + x))
      )
    })
    dp <- simulated.data.subset[wm, c("dead", "time_count")]
  } else {dp <- cbind(NA ,NA)}
  E.a <- apply(dp, 2, sum, na.rm = TRUE)
  return(c(simulated.data[1,3:7], E.a[1], E.a[2], a))
  gc()
}

# Function to implement the matching under each strata
par.match.noerrer <- function(a_i = a_i,
                              data.list,
                              GPS_mod = GPS_mod,
                              delta_n = delta_n,
                              scale = scale) {
  matching_noerror_level <- data.table(Reduce(rbind,mclapply(1:length(data.list),
                                                             function(i,
                                                                      a_i = a_i,
                                                                      GPS_mod = GPS_mod,
                                                                      delta_n = delta_n,
                                                                      scale = scale) {
                                                               return(matching.fun.dose.l1.caliper2(simulated.data = data.list[[i]],
                                                                                                    GPS_mod = GPS_mod_qd,
                                                                                                    a = a_i,
                                                                                                    delta_n = delta_n,
                                                                                                    scale = scale))
                                                             }, GPS_mod = GPS_mod, a_i = a_i, delta_n = delta_n, scale = scale, mc.cores = mc_cores)))
  colnames(matching_noerror_level) <- c("sex", "race", "dual", "entry_age_break", "followup_year", "dead", "time_count", "pm25_ensemble")
  return(matching_noerror_level)
  gc()
}

for( j in 1:length(aggregate_data.list)){
  for(i in 1:length(a.vals)){
    r<-matching.fun.dose.l1.caliper2(simulated.data = aggregate_data.list[[j]],
                                     GPS_mod = GPS_mod,
                                     a = a.vals[i],
                                     delta_n = delta_n,
                                     scale = 1)
  }
}

registerDoParallel(cores=5)
foreach(process = 0:49) %dopar%
  ( saveRDS(par.match.noerrer(a.vals[process + 1] + delta_n / 2,
                                               data.list = aggregate_data.list,
                                               GPS_mod = GPS_mod,
                                               delta_n = delta_n,
                                               scale = 1),paste0(dir_out, "Output/", process, ".rds")))

# Save the matched data
#saveRDS(match.noerrer, paste0(dir_out, "Output/", process, ".rds"))
#print(process)

# Load matched data
f <- list.files(paste0(dir_data, "Output/"), pattern = "\\.rds", full.names = TRUE)

matching <- rbindlist(lapply(f, readRDS))
matching2 <- subset(matching, time_count > 0)
matching2 <- as.data.frame(lapply(matching2, unlist))

# Fit univariate Poisson on matched data
matching_gnm2 <- summary(gnm(dead ~ pm25_ensemble+offset(log(time_count)),
                             eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                             data = matching2,
                             family = poisson(link = "log")))

matching_gnm <- summary(gnm(dead ~ pm25_ensemble+offset(log(time_count)),
                            eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data = subset(matching2, pm25_ensemble > quantile(matching2$pm25_ensemble, 0.01) &
                                            pm25_ensemble < quantile(matching2$pm25_ensemble, 0.99)),
                            family = poisson(link = "log")))
#exp(10*matching_gnm$coefficients[1])
#exp(10*matching_gnm2$coefficients[1])

save(matching2, matching_gnm, matching_gnm2, file = paste0(dir_out, "Matching.RData"))



