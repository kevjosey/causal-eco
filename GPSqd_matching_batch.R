##################GPS matching for the entire US data.
require(doParallel)
library(data.table)
library(fst)
library("parallel")
require(xgboost)

mc_cores = 20
process <- c(0:49)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]
###############aggregated level data, stratified by individual-level variables, each row is (zipcode * individual-level variables)
gc()
dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
mod_sd<-mod_sd_qd
feature_names<-feature_names_qd
GPS_mod<-GPS_mod_qd

# Matching on single exposure level a, a.vals selects the caliper
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])
rm(aggregate_data_qd, covariates_qd)

gc()

f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/FST_data_qd/",
                pattern = "\\.fst",
                full.names = TRUE)
mc_cores=4
aggregate_data.list <- mclapply(f, read_fst, mc.cores = mc_cores)

gc()
###########Matching on single exposure level a
matching.fun.dose.l1.caliper2 <- function(simulated.data,
                                          GPS_mod,
                                          a,
                                          delta_n=1,
                                          scale=1)
{
  ## cosmetic changes only
  simulated.data[["treat"]] <- simulated.data[["pm25_ensemble"]]
  simulated.data[["year_fac"]] <- as.factor(simulated.data[["year"]])
  simulated.data[["region"]] <- as.factor(simulated.data[["region"]])
  p.a <- dnorm(a,mean = predict(GPS_mod,data.matrix(simulated.data[feature_names])),
               sd=mod_sd)
  
  ## calculate min and max once, cache result
  treat.min <- min(simulated.data[["treat"]],na.rm=T)
  treat.max <- max(simulated.data[["treat"]],na.rm=T)
  GPS.min <- min(simulated.data[["GPS"]],na.rm=T)
  GPS.max <- max(simulated.data[["GPS"]],na.rm=T)

  if (nrow(simulated.data)>1){
    simulated.data <- transform(simulated.data,
                                std.treat = (treat - treat.min) / (treat.max - treat.min),
                                std.GPS = (GPS - GPS.min) / (GPS.max - GPS.min))
    std.a <- (a - treat.min) / (treat.max - treat.min)
    std.p.a <- (p.a - GPS.min) / (GPS.max - GPS.min)
  }else{
    simulated.data <- transform(simulated.data,
                                std.treat = treat ,
                                std.GPS = GPS )
    std.a <- a 
    std.p.a <- p.a 
  }
  
  simulated.data.subset <- simulated.data[abs(simulated.data[["treat"]] - a) <= (delta_n/2), ]
  ## doing the subtraction with `outer` is faster than looping over with sapply or parSapply
  if (nrow(simulated.data.subset)>=1){
    wm <- sapply(1:length(std.p.a),function(iter){
      return(apply(abs(outer(simulated.data.subset[["std.GPS"]], std.p.a[iter], `-`)) * scale,
                   2,
                   function(x) which.min(abs(simulated.data.subset[["std.treat"]] - std.a) * (1 - scale) + x))
      )
    })
    dp <- simulated.data.subset[wm, c("dead", "time_count")]
  }else{dp <-cbind(NA,NA)}
  
  E.a <- apply(dp, 2, sum, na.rm = T)
  return(c(simulated.data[1,3:7], E.a[1], E.a[2], a))
  gc()
}

########################function to implement the matching under each strata
par.match.noerrer<-function(a_i=a_i,data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=scale){
  matching_noerror_level <- data.table(Reduce(rbind,lapply(1:length(data.list),function(i,
                                                                                          a_i=a_i,
                                                                                          GPS_mod=GPS_mod,
                                                                                          delta_n=delta_n,
                                                                                          scale=scale){print("next");
    return(matching.fun.dose.l1.caliper2(simulated.data=data.list[[i]],
                                         GPS_mod=GPS_mod,
                                         a=a_i,
                                         delta_n=delta_n,
                                         scale=scale))
  },GPS_mod=GPS_mod,a_i=a_i,delta_n=delta_n,scale=scale)))
  colnames(matching_noerror_level) <-c("sex","race","dual","entry_age_break","followup_year","dead","time_count","pm25_ensemble")
  return(matching_noerror_level)
  gc()
}

delta_n<-a.vals[2]-a.vals[1]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=0.5)

saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/all/", process, ".rds"))

