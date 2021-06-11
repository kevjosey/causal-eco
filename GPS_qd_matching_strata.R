##################GPS matching for the entire US data.
require(doParallel)
library(data.table)
library(fst)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)

mc_cores = 20
process <- c(0:49)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]
###############aggregated level data, stratified by individual-level variables, each row is (zipcode * individual-level variables)
gc()
dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'


load(paste0(dir_data,"aggregate_data_qd.RData"))
# Matching on single exposure level a, a.vals selects the caliper
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])
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
  #p.a <- dnorm(a,mean = predict(GPS_mod,simulated.data),sd=summary(GPS_mod)[["sigma"]])
  simulated.data[["year_fac"]] <- as.factor(simulated.data[["year"]])
  simulated.data[["region"]] <- as.factor(simulated.data[["region"]])
  p.a <- dnorm(a,mean = predict(GPS_mod,data.matrix(simulated.data[feature_names])),
               sd=mod_sd)
  
  ## calculate min and max once, cache result
  treat.min <- min(simulated.data[["treat"]],na.rm=T)
  treat.max <- max(simulated.data[["treat"]],na.rm=T)
  GPS.min <- min(simulated.data[["GPS"]],na.rm=T)
  GPS.max <- max(simulated.data[["GPS"]],na.rm=T)
  ## using transform instead of $ is mostly cosmetic
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
  ## this subsetting doesn't depend on i, and therefore doesn't need to be done on each iteration
  simulated.data.subset <- simulated.data[abs(simulated.data[["treat"]] - a) <= (delta_n/2), ]
  ## doing the subtraction with `outer` is faster than looping over with sapply or parSapply
  if (nrow(simulated.data.subset)>=1){
    wm <- sapply(1:length(std.p.a),function(iter){
      return(apply(abs(outer(simulated.data.subset[["std.GPS"]], std.p.a[iter], `-`)) * scale,
                   2,
                   function(x) which.min(abs(simulated.data.subset[["std.treat"]] - std.a) * (1 - scale) + x))
      )
    })
    #dp<-simulated.data.subset[wm,]
    dp <- simulated.data.subset[wm, c("dead", "time_count")]
  }else{dp <-cbind(NA,NA)}
  ##else{dp<-data.frame(matrix(rep(NA, ncol(simulated.data)), nrow=1))}
  
  E.a <- apply(dp, 2, sum, na.rm = T)
  return(c(simulated.data[1,3:7], E.a[1], E.a[2], a))
  ##return(c(simulated.data, E.a[1], E.a[2], a))
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
  colnames(matching_noerror_level) <-c("sex", "race", "dual","entry_age_break","followup_year","dead","time_count","pm25_ensemble")
  return(matching_noerror_level)
  gc()
}

delta_n<-a.vals[2]-a.vals[1]

#Creating strata
#White female
#for(proces in 23:49){
white_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
load(paste0(dir_data,"balance_qd/covariates_white_female_qd.RData"))
covariates_qd<-covariates_white_female_qd
rm(covariates_white_female_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                     label = covariates_qd$pm25_ensemble,
                     nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
white_female_qd<-merge(white_female_qd,covariates_qd,
                                          by=c("zip","year"),all.x=T)
aggregate_data.list <- split(white_female_qd,
                                list(white_female_qd$dual,
                                     white_female_qd$entry_age_break,
                                     white_female_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/white_female/", process, ".rds"))
rm(match.noerrer, covariates_qd, white_female_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()
#}
#White Male
white_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
load(paste0(dir_data,"balance_qd/covariates_white_male_qd.RData"))
covariates_qd<-covariates_white_male_qd
rm(covariates_white_male_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
white_male_qd<-merge(white_male_qd,covariates_qd,
                       by=c("zip","year"),all.x=T)
aggregate_data.list <- split(white_male_qd,
                             list(white_male_qd$dual,
                                  white_male_qd$entry_age_break,
                                  white_male_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/white_male/", process, ".rds"))
rm(match.noerrer, covariates_qd, white_male_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()


#Black female
black_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
load(paste0(dir_data,"balance_qd/covariates_black_female_qd.RData"))
covariates_qd<-covariates_black_female_qd
rm(covariates_black_female_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
black_female_qd<-merge(black_female_qd,covariates_qd,
                     by=c("zip","year"),all.x=T)
aggregate_data.list <- split(black_female_qd,
                             list(black_female_qd$dual,
                                  black_female_qd$entry_age_break,
                                  black_female_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/black_female/", process, ".rds"))
rm(match.noerrer, covariates_qd, black_female_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()



#Black male
black_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
load(paste0(dir_data,"balance_qd/covariates_black_male_qd.RData"))
covariates_qd<-covariates_black_male_qd
rm(covariates_black_male_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
black_male_qd<-merge(black_male_qd,covariates_qd,
                       by=c("zip","year"),all.x=T)
aggregate_data.list <- split(black_male_qd,
                             list(black_male_qd$dual,
                                  black_male_qd$entry_age_break,
                                  black_male_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/black_male/", process, ".rds"))
rm(match.noerrer, covariates_qd, black_male_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()


#Hispanic female
hispanic_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
load(paste0(dir_data,"balance_qd/covariates_hispanic_female_qd.RData"))
covariates_qd<-covariates_hispanic_female_qd
rm(covariates_hispanic_female_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
hispanic_female_qd<-merge(hispanic_female_qd,covariates_qd,
                     by=c("zip","year"),all.x=T)
aggregate_data.list <- split(hispanic_female_qd,
                             list(hispanic_female_qd$dual,
                                  hispanic_female_qd$entry_age_break,
                                  hispanic_female_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/hispanic_female/", process, ".rds"))
rm(match.noerrer, covariates_qd, hispanic_female_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()


#Hispanic male
hispanic_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
load(paste0(dir_data,"balance_qd/covariates_hispanic_male_qd.RData"))
covariates_qd<-covariates_hispanic_male_qd
rm(covariates_hispanic_male_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
hispanic_male_qd<-merge(hispanic_male_qd,covariates_qd,
                          by=c("zip","year"),all.x=T)
aggregate_data.list <- split(hispanic_male_qd,
                             list(hispanic_male_qd$dual,
                                  hispanic_male_qd$entry_age_break,
                                  hispanic_male_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/hispanic_male/", process, ".rds"))
rm(match.noerrer, covariates_qd, hispanic_male_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()


#Asian female
asian_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
load(paste0(dir_data,"balance_qd/covariates_asian_female_qd.RData"))
covariates_qd<-covariates_asian_female_qd
rm(covariates_asian_female_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
asian_female_qd<-merge(asian_female_qd,covariates_qd,
                        by=c("zip","year"),all.x=T)
aggregate_data.list <- split(asian_female_qd,
                             list(asian_female_qd$dual,
                                  asian_female_qd$entry_age_break,
                                  asian_female_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/asian_female/", process, ".rds"))
rm(match.noerrer, covariates_qd, asian_female_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()


#Asian male
asian_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
load(paste0(dir_data,"balance_qd/covariates_asian_male_qd.RData"))
covariates_qd<-covariates_asian_male_qd
rm(covariates_asian_male_qd)
GPS_mod <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                  label = covariates_qd$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)])),
                         sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod,data.matrix(covariates_qd[,c(4:19)]))))
Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
asian_male_qd<-merge(asian_male_qd,covariates_qd,
                        by=c("zip","year"),all.x=T)
aggregate_data.list <- split(asian_male_qd,
                             list(asian_male_qd$dual,
                                  asian_male_qd$entry_age_break,
                                  asian_male_qd$followup_year))
aggregate_data.list<-aggregate_data.list[lapply(aggregate_data.list,nrow)>0]

match.noerrer<-par.match.noerrer(a.vals[process+1]+delta_n/2,
                                 data.list=aggregate_data.list,GPS_mod=GPS_mod,delta_n=delta_n,scale=1)
saveRDS(match.noerrer, paste0(dir_out, "Matching_Output_qd/asian_male/", process, ".rds"))
rm(match.noerrer, covariates_qd, asian_male_qd, GPS_mod, mod_sd, aggregate_data.list)
gc()




