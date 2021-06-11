##################GPS matching for the entire US data.
require(doParallel)
library(data.table)
library(fst)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)

mc_cores = 20
process<-c(0:499)[as.integer(as.character(commandArgs(trailingOnly=TRUE)))+1]

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
boots_id<-process
set.seed(boots_id)


#Creating strata
#White female
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_female_qd.list<-split(white_female_qd, list(white_female_qd$zip))
num_uniq_zip <- length(unique(white_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

white_female_qd_boots<-data.frame(Reduce(rbind, white_female_qd.list[zip_sample]))
covariates_boots<-aggregate(white_female_qd_boots[,c(10:25)], by=list(white_female_qd_boots$zip,
                                                                      white_female_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names
covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

white_female_qd_boots<-left_join(white_female_qd_boots,covariates_boots,by=c("zip","year"))

white_female_qd_boots.list <- split(white_female_qd_boots,list( 
                                                               white_female_qd_boots$dual,
                                                               white_female_qd_boots$entry_age_break,
                                                               white_female_qd_boots$followup_year))
white_female_qd_boots.list<-white_female_qd_boots.list[lapply(white_female_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=white_female_qd_boots.list,GPS_mod=GPS_mod,
                   delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/white_female/Matching_boot_", boots_id, ".rds") )
rm(matching2, white_female_qd, white_female_qd_boots.list, white_female_qd_boots, white_female_qd.list, GPS_mod,
   mod_sd, covariates_boots)


#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
white_male_qd.list<-split(white_male_qd, list(white_male_qd$zip))
num_uniq_zip <- length(unique(white_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

white_male_qd_boots<-data.frame(Reduce(rbind, white_male_qd.list[zip_sample]))
covariates_boots<-aggregate(white_male_qd_boots[,c(10:25)], by=list(white_male_qd_boots$zip,
                                                                      white_male_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

white_male_qd_boots<-left_join(white_male_qd_boots,covariates_boots,by=c("zip","year"))

white_male_qd_boots.list <- split(white_male_qd_boots,list( 
  white_male_qd_boots$dual,
  white_male_qd_boots$entry_age_break,
  white_male_qd_boots$followup_year))
white_male_qd_boots.list<-white_male_qd_boots.list[lapply(white_male_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=white_male_qd_boots.list,GPS_mod=GPS_mod,
                   delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/white_male/Matching_boot_", boots_id, ".rds") )

#Black female
black_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_female_qd.list<-split(black_female_qd, list(black_female_qd$zip))
num_uniq_zip <- length(unique(blakc_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

black_female_qd_boots<-data.frame(Reduce(rbind, black_female_qd.list[zip_sample]))
covariates_boots<-aggregate(black_female_qd_boots[,c(10:25)], by=list(black_female_qd_boots$zip,
                                                                    black_female_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

black_female_qd_boots<-left_join(black_female_qd_boots,covariates_boots,by=c("zip","year"))

black_female_qd_boots.list <- split(black_female_qd_boots,list( 
  black_female_qd_boots$dual,
  black_female_qd_boots$entry_age_break,
  black_female_qd_boots$followup_year))
black_female_qd_boots.list<-black_female_qd_boots.list[lapply(black_female_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=black_female_qd_boots.list,GPS_mod=GPS_mod,
                    delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/black_female/Matching_boot_", boots_id, ".rds") )


#Black male
black_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
black_male_qd.list<-split(black_male_qd, list(black_male_qd$zip))
num_uniq_zip <- length(unique(black_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

black_male_qd_boots<-data.frame(Reduce(rbind, black_male_qd.list[zip_sample]))
covariates_boots<-aggregate(black_male_qd_boots[,c(10:25)], by=list(black_male_qd_boots$zip,
                                                                      black_male_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

black_male_qd_boots<-left_join(black_male_qd_boots,covariates_boots,by=c("zip","year"))

black_male_qd_boots.list <- split(black_male_qd_boots,list( 
  black_male_qd_boots$dual,
  black_male_qd_boots$entry_age_break,
  black_male_qd_boots$followup_year))
black_male_qd_boots.list<-black_male_qd_boots.list[lapply(black_male_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=black_male_qd_boots.list,GPS_mod=GPS_mod,
                    delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/black_male/Matching_boot_", boots_id, ".rds") )


#Hispanic female
hispanic_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_female_qd.list<-split(hispanic_female_qd, list(hispanic_female_qd$zip))
num_uniq_zip <- length(unique(hispanic_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

hispanic_female_qd_boots<-data.frame(Reduce(rbind, hispanic_female_qd.list[zip_sample]))
covariates_boots<-aggregate(hispanic_female_qd_boots[,c(10:25)], by=list(hispanic_female_qd_boots$zip,
                                                                         hispanic_female_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

hispanic_female_qd_boots<-left_join(hispanic_female_qd_boots,covariates_boots,by=c("zip","year"))

hispanic_female_qd_boots.list <- split(hispanic_female_qd_boots,list( 
  hispanic_female_qd_boots$dual,
  hispanic_female_qd_boots$entry_age_break,
  hispanic_female_qd_boots$followup_year))
hispanic_female_qd_boots.list<-hispanic_female_qd_boots.list[lapply(hispanic_female_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=hispanic_female_qd_boots.list,GPS_mod=GPS_mod,
                    delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/hispanic_female/Matching_boot_", boots_id, ".rds") )


#Hispanic male
hispanic_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)

hispanic_male_qd.list<-split(hispanic_male_qd, list(hispanic_male_qd$zip))
num_uniq_zip <- length(unique(hispanic_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

hispanic_male_qd_boots<-data.frame(Reduce(rbind, hispanic_male_qd.list[zip_sample]))
covariates_boots<-aggregate(hispanic_male_qd_boots[,c(10:25)], by=list(hispanic_male_qd_boots$zip,
                                                                         hispanic_male_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

hispanic_male_qd_boots<-left_join(hispanic_male_qd_boots,covariates_boots,by=c("zip","year"))

hispanic_male_qd_boots.list <- split(hispanic_male_qd_boots,list( 
  hispanic_male_qd_boots$dual,
  hispanic_male_qd_boots$entry_age_break,
  hispanic_male_qd_boots$followup_year))
hispanic_male_qd_boots.list<-hispanic_male_qd_boots.list[lapply(hispanic_male_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=hispanic_male_qd_boots.list,GPS_mod=GPS_mod,
                   delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/hispanic_male/Matching_boot_", boots_id, ".rds") )


#Asian female
asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_female_qd.list<-split(asian_female_qd, list(asian_female_qd$zip))
num_uniq_zip <- length(unique(asian_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

asian_female_qd_boots<-data.frame(Reduce(rbind, asian_female_qd.list[zip_sample]))
covariates_boots<-aggregate(asian_female_qd_boots[,c(10:25)], by=list(asian_female_qd_boots$zip,
                                                                      asian_female_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])

GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names

covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

asian_female_qd_boots<-left_join(asian_female_qd_boots,covariates_boots,by=c("zip","year"))

asian_female_qd_boots.list <- split(asian_female_qd_boots,list( 
  asian_female_qd_boots$dual,
  asian_female_qd_boots$entry_age_break,
  asian_female_qd_boots$followup_year))
asian_female_qd_boots.list<-asian_female_qd_boots.list[lapply(asian_female_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=asian_female_qd_boots.list,GPS_mod=GPS_mod,
                   delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/asian_female/Matching_boot_", boots_id, ".rds") )


#Asian male
asian_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
asian_male_qd.list<-split(asian_female_qd, list(asian_male_qd$zip))
num_uniq_zip <- length(unique(asian_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 

asian_male_qd_boots<-data.frame(Reduce(rbind, asian_male_qd.list[zip_sample]))
covariates_boots<-aggregate(asian_male_qd_boots[,c(10:25)], by=list(asian_male_qd_boots$zip,
                                                                      asian_male_qd_boots$year), FUN=min)
colnames(covariates_boots)[1:2]<-c("zip","year")
covariates_boots$year_fac <- as.factor(covariates_boots$year)
covariates_boots$region <- as.factor(covariates_boots$region)
covariates_boots<-subset(covariates_boots[complete.cases(covariates_boots) ,])
GPS_mod <-xgboost(data = data.matrix(covariates_boots[,c(4:19)]), 
                  label = covariates_boots$pm25_ensemble,
                  nrounds=50)
mod_sd<- sd(covariates_boots$pm25_ensemble -predict(GPS_mod,data.matrix(covariates_boots[,c(4:19)])))
feature_names <- GPS_mod$feature_names
covariates_boots$GPS<-dnorm(covariates_boots$pm25_ensemble
                            ,mean = predict(GPS_mod,data.matrix(covariates_boots[,feature_names])),
                            sd=mod_sd)
Nm<-dnorm(covariates_boots$pm25_ensemble,mean=mean(covariates_boots$pm25_ensemble,na.rm=T),
          sd=sd(covariates_boots$pm25_ensemble,na.rm=T))
covariates_boots$IPW<-Nm/(covariates_boots$GPS)
covariates_boots<-covariates_boots[,c("zip","year","IPW","GPS")]

asian_male_qd_boots<-left_join(asian_male_qd_boots,covariates_boots,by=c("zip","year"))

asian_male_qd_boots.list <- split(asian_male_qd_boots,list( 
  asian_male_qd_boots$dual,
  asian_male_qd_boots$entry_age_break,
  asian_male_qd_boots$followup_year))
asian_male_qd_boots.list<-asian_male_qd_boots.list[lapply(asian_male_qd_boots.list,nrow)>0]

cl=makeCluster(20,outfile='')
registerDoParallel(cl)
matching <- rbindlist(lapply(a.vals+delta_n/2,function(a){print(a);
  par.match.noerrer(a_i=a,data.list=asian_male_qd_boots.list,GPS_mod=GPS_mod,
                    delta_n=delta_n,scale=1)
}))
stopCluster(cl)
# Obtain matched data
matching2<-subset(matching,time_count>0)
options(stringsAsFactors = FALSE)
matching2<-as.data.frame(lapply(matching2, unlist))
print(boots_id)
Sys.time()
saveRDS(matching2, 
        paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingqd/asian_male/Matching_boot_", boots_id, ".rds") )



