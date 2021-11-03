# Bootstrap on zip-code cluster to obtain robust CIs account for spatial correlation
library("mgcv")
library("parallel")
library("dplyr")
require(parallel)
library(data.table)
require(tidyr)
require(doParallel)
library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
library(fst)
library("parallel")
require(ggplot2)
require(cowplot)
require(ggExtra)
require(tidyr)
require(tidyverse)

set.seed(1234)
#Load Poisson model
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_qd/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,by=c("zip","year"),all.x=T)
rm(covariates_qd, GPS_mod_qd)
aggregate_data_qd <- subset(aggregate_data_qd,
                            pm25_ensemble < quantile(aggregate_data_qd$pm25_ensemble,0.95)&
                              pm25_ensemble > quantile(aggregate_data_qd$pm25_ensemble,0.05))

# rm
#all
cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(aggregate_data_qd, list(aggregate_data_qd$zip))
num_uniq_zip <- length(unique(aggregate_data_qd$zip))

loglinear_coefs_boots<-NULL
for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year),
  data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots<-c(loglinear_coefs_boots,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(loglinear_coefs_boots,file=paste0(dir_out,"loglinear_coefs_boots_qd_all_trimmed.RData"))
rm(gnm_raw, aggregate_data_qd, aggregate_data.list, cl)
#Race
#White female
require(dplyr)
load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,by=c("zip","year"),all.x=T)
rm(covariates_qd, GPS_mod_qd)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_female_qd <- subset(white_female_qd,
                          pm25_ensemble < quantile(white_female_qd$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(white_female_qd$pm25_ensemble,0.05))

cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_female_qd, list(white_female_qd$zip))
num_uniq_zip <- length(unique(white_female_qd$zip))
loglinear_coefs_boots_white_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                       mean_bmi + smoke_rate + hispanic + pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count))+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)
               , data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_female<-c(loglinear_coefs_boots_white_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_white_female_trimmed.RData"))
rm(white_female_qd)
gc()

#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
white_male_qd <- subset(white_male_qd,
                        pm25_ensemble < quantile(white_male_qd$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(white_male_qd$pm25_ensemble,0.05))
cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_male_qd, list(white_male_qd$zip))
num_uniq_zip <- length(unique(white_male_qd$zip))
loglinear_coefs_boots_white_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                              mean_bmi + smoke_rate + hispanic + pct_blk +
                              medhouseholdincome + medianhousevalue +
                              poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                              as.factor(year) + as.factor(region)
                            +offset(log(time_count))+
                 as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_male<-c(loglinear_coefs_boots_white_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_white_male_trimmed.RData"))
rm(white_male_qd)
gc()


#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_qd_female <- subset(black_qd_female,
                          pm25_ensemble < quantile(black_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(black_qd_female$pm25_ensemble,0.05))
cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_qd_female, list(black_qd_female$zip))
num_uniq_zip <- length(unique(black_qd_female$zip))
loglinear_coefs_boots_black_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                            mean_bmi + smoke_rate + hispanic + pct_blk +
                            medhouseholdincome + medianhousevalue +
                            poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                            as.factor(year) + as.factor(region)
                          +offset(log(time_count))+
                 as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_female<-c(loglinear_coefs_boots_black_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_black_female_trimmed.RData"))
rm(black_female_qd)

#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
black_qd_male <- subset(black_qd_male,
                        pm25_ensemble < quantile(black_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(black_qd_male$pm25_ensemble,0.05))

cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_qd_male, list(black_qd_male$zip))
num_uniq_zip <- length(unique(black_qd_male$zip))
loglinear_coefs_boots_black_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                              mean_bmi + smoke_rate + hispanic + pct_blk +
                              medhouseholdincome + medianhousevalue +
                              poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                              as.factor(year) + as.factor(region)
                            +offset(log(time_count))+
                 as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_male<-c(loglinear_coefs_boots_black_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_black_male_trimmed.RData"))


#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_qd_female <- subset(hispanic_qd_female,
                             pm25_ensemble< quantile(hispanic_qd_female$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(hispanic_qd_female$pm25_ensemble,0.05))

cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_qd_female, list(hispanic_qd_female$zip))
num_uniq_zip <- length(unique(hispanic_qd_female$zip))
loglinear_coefs_boots_hispanic_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)
                                         , data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_female<-c(loglinear_coefs_boots_hispanic_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_female_trimmed.RData"))
rm(hispanic_qd_female)


#Hispanic Male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
hispanic_qd_male <- subset(hispanic_qd_male,
                           pm25_ensemble < quantile(hispanic_qd_male$pm25_ensemble,0.95)&
                             pm25_ensemble > quantile(hispanic_qd_male$pm25_ensemble,0.05))

cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_qd_male, list(hispanic_qd_male$zip))
num_uniq_zip <- length(unique(hispanic_qd_male$zip))
loglinear_coefs_boots_hispanic_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_male<-c(loglinear_coefs_boots_hispanic_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_male_trimmed.RData"))
rm(hispanic_qd_male)


#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_qd_female <- subset(asian_qd_female,
                          pm25_ensemble < quantile(asian_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(asian_qd_female$pm25_ensemble,0.05))
cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_qd_female, list(asian_qd_female$zip))
num_uniq_zip <- length(unique(asian_qd_female$zip))
loglinear_coefs_boots_asian_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_female<-c(loglinear_coefs_boots_asian_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_asian_female_trimmed.RData"))
rm(asian_qd_female)

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
asian_qd_male <- subset(asian_qd_male,
                        pm25_ensemble < quantile(asian_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(asian_qd_male$pm2,0.05))
cl=makeCluster(10,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_qd_male, list(asian_qd_male$zip))
num_uniq_zip <- length(unique(asian_qd_male$zip))
loglinear_coefs_boots_asian_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25_ensemble + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_male<-c(loglinear_coefs_boots_asian_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_asian_male_trimmed.RData"))
rm(asian_qd_male)

#load(paste0(dir_data,"Poisson.RData"))
#exp(10*(Poisson$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
#exp(10*(Poisson$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

