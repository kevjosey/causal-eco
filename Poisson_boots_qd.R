library("parallel")
library("dplyr")
require(parallel)
library(data.table)
require(tidyr)
require(doParallel)
library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
library(fst)
require(xgboost)
require(parallel)
require(dplyr)

#Load Poisson model
dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/Poisson_qd/'

load(file= "/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main.RData")
load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)


# qd
#all
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots<-c(loglinear_coefs_boots,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(loglinear_coefs_boots,file=paste0(dir_out,"loglinear_coefs_boots_qd_all.RData"))



#Race

#White female
require(dplyr)
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_female<-c(loglinear_coefs_boots_white_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_white_female.RData"))
rm(white_female_qd)
gc()

#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_male<-c(loglinear_coefs_boots_white_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_white_male.RData"))
rm(white_male_qd)
gc()


#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_female<-c(loglinear_coefs_boots_black_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_black_female.RData"))
rm(black_qd_female)

#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_male<-c(loglinear_coefs_boots_black_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_black_male.RData"))


#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_female<-c(loglinear_coefs_boots_hispanic_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_female.RData"))
rm(hispanic_qd_female)


#Hispanic Male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_male<-c(loglinear_coefs_boots_hispanic_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_male.RData"))
rm(hispanic_qd_male)



#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_female<-c(loglinear_coefs_boots_asian_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_female,file=paste0(dir_out,"loglinear_coefs_boots_qd_asian_female.RData"))
rm(asian_qd_female)

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
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
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_qdax + winter_qdax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count)),eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_male<-c(loglinear_coefs_boots_asian_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_male,file=paste0(dir_out,"loglinear_coefs_boots_qd_asian_male.RData"))
rm(asian_qd_male)

