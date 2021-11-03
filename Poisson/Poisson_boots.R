# Bootstrap on zip-code cluster to obtain robust CIs account for spatial correlation
library("mgcv")
library("parallel")
library("dplyr")
require(parallel)
library(data.table)
require(tidyr)

set.seed(1234)
#Load Poisson model
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_rm/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)
rm(covariates_rm, GPS_mod_rm)


# RM
#all
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(aggregate_data_rm, list(aggregate_data_rm$zip))
num_uniq_zip <- length(unique(aggregate_data_rm$zip))

loglinear_coefs_boots<-NULL
for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))
               +(as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots<-c(loglinear_coefs_boots,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(loglinear_coefs_boots,file=paste0(dir_out,"loglinear_coefs_boots_rm_all.RData"))

#Race
#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_female_rm, list(white_female_rm$zip))
num_uniq_zip <- length(unique(white_female_rm$zip))
loglinear_coefs_boots_white_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                       mean_bmi + smoke_rate + hispanic + pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_female<-c(loglinear_coefs_boots_white_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_female,file=paste0(dir_out,"loglinear_coefs_boots_rm_white_female.RData"))
rm(white_female_rm)
gc()

#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_male_rm, list(white_male_rm$zip))
num_uniq_zip <- length(unique(white_male_rm$zip))
loglinear_coefs_boots_white_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                              mean_bmi + smoke_rate + hispanic + pct_blk +
                              medhouseholdincome + medianhousevalue +
                              poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                              as.factor(year) + as.factor(region)
                            +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_white_male<-c(loglinear_coefs_boots_white_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_white_male,file=paste0(dir_out,"loglinear_coefs_boots_rm_white_male.RData"))
rm(white_male_rm)
gc()


#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_rm_female, list(black_rm_female$zip))
num_uniq_zip <- length(unique(black_rm_female$zip))
loglinear_coefs_boots_black_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                            mean_bmi + smoke_rate + hispanic + pct_blk +
                            medhouseholdincome + medianhousevalue +
                            poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                            as.factor(year) + as.factor(region)
                          +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                  data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_female<-c(loglinear_coefs_boots_black_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_female,file=paste0(dir_out,"loglinear_coefs_boots_rm_black_female.RData"))
rm(black_female_rm)

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_rm_male, list(black_rm_male$zip))
num_uniq_zip <- length(unique(black_male_rm$zip))
loglinear_coefs_boots_black_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                              mean_bmi + smoke_rate + hispanic + pct_blk +
                              medhouseholdincome + medianhousevalue +
                              poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                              as.factor(year) + as.factor(region)
                            +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_black_male<-c(loglinear_coefs_boots_black_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_black_male,file=paste0(dir_out,"loglinear_coefs_boots_rm_black_male.RData"))



#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_rm_female, list(hispanic_rm_female$zip))
num_uniq_zip <- length(unique(hispanic_rm_female$zip))
loglinear_coefs_boots_hispanic_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_female<-c(loglinear_coefs_boots_hispanic_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_female,file=paste0(dir_out,"loglinear_coefs_boots_rm_hispanic_female.RData"))
rm(hispanic_rm_female)


#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_rm_male, list(hispanic_rm_male$zip))
num_uniq_zip <- length(unique(hispanic_rm_male$zip))
loglinear_coefs_boots_hispanic_male<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_hispanic_male<-c(loglinear_coefs_boots_hispanic_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_hispanic_male,file=paste0(dir_out,"loglinear_coefs_boots_rm_hispanic_male.RData"))
rm(hispanic_rm_male)



#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_rm_female, list(asian_rm_female$zip))
num_uniq_zip <- length(unique(asian_rm_female$zip))
loglinear_coefs_boots_asian_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_female<-c(loglinear_coefs_boots_asian_female,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_female,file=paste0(dir_out,"loglinear_coefs_boots_rm_asian_female.RData"))
rm(asian_rm_female)

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
cl=makeCluster(14,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_rm_male, list(asian_rm_male$zip))
num_uniq_zip <- length(unique(asian_rm_male$zip))
loglinear_coefs_boots_asian_female<-NULL

for (boots_id in 1:500){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_raw<-gnm(dead~  pm25 + 
                 mean_bmi + smoke_rate + hispanic + pct_blk +
                 medhouseholdincome + medianhousevalue +
                 poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                 as.factor(year) + as.factor(region)
               +offset(log(time_count))+
                 (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), data=aggregate_data_boots,family=poisson(link="log"))
  
  loglinear_coefs_boots_asian_male<-c(loglinear_coefs_boots_asian_male,summary(gnm_raw)$coefficients[1])
  rm(aggregate_data_boots)
}
stopCluster(cl)
save(loglinear_coefs_boots_asian_male,file=paste0(dir_out,"loglinear_coefs_boots_rm_asian_male.RData"))
rm(asian_rm_male)

#load(paste0(dir_data,"Poisson.RData"))
#exp(10*(Poisson$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
#exp(10*(Poisson$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

