# All types of statistical models implemented for the paper
#library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
#library(fst)
#require(xgboost)
require(parallel)
require(dplyr)
#require(mgcv)

#Load Poisson model
#load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main.RData")

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, covariates_rm)

#Randall Martin
# Cox-equvalent conditional Poisson Regression
gnm_raw_rm<-gnm(dead~  pm25 + 
               mean_bmi + smoke_rate + hispanic+ pct_blk +
                 medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year) + as.factor(region)
               +offset(log(time_count))+(as.factor(sex)+as.factor(race)+as.factor(dual)+
                                           as.factor(entry_age_break)+as.factor(followup_year)), 
               data=aggregate_data_rm,family=poisson(link="log"))
Poisson_rm<-(gnm_raw_rm$coefficients)
exp(10*Poisson_rm[2])


#Race
#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
gnm_raw_rm_white_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)+offset(log(time_count))+
                          (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=white_female_rm ,family=poisson(link="log"))
Poisson_rm_white_female<-(gnm_raw_rm_white_female$coefficients)
exp(10*Poisson_rm_white_female[2])

#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
gnm_raw_rm_white_male<-gnm(dead~  pm25 + 
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)+offset(log(time_count))
                           +(as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                             data=white_male_rm ,family=poisson(link="log"))
Poisson_rm_white_male<-(gnm_raw_rm_white_male$coefficients)
exp(10*Poisson_rm_white_male[2])

#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
gnm_raw_rm_black_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=black_rm_female,family=poisson(link="log"))
Poisson_rm_black_female<-(gnm_raw_rm_black_female$coefficients)
exp(10*Poisson_rm_black_female[2])

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
gnm_raw_rm_black_male<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=black_rm_male,family=poisson(link="log"))
Poisson_rm_black_male<-(gnm_raw_rm_black_male$coefficients)
exp(10*Poisson_rm_black_male[2])

#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
gnm_raw_rm_hispanic_female<-gnm(dead~  pm25 + 
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count))+
                           (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                         data=hispanic_rm_female,family=poisson(link="log"))
Poisson_rm_hispanic_female<-(gnm_raw_rm_hispanic_female$coefficients)
exp(10*Poisson_rm_hispanic_female[2])

#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
gnm_raw_rm_hispanic_male<-gnm(dead~  pm25 + 
                                  mean_bmi + smoke_rate + hispanic + pct_blk +
                                  medhouseholdincome + medianhousevalue +
                                  poverty + education + popdensity + pct_owner_occ +
                                  summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                  as.factor(year) + as.factor(region)
                                +offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                                data=hispanic_rm_male,family=poisson(link="log"))
Poisson_rm_hispanic_male<-(gnm_raw_rm_hispanic_male$coefficients)
exp(10*Poisson_rm_hispanic_male[2])


#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
gnm_raw_rm_asian_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=asian_rm_female,family=poisson(link="log"))
Poisson_rm_asian_female<-(gnm_raw_rm_asian_female$coefficients)
exp(10*Poisson_rm_asian_female[2])

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
gnm_raw_rm_asian_male<-gnm(dead~  pm25 + 
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count))+
                             (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                             data=asian_rm_male,family=poisson(link="log"))
Poisson_rm_asian_male<-(gnm_raw_rm_asian_male$coefficients)
exp(10*Poisson_rm_asian_male[2])

#QD
gnm_raw_qd<-gnm(dead~  pm25_ensemble + 
                  mean_bmi + smoke_rate + hispanic + pct_blk +
                  medhouseholdincome + medianhousevalue +
                  poverty + education + popdensity + pct_owner_occ +
                  summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                  as.factor(year) + as.factor(region)
                +offset(log(time_count))+
                  (as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                data=aggregate_data_qd,family=poisson(link="log"))
Poisson_qd<- (gnm_raw_qd$coefficients)
exp(10*Poisson_qd[2])

#By race
#White Female
white_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
gnm_raw_qd_white_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=white_qd_female,family=poisson(link="log"))
Poisson_qd_white_female<-(gnm_raw_qd_white_female$coefficients)
exp(10*Poisson_qd_white_female[2])

#White male
white_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
gnm_raw_qd_white_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count))+
                             (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                             data=white_qd_male,family=poisson(link="log"))
Poisson_qd_white_male<-(gnm_raw_qd_white_male$coefficients)
exp(10*Poisson_qd_white_male[2])

#Black Female
black_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
gnm_raw_qd_black_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=black_qd_female,family=poisson(link="log"))
Poisson_qd_black_female<-(gnm_raw_qd_black_female$coefficients)
exp(10*Poisson_qd_black_female[2])

#Black male
black_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
gnm_raw_qd_black_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count))+
                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year), 
                             data=black_qd_male,family=poisson(link="log"))
Poisson_qd_black_male<-(gnm_raw_qd_black_male$coefficients)
exp(10*Poisson_qd_black_male[2])

#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
gnm_raw_qd_hispanic_female<-gnm(dead~  pm25_ensemble+
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count))+
                           (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                         data=hispanic_qd_female,family=poisson(link="log"))
Poisson_qd_hispanic_female<-(gnm_raw_qd_hispanic_female$coefficients)
exp(10*Poisson_qd_hispanic_female[2])

#Hispanic male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
gnm_raw_qd_hispanic_male<-gnm(dead~  pm25_ensemble+
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count))+
                           (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                         data=hispanic_qd_male,family=poisson(link="log"))
Poisson_qd_hispanic_male<-(gnm_raw_qd_hispanic_male$coefficients)
exp(10*Poisson_qd_hispanic_male[2])

#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
gnm_raw_qd_asian_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count))+
                        (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                      data=asian_qd_female,family=poisson(link="log"))
Poisson_qd_asian_female<-(gnm_raw_qd_asian_female$coefficients)
exp(10*Poisson_qd_asian_female[2])

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
gnm_raw_qd_asian_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count))+
                             (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
                             data=asian_qd_male,family=poisson(link="log"))
Poisson_qd_asian_male<-(gnm_raw_qd_asian_male$coefficients)
exp(10*Poisson_qd_asian_male[2])

save(Poisson_qd, Poisson_qd_asian_female,Poisson_qd_asian_male, 
     Poisson_qd_black_female, Poisson_qd_black_male, 
     Poisson_qd_hispanic_female, Poisson_qd_hispanic_male,
     Poisson_qd_white_female, Poisson_qd_white_male, Poisson_rm, 
     Poisson_rm_asian_female, Poisson_rm_asian_male, 
    Poisson_rm_black_female, Poisson_rm_black_male, 
     Poisson_rm_hispanic_female, Poisson_rm_hispanic_male, 
     Poisson_rm_white_female, 
      Poisson_rm_white_male, file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main.RData")
