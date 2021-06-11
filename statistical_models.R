# All types of statistical models implemented for the paper
library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
library(fst)
require(xgboost)
require(parallel)
require(dplyr)
require(mgcv)

#Load Poisson model
#load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main.RData")

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, covariates_rm, GPS_qd, GPS_rm)

#Randall Martin
# Cox-equvalent conditional Poisson Regression
gnm_raw_rm<-gnm(dead~  pm25 + 
               mean_bmi + smoke_rate + hispanic+ pct_blk +
                 medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year) + as.factor(region)
               +offset(log(time_count)),
               eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
               data=aggregate_data_rm,family=poisson(link="log"))
Poisson_rm<-summary(gnm_raw_rm)
exp(10*Poisson_rm$coefficients[1])


#Race
#White
require(dplyr)
white_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1)
gnm_raw_rm_white<-gnm(dead~  pm25 + 
                   mean_bmi + smoke_rate + hispanic + pct_blk +
                   medhouseholdincome + medianhousevalue +
                   poverty + education + popdensity + pct_owner_occ +
                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                   as.factor(year) + as.factor(region)
                 +offset(log(time_count)),
                 eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                 data=white_rm,family=poisson(link="log"))
Poisson_rm_white<-summary(gnm_raw_rm_white)
exp(10*Poisson_rm_white$coefficients[1])


#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
gnm_raw_rm_white_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)+offset(log(time_count)),
                      eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=white_female_rm ,family=poisson(link="log"))
Poisson_rm_white_female<-summary(gnm_raw_rm_white_female)
exp(10*Poisson_rm_white_female$coefficients[1])

#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
gnm_raw_rm_white_male<-gnm(dead~  pm25 + 
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)+offset(log(time_count)),
                             eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=white_male_rm ,family=poisson(link="log"))
Poisson_rm_white_male<-summary(gnm_raw_rm_white_male)
exp(10*Poisson_rm_white_male$coefficients[1])

#Black
black_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2)
gnm_raw_rm_black<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=black_rm,family=poisson(link="log"))
Poisson_rm_black<-summary(gnm_raw_rm_black)
exp(10*Poisson_rm_black$coefficients[1])

#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
gnm_raw_rm_black_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=black_rm_female,family=poisson(link="log"))
Poisson_rm_black_female<-summary(gnm_raw_rm_black_female)
exp(10*Poisson_rm_black_female$coefficients[1])

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
gnm_raw_rm_black_male<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=black_rm_male,family=poisson(link="log"))
Poisson_rm_black_male<-summary(gnm_raw_rm_black_male)
exp(10*Poisson_rm_black_male$coefficients[1])

#Hispanic
hispanic_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5)
gnm_raw_rm_hispanic<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=hispanic_rm,family=poisson(link="log"))
Poisson_rm_hispanic<-summary(gnm_raw_rm_hispanic)
exp(10*Poisson_rm_hispanic$coefficients[1])

#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
gnm_raw_rm_hispanic_female<-gnm(dead~  pm25 + 
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count)),
                         eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                         data=hispanic_rm_female,family=poisson(link="log"))
Poisson_rm_hispanic_female<-summary(gnm_raw_rm_hispanic_female)
exp(10*Poisson_rm_hispanic_female$coefficients[1])

#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
gnm_raw_rm_hispanic_male<-gnm(dead~  pm25 + 
                                  mean_bmi + smoke_rate + hispanic + pct_blk +
                                  medhouseholdincome + medianhousevalue +
                                  poverty + education + popdensity + pct_owner_occ +
                                  summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                  as.factor(year) + as.factor(region)
                                +offset(log(time_count)),
                                eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                                data=hispanic_rm_male,family=poisson(link="log"))
Poisson_rm_hispanic_male<-summary(gnm_raw_rm_hispanic_male)
exp(10*Poisson_rm_hispanic_male$coefficients[1])

#Asian
asian_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4)
gnm_raw_rm_asian<-gnm(dead~  pm25 + 
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count)),
                         eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                         data=asian_rm,family=poisson(link="log"))
Poisson_rm_asian<-summary(gnm_raw_rm_asian)
exp(10*Poisson_rm_asian$coefficients[1])

#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
gnm_raw_rm_asian_female<-gnm(dead~  pm25 + 
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=asian_rm_female,family=poisson(link="log"))
Poisson_rm_asian_female<-summary(gnm_raw_rm_asian_female)
exp(10*Poisson_rm_asian_female$coefficients[1])

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
gnm_raw_rm_asian_male<-gnm(dead~  pm25 + 
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count)),
                             eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=asian_rm_male,family=poisson(link="log"))
Poisson_rm_asian_male<-summary(gnm_raw_rm_asian_male)
exp(10*Poisson_rm_asian_male$coefficients[1])

#QD
gnm_raw_qd<-gnm(dead~  pm25_ensemble + 
                  mean_bmi + smoke_rate + hispanic + pct_blk +
                  medhouseholdincome + medianhousevalue +
                  poverty + education + popdensity + pct_owner_occ +
                  summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                  as.factor(year) + as.factor(region)
                +offset(log(time_count)),
                eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                data=aggregate_data_qd,family=poisson(link="log"))
Poisson_qd<-summary(gnm_raw_qd)
exp(10*Poisson_qd$coefficients[1])

#By race
#White
white_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1)
gnm_raw_qd_white<-gnm(dead~  pm25_ensemble+
                   mean_bmi + smoke_rate + hispanic + pct_blk +
                   medhouseholdincome + medianhousevalue +
                   poverty + education + popdensity + pct_owner_occ +
                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                   as.factor(year) + as.factor(region)
                 +offset(log(time_count)),
                 eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                 data=white_qd,family=poisson(link="log"))
Poisson_qd_white<-summary(gnm_raw_qd_white)
exp(10*Poisson_qd_white$coefficients[1])

#White Female
white_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
gnm_raw_qd_white_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=white_qd_female,family=poisson(link="log"))
Poisson_qd_white_female<-summary(gnm_raw_qd_white_female)
exp(10*Poisson_qd_white_female$coefficients[1])

#White male
white_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
gnm_raw_qd_white_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count)),
                             eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=white_qd_male,family=poisson(link="log"))
Poisson_qd_white_male<-summary(gnm_raw_qd_white_male)
exp(10*Poisson_qd_white_male$coefficients[1])

#Black
black_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2)
gnm_raw_qd_black<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=black_qd,family=poisson(link="log"))
Poisson_qd_black<-summary(gnm_raw_qd_black)
exp(10*Poisson_qd_black$coefficients[1])

#Black Female
black_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
gnm_raw_qd_black_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=black_qd_female,family=poisson(link="log"))
Poisson_qd_black_female<-summary(gnm_raw_qd_black_female)
exp(10*Poisson_qd_black_female$coefficients[1])

#Black male
black_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
gnm_raw_qd_black_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count)),
                             eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=black_qd_male,family=poisson(link="log"))
Poisson_qd_black_male<-summary(gnm_raw_qd_black_male)
exp(10*Poisson_qd_black_male$coefficients[1])

#Hispanic
hispanic_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5)
gnm_raw_qd_hispanic<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=hispanic_qd,family=poisson(link="log"))
Poisson_qd_hispanic<-summary(gnm_raw_qd_hispanic)
exp(10*Poisson_qd_hispanic$coefficients[1])

#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
gnm_raw_qd_hispanic_female<-gnm(dead~  pm25_ensemble+
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count)),
                         eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                         data=hispanic_qd_female,family=poisson(link="log"))
Poisson_qd_hispanic_female<-summary(gnm_raw_qd_hispanic_female)
exp(10*Poisson_qd_hispanic_female$coefficients[1])

#Hispanic male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
gnm_raw_qd_hispanic_male<-gnm(dead~  pm25_ensemble+
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count)),
                         eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                         data=hispanic_qd_male,family=poisson(link="log"))
Poisson_qd_hispanic_male<-summary(gnm_raw_qd_hispanic_male)
exp(10*Poisson_qd_hispanic_male$coefficients[1])

#Asian
asian_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4)
gnm_raw_qd_asian<-gnm(dead~  pm25_ensemble+
                           mean_bmi + smoke_rate + hispanic + pct_blk +
                           medhouseholdincome + medianhousevalue +
                           poverty + education + popdensity + pct_owner_occ +
                           summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                           as.factor(year) + as.factor(region)
                         +offset(log(time_count)),
                         eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                         data=asian_qd,family=poisson(link="log"))
Poisson_qd_asian<-summary(gnm_raw_qd_asian)
exp(10*Poisson_qd_asian$coefficients[1])

#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
gnm_raw_qd_asian_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=asian_qd_female,family=poisson(link="log"))
Poisson_qd_asian_female<-summary(gnm_raw_qd_asian_female)
exp(10*Poisson_qd_asian_female$coefficients[1])

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
gnm_raw_qd_asian_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count)),
                             eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=asian_qd_male,family=poisson(link="log"))
Poisson_qd_asian_male<-summary(gnm_raw_qd_asian_male)
exp(10*Poisson_qd_asian_male$coefficients[1])

save(Poisson_qd, Poisson_qd_asian, Poisson_qd_asian_female,Poisson_qd_asian_male, 
     Poisson_qd_black, Poisson_qd_black_female, Poisson_qd_black_male, 
     Poisson_qd_hispanic, Poisson_qd_hispanic_female, Poisson_qd_hispanic_male,
     Poisson_qd_white, Poisson_qd_white_female, Poisson_qd_white_male, Poisson_rm, 
     Poisson_rm_asian, Poisson_rm_asian_female, Poisson_rm_asian_male, 
     Poisson_rm_black, Poisson_rm_black_female, Poisson_rm_black_male, 
     Poisson_rm_hispanic, Poisson_rm_hispanic_female, Poisson_rm_hispanic_male, 
     Poisson_rm_white, Poisson_rm_white_female, 
      Poisson_rm_white_male, file= "/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main.RData")
