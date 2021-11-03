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

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

rm( covariates_rm, GPS_rm)

#Randall Martin
# Cox-equvalent conditional Poisson Regression
aggregate_data_rm <- subset(aggregate_data_rm,
                        pm25 < quantile(aggregate_data_rm$pm25,0.95)&
                          pm25 > quantile(aggregate_data_rm$pm25,0.05))
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
rm(aggregate_data_rm)

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)
rm(covariates_rm, GPS_rm, GPS_mod_rm, gnm_raw_rm)

#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_female_rm <- subset(white_female_rm,
                          pm25 < quantile(white_female_rm$pm25,0.95)&
                            pm25 > quantile(white_female_rm$pm25,0.05))
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
rm(white_female_rm, gnm_raw_rm_white_female)

#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
white_male_rm <- subset(white_male_rm,
                        pm25 < quantile(white_male_rm$pm25,0.95)&
                          pm25 > quantile(white_male_rm$pm25,0.05))
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
rm(white_male_rm, gnm_raw_rm_white_male)

#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_rm_female <- subset(black_rm_female,
                          pm25 < quantile(black_rm_female$pm25,0.95)&
                            pm25 > quantile(black_rm_female$pm25,0.05))
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
rm(black_rm_female, gnm_raw_rm_black_female)

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
black_rm_male <- subset(black_rm_male,
                        pm25 < quantile(black_rm_male$pm25,0.95)&
                          pm25 > quantile(black_rm_male$pm25,0.05))
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
rm(black_rm_male, gnm_raw_rm_black_male)

#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_rm_female <- subset(hispanic_rm_female,
                             pm25< quantile(hispanic_rm_female$pm25,0.95)&
                               pm25 > quantile(hispanic_rm_female$pm25,0.05))

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
rm(hispanic_rm_female, gnm_raw_rm_hispanic_female)

#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
hispanic_rm_male <- subset(hispanic_rm_male,
                           pm25 < quantile(hispanic_rm_male$pm25,0.95)&
                             pm25 > quantile(hispanic_rm_male$pm25,0.05))
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
rm(hispanic_rm_male, gnm_raw_rm_hispanic_male)

#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_rm_female <- subset(asian_rm_female,
                          pm25 < quantile(asian_rm_female$pm25,0.95)&
                            pm25 > quantile(asian_rm_female$pm25,0.05))

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
rm(asian_rm_female, gnm_raw_rm_asian_female)

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
asian_rm_male <- subset(asian_rm_male,
                        pm25 < quantile(asian_rm_male$pm25,0.95)&
                          pm25 > quantile(asian_rm_male$pm2,0.05))

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
rm(asian_rm_male, gnm_raw_rm_asian_male)
rm(aggregate_data_rm, mod_sd_rm)

#QD
load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_rm, GPS_mod_qd)

aggregate_date_qd <- subset(aggregate_data_qd,
                            pm25_ensemble < quantile(aggregate_data_qd$pm25_ensemble,0.95)&
                              pm25_ensemble > quantile(aggregate_data_qd$pm25_ensemble,0.05))
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
rm(gnm_raw_qd)

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_qd)

#By race

#White Female
white_female_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_female_qd <- subset(white_female_qd,
                          pm25_ensemble < quantile(white_female_qd$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(white_female_qd$pm25_ensemble,0.05))
gnm_raw_qd_white_female<-gnm(dead~  pm25_ensemble+
                        mean_bmi + smoke_rate + hispanic + pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                      data=white_female_qd,family=poisson(link="log"))
Poisson_qd_white_female<-summary(gnm_raw_qd_white_female)
exp(10*Poisson_qd_white_female$coefficients[1])
rm(white_female_qd, gnm_raw_qd_white_female)

#White male
white_male_qd<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
white_male_qd <- subset(white_male_qd,
                        pm25_ensemble < quantile(white_male_qd$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(white_male_qd$pm25_ensemble,0.05))

gnm_raw_qd_white_male<-gnm(dead~  pm25_ensemble+
                               mean_bmi + smoke_rate + hispanic + pct_blk +
                               medhouseholdincome + medianhousevalue +
                               poverty + education + popdensity + pct_owner_occ +
                               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                               as.factor(year) + as.factor(region)
                             +offset(log(time_count)),
                             eliminate= (as.factor(sex):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
                             data=white_male_qd,family=poisson(link="log"))
Poisson_qd_white_male<-summary(gnm_raw_qd_white_male)
exp(10*Poisson_qd_white_male$coefficients[1])
rm(white_male_qd, gnm_raw_qd_white_male)


#Black Female
black_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_qd_female <- subset(black_qd_female,
                          pm25_ensemble < quantile(black_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(black_qd_female$pm25_ensemble,0.05))

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
rm(black_qd_female, gnm_raw_qd_black_female)

#Black male
black_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
black_qd_male <- subset(black_qd_male,
                        pm25_ensemble < quantile(black_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(black_qd_male$pm25_ensemble,0.05))

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
rm(black_qd_male, gnm_raw_qd_black_male)

#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_qd_female <- subset(hispanic_qd_female,
                             pm25_ensemble < quantile(hispanic_qd_female$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(hispanic_qd_female$pm25_ensemble,0.05))

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
rm(hispanic_qd_female, gnm_raw_qd_hispanic_female)

#Hispanic male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
hispanic_qd_male <- subset(hispanic_qd_male,
                           pm25_ensemble < quantile(hispanic_qd_male$pm25_ensemble,0.95)&
                             pm25_ensemble > quantile(hispanic_qd_male$pm25_ensemble,0.05))

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
rm(hispanic_qd_male, gnm_raw_qd_hispanic_male)

#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_qd_female <- subset(asian_qd_female,
                          pm25_ensemble < quantile(asian_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(asian_qd_female$pm25_ensemble,0.05))

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
rm(asian_qd_female, gnm_raw_qd_asian_female)

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
asian_qd_male <- subset(asian_qd_male,
                        pm25_ensemble < quantile(asian_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(asian_qd_male$pm25_ensemble,0.05))
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
rm(asian_qd_male, gnm_raw_qd_asian_male)
rm(aggregate_data_qd)

save(Poisson_qd, 
     Poisson_qd_asian_female,Poisson_qd_asian_male, 
     Poisson_qd_black_female, Poisson_qd_black_male, 
    Poisson_qd_hispanic_female, Poisson_qd_hispanic_male,
     Poisson_qd_white_female, Poisson_qd_white_male, 
    Poisson_rm, 
     Poisson_rm_asian_female, Poisson_rm_asian_male, 
    Poisson_rm_black_female, Poisson_rm_black_male, 
      Poisson_rm_hispanic_female, Poisson_rm_hispanic_male, 
      Poisson_rm_white_female, Poisson_rm_white_male, file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main_trimmed.RData")
