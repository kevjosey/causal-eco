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
#require(devtools)
require(mgcv)
require(gam)
require(KernSmooth)

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

#Randall Martin
# Cox-equvalent conditional gam_raw Regression
aggregate_data_rm <- subset(aggregate_data_rm,
                        pm25 < quantile(aggregate_data_rm$pm25,0.95)&
                          pm25 > quantile(aggregate_data_rm$pm25,0.05))

gam_raw_rm<-mgcv::bam(dead~  s(pm25, k=3) +
                        as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      data=aggregate_data_rm,family=poisson(link="log"))

save(gam_raw_rm, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_all_trimmed.RData")
rm(gam_raw_rm, aggregate_data_rm)
gc()


#White female
require(dplyr)
load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_female_rm <- subset(white_female_rm,
                          pm25 < quantile(white_female_rm$pm25,0.95)&
                            pm25 > quantile(white_female_rm$pm25,0.05))


gam_raw_rm_white_female<-mgcv::bam(dead~  s(pm25, k=3) +
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                                   data=white_female_rm ,family=poisson(link="log"))

rm(white_female_rm)
gc()
save(gam_raw_rm_white_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_female_trimmed.RData")
rm(gam_raw_rm_white_female)
gc()


#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
white_male_rm <- subset(white_male_rm,
                        pm25 < quantile(white_male_rm$pm25,0.95)&
                          pm25 > quantile(white_male_rm$pm25,0.05))


gam_raw_rm_white_male<-mgcv::bam(dead~  s(pm25, k=3) +
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)),
                                 data=white_male_rm ,family=poisson(link="log"))
rm(white_male_rm)
gc()
save(gam_raw_rm_white_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_male_trimmed.RData")
rm(gam_raw_rm_white_male)
gc()


#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_rm_female <- subset(black_rm_female,
                          pm25 < quantile(black_rm_female$pm25,0.95)&
                            pm25 > quantile(black_rm_female$pm25,0.05))


gam_raw_rm_black_female<-mgcv::bam(dead~  s(pm25, k=3) +
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                                   data=black_rm_female,family=poisson(link="log"))
rm(black_rm_female)
gc()
save(gam_raw_rm_black_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_female_trimmed.RData")
rm(gam_raw_rm_black_female)
gc()

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
black_rm_male <- subset(black_rm_male,
                        pm25 < quantile(black_rm_male$pm25,0.95)&
                          pm25 > quantile(black_rm_male$pm25,0.05))

gam_raw_rm_black_male<-mgcv::bam(dead~  s(pm25, k=3) +
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)),
                                 data=black_rm_male,family=poisson(link="log"))
rm(black_rm_male)
gc()
save(gam_raw_rm_black_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_male_trimmed.RData")
rm(gam_raw_rm_black_male)
gc()



#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)

hispanic_rm_female <- subset(hispanic_rm_female,
                             pm25< quantile(hispanic_rm_female$pm25,0.95)&
                               pm25 > quantile(hispanic_rm_female$pm25,0.05))

gam_raw_rm_hispanic_female<-mgcv::bam(dead~  s(pm25, k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                                        medhouseholdincome + medianhousevalue +
                                        poverty + education + popdensity + pct_owner_occ +
                                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                        as.factor(year) + as.factor(region)
                                      +offset(log(time_count)),
                                      data=hispanic_rm_female,family=poisson(link="log"))
rm(hispanic_rm_female)
gc()
save(gam_raw_rm_hispanic_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_female_trimmed.RData")
rm(gam_raw_rm_hispanic_female)
gc()


#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)

hispanic_rm_male <- subset(hispanic_rm_male,
                           pm25 < quantile(hispanic_rm_male$pm25,0.95)&
                             pm25 > quantile(hispanic_rm_male$pm25,0.05))

gam_raw_rm_hispanic_male<-mgcv::bam(dead~  s(pm25, k=3) +
                                      as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                      mean_bmi + smoke_rate + hispanic+ pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(year) + as.factor(region)
                                    +offset(log(time_count)),
                                    data=hispanic_rm_male,family=poisson(link="log"))
rm(hispanic_rm_male)
gc()
save(gam_raw_rm_hispanic_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_male_trimmed.RData")
rm(gam_raw_rm_hispanic_male)
gc()



#Asian female
asian_rm_female<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)

asian_rm_female <- subset(asian_rm_female,
                          pm25 < quantile(asian_rm_female$pm25,0.95)&
                            pm25 > quantile(asian_rm_female$pm25,0.05))

gam_raw_rm_asian_female<-mgcv::bam(dead~  s(pm25, k=3) +
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                                   data=asian_rm_female,family=poisson(link="log"))
rm(asian_rm_female)
gc()
save(gam_raw_rm_asian_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_female_trimmed.RData")
rm(gam_raw_rm_asian_female)
gc()

#Asian male
asian_rm_male<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)

asian_rm_male <- subset(asian_rm_male,
                        pm25 < quantile(asian_rm_male$pm25,0.95)&
                          pm25 > quantile(asian_rm_male$pm2,0.05))

gam_raw_rm_asian_male<-mgcv::bam(dead~  s(pm25, k=3) +
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)),
                                 data=asian_rm_male,family=poisson(link="log"))
rm(asian_rm_male)
gc()
save(gam_raw_rm_asian_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_male_trimmed.RData")
rm(gam_raw_rm_asian_male)
gc()

rm(aggregate_data_rm)
gc()

########QD
load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)
rm(covariates_qd, GPS_mod_qd)

aggregate_date_qd <- subset(aggregate_data_qd,
                            pm25_ensemble < quantile(aggregate_data_qd$pm25_ensemble,0.95)&
                              pm25_ensemble > quantile(aggregate_data_qd$pm25_ensemble,0.05))
gam_raw_qd<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                        as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      data=aggregate_data_qd,family=poisson(link="log"))

save(gam_raw_qd, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_all_trimmed.RData")
rm(gam_raw_qd, aggregate_data_qd)
gc()


#White female
require(dplyr)
load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_female_qd <- subset(white_female_qd,
                          pm25_ensemble < quantile(white_female_qd$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(white_female_qd$pm25_ensemble,0.05))


gam_raw_qd_white_female<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)), 
                                   data=white_female_qd ,family=poisson(link="log"))

rm(white_female_qd)
gc()
save(gam_raw_qd_white_female, 
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_female_trimmed.RData")
rm(gam_raw_qd_white_female)
gc()


#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
white_male_qd <- subset(white_male_qd,
                        pm25_ensemble < quantile(white_male_qd$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(white_male_qd$pm25_ensemble,0.05))


gam_raw_qd_white_male<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)), 
                                 data=white_male_qd ,family=poisson(link="log"))
rm(white_male_qd)
gc()
save(gam_raw_qd_white_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_male_trimmed.RData")
rm(gam_raw_qd_white_male)
gc()


#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_qd_female <- subset(black_qd_female,
                          pm25_ensemble < quantile(black_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(black_qd_female$pm25_ensemble,0.05))


gam_raw_qd_black_female<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                                   data=black_qd_female,family=poisson(link="log"))
rm(black_qd_female)
gc()
save(gam_raw_qd_black_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_female_trimmed.RData")
rm(gam_raw_qd_black_female)
gc()

#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
black_qd_male <- subset(black_qd_male,
                        pm25_ensemble < quantile(black_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(black_qd_male$pm25_ensemble,0.05))

gam_raw_qd_black_male<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)), 
                                 data=black_qd_male,family=poisson(link="log"))
rm(black_qd_male)
gc()
save(gam_raw_qd_black_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_male_trimmed.RData")
rm(gam_raw_qd_black_male)
gc()



#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)

hispanic_qd_female <- subset(hispanic_qd_female,
                             pm25_ensemble < quantile(hispanic_qd_female$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(hispanic_qd_female$pm25_ensemble,0.05))

gam_raw_qd_hispanic_female<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                                        medhouseholdincome + medianhousevalue +
                                        poverty + education + popdensity + pct_owner_occ +
                                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                        as.factor(year) + as.factor(region)
                                      +offset(log(time_count)),
                                      data=hispanic_qd_female,family=poisson(link="log"))
rm(hispanic_qd_female)
gc()
save(gam_raw_qd_hispanic_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_female_trimmed.RData")
rm(gam_raw_qd_hispanic_female)
gc()


#Hispanic Male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)

hispanic_qd_male <- subset(hispanic_qd_male,
                           pm25_ensemble < quantile(hispanic_qd_male$pm25_ensemble,0.95)&
                             pm25_ensemble > quantile(hispanic_qd_male$pm25_ensemble,0.05))

gam_raw_qd_hispanic_male<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                      as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                      mean_bmi + smoke_rate + hispanic+ pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(year) + as.factor(region)
                                    +offset(log(time_count)),
                                    data=hispanic_qd_male,family=poisson(link="log"))
rm(hispanic_qd_male)
gc()
save(gam_raw_qd_hispanic_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_male_trimmed.RData")
rm(gam_raw_qd_hispanic_male)
gc()



#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)

asian_qd_female <- subset(asian_qd_female,
                          pm25_ensemble < quantile(asian_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(asian_qd_female$pm25_ensemble,0.05))

gam_raw_qd_asian_female<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                                   data=asian_qd_female,family=poisson(link="log"))
rm(asian_qd_female)
gc()
save(gam_raw_qd_asian_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_female_trimmed.RData")
rm(gam_raw_qd_asian_female)
gc()

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

asian_qd_male <- subset(asian_qd_male,
                        pm25_ensemble < quantile(asian_qd_male$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(asian_qd_male$pm25_ensemble,0.05))

gam_raw_qd_asian_male<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)),
                                 data=asian_qd_male,family=poisson(link="log"))
rm(asian_qd_male)
gc()
save(gam_raw_qd_asian_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_male_trimmed.RData")
rm(gam_raw_qd_asian_male)
gc()

rm(aggregate_data_qd)
gc()

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_all_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_male_trimmed.RData")

require(schoenberg)

test.data.rm<-function(x, gm){
  test<-data.frame(pm25 = seq(4, 15,length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50),
                   mean_bmi=rep(mean(x$mean_bmi), 50),
                   smoke_rate= rep(mean(x$smoke_rate), 50),
                   hispanic=rep(mean(x$hispanic), 50),
                   pct_blk= rep(mean(x$pct_blk), 50),
                   medhouseholdincome=rep(mean(x$medhouseholdincome), 50),
                   medianhousevalue= rep(mean(x$medianhousevalue), 50),
                   poverty= rep(mean(x$poverty), 50),
                   education = rep(mean(x$education), 50),
                   popdensity = rep(mean(x$popdensity), 50),
                   pct_owner_occ= rep(mean(x$pct_owner_occ), 50),
                   summer_tmmx= rep(mean(x$summer_tmmx), 50),
                   winter_tmmx= rep(mean(x$winter_tmmx), 50),
                   summer_rmax = rep(mean(x$summer_rmax), 50),
                   winter_rmax= rep(mean(x$winter_rmax), 50),
                   year= rep(levels(as.factor(x$year))[1], 50),
                   region= rep(levels(as.factor(x$region))[1], 50),
                   time_count= rep(1, 50)) # This will give the appropriate RR estimate we desire
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  dly<- pred.vals$fit
  se.dly<-pred.vals$se.fit
  dly_low<-dly-qnorm(0.975)*se.dly
  dly_high<-dly+qnorm(0.975)*se.dly
  
  #logRR
  test$hr<-gm$family$linkinv(dly)
  test$hr_upper<-gm$family$linkinv(dly_high)
  test$hr_lower<-gm$family$linkinv(dly_low)
  
  test$y<-(dly - dly[1])
  test$upper<-(dly_high - dly_high[1])
  test$lower<-(dly_low - dly_low[1])
  return(test)
}

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
load(paste0(dir_data,"aggregate_data_rm_trimmed.RData"))
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
black_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
hispanic_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
asian_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
#Plotting
testallrm<-test.data.rm(aggregate_data_rm, gam_raw_rm)
test_white_female_rm<-test.data.rm(white_female_rm , gam_raw_rm_white_female)
test_white_male_rm<-test.data.rm(white_male_rm, gam_raw_rm_white_male)
test_black_female_rm<-test.data.rm(black_female_rm, gam_raw_rm_black_female)
test_black_male_rm<-test.data.rm(black_male_rm, gam_raw_rm_black_male)
test_hispanic_female_rm<-test.data.rm(hispanic_female_rm, gam_raw_rm_hispanic_female)
test_hispanic_male_rm<-test.data.rm(hispanic_male_rm, gam_raw_rm_hispanic_male)
test_asian_female_rm<-test.data.rm(asian_female_rm, gam_raw_rm_asian_female)
test_asian_male_rm<-test.data.rm(asian_male_rm, gam_raw_rm_asian_male)

pall_rm<-ggplot(data=testallrm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
 # geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="All HR")
p_white_female_rm<-ggplot(data=test_white_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="White female HR")
p_white_male_rm<-ggplot(data=test_white_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="White male HR")
p_black_female_rm<-ggplot(data=test_black_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Black female HR")
p_black_male_rm<-ggplot(data=test_black_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Black male HR")
p_hispanic_female_rm<-ggplot(data=test_hispanic_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Hispanic female HR")
p_hispanic_male_rm<-ggplot(data=test_hispanic_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Hispanic male HR")
p_asian_female_rm<-ggplot(data=test_asian_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Asian female HR")
p_asian_male_rm<-ggplot(data=test_asian_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(10*y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Asian male HR")

plot_grid(p_white_female_rm, p_white_male_rm,
          p_black_female_rm, p_black_male_rm,
          p_hispanic_female_rm, p_hispanic_male_rm,
          p_asian_female_rm, p_asian_male_rm, ncol=2)

plot(gam_raw_rm, all.terms = FALSE, trans=exp, shift=0.05, ylab="All RM HR")

par(mfrow=c(4,2))
plot(gam_raw_rm_white_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="White female HR")
plot(gam_raw_rm_white_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="White male HR")
plot(gam_raw_rm_black_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black female HR")
plot(gam_raw_rm_black_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black male HR")
plot(gam_raw_rm_hispanic_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic female HR")
plot(gam_raw_rm_hispanic_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic male HR")
plot(gam_raw_rm_asian_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian female HR")
plot(gam_raw_rm_asian_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian male HR")

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_all_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_male_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_female_trimmed.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_male_trimmed.RData")

plot(gam_raw_qd, all.teqds = FALSE, trans=exp, shift=0.05, ylab="All qd HR")
par(mfrow=c(4,2))
plot(gam_raw_qd_white_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="White female HR")
plot(gam_raw_qd_white_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="White male HR")
plot(gam_raw_qd_black_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black female HR")
plot(gam_raw_qd_black_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black male HR")
plot(gam_raw_qd_hispanic_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic female HR")
plot(gam_raw_qd_hispanic_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic male HR")
plot(gam_raw_qd_asian_female, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian female HR")
plot(gam_raw_qd_asian_male, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian male HR")

par(mfrow=c(1,2))
plot(gam_raw_qd, all.terms = FALSE, trans=exp, shift=0.05, ylab="All qd HR")
plot(gam_raw_rm, all.terms = FALSE, trans=exp, shift=0.05, ylab="All RM HR")


require(schoenberg)
test.data.qd<-function(x, gm){
  test<-data.frame(pm25_ensemble = seq(4, 15,length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50),
                   mean_bmi=rep(mean(x$mean_bmi), 50),
                   smoke_rate= rep(mean(x$smoke_rate), 50),
                   hispanic=rep(mean(x$hispanic), 50),
                   pct_blk= rep(mean(x$pct_blk), 50),
                   medhouseholdincome=rep(mean(x$medhouseholdincome), 50),
                   medianhousevalue= rep(mean(x$medianhousevalue), 50),
                   poverty= rep(mean(x$poverty), 50),
                   education = rep(mean(x$education), 50),
                   popdensity = rep(mean(x$popdensity), 50),
                   pct_owner_occ= rep(mean(x$pct_owner_occ), 50),
                   summer_tmmx= rep(mean(x$summer_tmmx), 50),
                   winter_tmmx= rep(mean(x$winter_tmmx), 50),
                   summer_rmax = rep(mean(x$summer_rmax), 50),
                   winter_rmax= rep(mean(x$winter_rmax), 50),
                   year= rep(levels(as.factor(x$year))[1], 50),
                   region= rep(levels(as.factor(x$region))[1], 50),
                   time_count= rep(1, 50)) # This will give the appropriate RR estimate we desire
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  dly<- pred.vals$fit
  se.dly<-pred.vals$se.fit
  dly_low<-dly-qnorm(0.975)*se.dly
  dly_high<-dly+qnorm(0.975)*se.dly
  
  #logRR
  test$hr<-gm$family$linkinv(dly)
  test$hr_upper<-gm$family$linkinv(dly_high)
  test$hr_lower<-gm$family$linkinv(dly_low)
  
  test$y<-(dly - dly[1])
  test$upper<-(dly_high - dly_high[1])
  test$lower<-(dly_low - dly_low[1])
  return(test)
}


dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
load(paste0(dir_data,"aggregate_data_qd_trimmed.RData"))
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
black_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
hispanic_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

#Plotting
testallqd<-test.data.qd(aggregate_data_qd, gam_raw_qd)
test_white_female_qd<-test.data.qd(white_female_qd , gam_raw_qd_white_female)
test_white_male_qd<-test.data.qd(white_male_qd, gam_raw_qd_white_male)
test_black_female_qd<-test.data.qd(black_female_qd, gam_raw_qd_black_female)
test_black_male_qd<-test.data.qd(black_male_qd, gam_raw_qd_black_male)
test_hispanic_female_qd<-test.data.qd(hispanic_female_qd, gam_raw_qd_hispanic_female)
test_hispanic_male_qd<-test.data.qd(hispanic_male_qd, gam_raw_qd_hispanic_male)
test_asian_female_qd<-test.data.qd(asian_female_qd, gam_raw_qd_asian_female)
test_asian_male_qd<-test.data.qd(asian_male_qd, gam_raw_qd_asian_male)


pall_qd<-ggplot(data=testallqd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="All HR")
p_white_female_qd<-ggplot(data=test_white_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White female HR")
p_white_male_qd<-ggplot(data=test_white_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White male HR")
p_black_female_qd<-ggplot(data=test_black_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black female HR")
p_black_male_qd<-ggplot(data=test_black_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black male HR")
p_hispanic_female_qd<-ggplot(data=test_hispanic_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic female HR")
p_hispanic_male_qd<-ggplot(data=test_hispanic_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic male HR")
p_asian_female_qd<-ggplot(data=test_asian_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian female HR")
p_asian_male_qd<-ggplot(data=test_asian_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
#  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian male HR")


plot_grid(p_white_female_qd, p_white_male_qd,
          p_black_female_qd, p_black_male_qd,
          p_hispanic_female_qd, p_hispanic_male_qd,
          p_asian_female_qd, p_asian_male_qd, ncol=2)



