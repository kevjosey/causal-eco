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
#pm<-aggregate_data_rm
#pm1<-pm %>% mutate(decile=ntile(pm25, 10))
#decile1<-subset(pm1, pm1$decile==1)
#decile10<-subset(pm1, pm1$decile==10)
#mean(decile1$pm25) #4.341009
#mean(decile10$pm25) #15.41103

aggregate_data_rm<-subset(aggregate_data_rm, aggregate_data_rm$pm25>=4.341009 & aggregate_data_rm$pm25<=15.41103)
#Randall Martin
# Cox-equvalent conditional gam_raw Regression
gam_raw_rm<-mgcv::bam(dead~  s(pm25, k=3) + 
                  as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                  mean_bmi + smoke_rate + hispanic+ pct_blk +
                  medhouseholdincome + medianhousevalue +
                  poverty + education + popdensity + pct_owner_occ +
                  summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                  as.factor(year) + as.factor(region)
                +offset(log(time_count)),
                data=aggregate_data_rm,family=poisson(link="log"))

save(gam_raw_rm, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_all.RData")
rm(gam_raw_rm)
gc()


#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_white_female.RData")
rm(gam_raw_rm_white_female)
gc()


#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_white_male.RData")
rm(gam_raw_rm_white_male)
gc()


#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_black_female.RData")
rm(gam_raw_rm_black_female)
gc()

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_black_male.RData")
rm(gam_raw_rm_black_male)
gc()



#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_hispanic_female.RData")
rm(gam_raw_rm_hispanic_female)
gc()


#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_hispanic_male.RData")
rm(gam_raw_rm_hispanic_male)
gc()



#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_asian_female.RData")
rm(gam_raw_rm_asian_female)
gc()

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_asian_male.RData")
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
#pm1<-aggregate_data_qd %>% mutate(decile=ntile(pm25_ensemble, 10))
#decile1<-subset(pm1, pm1$decile==1)
#decile10<-subset(pm1, pm1$decile==10)
#mean(decile1$pm25_ensemble) #4.581972
#mean(decile10$pm25_ensemble) #15.34556
aggregate_data_qd<-subset(aggregate_data_qd, aggregate_data_qd$pm25_ensemble>=4.581972 &
                            aggregate_data_qd$pm25_ensemble<=15.34556)

gam_raw_qd<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                        as.factor(race)+ as.factor(sex)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                data=aggregate_data_qd,family=poisson(link="log"))

save(gam_raw_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_all.RData")
rm(gam_raw_qd)
gc()



#White Female
white_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
gam_raw_qd_white_female<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                     as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                     mean_bmi + smoke_rate + hispanic+ pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(year) + as.factor(region)
                                   +offset(log(time_count)),
                             data=white_qd_female,family=poisson(link="log"))
rm(white_qd_female)
gc()
save(gam_raw_qd_white_female,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_white_female.RData")
rm(gam_raw_qd_white_female)
gc()

#White male
white_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
gam_raw_qd_white_male<-mgcv::bam(dead~  s(pm25_ensemble, k=3) + 
                                   as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                   mean_bmi + smoke_rate + hispanic+ pct_blk +
                                   medhouseholdincome + medianhousevalue +
                                   poverty + education + popdensity + pct_owner_occ +
                                   summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                   as.factor(year) + as.factor(region)
                                 +offset(log(time_count)),
                           data=white_qd_male,family=poisson(link="log"))
rm(white_qd_male)
gc()
save(gam_raw_qd_white_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_white_male.RData")
rm(gam_raw_qd_white_male)
gc()


#Black Female
black_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_black_female.RData")
rm(gam_raw_qd_black_female)
gc()

#Black male
black_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_black_male.RData")
rm(gam_raw_qd_black_male)
gc()


#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_hispanic_female.RData")
rm(gam_raw_qd_hispanic_female)
gc()

#Hispanic male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
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
save(gam_raw_qd_hispanic_male,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_hispanic_male.RData")
rm(gam_raw_qd_hispanic_male)
gc()



#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_asian_female.RData")
rm(gam_raw_qd_asian_female)
gc()

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_asian_male.RData")
rm(gam_raw_qd_asian_male)
gc()

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_all.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_white_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_white_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_black_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_black_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_hispanic_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_hispanic_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_asian_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_rm_asian_male.RData")

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

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_all.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_white_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_white_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_black_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_black_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_hispanic_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_hispanic_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_asian_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/GAM_decile/Main_gam_qd_asian_male.RData")

plot(gam_raw_qd, all.terms = FALSE, trans=exp, shift=0.05, ylab="All qd HR")
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

