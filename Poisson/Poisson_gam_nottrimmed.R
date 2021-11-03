# All types of statistical models implemented for the paper
library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
require(parallel)
require(dplyr)
#require(devtools)
require(mgcv)
require(gam)
require(KernSmooth)

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

#Randall Martin
# Cox-equvalent conditional gam_raw Regression

gam_raw_rm<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
                        as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      data=aggregate_data_rm,family=poisson(link="log"))

save(gam_raw_rm, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_all.RData")
rm(gam_raw_rm)
gc()


#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)


gam_raw_rm_white_female<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_white_female.RData")
rm(gam_raw_rm_white_female)
gc()


#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)


gam_raw_rm_white_male<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_white_male.RData")
rm(gam_raw_rm_white_male)
gc()


#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)


gam_raw_rm_black_female<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_black_female.RData")
rm(gam_raw_rm_black_female)
gc()

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)

gam_raw_rm_black_male<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_black_male.RData")
rm(gam_raw_rm_black_male)
gc()



#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)

gam_raw_rm_hispanic_female<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_hispanic_female.RData")
rm(gam_raw_rm_hispanic_female)
gc()


#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)


gam_raw_rm_hispanic_male<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_hispanic_male.RData")
rm(gam_raw_rm_hispanic_male)
gc()



#Asian female
asian_rm_female<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)

gam_raw_rm_asian_female<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_asian_female.RData")
rm(gam_raw_rm_asian_female)
gc()

#Asian male
asian_rm_male<-aggregate_data_rm %>%
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)


gam_raw_rm_asian_male<-mgcv::bam(dead~  s(pm25,bs='cr',  k=3) +
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_rm_asian_male.RData")
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

gam_raw_qd<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                        as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                        mean_bmi + smoke_rate + hispanic+ pct_blk +
                        medhouseholdincome + medianhousevalue +
                        poverty + education + popdensity + pct_owner_occ +
                        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                        as.factor(year) + as.factor(region)
                      +offset(log(time_count)),
                      data=aggregate_data_qd,family=poisson(link="log"))

save(gam_raw_qd, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_all.RData")
rm(gam_raw_qd)
gc()


#White female
require(dplyr)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)


gam_raw_qd_white_female<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_white_female.RData")
rm(gam_raw_qd_white_female)
gc()


#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)


gam_raw_qd_white_male<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_white_male.RData")
rm(gam_raw_qd_white_male)
gc()


#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)


gam_raw_qd_black_female<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_black_female.RData")
rm(gam_raw_qd_black_female)
gc()

#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)

gam_raw_qd_black_male<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_black_male.RData")
rm(gam_raw_qd_black_male)
gc()



#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)

gam_raw_qd_hispanic_female<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_hispanic_female.RData")
rm(gam_raw_qd_hispanic_female)
gc()


#Hispanic Male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)


gam_raw_qd_hispanic_male<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_hispanic_male.RData")
rm(gam_raw_qd_hispanic_male)
gc()



#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)

gam_raw_qd_asian_female<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_asian_female.RData")
rm(gam_raw_qd_asian_female)
gc()

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

gam_raw_qd_asian_male<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
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
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/Main_gam_qd_asian_male.RData")
rm(gam_raw_qd_asian_male)
gc()

rm(aggregate_data_qd)
gc()




