library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
require(doParallel)
library(data.table)
library(fst)
library("parallel")
require(xgboost)

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd, by=c("zip","year"),all.x=T)
rm(covariates_qd, covariates_rm, GPS_qd, GPS_rm)

#Randall Martin
#White female
require(dplyr)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
#Asian female
asian_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
#Asian male
asian_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)


#QD
#White Female
white_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
#White male
white_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
#Black Female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
#Hispanic male
hispanic_qd_male<-aggregate_data_qd %>%  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
#Asian female
asian_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
#Asian male
asian_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

test.data.qd<-function(x, gm){
  test<-data.frame(pm25_ensemble = seq(min(x$pm25_ensemble), max(x$pm25_ensemble),length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50),
                   year= rep(levels(as.factor(x$year))[1], 50),
                   region= rep(levels(as.factor(x$region))[1], 50),
                   time_count= rep(mean(x$time_count),50),
                   mean_bmi=rep(mean(x$mean_bmi), 50),
                   smoke_rate=rep(mean(x$smoke_rate), 50),
                   hispanic=rep(mean(x$hispanic), 50),
                   pct_blk=rep(mean(x$pct_blk), 50),
                   time_count=rep(mean(x$time_count), 50),
                   medhouseholdincome=rep(mean(x$medhouseholdincome), 50),
                   medianhousevalue=rep(mean(x$medianhousevalue), 50),
                   poverty=rep(mean(x$poverty), 50),
                   education=rep(mean(x$education), 50),
                   popdensity=rep(mean(x$popdensity), 50),
                   pct_owner_occ=rep(mean(x$pct_owner_occ), 50),
                   summer_tmmx=rep(mean(x$summer_tmmx), 50),
                   winter_tmmx=rep(mean(x$winter_tmmx), 50),
                   summer_rmax=rep(mean(x$summer_rmax), 50),
                   winter_rmax=rep(mean(x$winter_rmax), 50) 
  )
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  dly<- pred.vals$fit
  se.dly<-pred.vals$se.fit
  
  dly_low<-dly-qnorm(0.975)*se.dly
  dly_high<-dly+qnorm(0.975)*se.dly
  #logRR
  test$y<-gm$family$linkinv(dly)
  test$upper<-gm$family$linkinv(dly_high)
  test$lower<-gm$family$linkinv(dly_low)
  
  return(test)
}


test.data.rm<-function(x, gm){
  test<-data.frame(pm25 = seq(min(x$pm25), max(x$pm25),length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50),
                   year= rep(levels(as.factor(x$year))[1], 50),
                   region= rep(levels(as.factor(x$region))[1], 50),
                   time_count= rep(mean(x$time_count),50),
                   mean_bmi=rep(mean(x$mean_bmi), 50),
                   smoke_rate=rep(mean(x$smoke_rate), 50),
                   hispanic=rep(mean(x$hispanic), 50),
                   pct_blk=rep(mean(x$pct_blk), 50),
                   time_count=rep(mean(x$time_count), 50),
                   medhouseholdincome=rep(mean(x$medhouseholdincome), 50),
                   medianhousevalue=rep(mean(x$medianhousevalue), 50),
                   poverty=rep(mean(x$poverty), 50),
                   education=rep(mean(x$education), 50),
                   popdensity=rep(mean(x$popdensity), 50),
                   pct_owner_occ=rep(mean(x$pct_owner_occ), 50),
                   summer_tmmx=rep(mean(x$summer_tmmx), 50),
                   winter_tmmx=rep(mean(x$winter_tmmx), 50),
                   summer_rmax=rep(mean(x$summer_rmax), 50),
                   winter_rmax=rep(mean(x$winter_rmax), 50) 
  )
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  dly<- pred.vals$fit
  se.dly<-pred.vals$se.fit
  
  dly_low<-dly-qnorm(0.975)*se.dly
  dly_high<-dly+qnorm(0.975)*se.dly
  #logRR
  test$y<-gm$family$linkinv(dly)
  test$upper<-gm$family$linkinv(dly_high)
  test$lower<-gm$family$linkinv(dly_low)
  return(test)
}


dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/Poisson_qd/'

load(paste0(dir_out,"gam1_boots_qd_white_female.RData"))
testqdwhite_female<-test
load(paste0(dir_out,"gam1_boots_qd_white_male.RData"))
testqdwhite_male<-test
load(paste0(dir_out,"gam1_boots_qd_black_female.RData"))
testqdblack_female<-test
load(paste0(dir_out,"gam1_boots_qd_black_male.RData"))
testqdblack_male<-test
load(paste0(dir_out,"gam1_boots_qd_hispanic_female.RData"))
testqdhispanic_female<-test
load(paste0(dir_out,"gam1_boots_qd_hispanic_male.RData"))
testqdhispanic_male<-test
load(paste0(dir_out,"gam1_boots_qd_asian_female.RData"))
testqdasian_female<-test
load(paste0(dir_out,"gam1_boots_qd_asian_male.RData"))
testqdasian_male<-test
#Asian male needs to be redone


dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/Poisson_rm/'

load(paste0(dir_out,"gam1_boots_rm_white_female.RData"))
testrmwhite_female<-test
load(paste0(dir_out,"gam1_boots_rm_white_male.RData"))
testrmwhite_male<-test
load(paste0(dir_out,"gam1_boots_rm_black_female.RData"))
testrmblack_female<-test
load(paste0(dir_out,"gam1_boots_rm_black_male.RData"))
testrmblack_male<-test
load(paste0(dir_out,"gam1_boots_rm_hispanic_female.RData"))
testrmhispanic_female<-test
load(paste0(dir_out,"gam1_boots_rm_hispanic_male.RData"))
testrmhispanic_male<-test
load(paste0(dir_out,"gam1_boots_rm_asian_female.RData"))
testrmasian_female<-test
load(paste0(dir_out,"gam1_boots_rm_asian_male.RData"))
testrmasian_male<-test
#Asian male needs to be redone

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_all.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_white_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_black_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_hispanic_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_rm_asian_male.RData")

plot(gam_raw_rm, all.terms = FALSE, trans=exp, shift=0.05, ylab="All RM HR")

#RM
testallrm<-test.data.rm(aggregate_data_rm, gam_raw_rm)
test_white_female_rm<-test.data.rm(white_female_rm, gam_raw_rm_white_female)
test_white_male_rm<-test.data.rm(white_male_rm, gam_raw_rm_white_male)
test_black_female_rm<-test.data.rm(black_rm_female, gam_raw_rm_black_female)
test_black_male_rm<-test.data.rm(black_rm_male, gam_raw_rm_black_male)
test_hispanic_female_rm<-test.data.rm(hispanic_rm_female, gam_raw_rm_hispanic_female)
test_hispanic_male_rm<-test.data.rm(hispanic_rm_male, gam_raw_rm_hispanic_male)
test_asian_female_rm<-test.data.rm(asian_rm_female, gam_raw_rm_asian_female)
test_asian_male_rm<-test.data.rm(asian_rm_male, gam_raw_rm_asian_male)

require(ggplot2)
require(cowplot)
require(ggExtra)
pall_rm<-ggplot(data=testallrm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="All HR")
palll_rm1<-ggMarginal(pall_rm, type="histogram")
p_white_female_rm<-ggplot(data=test_white_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="White female HR")
p_white_female_rm1<-ggMarginal(p_white_female_rm, type="histogram")
p_white_male_rm<-ggplot(data=test_white_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="White male HR")
p_white_male_rm1<-ggMarginal(p_white_male_rm, type="histogram")
p_black_female_rm<-ggplot(data=test_black_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Black female HR")
p_black_female_rm1<-ggMarginal(p_black_female_rm, type="histogram")
p_black_male_rm<-ggplot(data=test_black_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Black male HR")
p_black_male_rm1<-ggMarginal(p_black_male_rm, type="histogram")
p_hispanic_female_rm<-ggplot(data=test_hispanic_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Hispanic female HR")
p_hispanic_female_rm1<-ggMarginal(p_hispanic_female_rm, type="histogram")
p_hispanic_male_rm<-ggplot(data=test_hispanic_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Hispanic male HR")
p_hispanic_male_rm1<-ggMarginal(p_hispanic_male_rm, type="histogram")
p_asian_female_rm<-ggplot(data=test_asian_female_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Asian female HR")
p_asian_female_rm1<-ggMarginal(p_asian_female_rm, type="histogram")
p_asian_male_rm<-ggplot(data=test_asian_male_rm, mapping=aes(x=pm25))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="rm PM2.5", y="Asian male HR")
p_asian_male_rm1<-ggMarginal(p_asian_male_rm, type="histogram")

p_hist_all<-ggplot(data=aggregate_data_rm, aes(x=pm25))+geom_histogram()+
  labs(x="All rm PM2.5", y="Freq")
p_white_female_rm2<-ggplot(data=white_female_rm, aes(x=pm25))+geom_histogram()+
  labs(x="White female rm PM2.5", y="Freq")
p_white_male_rm2<-ggplot(data=white_male_rm, aes(x=pm25))+geom_histogram()+
  labs(x="White male rm PM2.5", y="Freq")
p_black_female_rm2<-ggplot(data=black_rm_female, aes(x=pm25))+geom_histogram()+
  labs(x="Black female rm PM2.5", y="Freq")
p_black_male_rm2<-ggplot(data=black_rm_male, aes(x=pm25))+geom_histogram()+
  labs(x="Black male rm PM2.5", y="Freq")
p_hispanic_female_rm2<-ggplot(data=hispanic_rm_female, aes(x=pm25))+geom_histogram()+
  labs(x="Hispanic female rm PM2.5", y="Freq")
p_hispanic_male_rm2<-ggplot(data=hispanic_rm_male, aes(x=pm25))+geom_histogram()+
  labs(x="Hispanic male rm PM2.5", y="Freq")
p_asian_female_rm2<-ggplot(data=asian_rm_female, aes(x=pm25))+geom_histogram()+
  labs(x="Asian female rm PM2.5", y="Freq")
p_asian_male_rm2<-ggplot(data=asian_rm_female, aes(x=pm25))+geom_histogram()+
  labs(x="Asian male rm PM2.5", y="Freq")

plot_grid(
  p_white_female_rm,p_white_female_rm2,
  p_white_male_rm,p_white_male_rm2, ncol=2)

plot_grid(
  p_black_female_rm, p_black_female_rm2,
  p_black_male_rm, p_black_male_rm2, ncol=2)

plot_grid(
  p_hispanic_female_rm, p_hispanic_female_rm2,
  p_hispanic_male_rm, p_hispanic_male_rm2, ncol=2)

plot_grid(
  p_asian_female_rm, p_asian_female_rm2,
  p_asian_male_rm,p_asian_male_rm2, ncol=2)

plot_grid(p_white_female_rm, p_white_male_rm,
          p_black_female_rm, p_black_male_rm,
          p_hispanic_female_rm, p_hispanic_male_rm,
          p_asian_female_rm, p_asian_male_rm, ncol=2)


load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_all.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_white_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_black_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_hispanic_male.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_female.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam_qd_asian_male.RData")


#QD
testallqd<-test.data.qd(aggregate_data_qd, gam_raw_qd)
test_white_female_qd<-test.data.qd(white_qd_female, gam_raw_qd_white_female)
test_white_male_qd<-test.data.qd(white_qd_male, gam_raw_qd_white_male)
test_black_female_qd<-test.data.qd(black_qd_female, gam_raw_qd_black_female)
test_black_male_qd<-test.data.qd(black_qd_male, gam_raw_qd_black_male)
test_hispanic_female_qd<-test.data.qd(hispanic_qd_female, gam_raw_qd_hispanic_female)
test_hispanic_male_qd<-test.data.qd(hispanic_qd_male, gam_raw_qd_hispanic_male)
test_asian_female_qd<-test.data.qd(asian_qd_female, gam_raw_qd_asian_female)
test_asian_male_qd<-test.data.qd(asian_qd_male, gam_raw_qd_asian_male)

require(ggplot2)
require(cowplot)
require(ggExtra)
pall_qd<-ggplot(data=testallqd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="All HR")
palll_qd1<-ggMarginal(pall_qd, type="histogram")
p_white_female_qd<-ggplot(data=test_white_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White female HR")
p_white_female_qd1<-ggMarginal(p_white_female_qd, type="histogram")
p_white_male_qd<-ggplot(data=test_white_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White male HR")
p_white_male_qd1<-ggMarginal(p_white_male_qd, type="histogram")
p_black_female_qd<-ggplot(data=test_black_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black female HR")
p_black_female_qd1<-ggMarginal(p_black_female_qd, type="histogram")
p_black_male_qd<-ggplot(data=test_black_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black male HR")
p_black_male_qd1<-ggMarginal(p_black_male_qd, type="histogram")
p_hispanic_female_qd<-ggplot(data=test_hispanic_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic female HR")
p_hispanic_female_qd1<-ggMarginal(p_hispanic_female_qd, type="histogram")
p_hispanic_male_qd<-ggplot(data=test_hispanic_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic male HR")
p_hispanic_male_qd1<-ggMarginal(p_hispanic_male_qd, type="histogram")
p_asian_female_qd<-ggplot(data=test_asian_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian female HR")
p_asian_female_qd1<-ggMarginal(p_asian_female_qd, type="histogram")
p_asian_male_qd<-ggplot(data=test_asian_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian male HR")
p_asian_male_qd1<-ggMarginal(p_asian_male_qd, type="histogram")

p_hist_all<-ggplot(data=aggregate_data_qd, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="All QD PM2.5", y="Freq")
p_white_female_qd2<-ggplot(data=white_qd_female, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="White female QD PM2.5", y="Freq")
p_white_male_qd2<-ggplot(data=white_qd_male, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="White male QD PM2.5", y="Freq")
p_black_female_qd2<-ggplot(data=black_qd_female, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Black female QD PM2.5", y="Freq")
p_black_male_qd2<-ggplot(data=black_qd_male, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Black male QD PM2.5", y="Freq")
p_hispanic_female_qd2<-ggplot(data=hispanic_qd_female, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Hispanic female QD PM2.5", y="Freq")
p_hispanic_male_qd2<-ggplot(data=hispanic_qd_male, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Hispanic male QD PM2.5", y="Freq")
p_asian_female_qd2<-ggplot(data=asian_qd_female, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Asian female QD PM2.5", y="Freq")
p_asian_male_qd2<-ggplot(data=asian_qd_male, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Asian male QD PM2.5", y="Freq")

plot_grid(
  p_white_female_qd,p_white_female_qd2,
  p_white_male_qd,p_white_male_qd2, ncol=2)

plot_grid(
  p_black_female_qd, p_black_female_qd2,
  p_black_male_qd, p_black_male_qd2, ncol=2)

plot_grid(
  p_hispanic_female_qd, p_hispanic_female_qd2,
  p_hispanic_male_qd, p_hispanic_male_qd2, ncol=2)

plot_grid(
  p_asian_female_qd, p_asian_female_qd2,
  p_asian_male_qd,p_asian_male_qd2, ncol=2)

plot_grid(p_white_female_qd, p_white_male_qd,
          p_black_female_qd, p_black_male_qd,
          p_hispanic_female_qd, p_hispanic_male_qd,
          p_asian_female_qd, p_asian_male_qd, ncol=2)

