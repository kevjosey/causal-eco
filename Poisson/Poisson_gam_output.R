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

require(matrixStats)

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
hispanic_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
asian_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
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
                   mean_bmi=rep(mean(x$mean_bmi), 50),
                   smoke_rate=rep(mean(x$smoke_rate), 50),
                   hispanic=rep(mean(x$hispanic), 50),
                   pct_blk=rep(mean(x$pct_blk), 50),
                   time_count=rep(1, 50),
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
  return(gm$family$linkinv(dly))
  
}

dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Data/'
dir_out_main='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/Poisson_qd/'

#Not trimmed
#Main
load(paste0(dir_our_main, "Main_gam_qd_all.RData"))
main_all_qd<-test.data.qd(aggregate_data, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_white_female.RData"))
main_white_female_qd<-test.data.qd(white_female_qd, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_white_male.RData"))
main_white_male-qd<-test.data.qd(white_male_qd, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_black_female.RData"))
main_black_female_qd<-test.data.qd(black_qd_female, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_black_male.RData"))
main_black_male_qd<-test.data.qd(black_qd_male, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_hispanic_female.RData"))
main_hispanic_female_qd<-test.data.qd(hispanic_qd_female, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_hispanic_male.RData"))
main_hispanic_male_qd<-test.data.qd(hispanic_qd_male, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_asian_female.RData"))
main_asian_female_qd <-test.data.qd(asian_qd_female, gam_raw)

load(paste0(dir_our_main, "Main_gam_qd_asian_male.RData"))
main_asian_male_qd <-test.data.qd(asian_qd_male, gam_raw)

#Boostrap
#Main
load(paste0(dir_out,"gam1_boots_qd_all.RData"))
num_uniq_zip <- length(unique(aggregate_data_qd$zip))
main_all_qd$main<-exp(main_all_qd$main)
main_all_qd$se<- rowSds(test)
main_all_qd$ci_low<-exp(main_all_qd$main - 1.96*(main_all_qd$se*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
main_all_qd$ci_high<-exp(main_all_qd$main + 1.96*(main_all_qd$se*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


load(paste0(dir_out,"gam1_boots_qd_white_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_white_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_hispanic_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_hispanic_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_male.RData"))

