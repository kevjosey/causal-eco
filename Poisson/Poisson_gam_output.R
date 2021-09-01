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
require(ggplot2)
require(matrixStats)

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
black_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
hispanic_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

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

dir_out_main = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Trimmed/'
dir_out='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_qd/'

#Not trimmed
#Main
load(paste0(dir_out_main, "Main_gam_qd_all.RData"))
main_all_qd<-data.frame(test.data.qd(aggregate_data_qd, gam_raw_qd))
colnames(main_all_qd)<-"main"

#Whitefemale
load(paste0(dir_out_main, "Main_gam_qd_white_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_white_female.RData"))

main_white_female_qd<-data.frame(test.data.qd(white_female_qd, gam_raw_qd_white_female))
colnames(main_white_female_qd)<-"main"
main_white_female_qd$a.vals<-seq(min(white_female_qd$pm25_ensemble), max(white_female_qd$pm25_ensemble),length.out=50)

main_white_female_qd$main<-exp(main_white_female_qd$main)
main_white_female_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_white_female_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_white_female_qd<-ggplot(main_white_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("White Female HR")
p_white_female_qd

load(paste0(dir_out_main, "Main_gam_qd_white_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_white_male.RData"))

main_white_male_qd<-data.frame(test.data.qd(white_male_qd, gam_raw_qd_white_male))
colnames(main_white_male_qd)<-"main"
main_white_male_qd$a.vals<-seq(min(white_male_qd$pm25_ensemble), max(white_male_qd$pm25_ensemble),length.out=50)

main_white_male_qd$main<-exp(main_white_male_qd$main)
main_white_male_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_white_male_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_white_male_qd<-ggplot(main_white_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("White Male HR")
p_white_male_qd

load(paste0(dir_out_main, "Main_gam_qd_black_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_all.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_female.RData"))

main_black_female_qd<-data.frame(test.data.qd(black_female_qd, gam_raw_qd_black_female))
colnames(main_black_female_qd)<-"main"
main_black_female_qd$a.vals<-seq(min(black_female_qd$pm25_ensemble), max(black_female_qd$pm25_ensemble),length.out=50)

main_black_female_qd$main<-exp(main_black_female_qd$main)
main_black_female_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_black_female_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_black_female_qd<-ggplot(main_black_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Black Female HR")
p_black_female_qd

load(paste0(dir_out_main, "Main_gam_qd_black_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_male.RData"))

main_black_male_qd<-data.frame(test.data.qd(black_male_qd, gam_raw_qd_black_male))
colnames(main_black_male_qd)<-"main"
main_black_male_qd$a.vals<-seq(min(black_male_qd$pm25_ensemble), max(black_male_qd$pm25_ensemble),length.out=50)

main_black_male_qd$main<-exp(main_black_male_qd$main)
main_black_male_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_black_male_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_black_male_qd<-ggplot(main_black_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Black Male HR")
p_black_male_qd

load(paste0(dir_out_main, "Main_gam_qd_hispanic_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_hispanic_female.RData"))

main_hispanic_female_qd<-data.frame(test.data.qd(hispanic_female_qd, gam_raw_qd_hispanic_female))
colnames(main_hispanic_female_qd)<-"main"
main_hispanic_female_qd$a.vals<-seq(min(hispanic_female_qd$pm25_ensemble), max(hispanic_female_qd$pm25_ensemble),length.out=50)

main_hispanic_female_qd$main<-exp(main_hispanic_female_qd$main)
main_hispanic_female_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_hispanic_female_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_hispanic_female_qd<-ggplot(main_hispanic_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Female HR")
p_hispanic_female_qd

load(paste0(dir_out_main, "Main_gam_qd_hispanic_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_hispanic_male.RData"))

main_hispanic_male_qd<-data.frame(test.data.qd(hispanic_male_qd, gam_raw_qd_hispanic_male))
colnames(main_hispanic_male_qd)<-"main"
main_hispanic_male_qd$a.vals<-seq(min(hispanic_male_qd$pm25_ensemble), max(hispanic_male_qd$pm25_ensemble),length.out=50)

main_hispanic_male_qd$main<-exp(main_hispanic_male_qd$main)
main_hispanic_male_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_hispanic_male_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_hispanic_male_qd<-ggplot(main_hispanic_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Male HR")
p_hispanic_male_qd

load(paste0(dir_out_main, "Main_gam_qd_asian_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_female.RData"))

main_asian_female_qd<-data.frame(test.data.qd(asian_female_qd, gam_raw_qd_asian_female))
colnames(main_asian_female_qd)<-"main"
main_asian_female_qd$a.vals<-seq(min(asian_female_qd$pm25_ensemble), max(asian_female_qd$pm25_ensemble),length.out=50)

main_asian_female_qd$main<-exp(main_asian_female_qd$main)
main_asian_female_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_asian_female_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_asian_female_qd<-ggplot(main_asian_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Asian Female HR")
p_asian_female_qd

load(paste0(dir_out_main, "Main_gam_qd_asian_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_male.RData"))

main_asian_male_qd<-data.frame(test.data.qd(asian_male_qd, gam_raw_qd_asian_male))
colnames(main_asian_male_qd)<-"main"
main_asian_male_qd$a.vals<-seq(min(asian_male_qd$pm25_ensemble), max(asian_male_qd$pm25_ensemble),length.out=50)

main_asian_male_qd$main<-exp(main_asian_male_qd$main)
main_asian_male_qd$cilow<- exp(rowQuantiles(as.matrix(test), prob=0.025))
main_asian_male_qd$cihigh<-exp(rowQuantiles(as.matrix(test), prob=0.975))

p_asian_male_qd<-ggplot(main_asian_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=cilow, ymax=cihigh), alpha=0.1)+xlab("PM2.5")+ylab("Asian Male HR")
p_asian_male_qd

require(cowplot)
plot_grid(p_white_female_qd, p_white_male_qd,
          p_black_female_qd, p_black_male_qd,
          p_hispanic_female_qd, p_hispanic_male_qd,
          p_asian_female_qd)
