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
library("parallel")
require(ggplot2)
require(matrixStats)

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

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
  
  base<-data.frame(pm25_ensemble = rep(min(x$pm25_ensemble),50) ,
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
  base.vals<-predict(gm, newdata=base, type="link", se.fit = TRUE)
  dly<- pred.vals$fit - base.vals$fit
  return(gm$family$linkinv(dly))
  
}

dir_out_main = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Untrimmed/'
dir_out='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_qd/'

#Untrimmed
#Main
#Main
load(paste0(dir_out_main, "Main_gam_qd_all.RData"))
load(paste0(dir_out, "gam1_boots_qd_all.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)
num_uniq_zip <- length(unique(aggregate_data_qd$zip))

main_all_qd<-data.frame(test.data.qd(aggregate_data_qd, gam_raw_qd))
colnames(main_all_qd)<-"main"
main_all_qd$a.vals<-seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble),length.out=50)
main_all_qd$main<-(main_all_qd$main)

std<-apply((test), 1, sd, na.rm=TRUE)
main_all_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_all_qd$lower <-(main_all_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_all_qd$higher<-(main_all_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_all_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_all_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))

p_all_qd<-ggplot(main_all_qd, aes(x=a.vals)) +
  geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+
  xlab("PM2.5")+ylab("All HR")+theme(legend.position="none")
p_all_qd

#Whitefemale
load(paste0(dir_out_main, "Main_gam_qd_white_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_white_female.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
num_uniq_zip <- length(unique(white_female_qd$zip))
main_white_female_qd<-data.frame(test.data.qd(white_female_qd, gam_raw_qd_white_female))
colnames(main_white_female_qd)<-"main"
main_white_female_qd$a.vals<-seq(min(white_female_qd$pm25_ensemble), max(white_female_qd$pm25_ensemble),length.out=50)
main_white_female_qd$main<-(main_white_female_qd$main)

std<-apply((test), 1, sd, na.rm=TRUE)
main_white_female_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_white_female_qd$lower <-(main_white_female_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_white_female_qd$higher<-(main_white_female_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_white_female_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_white_female_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))

p_white_female_qd<-ggplot(main_white_female_qd, aes(x=a.vals)) +
  geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("White Female HR")+theme(legend.position="none")
p_white_female_qd

#White Male
#White male
load(paste0(dir_out_main, "Main_gam_qd_white_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_white_male.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
num_uniq_zip <- length(unique(white_male_qd$zip))

main_white_male_qd<-data.frame(test.data.qd(white_male_qd, gam_raw_qd_white_male))
colnames(main_white_male_qd)<-"main"
main_white_male_qd$a.vals<-seq(min(white_male_qd$pm25_ensemble), max(white_male_qd$pm25_ensemble),length.out=50)

main_white_male_qd$main<-(main_white_male_qd$main)
#main_white_male_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_white_male_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))
std<-apply((test), 1, sd, na.rm=TRUE)
main_white_male_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_white_male_qd$lower <-(main_white_male_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_white_male_qd$higher<-(main_white_male_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)


p_white_male_qd<-ggplot(main_white_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("White Male HR")+theme(legend.position="none")
p_white_male_qd

#Black female

load(paste0(dir_out_main, "Main_gam_qd_black_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_all.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_female.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)


black_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
num_uniq_zip <- length(unique(black_female_qd$zip))
main_black_female_qd<-data.frame(test.data.qd(black_female_qd, gam_raw_qd_black_female))
colnames(main_black_female_qd)<-"main"
main_black_female_qd$a.vals<-seq(min(black_female_qd$pm25_ensemble), max(black_female_qd$pm25_ensemble),length.out=50)

main_black_female_qd$main<-(main_black_female_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_black_female_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_black_female_qd$lower <-(main_black_female_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_black_female_qd$higher<-(main_black_female_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_black_female_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_black_female_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))

p_black_female_qd<-ggplot(main_black_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Black Female HR")+theme(legend.position="none")
p_black_female_qd


#Black male
load(paste0(dir_out_main, "Main_gam_qd_black_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_black_male.RData"))
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

black_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
num_uniq_zip <- length(unique(black_male_qd$zip))
main_black_male_qd<-data.frame(test.data.qd(black_male_qd, gam_raw_qd_black_male))
colnames(main_black_male_qd)<-"main"
main_black_male_qd$a.vals<-seq(min(black_male_qd$pm25_ensemble), max(black_male_qd$pm25_ensemble),length.out=50)

main_black_male_qd$main<-(main_black_male_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_black_male_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_black_male_qd$lower <-(main_black_male_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_black_male_qd$higher<-(main_black_male_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_black_male_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_black_male_qd$higher<- (rowQuantiles(as.matrix(test), prob=0.975))

p_black_male_qd<-ggplot(main_black_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Black Male HR")+theme(legend.position="none")
p_black_male_qd

#Hispanic female
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

hispanic_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
num_uniq_zip <- length(unique(hispanic_female_qd$zip))
main_hispanic_female_qd<-data.frame(test.data.qd(hispanic_female_qd, gam_raw_qd_hispanic_female))
colnames(main_hispanic_female_qd)<-"main"
main_hispanic_female_qd$a.vals<-seq(min(hispanic_female_qd$pm25_ensemble), max(hispanic_female_qd$pm25_ensemble),length.out=50)

main_hispanic_female_qd$main<-(main_hispanic_female_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_hispanic_female_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_hispanic_female_qd$lower <-(main_hispanic_female_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_hispanic_female_qd$higher<-(main_hispanic_female_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_hispanic_female_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_hispanic_female_qd$higher<- (rowQuantiles(as.matrix(test), prob=0.975))

p_hispanic_female_qd<-ggplot(main_hispanic_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Female HR")+theme(legend.position="none")
p_hispanic_female_qd

#Hispanic male
load(paste0(dir_out_main, "Main_gam_qd_hispanic_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_hispanic_male.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

hispanic_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
num_uniq_zip <- length(unique(hispanic_male_qd$zip))
main_hispanic_male_qd<-data.frame(test.data.qd(hispanic_male_qd, gam_raw_qd_hispanic_male))
colnames(main_hispanic_male_qd)<-"main"
main_hispanic_male_qd$a.vals<-seq(min(hispanic_male_qd$pm25_ensemble), max(hispanic_male_qd$pm25_ensemble),length.out=50)

main_hispanic_male_qd$main<-(main_hispanic_male_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_hispanic_male_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_hispanic_male_qd$lower <-(main_hispanic_male_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_hispanic_male_qd$higher<-(main_hispanic_male_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_hispanic_male_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_hispanic_male_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))

p_hispanic_male_qd<-ggplot(main_hispanic_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Male HR")+theme(legend.position="none")
p_hispanic_male_qd

#Asian female

load(paste0(dir_out_main, "Main_gam_qd_asian_female.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_female.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)


asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
num_uniq_zip <- length(unique(asian_female_qd$zip))

main_asian_female_qd<-data.frame(test.data.qd(asian_female_qd, gam_raw_qd_asian_female))
colnames(main_asian_female_qd)<-"main"
main_asian_female_qd$a.vals<-seq(min(asian_female_qd$pm25_ensemble), max(asian_female_qd$pm25_ensemble),length.out=50)

main_asian_female_qd$main<-(main_asian_female_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_asian_female_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_asian_female_qd$lower <-(main_asian_female_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_asian_female_qd$higher<-(main_asian_female_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_asian_female_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_asian_female_qd$higher<- (rowQuantiles(as.matrix(test), prob=0.975))

p_asian_female_qd<-ggplot(main_asian_female_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Asian Female HR")+theme(legend.position="none")
p_asian_female_qd

#Asian male
load(paste0(dir_out_main, "Main_gam_qd_asian_male.RData"))
load(paste0(dir_out,"gam1_boots_qd_asian_male.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)

asian_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
num_uniq_zip <- length(unique(asian_male_qd$zip))
main_asian_male_qd<-data.frame(test.data.qd(asian_male_qd, gam_raw_qd_asian_male))
colnames(main_asian_male_qd)<-"main"
main_asian_male_qd$a.vals<-seq(min(asian_male_qd$pm25_ensemble), max(asian_male_qd$pm25_ensemble),length.out=50)

main_asian_male_qd$main<-(main_asian_male_qd$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_asian_male_qd$median<-apply(test, 1, mean, na.rm=TRUE)
main_asian_male_qd$lower <-(main_asian_male_qd$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_asian_male_qd$higher<-(main_asian_male_qd$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_asian_male_qd$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_asian_male_qd$higher<-(rowQuantiles(as.matrix(test), prob=0.975))

p_asian_male_qd<-ggplot(main_asian_male_qd, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Asian Male HR")+theme(legend.position="none")
p_asian_male_qd

require(cowplot)
plot_grid(p_white_female_qd + ylim(0.5, 1.5)+theme_bw(), p_white_male_qd+ ylim(0.5, 1.5)+theme_bw(), 
          p_black_female_qd+ ylim(0.5, 1.5)+theme_bw(), p_black_male_qd+ ylim(0.5, 1.5)+theme_bw(), 
          p_hispanic_female_qd+ ylim(0.5, 1.5)+theme_bw(), p_hispanic_male_qd+ ylim(0.5, 1.5)+theme_bw(), 
          p_asian_female_qd+ ylim(0.5, 1.5)+theme_bw(), p_asian_male_qd+ ylim(0.5, 1.5)+theme_bw(), ncol=2, 
          labels="AUTO")

plot_grid(p_white_female_qd +theme_bw(), p_white_male_qd+theme_bw(), 
          p_black_female_qd+theme_bw(), p_black_male_qd+theme_bw(), 
          p_hispanic_female_qd+theme_bw(), p_hispanic_male_qd+theme_bw(), 
          p_asian_female_qd+theme_bw(), p_asian_male_qd+theme_bw(), ncol=2, 
          labels="AUTO")
