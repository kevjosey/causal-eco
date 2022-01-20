library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
require(cowplot)
require(ggplot2)
require(doParallel)
library(data.table)
library("parallel")
require(ggplot2)
require(matrixStats)

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)
white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_female_rm <- subset(white_female_rm,
                            pm25 < quantile(white_female_rm$pm25,0.95)&
                              pm25 > quantile(white_female_rm$pm25,0.05))

white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
white_male_rm <- subset(white_male_rm,
                          pm25 < quantile(white_male_rm$pm25,0.95)&
                            pm25 > quantile(white_male_rm$pm25,0.05))

black_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_female_rm <- subset(black_female_rm,
                          pm25 < quantile(black_female_rm$pm25,0.95)&
                            pm25 > quantile(black_female_rm$pm25,0.05))

black_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
black_male_rm <- subset(black_male_rm,
                          pm25 < quantile(black_male_rm$pm25,0.95)&
                            pm25 > quantile(black_male_rm$pm25,0.05))

hispanic_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_female_rm <- subset(hispanic_female_rm,
                          pm25 < quantile(hispanic_female_rm$pm25,0.95)&
                            pm25 > quantile(hispanic_female_rm$pm25,0.05))

hispanic_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
hispanic_male_rm <- subset(hispanic_male_rm,
                             pm25 < quantile(hispanic_male_rm$pm25,0.95)&
                               pm25 > quantile(hispanic_male_rm$pm25,0.05))

asian_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_female_rm <- subset(asian_female_rm,
                          pm25 < quantile(asian_female_rm$pm25,0.95)&
                            pm25 > quantile(asian_female_rm$pm25,0.05))

asian_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
asian_male_rm <- subset(asian_male_rm,
                          pm25 < quantile(asian_male_rm$pm25,0.95)&
                            pm25 > quantile(asian_male_rm$pm25,0.05))

aggregate_data_rm <- subset(aggregate_data_rm,
                        pm25 < quantile(aggregate_data_rm$pm25,0.95)&
                          pm25 > quantile(aggregate_data_rm$pm25,0.05))

test.data.rm<-function(x, gm){
  test<-data.frame(pm25 = seq(min(x$pm25), max(x$pm25),length.out=50) ,
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

base<-data.frame(pm25 = rep(min(x$pm25),50) ,
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

dir_out_main = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/GAM/Trimmed/'
dir_out='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_rm/'

#Trimmed
#Main
load(paste0(dir_out_main, "Main_gam_rm_all_trimmed.RData"))
load(paste0(dir_out, "gam1_boots_rm_all_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

aggregate_data_rm <- subset(aggregate_data_rm,
                            pm25 < quantile(aggregate_data_rm$pm25,0.95)&
                              pm25 > quantile(aggregate_data_rm$pm25,0.05))
num_uniq_zip <- length(unique(aggregate_data_rm$zip))

main_all_rm<-data.frame(test.data.rm(aggregate_data_rm, gam_raw_rm))
colnames(main_all_rm)<-"main"
main_all_rm$a.vals<-seq(min(aggregate_data_rm$pm25), max(aggregate_data_rm$pm25),length.out=50)
main_all_rm$main<-(main_all_rm$main)

std<-apply((test), 1, sd, na.rm=TRUE)
main_all_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_all_rm$lower <-(main_all_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_all_rm$higher<-(main_all_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_all_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_all_rm$higher<-(rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(aggregate_data_rm$pm25), max(aggregate_data_rm$pm25))
p_all_rm<-ggplot(main_all_rm, aes(x=a.vals)) +
  geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+
  xlab("PM2.5")+ylab("All HR")+
  theme(panel.grid=element_blank(), legend.position="none")+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_all_rm+theme(legend.position = "none")

p_all_hist<-ggplot(data=aggregate_data_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  coord_cartesian(xlim=range)+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_all_hist, p_all_rm+theme(legend.position="none"), align="hv", axis="tblr")
p_all_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Whitefemale
load(paste0(dir_out_main, "Main_gam_rm_white_female_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_white_female_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_female_rm <- subset(white_female_rm,
                          pm25 < quantile(white_female_rm$pm25,0.95)&
                            pm25 > quantile(white_female_rm$pm25,0.05))

num_uniq_zip <- length(unique(white_female_rm$zip))
main_white_female_rm<-data.frame(test.data.rm(white_female_rm, gam_raw_rm_white_female))
colnames(main_white_female_rm)<-"main"
main_white_female_rm$a.vals<-seq(min(white_female_rm$pm25), max(white_female_rm$pm25),length.out=50)
main_white_female_rm$main<-(main_white_female_rm$main)

std<-apply((test), 1, sd, na.rm=TRUE)
main_white_female_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_white_female_rm$lower <-(main_white_female_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_white_female_rm$higher<-(main_white_female_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

range=c(min(white_female_rm$pm25), max(white_female_rm$pm25))
p_white_female_rm<-ggplot(main_white_female_rm, aes(x=a.vals)) +
  geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("White Female HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_white_female_rm+theme(legend.position = "none")


p_white_female_hist<-ggplot(data=white_female_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_white_female_hist, p_white_female_rm+theme(legend.position="none"), 
                           align="hv", axis="tblr")
p_white_female_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#White Male
#White male
load(paste0(dir_out_main, "Main_gam_rm_white_male_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_white_male_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
white_male_rm <- subset(white_male_rm,
                        pm25 < quantile(white_male_rm$pm25,0.95)&
                          pm25 > quantile(white_male_rm$pm25,0.05))
num_uniq_zip <- length(unique(white_male_rm$zip))



main_white_male_rm<-data.frame(test.data.rm(white_male_rm, gam_raw_rm_white_male))
colnames(main_white_male_rm)<-"main"
main_white_male_rm$a.vals<-seq(min(white_male_rm$pm25), max(white_male_rm$pm25),length.out=50)

main_white_male_rm$main<-(main_white_male_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_white_male_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_white_male_rm$lower <-(main_white_male_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_white_male_rm$higher<-(main_white_male_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

range=c(min(white_male_rm$pm25), max(white_male_rm$pm25))
p_white_male_rm<-ggplot(main_white_male_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("White Male HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_white_male_rm+theme(legend.position = "none")


p_white_male_hist<-ggplot(data=white_male_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_white_male_hist, p_white_male_rm+theme(legend.position="none"), align="hv", axis="tblr")
p_white_male_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Black female

load(paste0(dir_out_main, "Main_gam_rm_black_female_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_all_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_black_female_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)


black_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_female_rm <- subset(black_female_rm,
                          pm25 < quantile(black_female_rm$pm25,0.95)&
                            pm25 > quantile(black_female_rm$pm25,0.05))

num_uniq_zip <- length(unique(black_female_rm$zip))
main_black_female_rm<-data.frame(test.data.rm(black_female_rm, gam_raw_rm_black_female))
colnames(main_black_female_rm)<-"main"
main_black_female_rm$a.vals<-seq(min(black_female_rm$pm25), max(black_female_rm$pm25),length.out=50)

main_black_female_rm$main<-(main_black_female_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_black_female_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_black_female_rm$lower <-(main_black_female_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_black_female_rm$higher<-(main_black_female_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_black_female_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_black_female_rm$higher<-(rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(black_female_rm$pm25), max(black_female_rm$pm25))
p_black_female_rm<-ggplot(main_black_female_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Black Female HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_black_female_rm+theme(legend.position = "none")


p_black_female_hist<-ggplot(data=black_female_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_black_female_hist, p_black_female_rm+theme(legend.position="none"), align="hv", axis="tblr")
p_black_female_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])


#Black male
load(paste0(dir_out_main, "Main_gam_rm_black_male_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_black_male_trimmed.RData"))
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

black_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
black_male_rm <- subset(black_male_rm,
                        pm25 < quantile(black_male_rm$pm25,0.95)&
                          pm25 > quantile(black_male_rm$pm25,0.05))

num_uniq_zip <- length(unique(black_male_rm$zip))
main_black_male_rm<-data.frame(test.data.rm(black_male_rm, gam_raw_rm_black_male))
colnames(main_black_male_rm)<-"main"
main_black_male_rm$a.vals<-seq(min(black_male_rm$pm25), max(black_male_rm$pm25),length.out=50)

main_black_male_rm$main<-(main_black_male_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_black_male_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_black_male_rm$lower <-(main_black_male_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_black_male_rm$higher<-(main_black_male_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_black_male_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_black_male_rm$higher<- (rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(black_male_rm$pm25), max(black_male_rm$pm25))
p_black_male_rm<-ggplot(main_black_male_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Black Male HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_black_male_rm+theme(legend.position = "none")

p_black_male_hist<-ggplot(data=black_male_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_black_male_hist, p_black_male_rm+theme(legend.position="none"), align="hv", axis="tblr")
p_black_male_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Hispanic female
load(paste0(dir_out_main, "Main_gam_rm_hispanic_female_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_hispanic_female_trimmed.RData"))
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

hispanic_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_female_rm <- subset(hispanic_female_rm,
                             pm25 < quantile(hispanic_female_rm$pm25,0.95)&
                               pm25 > quantile(hispanic_female_rm$pm25,0.05))

num_uniq_zip <- length(unique(hispanic_female_rm$zip))
main_hispanic_female_rm<-data.frame(test.data.rm(hispanic_female_rm, gam_raw_rm_hispanic_female))
colnames(main_hispanic_female_rm)<-"main"
main_hispanic_female_rm$a.vals<-seq(min(hispanic_female_rm$pm25), max(hispanic_female_rm$pm25),length.out=50)

main_hispanic_female_rm$main<-(main_hispanic_female_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_hispanic_female_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_hispanic_female_rm$lower <-(main_hispanic_female_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_hispanic_female_rm$higher<-(main_hispanic_female_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_hispanic_female_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_hispanic_female_rm$higher<- (rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(hispanic_female_rm$pm25), max(hispanic_female_rm$pm25))

p_hispanic_female_rm<-ggplot(main_hispanic_female_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Female HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_hispanic_female_rm+theme(legend.position = "none")


p_hispanic_female_hist<-ggplot(data=hispanic_female_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_hispanic_female_hist, p_hispanic_female_rm+theme(legend.position="none"), 
                           align="hv", axis="tblr")
p_hispanic_female_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Hispanic male
load(paste0(dir_out_main, "Main_gam_rm_hispanic_male_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_hispanic_male_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

hispanic_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
hispanic_male_rm <- subset(hispanic_male_rm,
                           pm25 < quantile(hispanic_male_rm$pm25,0.95)&
                             pm25 > quantile(hispanic_male_rm$pm25,0.05))

num_uniq_zip <- length(unique(hispanic_male_rm$zip))
main_hispanic_male_rm<-data.frame(test.data.rm(hispanic_male_rm, gam_raw_rm_hispanic_male))
colnames(main_hispanic_male_rm)<-"main"
main_hispanic_male_rm$a.vals<-seq(min(hispanic_male_rm$pm25), max(hispanic_male_rm$pm25),length.out=50)

main_hispanic_male_rm$main<-(main_hispanic_male_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_hispanic_male_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_hispanic_male_rm$lower <-(main_hispanic_male_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_hispanic_male_rm$higher<-(main_hispanic_male_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_hispanic_male_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_hispanic_male_rm$higher<-(rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(hispanic_male_rm$pm25), max(hispanic_male_rm$pm25))

p_hispanic_male_rm<-ggplot(main_hispanic_male_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Hispanic Male HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_hispanic_male_rm+theme(legend.position = "none")

p_hispanic_male_hist<-ggplot(data=hispanic_male_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_hispanic_male_hist, p_hispanic_male_rm+theme(legend.position="none"), 
                           align="hv", axis="tblr")
p_hispanic_male_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Asian female

load(paste0(dir_out_main, "Main_gam_rm_asian_female_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_asian_female_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)


asian_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_female_rm <- subset(asian_female_rm,
                          pm25 < quantile(asian_female_rm$pm25,0.95)&
                            pm25 > quantile(asian_female_rm$pm25,0.05))

num_uniq_zip <- length(unique(asian_female_rm$zip))

main_asian_female_rm<-data.frame(test.data.rm(asian_female_rm, gam_raw_rm_asian_female))
colnames(main_asian_female_rm)<-"main"
main_asian_female_rm$a.vals<-seq(min(asian_female_rm$pm25), max(asian_female_rm$pm25),length.out=50)

main_asian_female_rm$main<-(main_asian_female_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_asian_female_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_asian_female_rm$lower <-(main_asian_female_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_asian_female_rm$higher<-(main_asian_female_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_asian_female_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_asian_female_rm$higher<- (rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(asian_female_rm$pm25), max(asian_female_rm$pm25))
p_asian_female_rm<-ggplot(main_asian_female_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Asian Female HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_asian_female_rm+theme(legend.position = "none")

p_asian_female_hist<-ggplot(data=asian_female_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_asian_female_hist, p_asian_female_rm+theme(legend.position="none"), 
                           align="hv", axis="tblr")
p_asian_female_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

#Asian male
load(paste0(dir_out_main, "Main_gam_rm_asian_male_trimmed.RData"))
load(paste0(dir_out,"gam1_boots_rm_asian_male_trimmed.RData"))

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

asian_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
asian_male_rm <- subset(asian_male_rm,
                        pm25 < quantile(asian_male_rm$pm25,0.95)&
                          pm25 > quantile(asian_male_rm$pm25,0.05))
num_uniq_zip <- length(unique(asian_male_rm$zip))
main_asian_male_rm<-data.frame(test.data.rm(asian_male_rm, gam_raw_rm_asian_male))
colnames(main_asian_male_rm)<-"main"
main_asian_male_rm$a.vals<-seq(min(asian_male_rm$pm25), max(asian_male_rm$pm25),length.out=50)

main_asian_male_rm$main<-(main_asian_male_rm$main)
std<-apply((test), 1, sd, na.rm=TRUE)
main_asian_male_rm$median<-apply(test, 1, mean, na.rm=TRUE)
main_asian_male_rm$lower <-(main_asian_male_rm$main) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
main_asian_male_rm$higher<-(main_asian_male_rm$main)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#main_asian_male_rm$lower<- (rowQuantiles(as.matrix(test), prob=0.025))
#main_asian_male_rm$higher<-(rowQuantiles(as.matrix(test), prob=0.975))
range=c(min(asian_male_rm$pm25), max(asian_male_rm$pm25))

p_asian_male_rm<-ggplot(main_asian_male_rm, aes(x=a.vals)) + geom_line(aes( y = main, col="red"))+
  geom_ribbon(aes(ymin=lower, ymax=higher), alpha=0.1)+xlab("PM2.5")+ylab("Asian Male HR")+
  theme(panel.grid=element_blank(), legend.position="none")+ylim(0.5, 1.25)+
  coord_cartesian(xlim=range)+
  theme_cowplot()
p_asian_male_rm+theme(legend.position = "none")


p_asian_male_hist<-ggplot(data=asian_male_rm, aes(x=pm25))+
  geom_histogram(aes(y=..density..), color="grey", fill="white", bins=30, alpha=1)+
  coord_cartesian(xlim=range)+
  labs(x="", y="Density")+theme(panel.grid=element_blank())+
  scale_y_continuous(position="right")+theme_cowplot()

aligned_plots<-align_plots(p_asian_male_hist, p_asian_male_rm+theme(legend.position="none"), 
                           align="hv", axis="tblr")
p_asian_male_rm1<-ggdraw(aligned_plots[[1]])+
  draw_plot(aligned_plots[[2]])

require(cowplot)
plot_grid(p_white_female_rm + ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"),
          p_white_male_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"), 
          p_black_female_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"),
          p_black_male_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"), 
          p_hispanic_female_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"),
          p_hispanic_male_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"), 
          p_asian_female_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"),
          p_asian_male_rm+ ylim(0.5, 1.5)+theme_bw()+theme(legend.position="none"), 
          ncol=2, labels="AUTO")

plot_grid(p_white_female_rm +theme_bw() + theme(legend.position="none"),
          p_white_male_rm+theme_bw()+ theme(legend.position="none"), 
          p_black_female_rm+theme_bw()+ theme(legend.position="none"), 
          p_black_male_rm+theme_bw()+ theme(legend.position="none"), 
          p_hispanic_female_rm+theme_bw()+ theme(legend.position="none"),
          p_hispanic_male_rm+theme_bw()+ theme(legend.position="none"), 
          p_asian_female_rm+theme_bw()+ theme(legend.position="none"),
          p_asian_male_rm+theme_bw()+ theme(legend.position="none"), ncol=2, 
          labels="AUTO")

plot_grid(p_white_female_rm1, p_white_male_rm1,
          p_black_female_rm1, p_black_male_rm1,
          p_hispanic_female_rm1, p_hispanic_male_rm1,
          p_asian_female_rm1, p_asian_male_rm1, ncol=2)