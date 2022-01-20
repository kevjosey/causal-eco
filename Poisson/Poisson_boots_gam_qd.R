library("parallel")
library("dplyr")
require(parallel)
library(data.table)
require(tidyr)
require(doParallel)
library("survival")
library("gnm")
library("parallel")
require(doParallel)
library(data.table)
require(parallel)
require(dplyr)
require(mgcv)
require(zipcode)
require(ggplot2)


#Load Poisson model
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_qd/'

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,
                         by=c("zip","year"),all.x=T)

rm(covariates_qd, GPS_mod_qd)


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
process <- c(0:8)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]

# qd

if(process==1){
  print("All")
#all
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(aggregate_data_qd, list(aggregate_data_qd$zip))
num_uniq_zip <- length(unique(aggregate_data_qd$zip))


test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}


save(test,file=paste0(dir_out,"gam1_boots_qd_all.RData"))
rm(test, gam_raw)
gc()
stopCluster(cl)


#Race
} else if(process==2){
  print("White female")
#White female
require(dplyr)
white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_female_qd, list(white_female_qd$zip))
num_uniq_zip <- length(unique(white_female_qd$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw <-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_white_female.RData"))
rm(gam_raw, test)
gc()


}else if(process==3){
  print("White male")
#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(white_male_qd, list(white_male_qd$zip))
num_uniq_zip <- length(unique(white_male_qd$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_white_male.RData"))
gc()


#Black
} else if(process==4){
  print("Black female")
#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_qd_female, list(black_qd_female$zip))
num_uniq_zip <- length(unique(black_qd_female$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_black_female.RData"))
rm(test)
gc()


}else if(process==5){
#Black male
  print("Black male")
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(black_qd_male, list(black_qd_male$zip))
num_uniq_zip <- length(unique(black_qd_male$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)), 
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)


save(test,file=paste0(dir_out,"gam1_boots_qd_black_male.RData"))
rm(test)
gc()


} else if(process==6){
#Hispanic Female
  print("Hispanic female")
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_qd_female, list(hispanic_qd_female$zip))
num_uniq_zip <- length(unique(hispanic_qd_female$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~  s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}


save(test,file=paste0(dir_out,"gam1_boots_qd_hispanic_female.RData"))
gc()

stopCluster(cl)

} else if(process==7){
#Hispanic Male
  print("Hispanic male")
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(hispanic_qd_male, list(hispanic_qd_male$zip))
num_uniq_zip <- length(unique(hispanic_qd_male$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_hispanic_male.RData"))
gc()

} else if(process==8){
  print("Asian female")
#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_qd_female, list(asian_qd_female$zip))
num_uniq_zip <- length(unique(asian_qd_female$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gam_raw<-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',  k=3) + 
                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                       mean_bmi + smoke_rate + hispanic+ pct_blk +
                       medhouseholdincome + medianhousevalue +
                       poverty + education + popdensity + pct_owner_occ +
                       summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                       as.factor(year) + as.factor(region)
                     +offset(log(time_count)),
               data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_asian_female.RData"))
gc()

} else if(process==9){
#Asian male
  print("Asian male")
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
cl=makeCluster(20,outfile='')
registerDoParallel(cl)
aggregate_data.list<-split(asian_qd_male, list(asian_qd_male$zip))
num_uniq_zip <- length(unique(asian_qd_male$zip))

test<-data.frame(matrix(nrow=50, ncol=0))

for (boots_id in 1:500){
  set.seed(boots_id)
  print(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  #try({
    gam_raw<-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',  k=3) + 
                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                            mean_bmi + smoke_rate + hispanic+ pct_blk +
                            medhouseholdincome + medianhousevalue +
                            poverty + education + popdensity + pct_owner_occ +
                            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                            as.factor(year) + as.factor(region)
                          +offset(log(time_count)),
                          data=aggregate_data_boots,family=poisson(link="log"))
  
  test<-cbind(test, test.data.qd(aggregate_data_boots, gam_raw))
#}, silent=TRUE)
  print(boots_id)
  rm(aggregate_data_boots)
}
stopCluster(cl)

save(test,file=paste0(dir_out,"gam1_boots_qd_asian_male.RData"))
gc()
}
