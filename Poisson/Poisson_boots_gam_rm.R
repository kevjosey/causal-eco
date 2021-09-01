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
library(fst)
require(xgboost)
require(parallel)
require(dplyr)
require(mgcv)
require(zipcode)
require(ggplot2)


#Load Poisson model
dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/Poisson_rm/'

#load(file= "/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/Main_gam.RData")
load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,
                         by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm)

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
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  dly<- pred.vals$fit
  return(gm$family$linkinv(dly))
  
}

process <- c(0:8)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]

# rm

if(process==1){
  print("All")
  #all
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(aggregate_data_rm, list(aggregate_data_rm$zip))
  num_uniq_zip <- length(unique(aggregate_data_rm$zip))
  
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_all.RData"))
  rm(test, gam_raw)
  gc()
  stopCluster(cl)
  
  
  #Race
} else if(process==2){
  print("White female")
  #White female
  require(dplyr)
  white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(white_female_rm, list(white_female_rm$zip))
  num_uniq_zip <- length(unique(white_female_rm$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw <-mgcv::bam(dead~  s(pm25, k=3) + 
                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                          mean_bmi + smoke_rate + hispanic+ pct_blk +
                          medhouseholdincome + medianhousevalue +
                          poverty + education + popdensity + pct_owner_occ +
                          summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                          as.factor(year) + as.factor(region)
                        +offset(log(time_count)),
                        data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_white_female.RData"))
  rm(gam_raw, test)
  gc()
  
  
}else if(process==3){
  print("White male")
  #White Male
  white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(white_male_rm, list(white_male_rm$zip))
  num_uniq_zip <- length(unique(white_male_rm$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_white_male.RData"))
  gc()
  
  
  #Black
} else if(process==4){
  print("Black female")
  #Black female
  black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(black_rm_female, list(black_rm_female$zip))
  num_uniq_zip <- length(unique(black_rm_female$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_black_female.RData"))
  rm(test)
  gc()
  
  
}else if(process==5){
  #Black male
  print("Black male")
  black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(black_rm_male, list(black_rm_male$zip))
  num_uniq_zip <- length(unique(black_rm_male$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)), 
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_black_male.RData"))
  rm(test)
  gc()
  
  
} else if(process==6){
  #Hispanic Female
  print("Hispanic female")
  hispanic_rm_female<-aggregate_data_rm %>% 
    filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(hispanic_rm_female, list(hispanic_rm_female$zip))
  num_uniq_zip <- length(unique(hispanic_rm_female$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_hispanic_female.RData"))
  gc()
  
  stopCluster(cl)
  
} else if(process==7){
  #Hispanic Male
  print("Hispanic male")
  hispanic_rm_male<-aggregate_data_rm %>% 
    filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(hispanic_rm_male, list(hispanic_rm_male$zip))
  num_uniq_zip <- length(unique(hispanic_rm_male$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_hispanic_male.RData"))
  gc()
  
} else if(process==8){
  print("Asian female")
  #Asian female
  asian_rm_female<-aggregate_data_rm %>% 
    filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(asian_rm_female, list(asian_rm_female$zip))
  num_uniq_zip <- length(unique(asian_rm_female$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_asian_female.RData"))
  gc()
  
} else if(process==9){
  #Asian male
  print("Asian male")
  asian_rm_male<-aggregate_data_rm %>% 
    filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
  cl=makeCluster(20,outfile='')
  registerDoParallel(cl)
  aggregate_data.list<-split(asian_rm_male, list(asian_rm_male$zip))
  num_uniq_zip <- length(unique(asian_rm_male$zip))
  
  test<-data.frame(matrix(nrow=50, ncol=0))
  
  for (boots_id in 1:500){
    set.seed(boots_id)
    print(boots_id)
    zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
    aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
    
    gam_raw<-mgcv::bam(dead~  s(pm25, k=3) + 
                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                         mean_bmi + smoke_rate + hispanic+ pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region)
                       +offset(log(time_count)),
                       data=aggregate_data_boots,family=poisson(link="log"))
    
    test<-cbind(test, test.data.rm(aggregate_data_boots, gam_raw))
    print(boots_id)
    rm(aggregate_data_boots)
  }
  stopCluster(cl)
  
  save(test,file=paste0(dir_out,"gam1_boots_rm_asian_male.RData"))
  gc()
}