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
#Load Poisson model
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_qd/'

load(file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main_trimmed.RData")
#load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd <- subset(aggregate_data_qd,
                            pm25_ensemble < quantile(aggregate_data_qd$pm25_ensemble,0.95)&
                              pm25_ensemble > quantile(aggregate_data_qd$pm25_ensemble,0.05))




# qd
#all
aggregate_data.list<-split(aggregate_data_qd, list(aggregate_data_qd$zip))
num_uniq_zip <- length(unique(aggregate_data_qd$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_all_trimmed.RData"))
exp(10*Poisson_qd$coefficients[1])
exp(10*(Poisson_qd$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#White female
require(dplyr)
load(file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main_trimmed.RData")

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_female_qd <- subset(white_female_qd,
                          pm25_ensemble < quantile(white_female_qd$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(white_female_qd$pm25_ensemble,0.05))
aggregate_data.list<-split(white_female_qd, list(white_female_qd$zip))
num_uniq_zip <- length(unique(white_female_qd$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_white_female_trimmed.RData"))
rm(white_female_qd)
gc()
exp(10*Poisson_qd_white_female$coefficients[1])
exp(10*(Poisson_qd_white_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_white_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#White Male
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
white_male_qd <- subset(white_male_qd,
                        pm25_ensemble < quantile(white_male_qd$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(white_male_qd$pm25_ensemble,0.05))
aggregate_data.list<-split(white_male_qd, list(white_male_qd$zip))
num_uniq_zip <- length(unique(white_male_qd$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_white_male_trimmed.RData"))
rm(white_male_qd)
gc()
exp(10*Poisson_qd_white_male$coefficients[1])
exp(10*(Poisson_qd_white_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_white_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Black female
black_qd_female<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_qd_female <- subset(black_qd_female,
                        pm25_ensemble < quantile(black_qd_female$pm25_ensemble,0.95)&
                          pm25_ensemble > quantile(black_qd_female$pm25_ensemble,0.05))
aggregate_data.list<-split(black_qd_female, list(black_qd_female$zip))
num_uniq_zip <- length(unique(black_qd_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_black_female_trimmed.RData"))
rm(black_qd_female)
exp(10*Poisson_qd_black_female$coefficients[1])
exp(10*(Poisson_qd_black_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_black_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Black male
black_qd_male<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
black_qd_male <- subset(black_qd_male,
                          pm25_ensemble < quantile(black_qd_male$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(black_qd_male$pm25_ensemble,0.05))

aggregate_data.list<-split(black_qd_male, list(black_qd_male$zip))
num_uniq_zip <- length(unique(black_qd_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_black_male_trimmed.RData"))
exp(10*Poisson_qd_black_male$coefficients[1])
exp(10*(Poisson_qd_black_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_black_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Hispanic Female
hispanic_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_qd_female <- subset(hispanic_qd_female,
                          pm25_ensemble < quantile(hispanic_qd_female$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(hispanic_qd_female$pm25_ensemble,0.05))

aggregate_data.list<-split(hispanic_qd_female, list(hispanic_qd_female$zip))
num_uniq_zip <- length(unique(hispanic_qd_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_female_trimmed.RData"))
rm(hispanic_qd_female)
exp(10*Poisson_qd_hispanic_female$coefficients[1])
exp(10*(Poisson_qd_hispanic_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_hispanic_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Hispanic Male
hispanic_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
hispanic_qd_male <- subset(hispanic_qd_male,
                             pm25_ensemble < quantile(hispanic_qd_male$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(hispanic_qd_male$pm25_ensemble,0.05))

aggregate_data.list<-split(hispanic_qd_male, list(hispanic_qd_male$zip))
num_uniq_zip <- length(unique(hispanic_qd_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_hispanic_male_trimmed.RData"))
rm(hispanic_qd_male)
exp(10*Poisson_qd_hispanic_male$coefficients[1])
exp(10*(Poisson_qd_hispanic_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_hispanic_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Asian female
asian_qd_female<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_qd_female <- subset(asian_qd_female,
                             pm25_ensemble < quantile(asian_qd_female$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(asian_qd_female$pm25_ensemble,0.05))

aggregate_data.list<-split(asian_qd_female, list(asian_qd_female$zip))
num_uniq_zip <- length(unique(asian_qd_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_asian_female_trimmed.RData"))
rm(asian_qd_female)
exp(10*Poisson_qd_asian_female$coefficients[1])
exp(10*(Poisson_qd_asian_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_asian_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Asian male
asian_qd_male<-aggregate_data_qd %>% 
  filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
asian_qd_male <- subset(asian_qd_male,
                          pm25_ensemble < quantile(asian_qd_male$pm25_ensemble,0.95)&
                            pm25_ensemble > quantile(asian_qd_male$pm25_ensemble,0.05))
aggregate_data.list<-split(asian_qd_male, list(asian_qd_male$zip))
num_uniq_zip <- length(unique(asian_qd_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_qd_asian_male_trimmed.RData"))
rm(asian_qd_male)
exp(10*Poisson_qd_asian_male$coefficients[1])
exp(10*(Poisson_qd_asian_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_qd_asian_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#RM
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/Poisson_rm/'

load(file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main_trimmed.RData")
#load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm <- subset(aggregate_data_rm,
                            pm25 < quantile(aggregate_data_rm$pm25,0.95)&
                              pm25 > quantile(aggregate_data_rm$pm25,0.05))




# rm
#all
aggregate_data.list<-split(aggregate_data_rm, list(aggregate_data_rm$zip))
num_uniq_zip <- length(unique(aggregate_data_rm$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_all_trimmed.RData"))
exp(10*Poisson_rm$coefficients[1])
exp(10*(Poisson_rm$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#White female
require(dplyr)
load(file= "/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Poisson/Main_trimmed.RData")

white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_female_rm <- subset(white_female_rm,
                          pm25 < quantile(white_female_rm$pm25,0.95)&
                            pm25 > quantile(white_female_rm$pm25,0.05))
aggregate_data.list<-split(white_female_rm, list(white_female_rm$zip))
num_uniq_zip <- length(unique(white_female_rm$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_white_female_trimmed.RData"))
rm(white_female_rm)
gc()
exp(10*Poisson_rm_white_female$coefficients[1])
exp(10*(Poisson_rm_white_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_white_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#White Male
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
white_male_rm <- subset(white_male_rm,
                        pm25 < quantile(white_male_rm$pm25,0.95)&
                          pm25 > quantile(white_male_rm$pm25,0.05))
aggregate_data.list<-split(white_male_rm, list(white_male_rm$zip))
num_uniq_zip <- length(unique(white_male_rm$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_white_male_trimmed.RData"))
rm(white_male_rm)
gc()
exp(10*Poisson_rm_white_male$coefficients[1])
exp(10*(Poisson_rm_white_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_white_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Black female
black_rm_female<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_rm_female <- subset(black_rm_female,
                          pm25 < quantile(black_rm_female$pm25,0.95)&
                            pm25 > quantile(black_rm_female$pm25,0.05))
aggregate_data.list<-split(black_rm_female, list(black_rm_female$zip))
num_uniq_zip <- length(unique(black_rm_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_black_female_trimmed.RData"))
rm(black_rm_female)
exp(10*Poisson_rm_black_female$coefficients[1])
exp(10*(Poisson_rm_black_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_black_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Black male
black_rm_male<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
black_rm_male <- subset(black_rm_male,
                        pm25 < quantile(black_rm_male$pm25,0.95)&
                          pm25 > quantile(black_rm_male$pm25,0.05))

aggregate_data.list<-split(black_rm_male, list(black_rm_male$zip))
num_uniq_zip <- length(unique(black_rm_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_black_male_trimmed.RData"))
exp(10*Poisson_rm_black_male$coefficients[1])
exp(10*(Poisson_rm_black_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_black_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Hispanic Female
hispanic_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_rm_female <- subset(hispanic_rm_female,
                             pm25 < quantile(hispanic_rm_female$pm25,0.95)&
                               pm25 > quantile(hispanic_rm_female$pm25,0.05))

aggregate_data.list<-split(hispanic_rm_female, list(hispanic_rm_female$zip))
num_uniq_zip <- length(unique(hispanic_rm_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_hispanic_female_trimmed.RData"))
rm(hispanic_rm_female)
exp(10*Poisson_rm_hispanic_female$coefficients[1])
exp(10*(Poisson_rm_hispanic_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_hispanic_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Hispanic Male
hispanic_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
hispanic_rm_male <- subset(hispanic_rm_male,
                           pm25 < quantile(hispanic_rm_male$pm25,0.95)&
                             pm25 > quantile(hispanic_rm_male$pm25,0.05))

aggregate_data.list<-split(hispanic_rm_male, list(hispanic_rm_male$zip))
num_uniq_zip <- length(unique(hispanic_rm_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_hispanic_male_trimmed.RData"))
rm(hispanic_rm_male)
exp(10*Poisson_rm_hispanic_male$coefficients[1])
exp(10*(Poisson_rm_hispanic_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_hispanic_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))


#Asian female
asian_rm_female<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_rm_female <- subset(asian_rm_female,
                          pm25 < quantile(asian_rm_female$pm25,0.95)&
                            pm25 > quantile(asian_rm_female$pm25,0.05))

aggregate_data.list<-split(asian_rm_female, list(asian_rm_female$zip))
num_uniq_zip <- length(unique(asian_rm_female$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_asian_female_trimmed.RData"))
rm(asian_rm_female)
exp(10*Poisson_rm_asian_female$coefficients[1])
exp(10*(Poisson_rm_asian_female$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_asian_female$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#Asian male
asian_rm_male<-aggregate_data_rm %>% 
  filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
asian_rm_male <- subset(asian_rm_male,
                        pm25 < quantile(asian_rm_male$pm25,0.95)&
                          pm25 > quantile(asian_rm_male$pm25,0.05))
aggregate_data.list<-split(asian_rm_male, list(asian_rm_male$zip))
num_uniq_zip <- length(unique(asian_rm_male$zip))
load(paste0(dir_out,"loglinear_coefs_boots_rm_asian_male_trimmed.RData"))
rm(asian_rm_male)
exp(10*Poisson_rm_asian_male$coefficients[1])
exp(10*(Poisson_rm_asian_male$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(Poisson_rm_asian_male$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

