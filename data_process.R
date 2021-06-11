# Generate all required data from different data sources (see Table S1).
library(fst)
library(data.table)
library(parallel)
library("xgboost")
library("dplyr")
library("foreign")

f <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst",
                full.names = TRUE)

myvars <- c("qid", "year","zip","sex","race","age","dual","entry_age_break","statecode",
            "followup_year","followup_year_plus_one","dead","pm25_ensemble",
            "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
            "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")
national_merged2016_qd <- rbindlist(lapply(f,
                        read_fst,
                        columns = myvars,
                        as.data.table=TRUE))
#n<-national_merged2016_qd[!duplicated(national_merged2016_qd$qid),]

national_merged2016_qd<-as.data.frame(national_merged2016_qd)

national_merged2016_qd$zip <- sprintf("%05d", national_merged2016_qd$zip)

NORTHEAST=c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")  
SOUTH=c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST=c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST=c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

national_merged2016_qd$region=ifelse(national_merged2016_qd$state %in% NORTHEAST, "NORTHEAST", 
                                  ifelse(national_merged2016_qd$state %in% SOUTH, "SOUTH",
                                         ifelse(national_merged2016_qd$state %in% MIDWEST, "MIDWEST",
                                                ifelse(national_merged2016_qd$state %in% WEST, "WEST",
                                                       NA))))

national_merged2016_qd <- national_merged2016_qd[complete.cases(national_merged2016_qd[,c(1:27)]) ,]
#> dim(national_merged2016)
#[1] 573,370,257        27

#QD's analysis
save(national_merged2016_qd,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/national_merged2016_qd.RData")

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/national_merged2016_qd.RData")
# Main analysis with QD
covariates_qd<-aggregate(national_merged2016_qd[,c(12:27)], 
                         by=list(national_merged2016_qd$zip,national_merged2016_qd$year), 
                         FUN=min)
colnames(covariates_qd)[1:2]<-c("zip","year")

covariates_qd<-subset(covariates_qd[complete.cases(covariates_qd) ,])
covariates_qd$year_fac <- as.factor(covariates_qd$year)
covariates_qd$region <- as.factor(covariates_qd$region)

save(covariates_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")

#Strata
white_female_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==1 & national_merged2016_qd$sex==2)
covariates_white_female_qd<-aggregate(white_female_qd[,c(12:27)], 
                         by=list(white_female_qd$zip,white_female_qd$year), 
                         FUN=min)
colnames(covariates_white_female_qd)[1:2]<-c("zip","year")
covariates_white_female_qd<-subset(covariates_white_female_qd[complete.cases(covariates_white_female_qd) ,])
covariates_white_female_qd$year_fac <- as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region <- as.factor(covariates_white_female_qd$region)
save(covariates_white_female_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_female_qd.RData")
rm(white_female_qd, covariates_white_female_qd)
gc()

white_male_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==1 & national_merged2016_qd$sex==1)
covariates_white_male_qd<-aggregate(white_male_qd[,c(12:27)], 
                                      by=list(white_male_qd$zip,white_male_qd$year), 
                                      FUN=min)
colnames(covariates_white_male_qd)[1:2]<-c("zip","year")
covariates_white_male_qd<-subset(covariates_white_male_qd[complete.cases(covariates_white_male_qd) ,])
covariates_white_male_qd$year_fac <- as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region <- as.factor(covariates_white_male_qd$region)
save(covariates_white_male_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_male_qd.RData")
rm(white_male_qd, covariates_white_male_qd)

black_female_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==2 & national_merged2016_qd$sex==2)
covariates_black_female_qd<-aggregate(black_female_qd[,c(12:27)], 
                                    by=list(black_female_qd$zip, black_female_qd$year), 
                                    FUN=min)
colnames(covariates_black_female_qd)[1:2]<-c("zip","year")
covariates_black_female_qd<-subset(covariates_black_female_qd[complete.cases(covariates_black_female_qd) ,])
covariates_black_female_qd$year_fac <- as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region <- as.factor(covariates_black_female_qd$region)
save(covariates_black_female_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_female_qd.RData")
rm(black_female_qd, covariates_black_female_qd)
gc()

black_male_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==2 & national_merged2016_qd$sex==1)
covariates_black_male_qd<-aggregate(black_male_qd[,c(12:27)], 
                                      by=list(black_male_qd$zip, black_male_qd$year), 
                                      FUN=min)
colnames(covariates_black_male_qd)[1:2]<-c("zip","year")
covariates_black_male_qd<-subset(covariates_black_male_qd[complete.cases(covariates_black_male_qd) ,])
covariates_black_male_qd$year_fac <- as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region <- as.factor(covariates_black_male_qd$region)
save(covariates_black_male_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_male_qd.RData")
rm(black_male_qd, covariates_black_male_qd)
gc()

hispanic_female_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==5 & national_merged2016_qd$sex==2)
covariates_hispanic_female_qd<-aggregate(hispanic_female_qd[,c(12:27)], 
                                    by=list(hispanic_female_qd$zip, hispanic_female_qd$year), 
                                    FUN=min)
colnames(covariates_hispanic_female_qd)[1:2]<-c("zip","year")
covariates_hispanic_female_qd<-subset(covariates_hispanic_female_qd[complete.cases(covariates_hispanic_female_qd) ,])
covariates_hispanic_female_qd$year_fac <- as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region <- as.factor(covariates_hispanic_female_qd$region)
save(covariates_hispanic_female_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_female_qd.RData")
rm(hispanic_female_qd, covariates_hispanic_female_qd)
gc()

hispanic_male_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==5 & national_merged2016_qd$sex==1)
covariates_hispanic_male_qd<-aggregate(hispanic_male_qd[,c(12:27)], 
                                         by=list(hispanic_male_qd$zip, hispanic_male_qd$year), 
                                         FUN=min)
colnames(covariates_hispanic_male_qd)[1:2]<-c("zip","year")
covariates_hispanic_male_qd<-subset(covariates_hispanic_male_qd[complete.cases(covariates_hispanic_male_qd) ,])
covariates_hispanic_male_qd$year_fac <- as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region <- as.factor(covariates_hispanic_male_qd$region)
save(covariates_hispanic_male_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_male_qd.RData")
rm(covariates_hispanic_male_qd, hispanic_male_qd)

asian_female_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==4 & national_merged2016_qd$sex==2)
covariates_asian_female_qd<-aggregate(asian_female_qd[,c(12:27)], 
                                       by=list(asian_female_qd$zip, asian_female_qd$year), 
                                       FUN=min)
colnames(covariates_asian_female_qd)[1:2]<-c("zip","year")
covariates_asian_female_qd<-subset(covariates_asian_female_qd[complete.cases(covariates_asian_female_qd) ,])
covariates_asian_female_qd$year_fac <- as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region <- as.factor(covariates_asian_female_qd$region)
save(covariates_asian_female_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_female_qd.RData")
rm(asian_female_qd, covariates_asian_female_qd)
gc()

asian_male_qd<-national_merged2016_qd %>% filter(national_merged2016_qd$race==4 & national_merged2016_qd$sex==1)
covariates_asian_male_qd<-aggregate(asian_male_qd[,c(12:27)], 
                                      by=list(asian_male_qd$zip, asian_male_qd$year), 
                                      FUN=min)
colnames(covariates_asian_male_qd)[1:2]<-c("zip","year")
covariates_asian_male_qd<-subset(covariates_asian_male_qd[complete.cases(covariates_asian_male_qd) ,])
covariates_asian_male_qd$year_fac <- as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region <- as.factor(covariates_asian_male_qd$region)
save(covariates_asian_male_qd,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_male_qd.RData")
rm(asian_male_qd, covariates_asian_male_qd)
gc()

# Fit GPS model via Xgboost machine 
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")
GPS_mod_qd <-xgboost(data = data.matrix(covariates_qd[,c(4:19)]), 
                     label = covariates_qd$pm25_ensemble,
                     nrounds=50)
mod_sd_qd<- sd(covariates_qd$pm25_ensemble -predict(GPS_mod_qd,data.matrix(covariates_qd[,c(4:19)])))
feature_names_qd <- GPS_mod_qd$feature_names

covariates_qd$GPS<-dnorm(covariates_qd$pm25_ensemble,
                         mean = predict(GPS_mod_qd,data.matrix(covariates_qd[,c(4:19)])),
                      sd=sd(covariates_qd$pm25_ensemble-predict(GPS_mod_qd,data.matrix(covariates_qd[,c(4:19)]))))

Nm_qd<-dnorm(covariates_qd$pm25_ensemble,
             mean=mean(covariates_qd$pm25_ensemble,na.rm=T),
             sd=sd(covariates_qd$pm25_ensemble,na.rm=T))
covariates_qd$IPW<-Nm_qd/(covariates_qd$GPS)
covariates_qd<-covariates_qd[,c("zip","year","IPW","GPS")]
save(GPS_mod_qd,mod_sd_qd,feature_names_qd,covariates_qd,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/covariates_qd.RData")

# Generate count data for each individual characteristics and follow-up year
national_merged2016_qd$time_count<-national_merged2016_qd$followup_year_plus_one-national_merged2016_qd$followup_year
dead_personyear<-aggregate(cbind(national_merged2016_qd$dead,
                                 national_merged2016_qd$time_count), 
                           by=list(national_merged2016_qd$zip,
                                   national_merged2016_qd$year,
                                   national_merged2016_qd$sex,
                                   national_merged2016_qd$race,
                                   national_merged2016_qd$dual,
                                   national_merged2016_qd$entry_age_break,
                                   national_merged2016_qd$followup_year), 
                           FUN=sum)

confounders<-aggregate(national_merged2016_qd[,c(12:27)], 
                       by=list(national_merged2016_qd$zip,
                               national_merged2016_qd$year,
                                national_merged2016_qd$sex,
                               national_merged2016_qd$race,
                               national_merged2016_qd$dual,
                               national_merged2016_qd$entry_age_break,
                               national_merged2016_qd$followup_year), 
                       FUN=min)
aggregate_data_qd<-merge(dead_personyear,confounders
                      ,by=c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "Group.6", "Group.7"))
colnames(aggregate_data_qd)[8:9]<-c("dead","time_count")
colnames(aggregate_data_qd)[1:7]<-c("zip","year","sex","race","dual","entry_age_break","followup_year")
aggregate_data_qd<-subset(aggregate_data_qd[complete.cases(aggregate_data_qd) ,])

save(aggregate_data_qd,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/aggregate_data_qd.RData")


#load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/covariates_qd.RData")
#load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/aggregate_data_qd.RData")
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd,by=c("zip","year"),all.x=T)

aggregate_data.list_qd <- split(aggregate_data_qd,
                                list(aggregate_data_qd$sex,
                                     aggregate_data_qd$race, 
                                     aggregate_data_qd$dual,
                                     aggregate_data_qd$entry_age_break,
                                     aggregate_data_qd$followup_year))
aggregate_data.list_qd<-aggregate_data.list_qd[lapply(aggregate_data.list_qd,nrow)>0]

mclapply(aggregate_data.list_qd,function(data){
  write_fst(data,paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/FST_data_qd/",
                        data$sex[1],
                        data$race[1], 
                        data$dual[1],
                        data$entry_age_break[1],
                        data$followup_year[1],".fst"))
},mc.cores=16)


#Joining with Randall Martin's data
pm_rm<-NA
for(i in 2000:2016){
  temp<-read.csv(paste0("/nfs/home/P/prd789/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/rm_predictions/ben_2019_10_29/data_acag_pm25_zip-year/zip_pm25_",i, ".csv"))
  temp$YEAR<-rep(i, nrow(temp))
  pm_rm<-rbind(pm_rm, temp)
     
    }
rm(temp)
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/national_merged2016_qd.RData")
national_merged2016<-national_merged2016_qd
rm(national_merged2016_qd)

pm_rm<-subset(pm_rm, !is.na(pm_rm$ZIP))
pm_rm$ZIP <- sprintf("%05d", pm_rm$ZIP)
pm_rm<-subset(pm_rm, pm_rm$ZIP %in% unique(national_merged2016$zip))
pm_rm<-as.data.table(pm_rm)

national_merged2016<-as.data.table(national_merged2016)
setkey(national_merged2016, zip, year)
setkey(pm_rm, ZIP, YEAR)

national_merged2016[pm_rm, pm25:=i.pm25]
#Reordering to match Xiao's code
national_merged2016_rm<-national_merged2016[, c(1:11,28, 13:27, 12)]
rm(national_merged2016)
national_merged2016_rm<-subset(national_merged2016_rm, select=-pm25_ensemble)
save(national_merged2016_rm,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/national_merged2016_rm.RData")

# Main analysis
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/national_merged2016_rm.RData")
covariates_rm<-aggregate(national_merged2016_rm[,c(12:27)], by=list(national_merged2016_rm$zip,national_merged2016_rm$year), FUN=min)
colnames(covariates_rm)[1:2]<-c("zip","year")

covariates_rm<-subset(covariates_rm[complete.cases(covariates_rm) ,])
covariates_rm$year_fac <- as.factor(covariates_rm$year)
covariates_rm$region <- as.factor(covariates_rm$region)

save(covariates_rm,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_rm.RData")

#Strata
white_female_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==1 & national_merged2016_rm$sex==2)
covariates_white_female_rm<-aggregate(white_female_rm[,c(12:27)], 
                                      by=list(white_female_rm$zip,white_female_rm$year), 
                                      FUN=min)
colnames(covariates_white_female_rm)[1:2]<-c("zip","year")
covariates_white_female_rm<-subset(covariates_white_female_rm[complete.cases(covariates_white_female_rm) ,])
covariates_white_female_rm$year_fac <- as.factor(covariates_white_female_rm$year)
covariates_white_female_rm$region <- as.factor(covariates_white_female_rm$region)
save(covariates_white_female_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_white_female_rm.RData")
rm(white_female_rm, covariates_white_female_rm)
gc()

white_male_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==1 & national_merged2016_rm$sex==1)
covariates_white_male_rm<-aggregate(white_male_rm[,c(12:27)], 
                                    by=list(white_male_rm$zip,white_male_rm$year), 
                                    FUN=min)
colnames(covariates_white_male_rm)[1:2]<-c("zip","year")
covariates_white_male_rm<-subset(covariates_white_male_rm[complete.cases(covariates_white_male_rm) ,])
covariates_white_male_rm$year_fac <- as.factor(covariates_white_male_rm$year)
covariates_white_male_rm$region <- as.factor(covariates_white_male_rm$region)
save(covariates_white_male_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_white_male_rm.RData")
rm(white_male_rm, covariates_white_male_rm)

black_female_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==2 & national_merged2016_rm$sex==2)
covariates_black_female_rm<-aggregate(black_female_rm[,c(12:27)], 
                                      by=list(black_female_rm$zip, black_female_rm$year), 
                                      FUN=min)
colnames(covariates_black_female_rm)[1:2]<-c("zip","year")
covariates_black_female_rm<-subset(covariates_black_female_rm[complete.cases(covariates_black_female_rm) ,])
covariates_black_female_rm$year_fac <- as.factor(covariates_black_female_rm$year)
covariates_black_female_rm$region <- as.factor(covariates_black_female_rm$region)
save(covariates_black_female_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_black_female_rm.RData")
rm(black_female_rm, covariates_black_female_rm)
gc()

black_male_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==2 & national_merged2016_rm$sex==1)
covariates_black_male_rm<-aggregate(black_male_rm[,c(12:27)], 
                                    by=list(black_male_rm$zip, black_male_rm$year), 
                                    FUN=min)
colnames(covariates_black_male_rm)[1:2]<-c("zip","year")
covariates_black_male_rm<-subset(covariates_black_male_rm[complete.cases(covariates_black_male_rm) ,])
covariates_black_male_rm$year_fac <- as.factor(covariates_black_male_rm$year)
covariates_black_male_rm$region <- as.factor(covariates_black_male_rm$region)
save(covariates_black_male_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_black_male_rm.RData")
rm(black_male_rm, covariates_black_male_rm)
gc()

hispanic_female_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==5 & national_merged2016_rm$sex==2)
covariates_hispanic_female_rm<-aggregate(hispanic_female_rm[,c(12:27)], 
                                         by=list(hispanic_female_rm$zip, hispanic_female_rm$year), 
                                         FUN=min)
colnames(covariates_hispanic_female_rm)[1:2]<-c("zip","year")
covariates_hispanic_female_rm<-subset(covariates_hispanic_female_rm[complete.cases(covariates_hispanic_female_rm) ,])
covariates_hispanic_female_rm$year_fac <- as.factor(covariates_hispanic_female_rm$year)
covariates_hispanic_female_rm$region <- as.factor(covariates_hispanic_female_rm$region)
save(covariates_hispanic_female_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_hispanic_female_rm.RData")
rm(hispanic_female_rm, covariates_hispanic_female_rm)
gc()

hispanic_male_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==5 & national_merged2016_rm$sex==1)
covariates_hispanic_male_rm<-aggregate(hispanic_male_rm[,c(12:27)], 
                                       by=list(hispanic_male_rm$zip, hispanic_male_rm$year), 
                                       FUN=min)
colnames(covariates_hispanic_male_rm)[1:2]<-c("zip","year")
covariates_hispanic_male_rm<-subset(covariates_hispanic_male_rm[complete.cases(covariates_hispanic_male_rm) ,])
covariates_hispanic_male_rm$year_fac <- as.factor(covariates_hispanic_male_rm$year)
covariates_hispanic_male_rm$region <- as.factor(covariates_hispanic_male_rm$region)
save(covariates_hispanic_male_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_hispanic_male_rm.RData")
rm(covariates_hispanic_male_rm, hispanic_male_rm)

asian_female_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==4 & national_merged2016_rm$sex==2)
covariates_asian_female_rm<-aggregate(asian_female_rm[,c(12:27)], 
                                      by=list(asian_female_rm$zip, asian_female_rm$year), 
                                      FUN=min)
colnames(covariates_asian_female_rm)[1:2]<-c("zip","year")
covariates_asian_female_rm<-subset(covariates_asian_female_rm[complete.cases(covariates_asian_female_rm) ,])
covariates_asian_female_rm$year_fac <- as.factor(covariates_asian_female_rm$year)
covariates_asian_female_rm$region <- as.factor(covariates_asian_female_rm$region)
save(covariates_asian_female_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_asian_female_rm.RData")
rm(asian_female_rm, covariates_asian_female_rm)
gc()

asian_male_rm<-national_merged2016_rm %>% filter(national_merged2016_rm$race==4 & national_merged2016_rm$sex==1)
covariates_asian_male_rm<-aggregate(asian_male_rm[,c(12:27)], 
                                    by=list(asian_male_rm$zip, asian_male_rm$year), 
                                    FUN=min)
colnames(covariates_asian_male_rm)[1:2]<-c("zip","year")
covariates_asian_male_rm<-subset(covariates_asian_male_rm[complete.cases(covariates_asian_male_rm) ,])
covariates_asian_male_rm$year_fac <- as.factor(covariates_asian_male_rm$year)
covariates_asian_male_rm$region <- as.factor(covariates_asian_male_rm$region)
save(covariates_asian_male_rm,
     file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_asian_male_rm.RData")
rm(asian_male_rm, covariates_asian_male_rm)
gc()

# Fit GPS model via Xgboost machine
GPS_mod_rm <-xgboost(data = data.matrix(covariates_rm[,c(4:19)]), label = covariates_rm$pm25,nrounds=50)
mod_sd_rm<- sd(covariates_rm$pm25-predict(GPS_mod_rm,data.matrix(covariates_rm[,c(4:19)])))
feature_names_rm <- GPS_mod_rm$feature_names
covariates_rm$GPS<-dnorm(covariates_rm$pm25,mean = predict(GPS_mod_rm,data.matrix(covariates_rm[,c(4:19)])),
                      sd=sd(covariates_rm$pm25-predict(GPS_mod_rm,data.matrix(covariates_rm[,c(4:19)]))))
Nm_rm<-dnorm(covariates_rm$pm25,
             mean=mean(covariates_rm$pm25,na.rm=T),
             sd=sd(covariates_rm$pm25,na.rm=T))
covariates_rm$IPW<-Nm_rm/(covariates_rm$GPS)
covariates_rm<-covariates_rm[,c("zip","year","IPW","GPS")]

save(GPS_mod_rm,mod_sd_rm,feature_names_rm,covariates_rm,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/covariates_rm.RData")

# Generate count data for each individual characteristics and follow-up year
national_merged2016_rm$time_count<-national_merged2016_rm$followup_year_plus_one-national_merged2016_rm$followup_year
dead_personyear<-aggregate(cbind(national_merged2016_rm$dead,national_merged2016_rm$time_count),
                           by=list(national_merged2016_rm$zip,national_merged2016_rm$year,
                                   national_merged2016_rm$sex,
                                   national_merged2016_rm$race,
                                   national_merged2016_rm$dual,
                                   national_merged2016_rm$entry_age_break,
                                   national_merged2016_rm$followup_year), FUN=sum)

confounders<-aggregate(national_merged2016_rm[,c(12:27)], 
                       by=list(national_merged2016_rm$zip,
                               national_merged2016_rm$year,
                               national_merged2016_rm$sex,
                               national_merged2016_rm$race,
                               national_merged2016_rm$dual,
                               national_merged2016_rm$entry_age_break,
                               national_merged2016_rm$followup_year), 
                       FUN=min)

aggregate_data_rm<-merge(dead_personyear,confounders
                      ,by=c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "Group.6", "Group.7"))
colnames(aggregate_data_rm)[8:9]<-c("dead","time_count")
colnames(aggregate_data_rm)[1:7]<-c("zip","year","sex","race","dual","entry_age_break","followup_year")
aggregate_data_rm<-subset(aggregate_data_rm[complete.cases(aggregate_data_rm) ,])

save(aggregate_data_rm,file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/aggregate_data_rm.RData")

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/covariates_rm.RData")
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/aggregate_data_rm.RData")
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

aggregate_data.list_rm <- split(aggregate_data_rm,
                                list(aggregate_data_rm$sex,
                                     aggregate_data_rm$race, 
                                     aggregate_data_rm$dual,
                                     aggregate_data_rm$entry_age_break,
                                     aggregate_data_rm$followup_year))
aggregate_data.list_rm<-aggregate_data.list_rm[lapply(aggregate_data.list_rm,nrow)>0]

mclapply(aggregate_data.list_rm,function(data){
  write_fst(data,paste0("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/FST_data_rm/",
                        data$sex[1],
                        data$race[1], 
                        data$dual[1],
                        data$entry_age_break[1],
                        data$followup_year[1],".fst"))
},mc.cores=16)


