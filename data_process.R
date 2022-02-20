# Generate all required data from different data sources (see Table S1).
library(fst)
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(foreign)

### Big Data Cleaning (BEWARE!)

## Qian Di

# Condense Data
# f <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
#                 pattern = "\\.fst",
#                 full.names = TRUE)
# 
# myvars <- c("qid", "year","zip","sex","race","age","dual","entry_age_break","statecode",
#             "followup_year","followup_year_plus_one","dead","pm25_ensemble",
#             "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
#             "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")
# 
# national_merged2016_qd <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
# national_merged2016_qd$zip <- sprintf("%05d", national_merged2016_qd$zip)
# 
# NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")  
# SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
# MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
# WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")
# 
# # creates region
# national_merged2016_qd$region=ifelse(national_merged2016_qd$state %in% NORTHEAST, "NORTHEAST", 
#                                      ifelse(national_merged2016_qd$state %in% SOUTH, "SOUTH",
#                                             ifelse(national_merged2016_qd$state %in% MIDWEST, "MIDWEST",
#                                                    ifelse(national_merged2016_qd$state %in% WEST, "WEST", NA))))
# 
# national_merged2016_qd <- national_merged2016_qd[complete.cases(national_merged2016_qd[,c(1:27)]) ,]
# save(national_merged2016_qd, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
# 
# load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")

## Randall Martin

# pm_rm <- data.frame
# 
# for(i in 2000:2016){
#   temp <- read.csv(paste0("~/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/rm_predictions/ben_2019_10_29/data_acag_pm25_zip-year/zip_pm25_",i, ".csv"))
#   temp$YEAR <- rep(i, nrow(temp))
#   pm_rm<-rbind(pm_rm, temp)
# }
# 
# rm(temp); gc()
# 
# pm_rm <- subset(pm_rm, !is.na(pm_rm$ZIP))
# pm_rm$ZIP <- sprintf("%05d", pm_rm$ZIP)
# pm_rm <- subset(pm_rm, pm_rm$ZIP %in% unique(national_merged2016$zip))
# pm_rm <- as.data.table(pm_rm)
# 
# load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
# national_merged2016 <- national_merged2016_qd
# rm(national_merged2016_qd);gc()
# 
# national_merged2016 <- as.data.table(national_merged2016)
# setkey(national_merged2016, zip, year)
# setkey(pm_rm, ZIP, YEAR)
# 
# national_merged2016[pm_rm, pm25:=i.pm25]
# national_merged2016_rm <- national_merged2016[, c(1:11,28, 13:27, 12)] # reorder to match QD
# rm(national_merged2016)
# national_merged2016_rm <- subset(national_merged2016_rm, select=-pm25_ensemble)
# save(national_merged2016_rm, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_rm.RData")

### Build Clean Aggregate Data Set

# QD Outcomes, Exposures, Covariates, and Offsets
load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
national_merged2016_qd$time_count<-national_merged2016_qd$followup_year_plus_one-national_merged2016_qd$followup_year

dead_personyear <- aggregate(cbind(national_merged2016_qd$dead,national_merged2016_qd$time_count),
                           by=list(national_merged2016_qd$zip,
                                   national_merged2016_qd$year,
                                   national_merged2016_qd$sex,
                                   national_merged2016_qd$race),
                           FUN=sum)

confounders<-aggregate(national_merged2016_qd[,c(12:27)], 
                       by=list(national_merged2016_qd$zip,
                               national_merged2016_qd$year,
                               national_merged2016_qd$sex,
                               national_merged2016_qd$race), 
                       FUN=min)

zip_year <- paste0(national_merged2016_qd$zip, "-", national_merged2016_qd$year)
ind <- blp(w = with(national_merged2016_qd, cbind(dual, entry_age_break, followup_year)), w.id = zip_year)

aggregate_data_qd <- merge(dead_personyear, confounders, 
                           by = c("Group.1", "Group.2","Group.3", "Group.4"))

colnames(aggregate_data_qd)[8:9] <- c("dead","time_count")
colnames(aggregate_data_qd)[1:7] <- c("zip","year","sex","race")
aggregate_data_qd <- aggregate_data_qd[complete.cases(aggregate_data_qd),]

save(aggregate_data_qd, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_qd.RData")

# RM Outcomes, Exposures, Covariates, and Offsets
load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_rm.RData")
national_merged2016_rm$time_count<-national_merged2016_rm$followup_year_plus_one-national_merged2016_rm$followup_year

dead_personyear <- aggregate(cbind(national_merged2016_rm$dead,national_merged2016_rm$time_count),
                             by=list(national_merged2016_rm$zip,
                                     national_merged2016_rm$year,
                                     national_merged2016_rm$sex,
                                     national_merged2016_rm$race), 
                             FUN=sum)

confounders <- aggregate(national_merged2016_rm[,c(12:27)], 
                         by=list(national_merged2016_rm$zip,
                                 national_merged2016_rm$year,
                                 national_merged2016_rm$sex,
                                 national_merged2016_rm$race), 
                         FUN=min)

zip_year <- paste0(national_merged2016_rm$zip, "-", national_merged2016_rm$year)
ind <- blp(w = with(national_merged2016_qd, cbind(dual, entry_age_break, followup_year)), w.id = zip_year)

aggregate_data_rm<-merge(dead_personyear,confounders, 
                         by = c("Group.1", "Group.2", "Group.3", "Group.4"))

colnames(aggregate_data_rm)[8:9] <- c("dead","time_count")
colnames(aggregate_data_rm)[1:7] <- c("zip","year","sex","race")
aggregate_data_rm <- subset(aggregate_data_rm[complete.cases(aggregate_data_rm) ,])

save(aggregate_data_rm,file="~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_rm.RData")

### Create Strata Data

create_strata <- function(data, sex = c("both", "male", "female"), race = c("all", "white", "black", "hispanic", "asian")) {

  zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
               "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region") 
  ind_cov <- c("dual", "entry_age_break", "followup_year")
  
  if (sex == "male") {
    sex <- 1
  } else if (sex == "female") {
    sex <- 2
  } else {
    sex <- c(1,2)
  }
  
  if (race == "white") {
    race <- 1
  } else if (race == "black") {
    race <- 2
  } else if (race == "hispanic") {
    race <- 5
  } else if (race == "asian") {
    race <- 4
  } else {
    race <- c(0,1,2,3,4,5,6)
  }
  
  sub_data <- data %>% filter(data$race %in% race & data$sex %in% sex)
  
  # Covariates and Outcomes
  outcomes <- data.table(zip = sub_data$zip, year = sub_data$year,
                         dead = sub_data$dead, time_count = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year")]
  
  cmat <- data.frame(model.matrix(~ ., data = sub_data[,ind_cov])[,-1])
  ind_covariates <- data.frame(zip = outcomes$zip, year = outcomes$year)
  
  for (i in 1:ncol(cmat)){
    
    mat <- data.table(zip = sub_data$zip, year = sub_data$year, time_count = sub_data$time_count, val = cmat[,i])
    ind_tmp <- mat[,list(weighted.mean(val,time_count)), by = c("zip", "year")]
    colnames(ind_tmp)[3] <- colnames(cmat)[i]
    ind_covariates <- merge(ind_covariates, ind_tmp, by = c("zip", "year"), all.x = TRUE)
    
  }
  
  zip_covariates <- data.table(zip = sub_data$zip, year = sub_data$year,
                               model.matrix(~ ., data = sub_data[,zip_cov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
  
  new_data <- merge(outcomes, ind_covariates, by = c("zip", "year")) %>%
    merge(zip_covariates, by = c("zip", "year"))
  new_data$zip <- factor(new_data$zip)
  new_data$year <- factor(new_data$year)
  
  return(new_data)
  
}

# scenarios
scenarios <- expand.grid(sex = c("both", "male", "female"),
                         race = c("all", "white", "black", "hispanic", "asian"))

# Save Location
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd_strata/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm_strata/'

## QD Strata

load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd$pm25 <- aggregate_data_qd$pm25_ensemble
aggregate_data_qd$pm25_ensemble <- NULL
aggregate_data_qd$year <- factor(aggregate_data_qd$year)
aggregate_data_qd$sex <- factor(aggregate_data_qd$sex)
aggregate_data_qd$race <- factor(aggregate_data_qd$race)
aggregate_data_qd$entry_age_break <- factor(aggregate_data_qd$entry_age_break)
aggregate_data_qd$region <- factor(aggregate_data_qd$region)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data_qd, sex = scenario$sex, race = scenario$race)
  save(new_data, file = paste0(dir_out_qd, scenario$sex, "_", scenario$race, "_qd.RData"))
  
})

## RM Strata

load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm$year <- factor(aggregate_data_rm$year)
aggregate_data_rm$sex <- factor(aggregate_data_rm$sex)
aggregate_data_rm$race <- factor(aggregate_data_rm$race)
aggregate_data_rm$entry_age_break <- factor(aggregate_data_rm$entry_age_break)
aggregate_data_rm$region <- factor(aggregate_data_rm$region)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data_rm, sex = scenario$sex, race = scenario$race)
  save(new_data, file = paste0(dir_out_rm, scenario$sex, "_", scenario$race, "_rm.RData"))
  
})
