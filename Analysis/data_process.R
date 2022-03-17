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
# national_merged2016_rm <- national_merged2016[, c(1:11,28,13:27,12)] # reorder to match QD
# rm(national_merged2016)
# national_merged2016_rm <- subset(national_merged2016_rm, select=-pm25_ensemble)
# save(national_merged2016_rm, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_rm.RData")

### Build Clean Aggregate Data Set

## QD Outcomes, Exposures, Covariates, and Offsets

# load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
# national_merged2016_qd$time_count <- national_merged2016_qd$followup_year_plus_one-national_merged2016_qd$followup_year
# 
# dead_personyear <- aggregate(data.frame(dead = national_merged2016_qd$dead,
#                                         time_countnational_merged2016_qd$time_count),
#                              by=list(zip = national_merged2016_qd$zip,
#                                      year = national_merged2016_qd$year,
#                                      sex = national_merged2016_qd$sex,
#                                      race = national_merged2016_qd$race,
#                                      dual = national_merged2016_qd$dual,
#                                      entry_age_break = national_merged2016_qd$entry_age_break,
#                                      followup_year = national_merged2016_qd$followup_year),
#                              FUN=sum)
# 
# confounders <- aggregate(national_merged2016_qd[,c(12:27)],
#                          by=list(zip = national_merged2016_qd$zip,
#                                  year = national_merged2016_qd$year,
#                                  sex = national_merged2016_qd$sex,
#                                  race = national_merged2016_qd$race,
#                                  dual = national_merged2016_qd$dual,
#                                  entry_age_break = national_merged2016_qd$entry_age_break,
#                                  followup_year = national_merged2016_qd$followup_year),
#                          FUN=min)
# 
# aggregate_data_qd <- merge(dead_personyear, confounders,
#                            by = c("zip", "year","sex","race","dual",
#                                   "entry_age_break", "followup_year"))
# 
# aggregate_data_qd <- aggregate_data_qd[complete.cases(aggregate_data_qd),]
# 
# save(aggregate_data_qd, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_qd.RData")

## RM Outcomes, Exposures, Covariates, and Offsets

# load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_rm.RData")
# national_merged2016_rm$time_count <- national_merged2016_rm$followup_year_plus_one-national_merged2016_rm$followup_year
# 
# dead_personyear <- aggregate(data.frame(dead = national_merged2016_rm$dead,
#                                         time_countnational_merged2016_rm$time_count),
#                              by=list(zip = national_merged2016_rm$zip,
#                                      year = national_merged2016_rm$year,
#                                      sex = national_merged2016_rm$sex,
#                                      race = national_merged2016_rm$race,
#                                      dual = national_merged2016_rm$dual,
#                                      entry_age_break = national_merged2016_rm$entry_age_break,
#                                      followup_year = national_merged2016_rm$followup_year),
#                              FUN=sum)
# 
# confounders <- aggregate(national_merged2016_rm[,c(12:27)],
#                          by=list(zip = national_merged2016_rm$zip,
#                                  year = national_merged2016_rm$year,
#                                  sex = national_merged2016_rm$sex,
#                                  race = national_merged2016_rm$race,
#                                  dual = national_merged2016_rm$dual,
#                                  entry_age_break = national_merged2016_rm$entry_age_break,
#                                  followup_year = national_merged2016_rm$followup_year),
#                          FUN=min)
# 
# aggregate_data_rm <- merge(dead_personyear, confounders,
#                            by = c("zip", "year","sex","race","dual",
#                                   "entry_age_break", "followup_year"))
# 
# aggregate_data_rm <- aggregate_data_rm[complete.cases(aggregate_data_rm),]
# 
# save(aggregate_data_rm, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_rm.RData")

### Create Strata Data

create_strata <- function(data, dual = c(0,1,2), race = c("all", "white", "black")) {
  
  zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
               "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
  
  if (dual == 0) {
    dual0 <- 0
  } else if (dual == 1) {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  if (race == "white") {
    race0 <- 1
  } else if (race == "black") {
    race0 <- 2
  } else {
    race0 <- c(1,2,3,4,5)
  }
  
  sub_data <- data %>% filter(data$race %in% race0 & data$dual %in% dual0)
  
  # Covariates and Outcomes
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  sex = sub_data$sex, dual = sub_data$dual, age_break = sub_data$age_break,
                  dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "sex", "dual", "age_break")]
  
  x <- data.table(zip = sub_data$zip, year = sub_data$year,
                  model.matrix(~ ., data = sub_data[,zip_cov])[,-1])[
                    ,lapply(.SD, min), by = c("zip", "year")]
  
  return(list(w = w, x = x))
  
}

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2),
                         race = c("white", "black", "all"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/'
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'

## QD Strata

load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd$pm25 <- aggregate_data_qd$pm25_ensemble
aggregate_data_qd$pm25_ensemble <- NULL
aggregate_data_qd$zip <- factor(aggregate_data_qd$zip)
aggregate_data_qd$year <- factor(aggregate_data_qd$year)
aggregate_data_qd$sex <- factor(aggregate_data_qd$sex)
aggregate_data_qd$dual <- as.numeric(aggregate_data_qd$dual)
aggregate_data_qd$region <- factor(aggregate_data_qd$region)

aggregate_data_qd$age_break <- with(aggregate_data_qd, ifelse(entry_age_break %in% c(1,2), 1,
                                                              ifelse(entry_age_break %in% c(3,4), 2, 
                                                                     ifelse(entry_age_break %in% c(5,6), 3, 4))))
aggregate_data_qd$age_break <- factor(aggregate_data_qd$age_break)

aggregate_data_qd$race <- with(aggregate_data_qd, ifelse(!(race %in% c(1,2,3,4)), 5, race))
aggregate_data_qd$race <- factor(aggregate_data_qd$race)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data_qd, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
})

## RM Strata

load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm$zip <- factor(aggregate_data_rm$zip)
aggregate_data_rm$year <- factor(aggregate_data_rm$year)
aggregate_data_rm$sex <- factor(aggregate_data_rm$sex)
aggregate_data_rm$dual <- as.numeric(aggregate_data_rm$dual)
aggregate_data_rm$region <- factor(aggregate_data_rm$region)

aggregate_data_rm$age_break <- with(aggregate_data_rm, ifelse(entry_age_break %in% c(1,2), 1,
                                                              ifelse(entry_age_break %in% c(3,4), 2, 
                                                                     ifelse(entry_age_break %in% c(5,6), 3, 4))))
aggregate_data_rm$age_break <- factor(aggregate_data_rm$age_break)
aggregate_data_rm$race <- with(aggregate_data_rm, ifelse(!(race %in% c(1,2,3,4)), 5, race))
aggregate_data_rm$race <- factor(aggregate_data_rm$race)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data_rm, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
})