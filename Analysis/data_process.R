# Generate all required data from different data sources (see Table S1).
library(fst)
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(foreign)

### Big Data Cleaning (BEWARE!)

## Qian Di

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

## Randall Martin

# pm_rm <- data.frame()
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

load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_qd.RData")
national_merged2016_qd$time_count <- rep(1, nrow(national_merged2016_qd))
national_merged2016_qd$age_break <- cut(national_merged2016_qd$age, c(65,70,75,80,85,90,95,125), right = FALSE)
national_merged2016_qd$female <- national_merged2016_qd$sex - 1
colnames(national_merged2016_qd)[12] <- "pm25"

dead_personyear <- aggregate(data.frame(dead = national_merged2016_qd$dead,
                                        time_count = national_merged2016_qd$time_count),
                             by=list(zip = national_merged2016_qd$zip,
                                     year = national_merged2016_qd$year,
                                     female = national_merged2016_qd$female,
                                     race = national_merged2016_qd$race,
                                     dual = national_merged2016_qd$dual,
                                     age_break = national_merged2016_qd$age_break),
                             FUN=sum)

new_data <- national_merged2016_qd %>% distinct(zip, year, female, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016_qd, new_data); gc()

aggregate_data_qd <- merge(dead_personyear, confounders, by = c("zip", "year","female","race","dual","age_break"))
aggregate_data_qd <- aggregate_data_qd[complete.cases(aggregate_data_qd),]

save(aggregate_data_qd, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_qd.RData")

rm(dead_personyear, confounders, aggregate_data_qd); gc()

## RM Outcomes, Exposures, Covariates, and Offsets

load("~/shared_space/ci3_analysis/josey_erc_strata/Data/national_merged2016_rm.RData")
national_merged2016_rm$time_count <- rep(1, nrow(national_merged2016_rm))
national_merged2016_rm$age_break <- cut(national_merged2016_rm$age, c(65,70,75,80,85,90,95,125), right = FALSE)
national_merged2016_rm$female <- national_merged2016_rm$sex - 1

dead_personyear <- aggregate(data.frame(dead = national_merged2016_rm$dead,
                                        time_count = national_merged2016_rm$time_count),
                             by=list(zip = national_merged2016_rm$zip,
                                     year = national_merged2016_rm$year,
                                     female = national_merged2016_rm$female,
                                     race = national_merged2016_rm$race,
                                     dual = national_merged2016_rm$dual,
                                     age_break = national_merged2016_rm$age_break),
                             FUN=sum)

new_data <- national_merged2016_rm %>% distinct(zip, year, female, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016_rm, new_data); gc()

aggregate_data_rm <- merge(dead_personyear, confounders, by = c("zip", "year","female","race","dual","age_break"))
aggregate_data_rm <- aggregate_data_rm[complete.cases(aggregate_data_rm),]

save(aggregate_data_rm, file = "~/shared_space/ci3_analysis/josey_erc_strata/Data/aggregate_data_rm.RData")

rm(dead_personyear, confounders, aggregate_data_rm); gc()

### Create Strata Data

create_strata <- function(data, dual = c(0,1,2), race = c("all", "white", "black", "asian", "hispanic")) {
  
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
  } else if (race == "asian") {
    race0 <- 4
  } else if (race == "hispanic") {
    race0 <- 5
  } else {
    race0 <- c(1,2,3,4,5)
  }
  
  sub_data <- subset(data, race %in% race0 & dual %in% dual0)
  
  # Covariates and Outcomes
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  female = sub_data$female, dual = sub_data$dual, age_break = sub_data$age_break,
                  dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "age_break")]
  
  x <- data.table(zip = sub_data$zip, year = sub_data$year, model.matrix(~ ., data = sub_data[,zip_cov])[,-1])[
                    ,lapply(.SD, min), by = c("zip", "year")]
  
  return(list(w = w, x = x))
  
}

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2),
                         race = c("white", "black", "asian", "hispanic", "all"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/'
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'

## QD Strata

load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd$zip <- factor(aggregate_data_qd$zip)
aggregate_data_qd$year <- factor(aggregate_data_qd$year)
aggregate_data_qd$female <- as.numeric(aggregate_data_qd$female)
aggregate_data_qd$dual <- as.numeric(aggregate_data_qd$dual)
aggregate_data_qd$region <- factor(aggregate_data_qd$region)
aggregate_data_qd$age_break <- factor(aggregate_data_qd$age_break)
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
aggregate_data_rm$female <- as.numeric(aggregate_data_rm$female)
aggregate_data_rm$dual <- as.numeric(aggregate_data_rm$dual)
aggregate_data_rm$region <- factor(aggregate_data_rm$region)
aggregate_data_rm$age_break <- factor(aggregate_data_rm$age_break)
aggregate_data_rm$race <- factor(aggregate_data_rm$race)

lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(data = aggregate_data_rm, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  
})
