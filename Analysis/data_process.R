# Generate all required data from different data sources (see Table S1).
# library(fst)
# library(data.table)
# library(tidyr)
# library(dplyr)
# library(parallel)
# library(foreign)

### Big Data Cleaning (BEWARE!)

## Qian Di

# f <- list.files("/n/dominici_nsaph_l3/Lab/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
#                 pattern = "\\.fst",
#                 full.names = TRUE)

# myvars <- c("qid", "year","zip","sex","race","age","dual","entry_age_break","statecode",
#             "followup_year","followup_year_plus_one","dead","pm25_ensemble",
#             "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
#             "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")
# 
# national_merged2016 <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
# national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)

# NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")
# SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
# MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
# WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

# creates region
# national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
#                                      ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
#                                             ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
#                                                    ifelse(national_merged2016$state %in% WEST, "WEST", NA))))
# 
# national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:27)]) ,]
# save(national_merged2016, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")

## Randall Martin

# pm_rm <- data.frame()

# for(i in 2000:2016){
#   temp <- read.csv(paste0("/n/dominici_nsaph_l3/Lab/ci3_exposure/pm25/whole_us/annual/zipcode/rm_predictions/ben_2019_10_29/data_acag_pm25_zip-year/zip_pm25_",i, ".csv"))
#   temp$YEAR <- rep(i, nrow(temp))
#   pm_rm<-rbind(pm_rm, temp)
# }

# rm(temp); gc()

# pm_rm <- subset(pm_rm, !is.na(pm_rm$ZIP))
# pm_rm$ZIP <- sprintf("%05d", pm_rm$ZIP)
# pm_rm <- subset(pm_rm, pm_rm$ZIP %in% unique(national_merged2016$zip))
# pm_rm <- as.data.table(pm_rm)

# national_merged2016 <- as.data.table(national_merged2016)
# setkey(national_merged2016, zip, year)
# setkey(pm_rm, ZIP, YEAR)

# national_merged2016[pm_rm, pm25:=i.pm25]
# national_merged2016 <- national_merged2016[, c(1:11,28,13:27,12)] # reorder to match QD
# rm(national_merged2016)
# national_merged2016 <- subset(national_merged2016, select=-pm25_ensemble)
# save(national_merged2016, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016_rm.RData")

### Build Clean Aggregate Data Set

## QD Outcomes, Exposures, Covariates, and Offsets

load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016.RData")
national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,75,85,95,125), right = FALSE)
national_merged2016$sex <- national_merged2016$sex - 1
colnames(national_merged2016)[12] <- "pm25"

dead_personyear <- aggregate(data.frame(dead = national_merged2016$dead,
                                        time_count = national_merged2016$time_count),
                             by=list(zip = national_merged2016$zip,
                                     year = national_merged2016$year,
                                     sex = national_merged2016$sex,
                                     race = national_merged2016$race,
                                     dual = national_merged2016$dual,
                                     age_break = national_merged2016$age_break),
                             FUN=sum)

new_data <- national_merged2016 %>% distinct(zip, year, sex, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016, new_data); gc()

aggregate_data <- merge(dead_personyear, confounders, by = c("zip", "year","sex","race","dual","age_break"))
aggregate_data <- aggregate_data[complete.cases(aggregate_data),]

aggregate_data$zip <- factor(aggregate_data$zip)
aggregate_data$year <- factor(aggregate_data$year)
aggregate_data$region <- factor(aggregate_data$region)
aggregate_data$sex <- as.numeric(aggregate_data$sex)
aggregate_data$dual <- as.numeric(aggregate_data$dual)
aggregate_data$age_break <- factor(aggregate_data$age_break)
aggregate_data$race <- factor(aggregate_data$race)

save(aggregate_data, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data.RData")

rm(dead_personyear, confounders, aggregate_data); gc()

## RM Outcomes, Exposures, Covariates, and Offsets

load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/national_merged2016_rm.RData")
national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,75,85,95,125), right = FALSE)
national_merged2016$sex <- national_merged2016$sex - 1

dead_personyear <- aggregate(data.frame(dead = national_merged2016$dead,
                                        time_count = national_merged2016$time_count),
                             by=list(zip = national_merged2016$zip,
                                     year = national_merged2016$year,
                                     sex = national_merged2016$sex,
                                     race = national_merged2016$race,
                                     dual = national_merged2016$dual,
                                     age_break = national_merged2016$age_break),
                             FUN=sum)

new_data <- national_merged2016 %>% distinct(zip, year, sex, race, dual, age_break, .keep_all = TRUE)
confounders <- new_data[,c(1,2,4,6,12:27,29,30)]

rm(national_merged2016, new_data); gc()

aggregate_data <- merge(dead_personyear, confounders, by = c("zip", "year","sex","race","dual","age_break"))
aggregate_data <- aggregate_data[complete.cases(aggregate_data),]

aggregate_data$zip <- factor(aggregate_data$zip)
aggregate_data$year <- factor(aggregate_data$year)
aggregate_data$region <- factor(aggregate_data$region)
aggregate_data$sex <- as.numeric(aggregate_data$sex)
aggregate_data$dual <- as.numeric(aggregate_data$dual)
aggregate_data$age_break <- factor(aggregate_data$age_break)
aggregate_data$race <- factor(aggregate_data$race)

save(aggregate_data, file = "/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data_rm.RData")

rm(dead_personyear, confounders, aggregate_data); gc()
