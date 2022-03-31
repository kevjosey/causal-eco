library(matrixStats)

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2),
                         race = c("white", "black", "all"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/'
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'

# ZIP Code Statistics

load(paste0(dir_data_qd, "2_all_qd.RData"))
sum(duplicated(paste0(new_data$x$zip, new_data$x$year))) # check

pop_dat <- aggregate(new_data$w$time_count, 
                     by = list(zip = new_data$w$zip, year = new_data$w$year), FUN = sum)
colnames(pop_dat)[3] <- "population"
zip_data <- merge(new_data$x, pop_dat, by = c("year", "zip"))

zip_mean <- apply(zip_data[,3:(ncol(zip_data) - 1)], 2, weighted.mean, w = zip_data$population)
zip_sd <- colWeightedSds(zip_data[,3:(ncol(zip_data) - 1)], w = zip_data$population)

out1 <- cbind(mean = round(zip_mean, 3), sd = round(zip_sd, 3))

# Individual Level Statistics

load(paste0(dir_data,"aggregate_data_qd.RData"))
ind_data <- select(aggregate_data_qd, c(sex, race, dual, dual, entry_age_break,
                                        followup_year, dead, time_count, region))

count(x = data.frame(ind_data), sex, wt = time_count)

