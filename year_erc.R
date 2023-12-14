library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(scam)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/model_erc.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/bootstrap.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white","black","hispanic","asian"),
                         years = c("2001-2004", "2005-2008", "2009-2012", "2013-2016"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$years <- as.character(scenarios$years)
a.vals <- seq(4, 16, length.out = 121)
nboot <- 1000

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights

# run it!
mclapply(seq(1,60,1), function(i, ...) {
  
  scenario <- scenarios[i,]
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  if (scenario$race == "white") {
    race0 <- 1
  } else if (scenario$race == "black") {
    race0 <- 2
  } else if (scenario$race == "asian") {
    race0 <- 4
  } else if (scenario$race == "hispanic") {
    race0 <- 5
  } else if (scenario$race == "other") {
    race0 <- 3
  } else {
    race0 <- c(0,1,2,3,4,5,6)
  }
  
  if (scenario$years == "2001-2004") {
    year0 <- c("2001", "2002", "2003", "2004")
  } else if (scenario$years == "2005-2008") {
    year0 <- c("2005", "2006", "2007", "2008")
  } else if (scenario$years == "2009-2012") {
    year0 <- c("2009", "2010", "2011", "2012")
  } else if (scenario$years == "2013-2016") {
    year0 <- c("2013", "2014", "2015", "2016")
  }
  
  aggregate_data$zip <- as.character(aggregate_data$zip)
  aggregate_data$year <- as.character(aggregate_data$year)
  aggregate_data$region <- as.character(aggregate_data$region)
  
  ## ZIP Code Covariates
  z <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
         "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  # zip-code-specific data
  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0 & year %in% year0)
  x0 <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                   model.matrix(~ ., data = sub_data[,z])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x0 <- data.table(setDF(x0)[,-which(colnames(x0) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x0$pm25)
  x <- merge(x0, data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                            m = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")], by = c("zip", "year", "region"))
  
  # individual-level predictors
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  # w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year),region = factor(sub_data$region),
  #                 dual = factor(sub_data$dual), race = factor(sub_data$race),
  #                 sex = factor(sub_data$sex), age_break = factor(sub_data$age_break),
  #                 y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region", "dual", "race", "sex", "age_break")]
  
  # issues with factor variables...
  x$zip <- as.character(x$zip)
  x$year <- as.character(x$year)
  x$region <- as.character(x$region)
  
  w$zip <- as.character(w$zip)
  w$year <- as.character(w$year)
  w$region <- as.character(w$region)
  
  # create id variable necessary for bootstrap
  x$id <- paste(x$zip, x$year, sep = "-")
  w$id <- paste(w$zip, w$year, sep = "-")
  
  # fit model on full data
  full_data <- model_erc(x = x, w = w, z = z, a.vals = a.vals, se.fit = TRUE)
  
  # fit model on bootstrap data
  # boot_data <- replicate(nboot, bootable(w = w, x = x, z = z, a.vals = a.vals))
  # 
  # boot_erc <- do.call(rbind, lapply(boot_data, function(iter, ...) iter$erc))
  # boot_ed <- do.call(rbind, lapply(boot_data, function(iter, ...) iter$ed))
  # colnames(boot_erc) <- colnames(boot_ed) <- a.vals
  
  new_data <- list(est_data = full_data$est_data, 
                   excess_death = full_data$excess_death, 
                   # boot_erc = boot_erc, boot_ed = boot_ed, 
                   wx = full_data$wx)
  
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_", scenario$years, ".RData"))
  
}, mc.cores = 15)
