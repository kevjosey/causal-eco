library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(scam)
library(SuperLearner)
library(earth)
library(glmnet)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/erc_fun.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/bootstrap.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"),
                         region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST", "US"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$region <- as.character(scenarios$region)
a.vals <- seq(4, 16, length.out = 121)
nboot <- 200

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'
load(paste0(dir_data,"aggregate_data.RData"))
# aggregate_data <- aggregate_data[aggregate_data$pm25 <= 20,]
# aggregate_data <- subset(aggregate_data, year %in% c(2009,2010,2011,2012,2013,2014))

# run it!
lapply(seq(1,15,by = 2), function(i, ...) {
  
  scenario <- scenarios[i,]
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  if (scenario$region == "US") {
    region0 <- unique(aggregate_data$region)
  } else {
    region0 <- scenario$region
  }
  
  ## ZIP Code Covariates
  z <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", 
         "medhouseholdincome", "medianhousevalue", "poverty", "education", "pct_owner_occ", 
         "popdensity", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  # zip-code-specific data
  sub_data <- subset(aggregate_data, region %in% region0 & dual %in% dual0)
  
  x0 <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                   model.matrix(~ ., data = sub_data[,z])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x0 <- data.table(setDF(x0)[,-which(colnames(x0) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x0$pm25)
  x <- merge(x0, data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                            m = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")], by = c("zip", "year", "region"))
  
  # individual-level predictors
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # create id variable necessary for bootstrap
  x$id <- paste(x$zip, x$year, sep = "-")
  w$id <- paste(w$zip, w$year, sep = "-")
  
  # fit model on full data
  full_data <- erc_implement(x = x, w = w, z = z, a.vals = a.vals, 
                             se.fit = TRUE, boot = FALSE, 
                             state = scenario$region, dual = scenario$dual)
  
  # fit model on bootstrap data
  # boot_data <- replicate(nboot, bootable(w = w, x = x, z = z, a.vals = a.vals,
  #                       state = scenario$region, dual = scenario$dual), simplify = FALSE)
  # boot_erc <- do.call(rbind, lapply(boot_data, function(iter, ...) iter$erc))
  # colnames(boot_erc) <- colnames(boot_ed) <- a.vals
  
  new_data <- list(est_data = full_data$est_data, wx = full_data$wx)
  save(new_data, file = paste0(dir_out, scenario$region, "_", scenario$dual, ".RData"))
  
}, mc.cores = 8)
