library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(scam)
library(SuperLearner)
library(earth)
library(glmnet)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/srf_fun.R')
set.seed(42)

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_SRF/'
load(paste0(dir_data,"aggregate_data.RData"))
# aggregate_data <- subset(aggregate_data, pm25 <= 20)
# aggregate_data <- subset(aggregate_data, year %in% c(2009,2010,2011,2012,2013,2014))

# scenarios
states <- as.character(unique(aggregate_data$statecode))
scenarios <- expand.grid(state = states, dual = c("both"))
scenarios$state <- as.character(scenarios$state)
scenarios$dual <- as.character(scenarios$dual)

scenarios <- subset(scenarios, state %in% c("AR", "GA", "FL", "OK"))

# run it!
mclapply(seq(1, nrow(scenarios), by = 2), function(i, ...) {
  
  scenario <- scenarios[i,]
  
  if (scenario$state == "US") {
    state0 <- states
  } else {
    state0 <- scenario$state
  }
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  ## ZIP Code Covariates
  z <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
         "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  # zip-code-specific data
  sub_data <- subset(aggregate_data, statecode %in% state0 & dual %in% dual0)
  x0 <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                   model.matrix(~ ., data = sub_data[,z])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x0 <- data.table(setDF(x0)[,-which(colnames(x0) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x0$pm25)
  x <- merge(x0, data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                            m = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")], by = c("zip", "year", "region"))
  
  # individual-level predictors
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  dual = factor(sub_data$dual), race = factor(sub_data$race),
                  sex = factor(sub_data$sex), age_break = factor(sub_data$age_break),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region", "dual", "race", "sex", "age_break")]
  
  # create id variable necessary for bootstrap
  x$id <- paste(x$zip, x$year, sep = "-")
  w$id <- paste(w$zip, w$year, sep = "-")
  delta <- c(6,7,8,9,10,11,12)
  
  # fit model on full data
  full_data <- lapply(delta, srf_implement, x = x, w = w, z = z,
                        state = scenario$state, dual = scenario$dual)
  
  est_data <- data.frame(delta = delta, do.call(rbind, lapply(full_data, function(arg) cbind(est = arg$theta, se = sqrt(arg$omega2)))))
  save(est_data, file = paste0(dir_out, scenario$state, "_", scenario$dual, ".RData"))
  
}, mc.cores = 5)
