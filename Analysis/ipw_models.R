library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(xgboost)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/pm-risk/Functions/calibrate.R')
set.seed(42)

# scenarios
# scenarios
scenarios <- expand.grid(dual = c("high", "low"), race = c("white","black","hispanic","asian","other"),
                         sex = c("female","male"), age_break = c("[65,75)","[75,85)","[85,95)","[95,125)"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$sex <- as.numeric(scenarios$sex)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(2, 31, length.out = 146)

### Fit Balancing Weights

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
load(paste0(dir_data,"aggregate_data_qd.RData"))

# Function for Fitting Weights
create_strata <- function(aggregate_data,
                          dual = c("high","low"),
                          race = c("white", "black", "asian", "hispanic", "other"),
                          sex = c("male", "female"),
                          age_break = c("[65,75)","[75,85)","[85,95)","[95,125)")) {
  
  if (dual == "high") {
    dual0 <- 0
  } else if (dual == "low") {
    dual0 <- 1
  }
  
  if (race == "white") {
    race0 <- 1
  } else if (race == "black") {
    race0 <- 2
  } else if (race == "asian") {
    race0 <- 4
  } else if (race == "hispanic") {
    race0 <- 5
  } else if (race == "other") {
    race0 <- 3
  }
  
  if (sex == "male") {
    dual0 <- 0
  } else if (dual == "female") {
    dual0 <- 1
  }
    
  age_break0 <- age_break
  
  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0 & sex %in% sex0 & age_break %in% age_break0 )
  
  # Outcome and Person-Years At-Risk
  w <- data.table(zip = sub_data$zip, year = sub_data$year, race = sub_data$race,
                  female = sub_data$female, dual = sub_data$dual, age_break = sub_data$age_break,
                  dead = sub_data$dead, time_count = sub_data$time_count)[
                    ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "entry_age_break", "followup_year")]
  
  # ZIP Code Covariates
  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")
  x <- data.table(zip = sub_data$zip, year = sub_data$year,
                  model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
  x.all <- data.table(zip = aggregate_data$zip, year = sub_data$year,
                      model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
  
  x.tmp <- subset(x, select = -c(zip, pm25))
  x.tmp$year <- factor(x.tmp$year)
  x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
  x.all.tmp <- subset(x.all, select = -c(zip, pm25))
  x.all.tmp$year <- factor(x.all.tmp$year)
  x.all.tmp <- x.all.tmp %>% mutate_if(is.numeric, scale)
  
  # target margin
  tm <- colMeans(model.matrix(~ ., data = data.frame(x.all.tmp)))
  
  ## LM GPS
  
  pimod <- lm(a ~ ., data = data.frame(a = x$pm25, x.tmp))
  pimod.vals <- c(pimod$fitted.values)
  pimod.sd <- sigma(pimod)
  
  # nonparametric density
  a.std <- c(x$pm25 - pimod.vals) / pimod.sd
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # ipw numerator
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / pimod.sd
    approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  })
  
  phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  phat <- predict(smooth.spline(a.vals, phat.vals), x = x$pm25)$y
  phat[phat < 0] <- .Machine$double.eps
  
  x$ipw <- phat/pihat # LM GPS
  
  ## Calibration Weights
  
  n <- nrow(x.tmp)
  x.mat <- model.matrix(~ ., data = data.frame(x.tmp))
  astar <- c(x$pm25 - mean(x$pm25))/var(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  mod <- calibrate(cmat = cbind(x.mat*astar, astar2, x.mat), 
                   target = c(rep(0, ncol(x.mat) + 1), n*tm))
  
  x$cal <- mod$weights # CALIBRATION
  
  # truncation
  trunc0 <- quantile(x$ipw, 0.005)
  trunc1 <- quantile(x$ipw, 0.995)
  x$ipw[x$ipw < trunc0] <- trunc0
  x$ipw[x$ipw > trunc1] <- trunc1
  
  x$cal_trunc <- x$cal # TRUNCATED CALIBRATION
  
  trunc0 <- quantile(x$cal, 0.005)
  trunc1 <- quantile(x$cal, 0.995)
  x$cal_trunc[x$cal < trunc0] <- trunc0
  x$cal_trunc[x$cal > trunc1] <- trunc1
  
  # format variables
  w$zip <- factor(w$zip)
  w$year <- factor(w$year)
  w$female <- as.numeric(w$female)
  w$race <- factor(w$race)
  w$dual <- as.numeric(w$dual)
  w$age_break <- factor(w$age_break)
  
  x$zip <- factor(x$zip)
  x$year <- factor(x$year)
  x$id <- paste(x$zip, x$year, sep = "-")
  
  return(list(w = w, x = x, phat.vals = phat.vals))
  
}

# collate
lapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(aggregate_data = aggregate_data, dual = scenario$dual, race = scenario$race,
                            sex = scenario$sex, age_break = scenario$age_break)
  
  print(i)
  save(new_data, file = paste0(dir_data, "qd2/", scenario$dual, "_", scenario$race, "_", 
                               scenario$sex, "_", scenario$age_break, ".RData"))
  
})
