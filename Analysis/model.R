library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(KernSmooth)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/kwls.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("high", "low","both"), race = c("white","black","hispanic","asian","all"),
                         sex = c("both"), age_break = c("[65,75)","[75,85)","[85,95)","[95,125)","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$sex <- as.character(scenarios$sex)
scenarios$age_break <- as.character(scenarios$age_break)

### Fit Balancing Weights

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights
create_strata <- function(aggregate_data,
                          dual = c("high","low","both"),
                          race = c("white","black","asian","hispanic","other","all"),
                          sex = c("male","female","both"),
                          age_break = c("[65,75)","[75,85)","[85,95)","[95,125)","all"),
                          a.vals = seq(2, 31, length.out = 146),
                          bw.seq = seq(0.1, 5, length.out = 25)) {
  
  if (age_break != "all") {
    age_break0 <- age_break
  } else {
    age_break0 <- c("[65,75)","[75,85)","[85,95)","[95,125)")
  }
  
  if (dual == "high") {
    dual0 <- 0
  } else if (dual == "low") {
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
  } else if (race == "other") {
    race0 <- 3
  } else {
    race0 <- c(0,1,2,3,4,5,6)
  }
  
  if (sex == "male") {
    sex0 <- 0
  } else if (sex == "female") {
    sex0 <- 1
  } else {
    sex0 <- c(0,1)
  }
    
  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0 & sex %in% sex0 & age_break %in% age_break0 )
  
  ## ZIP Code Covariates
  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  x <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                  model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  w <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # data format
  x$zip <- factor(x$zip)
  w$zip <- factor(w$zip)
  x$year <- factor(x$year)
  w$year <- factor(w$year)
  x$region <- factor(x$region)
  w$region <- factor(w$region)
  x$id <- paste(x$zip, x$year, sep = "-")
  w$id <- paste(w$zip, w$year, sep = "-")
  
  ## Strata-specific design matrix
  x.tmp <- subset(x, select = -c(zip, pm25, id))
  x.tmp$year <- factor(x.tmp$year)
  x.tmp$region <- factor(x.tmp$region)
  x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
  
  ## LM GPS
  # pimod <- lm(a ~ ., data = data.frame(a = x$pm25, x.tmp))
  # pimod.vals <- c(pimod$fitted.values)
  # pimod.sd <- sigma(pimod)
  # 
  # # nonparametric density
  # a.std <- c(x$pm25 - pimod.vals) / pimod.sd
  # dens <- density(a.std)
  # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  # 
  # # ipw numerator
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   std <- c(a.tmp - pimod.vals) / pimod.sd
  #   approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  # })
  # 
  # phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  # phat <- predict(smooth.spline(a.vals, phat.vals), x = x$pm25)$y
  # phat[phat < 0] <- .Machine$double.eps
  # 
  # x$ipw <- phat/pihat # LM GPS
  
  ## Strata-specific Calibration Weights
  x.mat <- model.matrix(~ ., data = data.frame(x.tmp))
  astar <- c(x$pm25 - mean(x$pm25))/var(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  
  # components for later
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  
  # fit calibration weights
  mod <- calibrate(cmat = cmat, target = tm)
  x$cal <- mod$weights
  
  # truncation
  x$trunc <- x$cal
  trunc0 <- quantile(x$cal, 0.001)
  trunc1 <- quantile(x$cal, 0.999)
  x$trunc[x$cal < trunc0] <- trunc0
  x$trunc[x$cal > trunc1] <- trunc1
    
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("zip", "year", "region", "id"))
    
  wx$psi <- wx$cal*wx$y/wx$n
  wx$psi_trunc <- wx$trunc*wx$y/wx$n
  
  risk.est <- sapply(bw.seq, risk.fn, a.vals = a.vals, psi = wx$psi_trunc, a = wx$pm25)
  bw <- c(bw.seq[which.min(risk.est)])

  target <- sapply(a.vals, kern_est_eco, a = wx$pm25, psi = wx$psi_trunc, weights = wx$n, bw = bw, se.fit = TRUE,
                   x = x.mat, astar = astar, astar2 = astar2, cmat = cmat, ipw = wx$trunc, eco = TRUE, sandwich = FALSE)
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals, estimate = target[1,], se = sqrt(target[2,]))
  
  return(list(est_data = est_data, individual_data = w, zip_data = x))
  
}

# run it all
mclapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(aggregate_data = aggregate_data, dual = scenario$dual, race = scenario$race,
                            sex = scenario$sex, age_break = scenario$age_break)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_", 
                               scenario$sex, "_", scenario$age_break, ".RData"))
  
}, mc.cores = 25)
