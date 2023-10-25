library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(splines)
library(splines2)
library(sandwich)
library(haldensify)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_std.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("high", "low","both"), race = c("white","black","hispanic","asian","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_SRF/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights
model_si <- function(aggregate_data,
                     dual = c("high","low","both"),
                     race = c("white","black","asian","hispanic","other","all"),
                     d = c(8,9,10,11,12)) {
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
  
  ## ZIP Code Covariates
  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  # Strata-specific outcomes and subset
  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0)
  x <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x <-  data.table(setDF(x)[,-which(colnames(x) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x$pm25)
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("zip", "year", "region"))
  
  # data format
  wx$ybar <- wx$y/wx$n
  wx$id <- paste(wx$zip, wx$year, sep = "-")
  
  ## Outcome models
  
  # estimate nuisance outcome model with splines
  covar <- subset(wx, select = c("year","region",zcov[-1]))
  inner <- paste(colnames(covar), collapse = " + ")
  nsa <- ns(wx$pm25, intercept = TRUE, df = 5)
  
  w.mat <- cbind(nsa, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                   data = data.frame(aa = wx$pm25, covar)))
  w.mat <- w.mat[,-which(colnames(w.mat) %in% c("year2000", "year2000:aa"))]
  
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = wx$ybar, w.mat),
               weights = wx$n, family = quasipoisson())
  muhat <- mumod$fitted.value
  
  # Fit baseline values
  x.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n)))
  astar <- c(wx$pm25 - mean(wx$pm25))/var(wx$pm25)
  astar2 <- c((wx$pm25 - mean(wx$pm25))^2/var(wx$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  
  psi <- weighted.mean(wx$ybar, w = wx$n)
  eif <- weighted.mean((wx$ybar - psi)^2, w = wx$n)/nrow(wx)
  
  # variance estimation
  vals <- sapply(d, function(d.tmp, ...) {
    
    shift <- ifelse(wx$pm25 > d.tmp, d.tmp, wx$pm25)
    
    ## Strata-specific Calibration Weights
    bstar <- c(shift - mean(shift))/var(shift)
    bstar2 <- c((shift - mean(shift))^2/var(shift) - 1)
    dmat <- cbind(x.mat*bstar, bstar2, x.mat)
    tm <- colSums(dmat*wx$n)
    
    # fit calibration model
    ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = wx$n)
    cal.tmp <- ipwmod$weights/ipwmod$base_weights
    
    # truncation
    trunc <- cal.tmp
    trunc0 <- quantile(cal.tmp, 0.001)
    trunc1 <- quantile(cal.tmp, 0.999)
    trunc[cal.tmp < trunc0] <- trunc0
    trunc[cal.tmp > trunc1] <- trunc1
    
    # index from target value
    nsa.tmp <- predict(nsa, newx = shift)
    
    # Target Stochastic Interventions
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                         data = data.frame(aa = shift, covar)))
    w.tmp <- w.tmp[,-which(colnames(w.tmp) %in% c("year2000", "year2000:aa"))]
    muhat.tmp <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    
    mu <- weighted.mean(cal.tmp*(wx$ybar - muhat) + muhat.tmp, w = wx$n)
    sig2 <- weighted.mean((cal.tmp*(wx$ybar - muhat) + muhat.tmp - psi)^2, w = wx$n)/nrow(wx)
    
    # Excess Deaths
    lambda <- sum(wx$y - wx$n*(cal.tmp*(wx$ybar - muhat) + muhat.tmp))
    omega2 <- sum(c(wx$n*(cal.tmp*(wx$ybar - muhat) + muhat.tmp - psi)^2))
    
    return(c(mu = mu, sig2 = sig2, lambda = lambda, omega2 = omega2))
    
  })
  
  # extract estimates
  est_data <- data.frame(d = c(0, d), estimate = c(psi, vals[1,]), se = c(sqrt(eif), sqrt(vals[2,])))
  excess_death <- data.frame(d = d, estimate = vals[3,], se = sqrt(vals[4,])) 
  
  return(list(est_data = est_data, excess_death = excess_death, wx = wx))
  
}

# run it all
mclapply(seq(1,15,2), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- model_si(aggregate_data = aggregate_data, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
}, mc.cores = 8)