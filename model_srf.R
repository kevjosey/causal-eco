library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(splines)
library(splines2)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_std.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("high", "low","both"), race = c("white","black","hispanic","asian","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

### Fit Balancing Weights

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_DR/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights
create_strata <- function(aggregate_data,
                          dual = c("high","low","both"),
                          race = c("white","black","asian","hispanic","other","all"),
                          d = c(8,9,10,11,12)) {
  
  ## ZIP Code Covariates
  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  x <- data.table(zip = factor(aggregate_data$zip), year = factor(aggregate_data$year), region = factor(aggregate_data$region),
                  model.matrix(~ ., data = aggregate_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x <-  data.table(setDF(x)[,-which(colnames(x0) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x0$pm25)
  w <- data.table(zip = factor(aggregate_data$zip), year = factor(aggregate_data$year), region = factor(aggregate_data$region),
                  y = aggregate_data$dead, n = aggregate_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("zip", "year", "region"))
  
  # data format
  wx$ybar <- wx$y/wx$n
  wx$id <- paste(wx$zip, wx$year, sep = "-")
  
  ## Outcome models
  
  # estimate nuisance outcome model with splines
  covar <- subset(wx, select = c("year","region",zcov[-1]))
  inner <- paste(colnames(covar), collapse = " + ")
  nsa <- nsa(wx$pm25, intercept = TRUE, df = 7)
  
  w.mat <- cbind(nsa, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                   data = data.frame(aa = wx$pm25, covar)))
  w.mat <- w.mat[,-which(colnames(w.mat) %in% c("year2000", "year2000:aa"))]
  
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = wx$ybar, w.mat),
               weights = wx$n, family = quasipoisson())
  muhat <- predict

  # variance estimation
  vals <- sapply(d, function(d.tmp, ...) {
    
    shift <- ifelse(wx$pm25 > d.tmp, d.tmp, wx$pm25)
    
    ## Strata-specific Calibration Weights
    x.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n)))
    astar <- c(wx$pm25 - mean(wx$pm25))/var(wx$pm25)
    astar2 <- c((wx$pm25 - mean(wx$pm25))^2/var(wx$pm25) - 1)
    bstar <- c(shift - mean(wx$pm25))/var(wx$pm25)
    bstar2 <- c((shift - mean(wx$pm25))^2/var(wx$pm25) - 1)
    cmat <- cbind(x.mat*astar, astar2, x.mat)
    
    tm <- sum(wx$n)*c(apply(bstar*x.mat, 2, weighted.mean, w = wx$n),
                      weighted.mean(bstar2, w = wx$n),
                      apply(x.mat, 2, weighted.mean, w = wx$n))
    
    # fit calibration model
    ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = wx$n)
    wx$cal <- ipwmod$weights/ipwmod$base_weights
    
    # truncation
    wx$trunc <- wx$cal
    trunc0 <- quantile(wx$cal, 0.001)
    trunc1 <- quantile(wx$cal, 0.999)
    wx$trunc[wx$cal < trunc0] <- trunc0
    wx$trunc[wx$cal > trunc1] <- trunc1
    
    # index from target value
    nsa.tmp <- predict(nsa, newx = shift)
    
    # target sample values
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                         data = data.frame(aa = shift, covar)))
    w.tmp <- w.tmp[,-which(colnames(w.tmp) %in% c("year2000", "year2000:aa"))]
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))

    psi <- weighted.mean(wx$trunc*(wx$ybar - muhat) + mhat, w = wx$n)
    eif <- weighted.mean((wx$trunc*(wx$ybar - muhat) + mhat - psi)^2, w = wx$n)/nrow(wx)
    
    return(c(estimate = psi, variance = eif))
    
  })
  
  # extract estimates
  est_data <- data.frame(d = d, estimate = vals[1,], se = sqrt(vals[2,]))
  
  return(list(est_data = est_data, wx = wx))
  
}

# run it all
mclapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(aggregate_data = aggregate_data, dual = scenario$dual, race = scenario$race,
                            sex = scenario$sex, age_break = scenario$age_break)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_", 
                               scenario$sex, "_", scenario$age_break, ".RData"))
  
}, mc.cores = 25)