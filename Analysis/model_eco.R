library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(splines)
library(splines2)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_eco.R')
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
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_DR/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights
create_strata <- function(aggregate_data,
                          dual = c("high","low","both"),
                          race = c("white","black","asian","hispanic","other","all"),
                          sex = c("male","female","both"),
                          age_break = c("[65,75)","[75,85)","[85,95)","[95,125)","all"),
                          a.vals = seq(4, 16, length.out = 121)) {
  
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
  
  x <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("zip", "year", "region"))
  
  # data format
  wx$ybar <- wx$y/wx$n
  wx$id <- paste(wx$zip, wx$year, sep = "-")
  
  ## Strata-specific design matrix
  x.tmp <- subset(wx, select = -c(zip, id, pm25, y, ybar, n))
  x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)
  
  ## Strata-specific Calibration Weights
  x.mat <- model.matrix(~ ., data = data.frame(x.tmp))
  astar <- c(wx$pm25 - mean(wx$pm25))/var(wx$pm25)
  astar2 <- c((wx$pm25 - mean(wx$pm25))^2/var(wx$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  
  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm)
  wx$cal <- ipwmod$weights/ipwmod$base_weights
  
  # truncation
  wx$trunc <- wx$cal
  trunc0 <- quantile(wx$cal, 0.001)
  trunc1 <- quantile(wx$cal, 0.999)
  wx$trunc[wx$cal < trunc0] <- trunc0
  wx$trunc[wx$cal > trunc1] <- trunc1
  
  ## Outcome models
  
  # estimate nuisance outcome model with splines
  covar <- subset(wx, select = c("year","region",zcov[-1])) %>% 
    mutate_if(is.numeric, scale)
  inner <- paste(colnames(covar), collapse = " + ")
  nsa <- nsa(wx$pm25, intercept = TRUE, df = 7)
  w.mat <- cbind(nsa, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                    data = data.frame(aa = wx$pm25, covar)))
  w.mat <- w.mat[,-which(colnames(w.mat) %in% c("year2000", "year2000:aa"))]
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = wx$ybar, w.mat),
               weights = wx$n, family = quasipoisson())
  
  target <- gam_eco(a = wx$pm25, y = wx$ybar, family = mumod$family, weights = wx$n, 
                    se.fit = TRUE, a.vals = a.vals, x = x.mat, w = w.mat,
                    ipw = wx$trunc, muhat = mumod$fitted.values, 
                    astar = astar, astar2 = astar2, cmat = cmat)
  
  # variance estimation
  vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # index from target value
    idx <- which.min(abs(a.vals - a.tmp))
    
    nsa.tmp <- predict(nsa, newx = rep(a.tmp, nrow(wx)))
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                          data = data.frame(aa = rep(a.tmp, nrow(wx)), covar)))
    w.tmp <- w.tmp[,-which(colnames(w.tmp) %in% c("year2000", "year2000:aa"))]
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    delta <- c(mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    
    # Naive Variance
    # Sig <- vcovHC(mumod, type = "HC3", sandwich = TRUE)
    # sig2 <- c(t(delta) %*% w.tmp %*% Sig %*% t(w.tmp) %*% delta)/(sum(wx$n)^2) + target$sig2[idx]
    
    # Robust Variance
    l <- ncol(w.tmp)
    o <- ncol(target$g.vals)
    g.val <- c(target$g.vals[idx,])
    Sig <- as.matrix(target$Sig)
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(nrow(wx)^2) +
      2*c(t(delta) %*% w.tmp %*% Sig[1:l, (l + 1):(l + o)] %*% g.val)/nrow(wx)
    sig2 <- first + c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    
    mu <- mean(mhat) + target$eta.vals[idx]
    
    return(c(mu = mu, sig2 = sig2))
    
  })
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals, estimate = vals[1,], se = sqrt(vals[2,]))
  
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