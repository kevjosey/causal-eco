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

u.zip <- unique(aggregate_data$zip)
s.zip <- sample(u.zip, size = 1000)
sample_data <- subset(aggregate_data, zip %in% s.zip)

# Function for Fitting Weights
create_strata <- function(sample_data,
                          dual = c("high","low","both"),
                          race = c("white","black","asian","hispanic","other","all"),
                          a.vals = seq(4, 16, length.out = 121)) {
  
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
  sub_data <- subset(sample_data, race %in% race0 & dual %in% dual0)
  x <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                   model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  dual = factor(sub_data$dual), race = factor(sub_data$race),
                  sex = factor(sub_data$sex), age_break = factor(sub_data$age_break),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region", "dual", "race", "sex", "age_break")]
  
  # merge data components such as outcomes and exposures
  wx <- rbind(merge(w, x, by = c("zip", "year", "region")))

  # data format
  wx$ybar <- wx$y/wx$n
  wx$id <- paste(wx$zip, wx$year, sep = "-")
  
  ## Strata-specific design matrix
  
  if (race == "all" & dual == "both") {
    x.mat <- cbind(model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n))))
  } else if (race != "all" & dual == "both") {
    x.mat <- cbind(model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, race))))
  } else if (race == "all" & dual != "both") {
    x.mat <- cbind(model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, dual))))
  } else {
    x.mat <- cbind(model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, race, dual))))
  }
  
  astar <- c(wx$pm25 - mean(wx$pm25))/var(wx$pm25)
  astar2 <- c((wx$pm25 - mean(wx$pm25))^2/var(wx$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat*wx$n))
  
  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = wx$n)
  wx$cal <- ipwmod$weights/ipwmod$base_weights
  
  # truncation
  wx$trunc <- wx$cal
  trunc0 <- quantile(wx$cal, 0.001)
  trunc1 <- quantile(wx$cal, 0.999)
  wx$trunc[wx$cal < trunc0] <- trunc0
  wx$trunc[wx$cal > trunc1] <- trunc1

  ## Outcome models
  
  # estimate nuisance outcome model with gam
  # covar <- subset(wx, select = c("year","region",zcov[-1])) %>%
  #   mutate_if(is.numeric, scale)
  # inner <- paste(colnames(covar), collapse = " + ")
  # fmla <- as.formula(paste0("ybar ~ s(aa) + ", inner)) # , " + aa:(year + region)"))
  # mumod <- bam(fmla, data = data.frame(ybar = wx$ybar, aa = wx$pm25, wx),
  #              weights = wx$n, family = quasipoisson())
  # w.mat <- predict(mumod, type = "lpmatrix")
  
  # estimate nuisance outcome model with splines
  covar <- subset(wx, select = c("year", "region", zcov[-1]))
  
  if (race == "all" & dual == "both") {
    covar <- subset(wx, select = c("year", "region", "race", "dual", zcov[-1]))
  } else if (race != "all" & dual == "both") {
    covar <- subset(wx, select = c("year", "region", "dual", zcov[-1]))
  } else if (race == "all" & dual != "both") {
    covar <- subset(wx, select = c("year", "region", "race", zcov[-1]))
  } else {
    covar <- subset(wx, select = c("year", "region", zcov[-1]))
  }
  
  inner <- paste(colnames(covar), collapse = " + ")
  nsa <- ns(wx$pm25, intercept = TRUE, df = 7)
  
  w.mat <- cbind(predict(nsa, newx = wx$pm25), 
                 model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                              data = data.frame(aa = wx$pm25, covar)))
  w.mat <- w.mat[,-which(colnames(w.mat) %in% c("year2000", "year2000:aa"))]
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = wx$ybar, w.mat),
               weights = wx$n, family = quasipoisson())
  
  muhat <- predict(mumod, newdata = data.frame(w.mat), type = "response")
  target <- gam_std(a = wx$pm25, y = wx$ybar, family = mumod$family, weights = wx$n, 
                    se.fit = TRUE, a.vals = a.vals, s = s, x = x.mat, w = w.mat,
                    ipw = wx$trunc, muhat = muhat, astar = astar, astar2 = astar2, cmat = cmat)
  
  # variance estimation
  dr.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # index from target value
    idx <- which.min(abs(a.vals - a.tmp))
    nsa.tmp <- predict(nsa, newx = rep(a.tmp, nrow(wx)))
    g.val <- c(target$g.vals[idx,])
    
    # Outcome Prediction
    # w.tmp <- predict(mumod, type = "lpmatrix", newdata = data.frame(aa = a.tmp, covar),
    #                  newdata.guaranteed = TRUE, block.size = nrow(wx))
    
    # target sample values
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                data = data.frame(aa = rep(a.tmp, nrow(wx)), covar)))
    w.tmp <- w.tmp[,-which(colnames(w.tmp) %in% c("year2000", "year2000:aa"))]
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    
    # Naive Variance
    # Sig <- vcovHC(mumod, type = "HC3", sandwich = TRUE)
    # sig2 <- c(t(delta0) %*% w.tmp %*% Sig %*% t(w.tmp) %*% delta0)/(sum(wx$n)^2) + target$sig2[idx]

    # Robust Variance
    l <- ncol(w.tmp0)
    o <- length(g.val)
    
    Sig <- as.matrix(target$Sig)
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(wx$n)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2 <- first + second
    
    mu <- weighted.mean(mhat, w = wx$n) + target$eta.vals[idx]
    
    # Excess Deaths
    cut <- as.numeric(I(wx$pm25 > a.tmp))
    lambda <- sum(cut*(wx$y - wx$n*(mhat + wx$trunc(wx$ybar - muhat))))
    var.tmp <- diag((delta*w.tmp) %*% Sig[1:l,1:l] %*% t(delta*w.tmp))
    omega2 <- sum(cut*c(wx$n)^2*wx$trunc*(wx$ybar - muhat)^2) + c(t(cut) %*% var.tmp %*% cut)
      
    return(c(mu = mu, sig2 = sig2, lambda = lambda, omega2 = omega2))
    # return(c(mu = mu, sig2 = sig2))
    
  })
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals, estimate = dr.vals[1,], se = sqrt(dr.vals[2,]))
  excess_death <- data.frame(a.vals = a.vals, estimate = dr.vals[3,], se = sqrt(dr.vals[4,])) 
  
  return(list(est_data = est_data, excess_death = excess_death, wx = wx))
  
}

# run it all
mclapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(sample_data = sample_data, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
}, mc.cores = 15)
