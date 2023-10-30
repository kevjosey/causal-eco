library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(scam)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_std.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(dual = c("high", "low","both"), race = c("white","black","hispanic","asian","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'
load(paste0(dir_data,"aggregate_data.RData"))

# Function for Fitting Weights
model_erc <- function(aggregate_data,
                      dual = c("high","low","both"),
                      race = c("white","black","asian","hispanic","other","all"),
                      a.vals = seq(4, 16, length.out = 121),
                      se.fit = TRUE) {
  
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
  x0 <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                   model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x0 <- data.table(setDF(x0)[,-which(colnames(x0) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x0$pm25)
  x <- merge(x0, data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                            m = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")], by = c("zip", "year", "region"))
  # w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
  #                 dual = factor(sub_data$dual), race = factor(sub_data$race),
  #                 sex = factor(sub_data$sex), age_break = factor(sub_data$age_break),
  #                 y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region", "dual", "race", "sex", "age_break")]
  w <- data.table(zip = factor(sub_data$zip), year = factor(sub_data$year), region = factor(sub_data$region),
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # data format
  w$ybar <- w$y/w$n
  x$id <- paste(x$zip, x$year, sep = "-")
  
  ## IPW Models
  
  # design matrix
  x.mat <- model.matrix(~ ., data = subset(setDF(x), select = -c(zip, id, pm25, m)))
  astar <- c(x$pm25 - mean(x$pm25))/var(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat*x$m))
  
  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = x$m)
  x$cal <- ipwmod$weights/ipwmod$base_weights
  
  # truncation
  x$trunc <- x$cal
  trunc0 <- quantile(x$cal, 0.01)
  trunc1 <- quantile(x$cal, 0.99)
  x$trunc[x$cal < trunc0] <- trunc0
  x$trunc[x$cal > trunc1] <- trunc1
  
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("zip", "year", "region"))
  # w.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, n, m, ybar, dual, race, sex, age_break, cal, trunc)))
  w.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, n, m, ybar, cal, trunc)))
  bstar <- c(wx$pm25 - mean(x$pm25))/var(x$pm25)
  bstar2 <- c((wx$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  dmat <- cbind(w.mat*bstar, bstar2, w.mat)

  ## Outcome models
  
  # estimate nuisance outcome model with splines
  # if (race == "all" & dual == "both") {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", "race", "dual", zcov[-1]))
  # } else if (race != "all" & dual == "both") {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", "dual", zcov[-1]))
  # } else if (race == "all" & dual != "both") {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", "race", zcov[-1]))
  # } else {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", zcov[-1]))
  # }
  
  covar <- subset(wx, select = c("year", "region", zcov[-1]))
  inner <- paste(colnames(covar), collapse = " + ")
  
  fmla <- formula(paste0("ybar ~ s(a, bs = 'cr', k = 7) + ", inner))
  mumod <- scam(fmla, data = data.frame(ybar = wx$ybar, a = wx$pm25, covar),
                weights = wx$n, family = quasipoisson())
  muhat <- predict(mumod, type = "response")
  wx.mat <- predict(mumod, type = "lpmatrix")
  
  # marginal values
  mhat.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    muhat.tmp <- predict(mumod, newdata = data.frame(a = rep(a.tmp, nrow(wx)), covar), type = "response")
    return(weighted.mean(muhat.tmp, w = wx$n))

  })

  mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = wx$pm25)$y
  
  target <- gam_std(a = wx$pm25, y = wx$ybar, family = mumod$family, weights = wx$n, 
                    a.vals = a.vals, x = w.mat, w = wx.mat, se.fit = TRUE,
                    ipw = wx$trunc, muhat = muhat, mhat = mhat, 
                    astar = bstar, astar2 = bstar2, cmat = dmat)
  
  # excess death estimation and variance
  dr.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Index of Target Value
    idx <- which.min(abs(a.vals - a.tmp))
    
    # Target Means and Design
    muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar), type = "response")
    wx.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar), type = "lpmatrix")
    
    # Excess Deaths
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))
    cut <- as.numeric(I(wx$pm25 > a.tmp))
    lambda <- sum(cut*(wx$y - wx$n*(muhat.tmp + wx$trunc*(wx$ybar - muhat))))
    
    # Naive Variance
    # Sig <- vcovHC(mumod, type = "HC3", sandwich = TRUE)
    # sig2 <- c(t(delta0) %*% wx.tmp %*% Sig %*% t(wx.tmp) %*% delta0)/(sum(wx$n)^2) + target$sig2[idx]
    
    # Robust Variance
    if (se.fit) {

      # extract estimates from target
      mu <- target$eta.vals[idx]
      Sig <- as.matrix(target$Sig)
      g.val <- c(target$g.vals[idx,])
      
      # Dimensions
      l <- ncol(wx.tmp)
      o <- length(g.val)
      
      # ERC Variance
      first <- c(t(delta) %*% wx.tmp %*% Sig[1:l,1:l] %*% t(wx.tmp) %*% delta)/(sum(wx$n)^2)
      second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
      sig2 <- first + second
      
      # Excess Death Variance
      var.tmp <- c(t(cut) %*% (delta*wx.tmp) %*% Sig[1:l,1:l] %*% c(t(delta*wx.tmp) %*% cut))
      omega2 <- sum(cut*c(wx$n*wx$trunc*(wx$ybar - muhat))^2) + var.tmp
      return(c(mu = mu, sig2 = sig2, lambda = lambda, omega2 = omega2))
      
    } else {
      
      mu <- target[idx]
      return(c(mu = mu, lambda = lambda))
      
    }
    
  })
  
  # extract estimates
  if (se.fit) {
    est_data <- data.frame(a.vals = a.vals, estimate = dr.vals[1,], se = sqrt(dr.vals[2,]))
    excess_death <- data.frame(a.vals = a.vals, estimate = dr.vals[3,], se = sqrt(dr.vals[4,])) 
  } else {
    est_data <- data.frame(a.vals = a.vals, estimate = dr.vals[1,])
    excess_death <- data.frame(a.vals = a.vals, estimate = dr.vals[2,]) 
  }
  
  return(list(est_data = est_data, excess_death = excess_death, wx = wx))
  
}

# run it all
mclapply(seq(2,15,2), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- model_erc(aggregate_data = aggregate_data, dual = scenario$dual, race = scenario$race)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
}, mc.cores = 8)
