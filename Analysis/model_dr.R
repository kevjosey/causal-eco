library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(splines)
library(splines2)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_dr.R')
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
  
  x <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                   model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  w <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  x0 <- data.table(zip = aggregate_data$zip, year = aggregate_data$year, region = aggregate_data$region,
                   model.matrix(~ ., data = aggregate_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  w0 <- data.table(zip = aggregate_data$zip, year = aggregate_data$year, region = aggregate_data$region,
                   y = aggregate_data$dead, n = aggregate_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]
  
  # merge data components such as outcomes and exposures
  wx <- rbind(merge(w, x, by = c("zip", "year", "region")),
              merge(w0, x0, by = c("zip", "year", "region")))
  s <- rep(c(1,0), times = c(nrow(w), nrow(w0)))

  # data format
  wx$ybar <- wx$y/wx$n
  wx$zip <- factor(wx$zip)
  wx$year <- factor(wx$year)
  wx$region <- factor(wx$region)
  wx$id <- paste(wx$zip, wx$year, sep = "-")
  
  ## Strata-specific design matrix
  x.tmp <- subset(wx, select = -c(zip, id, pm25, y, ybar, n)) %>%
    mutate_if(is.numeric, scale)
  
  ## Strata-specific Calibration Weights
  x.mat <- cbind(model.matrix(~ ., data = data.frame(x.tmp)))
  astar <- c(wx$pm25 - mean(wx$pm25[s == 1]))/var(wx$pm25[s == 1])
  astar2 <- c((wx$pm25 - mean(wx$pm25[s == 1]))^2/var(wx$pm25[s == 1]) - 1)
  cmat <- cbind(s*x.mat*astar, s*astar2, s*x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), c(t(x.mat) %*% c((1 - s)*wx$n))*(sum(s)/sum(1 - s)))
  
  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = wx$n)
  wx$cal <- ipwmod$weights/ipwmod$base_weights
  
  # truncation
  wx$trunc <- wx$cal
  trunc0 <- quantile(wx$cal[s == 1], 0.001)
  trunc1 <- quantile(wx$cal[s == 1], 0.999)
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
  covar <- subset(wx, select = c("year", "region", zcov[-1])) %>% 
    mutate_if(is.numeric, scale)
  inner <- paste(colnames(covar), collapse = " + ")
  nsa <- ns(wx$pm25[s == 1], intercept = TRUE, df = 7)
  w.mat <- cbind(predict(nsa, newx = wx$pm25), 
                 model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                              data = data.frame(aa = wx$pm25, covar)))
  w.mat <- w.mat[,-which(colnames(w.mat) %in% c("year2000", "year2000:aa"))]
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = wx$ybar, w.mat),
               weights = wx$n, family = quasipoisson(), subset = c(s == 1))
  muhat <- predict(mumod, newdata = data.frame(w.mat), type = "response")
  
  target <- gam_dr(a = wx$pm25, y = wx$ybar, family = mumod$family, weights = wx$n, 
                    se.fit = TRUE, a.vals = a.vals, s = s, x = x.mat, w = w.mat,
                    ipw = wx$trunc, muhat = muhat, astar = astar, astar2 = astar2, cmat = cmat)
  
  # variance estimation
  vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Outcome Prediction
    # w.tmp <- predict(mumod, type = "lpmatrix", newdata = data.frame(aa = a.tmp, covar),
    #                  newdata.guaranteed = TRUE, block.size = nrow(wx))
    
    nsa.tmp <- predict(nsa, newx = rep(a.tmp, sum(1 - s)))
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ 0 +", inner, "+ aa:(year + region)")), 
                                         data = subset(data.frame(aa = rep(a.tmp, nrow(wx)), covar), c(s == 0))))
    w.tmp <- w.tmp[,-which(colnames(w.tmp) %in% c("year2000", "year2000:aa"))]
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    delta <- c(wx$n[s == 0]*mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    
    # index from target value
    idx <- which.min(abs(a.vals - a.tmp))
    
    # Naive Variance
    # Sig <- vcovHC(mumod, type = "HC3", sandwich = TRUE)
    # sig2 <- c(t(delta) %*% w.tmp %*% Sig %*% t(w.tmp) %*% delta)/(sum(wx$n)^2) + target$sig2[idx]

    # Robust Variance
    l <- ncol(w.tmp)
    o <- ncol(target$g.vals)
    g.val <- c(target$g.vals[idx,])
    Sig <- as.matrix(target$Sig)
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(wx$n[s == 0])^2) +
      2*c(t(delta) %*% w.tmp %*% Sig[1:l, (l + 1):(l + o)] %*% g.val)/sum(wx$n[s == 0])
    sig2 <- first + c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    
    mu <- weighted.mean(mhat, w = wx$n[s == 0]) + target$eta.vals[idx]
    
    return(c(mu = mu, sig2 = sig2))
    
  })
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals, estimate = vals[1,], se = sqrt(vals[2,]))
  
  return(list(est_data = est_data, wx1 = wx[s == 1,], wx0 = wx[s == 0,]))
  
}

# run it all
mclapply(1:nrow(scenarios), function(i, ...) {
  
  scenario <- scenarios[i,]
  new_data <- create_strata(aggregate_data = aggregate_data,
                            dual = scenario$dual, race = scenario$race,
                            sex = scenario$sex, age_break = scenario$age_break)
  save(new_data, file = paste0(dir_out, scenario$dual, "_", scenario$race, "_", 
                               scenario$sex, "_", scenario$age_break, ".RData"))
  
}, mc.cores = 25)
