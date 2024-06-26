library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(scam)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_ipw.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
set.seed(42)

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Eco/'
load(paste0(dir_data,"aggregate_data_rti.RData"))
a.vals <- seq(4, 16, length.out = 121)

## ZIP Code Covariates
z <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
          "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")

x <- data.table(zip = factor(aggregate_data$zip), year = factor(aggregate_data$year), region = factor(aggregate_data$region),
                model.matrix(~ ., data = aggregate_data[,z])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
x <-  data.table(setDF(x)[,-which(colnames(x) == "pm25")] %>% mutate_if(is.numeric, scale), pm25 = x$pm25)
w <- data.table(zip = factor(aggregate_data$zip), year = factor(aggregate_data$year), region = factor(aggregate_data$region),
                y = aggregate_data$dead, n = aggregate_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")]

# merge data components such as outcomes and exposures
wx <- merge(w, x, by = c("zip", "year", "region"))

# data format
wx$ybar <- wx$y/wx$n
wx$id <- paste(wx$zip, wx$year, sep = "-")

## Strata-specific Calibration Weights

# constraints and target margin
x.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n)))
astar <- c(wx$pm25 - mean(wx$pm25))/var(wx$pm25)
astar2 <- c((wx$pm25 - mean(wx$pm25))^2/var(wx$pm25) - 1)
cmat <- cbind(x.mat*astar, astar2, x.mat)
tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))

# fit calibration model
ipwmod <- calibrate(cmat = cmat, target = tm)
wx$cal <- ipwmod$weights/ipwmod$base_weights

# truncation
wx$trunc <- wx$cal
trunc0 <- quantile(wx$cal, 0.005)
trunc1 <- quantile(wx$cal, 0.995)
wx$trunc[wx$cal < trunc0] <- trunc0
wx$trunc[wx$cal > trunc1] <- trunc1

## Outcome models

# estimate nuisance model with gam
covar <- subset(wx, select = c("year","region",z[-1]))
inner <- paste(colnames(covar), collapse = " + ")
fmla <- formula(paste0("ybar ~ s(a, bs = 'tp') + ", inner))
mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
              data = data.frame(a = wx$pm25, ybar = wx$ybar, covar))
muhat <- predict(mumod, type = "response")
w.mat <- predict(mumod, type = "lpmatrix")

target <- gam_dr(a = wx$pm25, y = wx$ybar, family = mumod$family,
                  se.fit = TRUE, a.vals = a.vals, x = x.mat, w = w.mat,
                  ipw = wx$trunc, muhat = mumod$fitted.values, eco = TRUE,
                  astar = astar, astar2 = astar2)

## Variance Estimation

vals <- sapply(a.vals, function(a.tmp, ...) {

  # Index of Target Value
  idx <- which.min(abs(a.vals - a.tmp))
  
  # Target Means and Design
  w.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar), type = "lpmatrix")
  muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar), type = "response")
  
  # Excess Deaths
  cut <- as.numeric(I(wx$pm25 > a.tmp))
  delta <- c(mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))

  # Extract Target Values
  mu <- target$mu.vals[idx] + mean(muhat.tmp)
  Sig <- as.matrix(target$Sig)
  g.val <- c(target$g.vals[idx,])
  
  # Dimensions
  l <- ncol(w.tmp)
  o <- length(g.val)
  
  # ERC Variance
  first <-  c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(nrow(wx)^2)
  second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
  third <- c(t(g.val) %*% Sig[(l + 1):(l + o), 1:l] %*% t(w.tmp) %*% delta)/nrow(wx)
  sig2 <- first + second + 2*third
  
  return(c(mu = mu, sig2 = sig2))
  
})

# extract estimates
est_data <- data.frame(a.vals = a.vals, estimate = vals[1,], se = sqrt(vals[2,]))

# save estimates
new_data <- list(est_data = est_data, wx = wx)
save(new_data, file = paste0(dir_out, "Ecological.RData"))
