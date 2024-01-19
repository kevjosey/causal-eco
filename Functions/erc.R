model_erc <- function(x, w, z, se.fit = TRUE, boot = FALSE, 
                      a.vals = seq(4, 16, length.out = 121), 
                      region = "US", race = "all", dual = "both") {
  
  # data format
  w$ybar <- w$y/w$n
  
  # merge data components such as outcomes and exposures
  
  if (boot) {
    wx <- merge(w, x, by = c("boot.id", "id", "zip", "year", "region"))
    df.tmp <- subset(setDF(wx), select = -c(zip, id, boot.id, pm25, y, ybar, n, m,
                                            dual, race, sex, age_break))
    # df.tmp <- subset(setDF(wx), select = -c(zip, id, boot.id, pm25, y, ybar, n, m, ipw, trunc))
  } else {
    wx <- merge(w, x, by = c("id", "zip", "year", "region"))
    df.tmp <- subset(setDF(wx), select = -c(zip, id, pm25, y, ybar,
                                            n, m, dual, race, sex, age_break))
    # df.tmp <- subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, m, ipw, trunc))
  }

  ## Outcome Models
  
  # estimate nuisance outcome model with GAM
  if (race == "all" & dual == "both") {
    covar <- subset(wx, select = c("sex", "age_break", "race", "dual", z[-1]))
  } else if (race != "all" & dual == "both") {
    covar <- subset(wx, select = c("sex", "age_break", "dual", z[-1]))
  } else if (race == "all" & dual != "both") {
    covar <- subset(wx, select = c("sex", "age_break", "race", z[-1]))
  } else {
    covar <- subset(wx, select = c("sex", "age_break", z[-1]))
  }
  
  # covar <- subset(wx, select = z[-1])
  inner <- paste(colnames(covar), collapse = " + ")
  
  if (region == "US") {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'mpi') + year*region +", inner))
  } else if (region != "US") {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'mpi') + year +", inner))
  }

  mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
                data = data.frame(a = wx$pm25, ybar = wx$ybar, covar,
                                  year = wx$year, region = wx$region))

  # excess death estimation and variance
  out.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Index of Target Value
    idx <- which.min(abs(a.vals - a.tmp))
    
    muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar, year = wx$year, region = wx$region), type = "response")
    wx.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, year = wx$year, region = wx$region, covar), type = "lpmatrix")
    
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))
    cut <- as.numeric(I(wx$pm25 > a.tmp))
  
    mu <- weighted.mean(muhat.tmp, w = wx$n)
    lambda <- sum(cut*(wx$y - wx$n*muhat.tmp))
    
    if (se.fit) {
      sig2 <- c(t(delta) %*% wx.tmp %*% mumod$Ve %*% t(wx.tmp) %*% delta)/(sum(wx$n)^2)
      omega2 <- c(t(cut) %*% (delta*wx.tmp) %*% mumod$Ve %*% c(t(delta*wx.tmp) %*% cut))
      return(c(mu = mu, sig2 = sig2, lambda = lambda, omega2 = omega2))
    } else {
      return(c(mu = mu, lambda = lambda))
    }
    
  })
  
  # extract estimates
  if (se.fit) {
    est_data <- data.frame(a.vals = a.vals, estimate = out.vals[1,], se = sqrt(out.vals[2,]))
    excess_death <- data.frame(a.vals = a.vals, estimate = out.vals[3,], se = sqrt(out.vals[4,])) 
  } else {
    est_data <- data.frame(a.vals = a.vals, estimate = out.vals[1,])
    excess_death <- data.frame(a.vals = a.vals, estimate = out.vals[2,]) 
  }
  
  return(list(est_data = est_data, excess_death = excess_death, wx = wx))
  
}