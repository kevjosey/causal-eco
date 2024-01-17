model_erc <- function(x, w, z, se.fit = TRUE, boot = FALSE, 
                      a.vals = seq(4, 16, length.out = 121), 
                      region = "US", race = "all") {
  
  # data format
  w$ybar <- w$y/w$n
  
  # exposure design matrix
  if (boot) {
    
    if (region != "US") {
      x.mat <- model.matrix(~ ., data = subset(setDF(x), select = -c(region, zip, id, boot.id, pm25, m)))
    } else {
      x.mat <- model.matrix(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, boot.id, pm25, m)))
    }
    
  } else {
    
    if (region != "US") {
      x.mat <- model.matrix(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
    } else {
      x.mat <- model.matrix(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
    }

  }
  
  ## Linear IPW Model
  # pimod <- lm(a ~ 0 + ., data = data.frame(a = x$pm25, x.mat))
  # pimod.vals <- c(pimod$fitted.values)
  # pimod.sd <- sigma(pimod)

  # nonparametric density / denominator
  # a.std <- c(x$pm25 - pimod.vals) / pimod.sd
  # dens <- density(a.std)
  # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd

  # ipw numerator
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   std <- c(a.tmp - pimod.vals) / pimod.sd
  #   approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  # })
  # 
  # phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  # phat <- predict(smooth.spline(a.vals, phat.vals), x = x$pm25)$y
  # phat[phat < 0] <- .Machine$double.eps

  # combine numerator and denominator
  # x$ipw <- phat/pihat # LM GPS
  # x$ipw <- ifelse(is.na(x$ipw), 1, x$ipw)
  
  ## Calibration Models
  astar <- c(x$pm25 - mean(x$pm25))/sd(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))

  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm)
  x$ipw <- ipwmod$weights/ipwmod$base_weights
  if (!ipwmod$converged)
    x$ipw <- 1
  
  # truncation
  x$trunc <- x$ipw
  trunc0 <- quantile(x$ipw, 0.005)
  trunc1 <- quantile(x$ipw, 0.995)
  x$trunc[x$ipw < trunc0] <- trunc0
  x$trunc[x$ipw > trunc1] <- trunc1
  
  # merge data components such as outcomes and exposures
  
  if (boot) {
    wx <- merge(w, x, by = c("boot.id", "id", "zip", "year", "region"))
    # df.tmp <- subset(setDF(wx), select = -c(zip, id, boot.id, pm25, y, ybar,
    #                                         n, m, ipw, trunc, 
    #                                         dual, race, sex, age_break))
    df.tmp <- subset(setDF(wx), select = -c(zip, id, boot.id, pm25, y, ybar, n, m, ipw, trunc))
  } else {
    wx <- merge(w, x, by = c("id", "zip", "year", "region"))
    # df.tmp <- subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, 
    #                                         n, m, ipw, trunc,
    #                                         dual, race, sex, age_break))
    df.tmp <- subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, m, ipw, trunc))
  }
  
  if (region != "US") 
    df.tmp <- subset(df.tmp, select = -c(region))
  
  w.mat <- model.matrix(~ ., data = df.tmp)
  
  # calibration components
  bstar <- c(wx$pm25 - mean(x$pm25))/sd(x$pm25)
  bstar2 <- c((wx$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  
  ## Outcome Models
  
  # estimate nuisance outcome model with GAM
  # if (race == "all") {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", "race", "dual", z[-1]))
  # } else if (race != "all") {
  #   covar <- subset(wx, select = c("year", "region", "sex", "age_break", "dual", z[-1]))
  # }
  # 
  # if (region != "US")
  #   covar <- subset(covar, select = -c(region))
  
  covar <- subset(wx, select = z[-1])
  inner <- paste(colnames(covar), collapse = " + ")
  
  if (region == "US") {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'mpi') + year + region + year:region +", inner))
  } else if (region != "US") {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'mpi') + year +", inner))
  }

  mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
                data = data.frame(a = wx$pm25, ybar = wx$ybar, covar,
                                  year = wx$year, region = wx$region))
  
  muhat <- predict(mumod, type = "response")
  wx.mat <- predict(mumod, type = "lpmatrix")
  
  # marginal values
  # mhat.vals <- sapply(a.vals, function(a.tmp, ...) {
  #   muhat.tmp <- predict(mumod, newdata = data.frame(a = rep(a.tmp, nrow(wx)), covar), type = "response")
  #   return(weighted.mean(muhat.tmp, w = wx$n))
  # })
  # 
  # mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), wx$pm25)$y
  
  target <- gam_std(a = wx$pm25, y = wx$ybar, family = mumod$family, weights = wx$n, 
                    a.vals = a.vals, x = w.mat, w = wx.mat, se.fit = se.fit,
                    ipw = wx$trunc, muhat = muhat, astar = bstar, astar2 = bstar2)
  
  # excess death estimation and variance
  dr.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Index of Target Value
    idx <- which.min(abs(a.vals - a.tmp))
    
    # Target Means and Design
    muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, year = wx$year, region = wx$region, covar), type = "response")
    wx.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, year = wx$year, region = wx$region, covar), type = "lpmatrix")
    
    # Excess Deaths
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))
    cut <- as.numeric(I(wx$pm25 > a.tmp))
    lambda <- sum(cut*(wx$y - wx$n*(muhat.tmp + wx$trunc*(wx$ybar - muhat))))
    
    # Naive Variance
    # Sig <- vcovHC(mumod, type = "HC3", sandwich = TRUE)
    # sig2 <- c(t(delta0) %*% wx.tmp %*% Sig %*% t(wx.tmp) %*% delta0)/(sum(wx$n)^2) + target$sig2[idx]
    
    # Robust Variance
    if (se.fit) {
      
      mu <- target$mu.vals[idx] + weighted.mean(muhat.tmp, w = wx$n)
      # mu <- target$mu.vals[idx]
      
      # Extract Values from Target
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
      
      mu <- target[idx] + weighted.mean(muhat.tmp, w = wx$n)
      # mu <- target[idx]
      
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