erc_implement <- function(x, w, z, se.fit = TRUE, boot = FALSE, 
                          a.vals = seq(4, 16, length.out = 121), 
                          sl.lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth"),
                          state = "US", dual = "both") {
  
  # data format
  w$ybar <- w$y/w$n
  
  if (boot & state != "US") {
    x.mat <- model.frame(~ ., data = subset(setDF(x), select = -c(region, zip, id, boot.id, pm25, m)))
  } else if (boot & state == "US") {
    x.mat <- model.frame(~ ., data = subset(setDF(x), select = -c(zip, id, boot.id, pm25, m)))
  } else if (!boot & state != "US") {
    x.mat <- model.frame(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
  } else {
    x.mat <- model.frame(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
  }
  
  ## Linear IPW Model
  pimod <- SuperLearner(Y = x$pm25, X = x.mat, SL.lib = sl.lib)
  pimod.vals <- c(pimod$SL.predict)
  pimod.sd <- sd(x$pm25 - pimod.vals)
  
  # nonparametric density / denominator
  a.std <- c(x$pm25 - pimod.vals) / pimod.sd
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # ipw numerator
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / pimod.sd
    approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  })
  
  phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  phat <- predict(smooth.spline(a.vals, phat.vals), x = x$pm25)$y
  phat[phat < 0] <- .Machine$double.eps
  
  # combine numerator and denominator
  x$ipw <- phat/pihat # LM GPS
  x$ipw <- ifelse(is.na(x$ipw), 1, x$ipw)
  
  # truncation
  x$trunc <- x$ipw
  trunc0 <- quantile(x$ipw, 0.005)
  trunc1 <- quantile(x$ipw, 0.995)
  x$trunc[x$ipw < trunc0] <- trunc0
  x$trunc[x$ipw > trunc1] <- trunc1
  
  # merge data components such as outcomes and exposures
  if (boot) {
    wx <- merge(w, x, by = c("boot.id", "id", "zip", "year", "region"))
  } else {
    wx <- merge(w, x, by = c("id", "zip", "year", "region"))
  }
  
  if (state != "US") {
    w.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(region, zip, id, pm25, y, ybar, n, m, ipw, trunc)))
  } else {
    w.mat <- model.matrix(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, y, ybar, n, m, ipw, trunc)))
  }
  
  # calibration components
  bstar <- c(wx$pm25 - mean(x$pm25))/sd(x$pm25)
  bstar2 <- c((wx$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  
  ## Outcome Models
  # estimate nuisance outcome model with GAM
  if (dual == "both") {
    covar <- subset(wx, select = c("sex", "age_break", "race", "dual", z[-1]))
  } else {
    covar <- subset(wx, select = c("sex", "age_break", "race", "dual", z[-1]))
  }

  inner <- paste(colnames(covar), collapse = " + ")
  
  if (state == "US") {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'tp') + year + region + year:region +", inner))
  } else {
    fmla <- formula(paste0("ybar ~ s(a, bs = 'tp') + year +", inner))
  }
  
  mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
                data = data.frame(a = wx$pm25, ybar = wx$ybar, covar,
                                  year = wx$year, region = wx$region))
  
  muhat <- predict(mumod, type = "response")
  wx.mat <- predict(mumod, type = "lpmatrix")
  
  # fit model on full data
  target <- gam_dr(a = wx$pm25, y = wx$ybar, weights = wx$n, 
                   ipw = wx$trunc, muhat = muhat, a.vals = a.vals,
                   x = w.mat, w = wx.mat, astar = bstar, astar2 = bstar2,
                   k = -1, se.fit = se.fit, eco = FALSE)
  
  vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Index of Target Value
    idx <- which.min(abs(a.vals - a.tmp))
    
    # Target Means and Design
    w.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar, year = wx$year, region = wx$region), type = "lpmatrix")  
    muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar, year = wx$year, region = wx$region), type = "response")
    
    # Excess Deaths
    cut <- as.numeric(I(wx$pm25 > a.tmp))
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))
    lambda <- sum(cut*(wx$y - wx$n*(muhat.tmp + wx$trunc*(wx$ybar - mumod$fitted.values))))
    
    # ERC Estimate
    mu <- target$mu.vals[idx] + weighted.mean(muhat.tmp, w = wx$n)
    
    if (se.fit) {
    
      # Extract Target Values
      Sig <- as.matrix(target$Sig)
      g.val <- c(target$g.vals[idx,])
      
      # Dimensions
      l <- ncol(w.tmp)
      o <- length(g.val)
      
      # ERC Variance
      first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(wx$n)^2)
      second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
      sig2 <- first + second
    
      # Excess Death Variance
      var.tmp <- c(t(cut) %*% (delta*w.tmp) %*% Sig[1:l,1:l] %*% c(t(delta*w.tmp) %*% cut))
      omega2 <- sum(cut*c(wx$n*wx$trunc*(wx$ybar - mumod$fitted.values))^2) + var.tmp
      
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
                             
}