erc_implement <- function(x, w, z, se.fit = TRUE, a.vals = seq(4, 16, length.out = 121), 
                          sl.lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth"),
                          state = "US", dual = "both") {
  
  # data format
  w$ybar <- w$y/w$n
  
 if (state != "US") {
    x.mat <- model.matrix(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
    x.df <- model.frame(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
  } else {
    x.mat <- model.matrix(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
    x.df <- model.frame(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
  }
  
  ## Linear IPW Model
  # pimod <- SuperLearner(Y = x$pm25, X = x.df, SL.lib = sl.lib)
  # pimod.vals <- c(pimod$SL.predict)
  # pimod.sd <- sd(x$pm25 - pimod.vals)
  
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
  
  # Calibration Weights
  astar <- c(x$pm25 - mean(x$pm25))/var(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  
  # fit calibration model
  ipwmod <- calibrate(cmat = cmat, target = tm)
  x$ipw <- ipwmod$weights/ipwmod$base_weights
  
  # truncation
  x$trunc <- x$ipw
  trunc0 <- quantile(x$ipw, 0.005)
  trunc1 <- quantile(x$ipw, 0.995)
  x$trunc[x$ipw < trunc0] <- trunc0
  x$trunc[x$ipw > trunc1] <- trunc1
  
  # merge data components such as outcomes and exposures
    wx <- merge(w, x, by = c("id", "zip", "year", "region"))
  
  if (state != "US") {
    w.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(region, zip, id, pm25, y, ybar, n, m, ipw, trunc))) 
  } else {
    w.mat <- model.matrix(~ ., data = subset(setDF(wx), select = -c(zip, id, pm25, y, ybar, n, m, ipw, trunc))) 
  }
  
  # calibration components
  bstar <- c(wx$pm25 - mean(x$pm25))/var(x$pm25)
  bstar2 <- c((wx$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  
  ## Outcome Models
  covar <- subset(wx, select = c(z[-1]))
  inner <- paste(colnames(covar), collapse = " + ")
  
  if (state == "US") {
    fmla <- formula(paste0("ybar ~ s(a, k = 6, bs = 'mpi') + year + region + year:region +", inner))
  } else {
    fmla <- formula(paste0("ybar ~ s(a, k = 6, bs = 'mpi') + year +", inner))
  }
  
  # estimate nuisance outcome model with GAM
  mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
                data = data.frame(a = wx$pm25, ybar = wx$ybar, covar,
                                  year = wx$year, region = wx$region))
  
  muhat <- predict(mumod, type = "response")
  wx.mat <- predict(mumod, type = "lpmatrix")
  
  # fit model on full data
  target <- gam_dr(a = wx$pm25, y = wx$ybar, weights = wx$n, 
                   ipw = wx$trunc, muhat = muhat, a.vals = a.vals,
                   x = w.mat, w = wx.mat, astar = bstar, astar2 = bstar2,
                   se.fit = se.fit, family = quasipoisson(), eco = TRUE)
  
  out.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    # Index of Target Value
    idx <- which.min(abs(a.vals - a.tmp))
    
    # Target Means and Design
    w.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar, year = wx$year, region = wx$region), type = "lpmatrix")  
    muhat.tmp <- predict(mumod, newdata = data.frame(a = a.tmp, covar, year = wx$year, region = wx$region), type = "response")
    
    # Excess Deaths
    delta <- c(wx$n*mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp)))

    # ERC Estimate
    mu <- target$mu.vals[idx] + weighted.mean(muhat.tmp, w = wx$n)
    
    if (se.fit) {
    
      # Extract Target Values
      Sig <- as.matrix(target$Sig)
      g.val <- c(target$g.vals[idx,])
      
      # Dimensions
      l <- ncol(w.tmp)
      o <- length(g.val)
      m <- nrow(wx)
      
      # ERC Variance
      first <-  c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(wx$n)^2)
      second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
      third <- c(t(g.val) %*% Sig[(l + 1):(l + o), 1:l] %*% t(w.tmp) %*% delta)/sum(wx$n)
      sig2 <- first + second + 2*third
      
      return(c(mu = mu, sig2 = sig2))
      
    } else
      return(mu)
    
  })
  
  # extract estimates
  if (se.fit) {
    est_data <- data.frame(a.vals = a.vals, estimate = out.vals[1,], se = sqrt(out.vals[2,]))
  } else {
    est_data <- data.frame(a.vals = a.vals, estimate = out.vals)
  }
  
  return(list(est_data = est_data, wx = wx))
                             
}