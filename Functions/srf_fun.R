srf_implement <- function(delta, x, w, z, 
                          sl.lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth"),
                          state = "US", dual = "both", ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  if (state != "US") {
    x.df <- model.frame(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
    x.mat <- model.matrix(~ ., data = subset(setDF(x), select = -c(region, zip, id, pm25, m)))
  } else {
    x.df <- model.frame(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
    x.mat <- model.matrix(~ . + year:region, data = subset(setDF(x), select = -c(zip, id, pm25, m)))
  }

  w$ybar <- w$y/w$n
  x$pm25.new <- ifelse(x$pm25 > delta, delta, x$pm25)
  
  ## Linear IPW Model
  pimod <- SuperLearner(Y = x$pm25, X = x.df, SL.lib = sl.lib, family = gaussian())
  pimod.vals <- c(pimod$SL.predict)
  pimod.sd <- sd(x$pm25 - pimod.vals)
  
  # nonparametric density / denominator
  PrA <- mean(I(x$pm25 < delta))
  a.std <- c(x$pm25 - pimod.vals) / pimod.sd
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # ipw numerator
  pitilde <- PrA*pihat + (1 - PrA)*mean(pihat)
  
  # combine numerator and denominator
  x$ipw <- pitilde/pihat # LM GPS
  x$ipw <- ifelse(is.na(x$ipw), 1, x$ipw)
  
  # Truncate IPW
  x$ipw_trunc <- x$ipw
  trunc0 <- quantile(x$ipw, 0.005)
  trunc1 <- quantile(x$ipw, 0.995)
  x$ipw_trunc[x$ipw < trunc0] <- trunc0
  x$ipw_trunc[x$ipw > trunc1] <- trunc1
  
  ## Calibratio Weights
  astar <- c(x$pm25 - mean(x$pm25))/sd(x$pm25)
  astar2 <- c((x$pm25 - mean(x$pm25))^2/var(x$pm25) - 1)
  cmat1 <- cbind(x.mat*astar, astar2, x.mat)
  atilde <- c(x$pm25.new - mean(x$pm25))/sd(x$pm25)
  atilde2 <- c((x$pm25.new - mean(x$pm25))^2/var(x$pm25) - 1)
  cmat0 <- cbind(x.mat*atilde, atilde2, x.mat)
  
  # construct linear term vector
  ipwmod <- calibrate(cmat = cmat1, target = colSums(cmat0))
  x$cal <- ipwmod$weights
  
  # Truncate Calibration
  x$cal_trunc <- x$cal
  trunc0 <- quantile(x$cal, 0.005)
  trunc1 <- quantile(x$cal, 0.995)
  x$cal_trunc[x$cal < trunc0] <- trunc0
  x$cal_trunc[x$cal > trunc1] <- trunc1
  
  # compute imbalances
  imbalance <- c(colMeans(cmat0)) - c(t(cmat1) %*% x$cal)/nrow(cmat1)
  names(imbalance) <- colnames(cmat1)
  
  # merge data components such as outcomes and exposures
  wx <- merge(w, x, by = c("id", "zip", "year", "region"))
  
  ## Outcome Models
  covar <- subset(wx, select = c(z[-1]))
  inner <- paste(colnames(covar), collapse = " + ")
  
  if (state == "US") {
    fmla <- formula(paste0("ybar ~ s(a, k = 6, bs = 'tp') + year + region + year:region +", inner))
  } else {
    fmla <- formula(paste0("ybar ~ s(a, k = 6, bs = 'tp') + year +", inner))
  }
  
  # estimate nuisance outcome model with GAM
  mumod <- scam(fmla, weights = wx$n, family = quasipoisson(),
                data = data.frame(a = wx$pm25, ybar = wx$ybar, covar,
                                  year = wx$year, region = wx$region))
  
  muhat <- predict(mumod, type = "response")
  mutilde <- predict(mumod, type = "response", newdata = data.frame(a = wx$pm25.new, wx))

  # Point Estimates
  wx$psi_ipw <- c(wx$ybar - muhat)*wx$ipw_trunc + mutilde
  wx$psi_cal <- c(wx$ybar - muhat)*wx$cal_trunc + mutilde
  mu_ipw <- weighted.mean(wx$psi_ipw, w = wx$n)
  theta_ipw <- weighted.mean(wx$psi_ipw - wx$ybar, w = wx$n)
  mu_cal <- weighted.mean(wx$psi_cal, w = wx$n)
  theta_cal <- weighted.mean(wx$psi_cal - wx$ybar, w = wx$n)
  
  # variances (needs work)
  # var.tmp <- wx %>% group_by(zip) %>% 
  #   summarize(mu_eif = weighted.mean(psi, w = n) - mu,
  #             theta_eif = weighted.mean(psi - ybar, w = n) - theta,
  #             m = sum(n))
  
  var.tmp.wi <- wx %>% group_by(zip) %>% 
    summarize(mu_eif = psi_ipw - weighted.mean(psi_ipw, w = n),
              theta_eif = psi_ipw - ybar - weighted.mean(psi_ipw - ybar, w = n),
              m = sum(n))
  var.tmp.bw <- wx %>% group_by(zip) %>% 
    summarize(mu_eif = weighted.mean(psi_ipw, w = n) - mu_ipw,
              theta_eif = weighted.mean(psi_ipw - ybar, w = n) - theta_ipw,
              m = sum(n))
  
  sig2 <- sum(var.tmp.wi$m*(var.tmp.wi$mu_eif)^2/sum(var.tmp.wi$m))/(nrow(var.tmp.wi) - nrow(var.tmp.bw)) + 
    sum(var.tmp.bw$m*(var.tmp.bw$mu_eif)^2/sum(var.tmp.bw$m))/(nrow(var.tmp.bw) - 1)
  omega2 <- sum(var.tmp.wi$m*(var.tmp.wi$theta_eif)^2/sum(var.tmp.wi$m))/(nrow(var.tmp.wi) - nrow(var.tmp.bw)) + 
    sum(var.tmp.bw$m*(var.tmp.bw$theta_eif)^2/sum(var.tmp.bw$m))/(nrow(var.tmp.bw) - 1)
  
  # return output
  return(list(mu = mu_cal, theta = theta_cal, sig2 = sig2, omega2 = omega2, ipw = wx$ipw))
  
}