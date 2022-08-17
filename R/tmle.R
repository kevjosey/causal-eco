
tmle_glm <- function(a_w, y, w, a_x = a_w, x = w, a.vals,
                     log.pop = NULL, trunc = 0.01, df = 4) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  # number of clusters
  n <- nrow(x)
  m <- nrow(w)
  x <- data.frame(x)
  w <- data.frame(w)  

  # error checks
  if (length(a_w) != nrow(w))
    stop("length(a_w) != nrow(w)")
  
  if (length(a_x) != nrow(x))
    stop("length(a_x) != nrow(x)")
  
  if (length(log.pop) != length(y))
    stop("length(log.pop) != length(y)")
  
  if (nrow(w) != length(y))
    stop("nrow(w) != length(y)")
  
  # estimate nuisance outcome model with glm
  mumod <- gam(y ~ s(a, 5) + . - a, weights = exp(log.pop),
               data = data.frame(y = y, a = a_w, w), 
               family = quasibinomial())
  muhat <- mumod$fitted.values
  
  # estimate nuisance GPS parameters with SuperLearner
  pimod <- lm(a ~ ., data = data.frame(a = a_x, x))
  pimod.vals <- c(pimod$fitted.values, predict(pimod, newdata = data.frame(w)))
  pimod.sd <- sigma(pimod)
  
  # parametric density
  # pihat <- dnorm(c(a_x, a_w), pimod.vals, pimod.sd)
  #
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a.tmp, pimod.vals, pimod.sd)
  # })
  #
  # phat.vals <- colMeans(pihat.mat[1:n,], na.rm = TRUE)
  # phat <- predict(smooth.spline(a.vals, phat.vals), x = c(a_x, a_w))$y
  # phat[phat<0] <- 0
  
  # nonparametric density - SL
  a.std <- c(c(a_x, a_w) - pimod.vals) / pimod.sd
  dens <- density(a.std[1:n])
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / pimod.sd
    approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  })
  
  phat.vals <- colMeans(pihat.mat[1:n,], na.rm = TRUE)
  phat <- predict(smooth.spline(a.vals, phat.vals), x = c(a_x, a_w))$y
  phat[phat < 0] <- 1e-6
  
  # TMLE update
  ipw <- phat/pihat
  trunc0 <- quantile(ipw[1:n], trunc)
  trunc1 <- quantile(ipw[1:n], 1 - trunc)
  ipw[ipw < trunc0] <- trunc0
  ipw[ipw > trunc1] <- trunc1
  
  nsa_x <- ns(a_x, df = df, intercept = TRUE)
  nsa_w <- predict(nsa_x, newx = a_w)
  base <- nsa_w*ipw[-(1:n)]
  new_mod <- glm(y ~ 0 + base, weights = exp(log.pop), family = quasibinomial())
  param <- coef(new_mod)
  a.mat <- predict(nsa_x, newx = a.vals)
  
  # predict potential outcomes and aggregate by person years
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    # subset
    mat.tmp <- a.mat[k,,drop = FALSE]
    pihat.tmp <- pihat.mat[,k]
    muhat.tmp <- predict(mumod, newdata = data.frame(a = rep(a.vals[k], m), w), type = "response")
    
    # get weights
    wts <- c(mean(pihat.tmp[1:n], na.rm = TRUE)/pihat.tmp[-(1:n)])
    wts[wts < trunc0] <- trunc0
    wts[wts > trunc1] <- trunc1
    mat <- mat.tmp[rep(1,m),]*wts
    
    # update
    muhat.val <- plogis(qlogis(muhat.tmp) + c(mat%*%param))
    return(sum(exp(log.pop)*muhat.val, na.rm = TRUE)/sum(exp(log.pop)))
    
  })
  
  return(list(estimate = estimate, weights_w = ipw[-(1:n)], weights_x = ipw[1:n]))
  
}
