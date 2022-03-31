
tmle_glm <- function(a_w, y, w, a_x = a_w, x = w, a.vals, offset = rep(0, length(y)),
                     nsa = NULL, family = gaussian(), df = 4, trunc = 0.01) {
  
  # number of clusters
  n <- nrow(x)
  
  # error checks
  if (length(a_w) != nrow(w))
    stop("length(a_w) != nrow(w)")
  
  if (length(a_x) != nrow(x))
    stop("length(a_x) != nrow(x)")
  
  if (length(offset) != length(y))
    stop("length(offset) != length(y)")
  
  if (nrow(w) != length(y))
    stop("nrow(w) != length(y)")
  
  if (!is.null(nsa) & !inherits(nsa, c("bs", "ns")))
    stop("nsa must be a spline basis from the splines package.")

  if (is.null(nsa)) {
    nsa <- ns(a_x, df = df)
  } else{
    df <- ncol(nsa)
  }
  
  nsa_w <- predict(nsa, newx = a_w)
  
  # estimate nuisance outcome model with glm
  mumod <- glm(y ~ ., data = data.frame(nsa_w, w), offset = offset, family = family)
  muhat <- predict(mumod, newdata = data.frame(nsa_w, w), type = "response")
  
  # estimate nuisance GPS parameters with SuperLearner
  pimod <- lm(a ~ ., data = data.frame(x, a = a_x))
  pimod.vals <- c(pimod$fitted.values, predict(pimod, newdata = data.frame(w)))
  pi2mod.vals <- sigma(pimod)^2
  
  # parametric density
  # pihat <- dnorm(c(a_x, a_w), pimod.vals, sqrt(pi2mod.vals))
  #
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))
  # })
  #
  # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n,], na.rm = T)), x = c(a_x, a_w))$y
  # phat[phat<0] <- 0
  
  # nonparametric denisty
  a.std <- c(c(a_x, a_w) - pimod.vals) / sqrt(pi2mod.vals)
  dens <- density(a.std[1:n])
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
  
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / sqrt(pi2mod.vals)
    approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
  })
  
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n,], na.rm = T)), x = c(a_x, a_w))$y
  phat[phat<0] <- 1e-6
  
  # TMLE update
  ipw <- phat/pihat
  trunc0 <- quantile(ipw[1:n], trunc)
  trunc1 <- quantile(ipw[1:n], 1 - trunc)
  ipw[ipw < trunc0] <- trunc0
  ipw[ipw > trunc1] <- trunc1
  
  base <- cbind(1, nsa_w)*ipw[-(1:n)]
  new_mod <- glm(y ~ 0 + base, offset = family$linkfun(muhat) + offset, family = family)
  param <- coef(new_mod)
  wa.tmp <- model.matrix(~ ., data = data.frame(nsa_w, w))
  a.mat <- predict(nsa, newx = a.vals)
  
  # predict potential outcomes and aggregate by person years
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    # subset
    mat.tmp <- a.mat[k,,drop = FALSE]
    wa.tmp[,2:(df + 1)] <- mat.tmp[rep(1,nrow(wa.tmp)),]
    pihat.tmp <- pihat.mat[,k]
    muhat.tmp <- family$linkinv(c(wa.tmp%*%(mumod$coefficients)))
    
    # get weights
    wts <- c(mean(pihat.tmp[1:n], na.rm = TRUE)/pihat.tmp[-(1:n)])
    wts[wts < trunc0] <- trunc0
    wts[wts > trunc1] <- trunc1
    mat <- cbind(1, mat.tmp[rep(1,length(wts)),])*wts
    
    # update
    muhat.val <- family$linkinv(family$linkfun(muhat.tmp) + c(mat%*%param))
    return(weighted.mean(muhat.val, w = family$linkinv(offset), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, weights_w = ipw[-(1:n)], weights_x = ipw[1:n]))
  
}
