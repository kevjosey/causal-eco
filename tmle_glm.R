
tmle_glm <- function(a_w, a_x = a_w, w, x = w, y, offset, a.vals, 
                     family = gaussian(), df = 4, trunc = 0.01) {
  
  # number of clusters
  n <- nrow(x)
  
  if (length(a_w) != nrow(w))
    stop("length(a_w) != nrow(w)")
  
  if (length(a_x) != nrow(x))
    stop("length(a_x) != nrow(x)")
  
  if (length(offset) != length(y))
    stop("length(offset) != length(y)")
  
  if (nrow(w) != length(y))
    stop("nrow(w) != length(y)")
  
  # estimate nuisance outcome model with splines
  fmla <- formula(paste0("y ~ ns(a,", df,") +", paste0(colnames(w), collapse = "+")))
  mumod <- glm(fmla, data = data.frame(w, a = a_w), offset = offset, family = family)
  muhat <- predict(mumod, newdata = data.frame(w, a = a_w))
  
  # estimate nuisance GPS parameters with lm
  pimod <- lm(a ~ ., data = data.frame(x, a = a_x))
  pimod.vals <- c(pimod$fitted.values, predict(pimod, newdata = data.frame(w)))
  pi2mod.vals <- sigma(pimod)^2
  
  # parametric density
  # pihat <- dnorm(c(a_x, a_w), pimod.vals, sqrt(pi2mod.vals))
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))
  # })
  # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n,], na.rm = T)), x = c(a_x, a_w))$y
  # phat[phat<0] <- 1e-6
  
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
  nsa <- ns(a_x, df = df + 1, intercept = TRUE)
  ipw <- phat/pihat
  trunc0 <- quantile(ipw[1:n], trunc)
  trunc1 <- quantile(ipw[1:n], 1 - trunc)
  ipw[ipw < trunc0] <- trim0
  ipw[ipw > trunc1] <- trim1
  base <- predict(nsa, newx = a_w)*ipw[-(1:n)]
  new_mod <- glm(y ~ 0 + base, offset = family$linkfun(muhat) + offset, family = family)
  param <- coef(new_mod)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    muhat.tmp <- predict(mumod, newdata = data.frame(w, a = a.vals[k]))
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp[1:n], na.rm = TRUE)/pihat.tmp[-(1:n)])
    wts[wts < trunc0] <- trunc0
    wts[wts > trunc1] <- trunc1
    mat <- predict(nsa, newx = rep(a.tmp, length(wts)))*wts
    return(weighted.mean(family$linkinc(log(muhat.tmp) + c(mat%*%param)), 
                         w = family$linkinv(offset), na.rm = TRUE))
  })
  
  return(list(estimate = estimate, weights = ipw[-(1:n)]))
  
}
