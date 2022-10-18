# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(psi.lm, psi.sl, psi.cal, w.id, log.pop, int.mat, x.id, a_x,
                      loess = FALSE, bw = NULL,  a.vals = seq(min(a), max(a), length.out = 100),
                      bw.seq = seq(0.1, 2, by = 0.1), folds = 5, se.fit = TRUE) {	
  
  # select bw if null
  if (is.null(bw))
    bw <- cv_bw(a = a_w, psi = psi.lm, weights = exp(log.pop), folds = folds, bw.seq = bw.seq)
  
  # marginalize psi within zip-year
  wts <- do.call(c, lapply(split(exp(log.pop), w.id), sum))
  
  list.lm <- split(data.frame(psi = psi.lm, wts = exp(log.pop)), w.id)
  psi.lm.new <- data.frame(psi = do.call(c, lapply(list.lm, function(df) sum(df$psi*df$wts)/sum(df$wts))), 
                            wts = wts, id = names(list.lm))
  lm.dat <- inner_join(psi.lm.new, data.frame(a = a_x, id = x.id), by = "id")

  list.sl <- split(data.frame(psi = psi.sl, wts = exp(log.pop)) , w.id)
  psi.sl.new <- data.frame(psi = do.call(c, lapply(list.sl, function(df) sum(df$psi*df$wts)/sum(df$wts))),
                           wts = wts, id = names(list.sl))
  sl.dat <- inner_join(psi.sl.new, data.frame(a = a_x, id = x.id), by = "id")
  
  list.cal <- split(data.frame(psi = psi.cal, wts = exp(log.pop)) , w.id)
  psi.cal.new <- data.frame(psi = do.call(c, lapply(list.cal, function(df) sum(df$psi*df$wts)/sum(df$wts))), 
                            wts = wts, id = names(list.cal))
  cal.dat <- inner_join(psi.cal.new, data.frame(a = a_x, id = x.id), by = "id")
  
  mat.list <- split(cbind(exp(log.pop), int.mat), w.id)
  int.mat.new <- do.call(rbind, lapply(split(cbind(exp(log.pop), int.mat), w.id), 
                      function(vec) { mat <- matrix(vec, ncol = length(a.vals) + 1);
                       colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1]) } ))
  
  # KWLS regression
  out.lm <- sapply(a.vals, kern_est, psi = lm.dat$psi, a = lm.dat$a, weights = lm.dat$wts,
                   loess = loess, bw = bw, se.fit = se.fit, int.mat = int.mat.new, a.vals = a.vals)
  
  out.sl <- sapply(a.vals, kern_est, psi = sl.dat$psi, a = sl.dat$a, weights = sl.dat$wts,
                   loess = loess, bw = bw, se.fit = se.fit, int.mat = int.mat.new, a.vals = a.vals)
  
  out.cal <- sapply(a.vals, kern_est, psi = cal.dat$psi, a = cal.dat$a, weights = cal.dat$wts,
                    loess = loess, bw = bw, se.fit = se.fit, int.mat = int.mat.new, a.vals = a.vals)
  
  # linear model approx
  fit.lm <- lm(psi ~ a, weights = lm.dat$wts, data = lm.dat) 
  fit.sl <- lm(psi ~ a, weights = sl.dat$wts, data = sl.dat)
  fit.cal <- lm(psi ~ a, weights = cal.dat$wts, data = cal.dat)
  
  if (se.fit) {
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    estimate.sl <- out.sl[1,]
    variance.sl <- out.sl[2,]
    estimate.cal <- out.cal[1,]
    variance.cal <- out.cal[2,]
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm, fit.lm = fit.lm, 
                estimate.sl = estimate.sl, variance.sl = variance.sl, fit.sl = fit.sl,
                estimate.cal = estimate.cal, variance.cal = variance.cal, fit.cal = fit.cal))
    
  } else {
    
    return(list( estimate.lm = out.lm, fit.lm = fit.lm, 
                 estimate.sl = out.sl, fit.sl = fit.sl,
                 estimate.cal = out.cal, fit.cal = fit.cal))
                
  }
  
}

# estimate glm outcome model (Keith - adapt code to fit a Poisson outcome as an excercise)
gam_est <- function(a_w, y, w, w.id, a_x, x, x.id, 
                    a.vals, log.pop = NULL, trunc = 0.01,
                    sl.lib = c("SL.mean","SL.glm"), ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  m <- nrow(w)
  x <- data.frame(x)
  w <- data.frame(w)
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - 1e-6
  
  # estimate nuisance outcome model with glm
  mumod <- gam(ybar ~ s(a, 5) + . - a + 
                 a:(regionWEST + regionNORTHEAST + regionSOUTH), # need better coding to generalize
               weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a_w, w), 
               family = quasipoisson())
  muhat <- mumod$fitted.values
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = exp(log.pop))
  
  # LM
  pimod.lm <- lm(a ~ ., data = data.frame(a = a_x, x))
  pimod.vals.lm <- c(pimod.lm$fitted.values, predict(pimod.lm, newdata = data.frame(w)))
  pimod.sd.lm <- sigma(pimod.lm)
  
  # nonparametric denisty - LM
  a.std.lm <- c(c(a_x, a_w) - pimod.vals.lm) / pimod.sd.lm
  dens.lm <- density(a.std.lm[1:n])
  pihat.lm <- approx(x = dens.lm$x, y = dens.lm$y, xout = a.std.lm)$y / pimod.sd.lm
  
  pihat.mat.lm <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals.lm) / pimod.sd.lm
    approx(x = dens.lm$x, y = dens.lm$y, xout = std)$y / pimod.sd.lm
  })
  
  phat.vals.lm <- colMeans(pihat.mat.lm[1:n,], na.rm = TRUE)
  phat.lm <- predict(smooth.spline(a.vals, phat.vals.lm), x = c(a_x, a_w))$y
  phat.lm[phat.lm < 0] <- .Machine$double.eps
  
  # SuperLearner
  pimod.sl <- SuperLearner(Y = a_x, X = x, SL.library = sl.lib, family = gaussian())
  pimod.vals.sl <- c(c(pimod.sl$SL.predict), c(predict(pimod.sl, newdata = data.frame(w))$pred))
  pimod.sd.sl <- sd(a_x - pimod.vals.sl[1:n])
  
  # nonparametric density - SL
  a.std.sl <- c(c(a_x, a_w) - pimod.vals.sl) / pimod.sd.sl
  dens.sl <- density(a.std.sl[1:n])
  pihat.sl <- approx(x = dens.sl$x, y = dens.sl$y, xout = a.std.sl)$y / pimod.sd.sl

  pihat.mat.sl <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals.sl) / pimod.sd.sl
    approx(x = dens.sl$x, y = dens.sl$y, xout = std)$y / pimod.sd.sl
  })

  phat.vals.sl <- colMeans(pihat.mat.sl[1:n,], na.rm = TRUE)
  phat.sl <- predict(smooth.spline(a.vals, phat.vals.sl), x = c(a_x, a_w))$y
  phat.sl[phat.sl < 0] <- .Machine$double.eps
  
  # truncation
  ipw.lm <- phat.lm/pihat.lm
  trunc0.lm <- quantile(ipw.lm[1:n], trunc)
  trunc1.lm <- quantile(ipw.lm[1:n], 1 - trunc)
  ipw.lm[ipw.lm < trunc0.lm] <- trunc0.lm
  ipw.lm[ipw.lm > trunc1.lm] <- trunc1.lm

  ipw.sl <- phat.sl/pihat.sl
  trunc0.sl <- quantile(ipw.sl[1:n], trunc)
  trunc1.sl <- quantile(ipw.sl[1:n], 1 - trunc)
  ipw.sl[ipw.sl < trunc0.sl] <- trunc0.sl
  ipw.sl[ipw.sl > trunc1.sl] <- trunc1.sl
  
  # full calibration weights
  x <- x %>% mutate_if(is.numeric, scale)
  x.mat <- model.matrix(~ ., data = data.frame(x))
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a_x - mean(a_x))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(n, rep(0, ncol(x.mat) + 1)))
  
  ipw.cal_x <- mod$weights
  ipw.cal_mat <- inner_join(x = data.frame(id = w.id), 
                            y = data.frame(id = x.id, wts = ipw.cal_x), 
                            by = "id")
  ipw.cal_w <- ipw.cal_mat$wts
  ipw.cal <- c(ipw.cal_x, ipw.cal_w)
  
  if (any(ipw.cal_mat$id != w.id))
    stop("id is getting scrambled!")
  
  # pseudo outcome
  resid.lm <- c(ybar - muhat)*ipw.lm[-(1:n)]
  resid.sl <- c(ybar - muhat)*ipw.sl[-(1:n)]
  resid.cal <- c(ybar - muhat)*ipw.cal[-(1:n)]
  
  out <- list(resid.lm = resid.lm, weights.lm_x = ipw.lm[1:n], weights.lm_w = ipw.lm[-(1:n)],
              resid.sl = resid.sl, weights.sl_x = ipw.sl[1:n], weights.sl_w = ipw.sl[-(1:n)],
              resid.cal = resid.cal, weights.cal_x = ipw.cal_x, weights.cal_w = ipw.cal_w, 
              muhat.mat = muhat.mat, phat.vals = phat.vals.lm, w.id = w.id, log.pop = log.pop)
  
  return(out)
  
}

# Kernel weighted least squares
kern_est <- function(a.new, a, psi, bw, weights, se.fit = FALSE, 
                     int.mat = NULL, a.vals = NULL, loess = FALSE) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Tricube Weights
  if (loess) {
  
    a.std <- a - a.new
    k <- floor(min(bw, 1)*length(a))
    idx <- order(abs(a.std))[1:k]
    
    # subset
    a.std <- a.std[idx]
    psi <- psi[idx]
    weights <- weights[idx]
    int.mat <- int.mat[idx,]
    max.a.std <- max(abs(a.std))
    
    # construct kernel weight
    k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
    g.std <- cbind(1, a.std)

  } else {   # Gaussian Kernel
    
    a.std <- (a - a.new) / bw
    k.std <- dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  }
    
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    # LOESS
    if (loess) {
      
      kern.mat <- matrix(rep(c((1 - abs((a.vals - a.new)/max.a.std)^3)^3), k), byrow = T, nrow = k)
      kern.mat[matrix(rep(abs(a.vals - a.new)/max.a.std, k), byrow = T, nrow = k) > 1] <- 0
      g.vals <- matrix(rep(c(a.vals - a.new), k), byrow = T, nrow = k)
      
      intfn1.mat <- kern.mat * int.mat
      intfn2.mat <- g.vals * kern.mat * int.mat
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
      
    } else {  # Gaussian Kernel
      
      kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
      g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
      
      intfn1.mat <- kern.mat * int.mat
      intfn2.mat <- g.vals * kern.mat * int.mat
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
      
    }
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select bw
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.1, 2, by = 0.1)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  idx <- sample(x = folds, size = n, replace = TRUE)
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a[idx == k], kern_est, psi = psi[idx != k], a = a[idx != k], 
                      weights = weights[idx != k], bw = h, se.fit = FALSE)
      cv.vec[k] <- mean((psi[idx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}

opt_fun <- function(param, psi, g.std, k.std) {
  
  sum(k.std*(psi - exp(c(g.std %*% param)))^2)
  
}
