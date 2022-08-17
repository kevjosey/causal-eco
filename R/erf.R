# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
count_erf <- function(a_w, y, w, a_x, x, w.id, x.id,
                      log.pop = NULL, trunc = 0.01,
                      a.vals = seq(min(a), max(a), length.out = 100),
                      bw = NULL, bw.seq = seq(0.1, 2, by = 0.1), folds = 5,
                      sl.lib = c("SL.mean", "SL.glm"), se.fit = TRUE) {	
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  n <- length(a_x)
  m <- length(a_w)

  wrap <- gam_est(a_w = a_w, y = y, w = w, w.id = w.id, 
                  a_x = a_x, x = x, x.id = x.id, 
                  a.vals = a.vals, log.pop = log.pop, 
                  trunc = trunc, sl.lib = sl.lib)
  
  psi.lm <- wrap$psi.lm
  psi.sl <- wrap$psi.sl
  psi.cal <- wrap$psi.cal
  int.mat <- wrap$int.mat
  
  # select bw if null
  if (is.null(bw))
    bw <- cv_bw(a = a_w, psi = psi.sl, weights = exp(log.pop), folds = folds, bw.seq = bw.seq)
  
  # asymptotics
  out.lm <- sapply(a.vals, kern_est, psi = psi.lm, a = a_w, weights = exp(log.pop),
                   bw = bw, se.fit = se.fit, int.mat = int.mat, a.vals = a.vals)
  
  out.sl <- sapply(a.vals, kern_est, psi = psi.sl, a = a_w, weights = exp(log.pop),
                   bw = bw, se.fit = se.fit, int.mat = int.mat, a.vals = a.vals)
  
  out.cal <- sapply(a.vals, kern_est, psi = psi.cal, a = a_w, weights = exp(log.pop),
                    bw = bw, se.fit = se.fit, int.mat = int.mat, a.vals = a.vals)
  
  if (se.fit) {
    
    estimate.lm <- out.lm[1,]
    variance.lm <- out.lm[2,]
    estimate.sl <- out.sl[1,]
    variance.sl <- out.sl[2,]
    estimate.cal <- out.cal[1,]
    variance.cal <- out.cal[2,]
    
    names(estimate.lm) <- names(estimate.sl) <- names(estimate.cal) <- a.vals
    names(variance.lm) <- names(variance.sl) <- names(variance.cal) <- a.vals
    
    return(list(estimate.lm = estimate.lm, variance.lm = variance.lm,
                estimate.sl = estimate.sl, variance.sl = variance.sl,
                estimate.cal = estimate.cal, variance.cal = variance.cal,
                weights.lm_x = wrap$weights.lm_x, weights.lm_w = wrap$weights.lm_w,
                weights.sl_x = wrap$weights.sl_x, weights.sl_w = wrap$weights.sl_w,
                weights.cal_x = wrap$weights.cal_x, weights.cal_w = wrap$weights.cal_w))
    
  } else {
    
    names(out.lm) <- names(out.sl) <- names(out.cal) <- a.vals
    return(list(estimate.lm = out.lm, estimate.sl = out.sl, estimate.cal = out.cal,
                weights.lm_x = wrap$weights.lm_x, weights.lm_w = wrap$weights.lm_w,
                weights.sl_x = wrap$weights.sl_x, weights.sl_w = wrap$weights.sl_w,
                weights.cal_x = wrap$weights.cal_x, weights.cal_w = wrap$weights.cal_w))
    
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
  mumod <- gam(ybar ~ s(a, 5) + . - a, weights = exp(log.pop),
               data = data.frame(ybar = ybar, a = a_w, w), 
               family = quasibinomial())
  muhat <- mumod$fitted.values
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  mhat.vals <- colSums(exp(log.pop)*muhat.mat, na.rm = TRUE)/sum(exp(log.pop))
  mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a_w)$y
  
  # LM
  pimod.lm <- lm(a ~ ., data = data.frame(a = a_x, x))
  pimod.vals.lm <- c(pimod.lm$fitted.values, predict(pimod.lm, newdata = data.frame(w)))
  pimod.sd.lm <- sigma(pimod.lm)
  
  # SuperLearner
  pimod.sl <- SuperLearner(Y = a_x, X = x, SL.library = sl.lib, family = gaussian())
  pimod.vals.sl <- c(c(pimod.sl$SL.predict), c(predict(pimod.sl, newdata = data.frame(w))$pred))
  pimod.sd.sl <- sd(a_x - pimod.vals.sl[1:n])
  
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
  phat.lm[phat.lm < 0] <- 1e-6
  
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
  phat.sl[phat.sl < 0] <- 1e-6
  
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
  m <- ncol(x.mat) + 1
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a_x - mean(a_x))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), target = c(n, rep(0, m)))
  ipw.cal_x <- mod$weights
  ipw.cal_mat <- inner_join(x = data.frame(id = w.id), 
                            y = data.frame(id = x.id, wts = ipw.cal_x), 
                            by = "id")
  ipw.cal_w <- ipw.cal_mat$wts
  ipw.cal <- c(ipw.cal_x, ipw.cal_w)
  
  if (any(ipw.cal_mat$id != w.id))
    stop("id is getting scrambled!")
  
  # pseudo outcome
  psi.lm <- c(ybar - muhat)*ipw.lm[-(1:n)] + mhat
  psi.sl <- c(ybar - muhat)*ipw.sl[-(1:n)] + mhat
  psi.cal <- c(ybar - muhat)*ipw.cal[-(1:n)] + mhat
  
  # integration matrix
  mhat.mat <- matrix(rep(mhat.vals, m), byrow = TRUE, nrow = m)
  phat.mat <- matrix(rep(phat.vals.sl, m), byrow = TRUE, nrow = m)  
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  out <- list(psi.lm = psi.lm, psi.sl = psi.sl, psi.cal = psi.cal, int.mat = int.mat, 
              weights.lm_x = ipw.lm[1:n], weights.lm_w = ipw.lm[-(1:n)],
              weights.sl_x = ipw.sl[1:n], weights.sl_w = ipw.sl[-(1:n)],
              weights.cal_x = ipw.cal_x, weights.cal_w = ipw.cal_w)
  
  return(out)
  
}

# Kernel weighted least squares
kern_est <- function(a.new, a, psi, bw, weights, se.fit = FALSE, int.mat = NULL, a.vals = NULL) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    # Gaussian Kernel Matrix
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.vals * kern.mat * int.mat
    
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
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
