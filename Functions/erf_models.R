
# count_erf is a wrapper for the KWLS algorithms
count_erf <- function(psi, log.pop, w.id, a, x.id,
                      a.vals = seq(min(a), max(a), length.out = 100), 
                      bw.seq = seq(0.1, 3, length.out = 20), bw = NULL) {	
  
  # Separate Data into List
  mat.list <- split(cbind(exp(log.pop), psi), w.id)
  
  # Aggregate by ZIP-code-year
  psi.pool <- do.call(c, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = 2)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.pool <- data.frame(id = names(mat.list), psi = psi.pool)
  dat <- inner_join(mat.pool, data.frame(a = a, id = x.id), by = "id")
  
  if (is.null(phat.vals)) {
    warning("Setting phat.vals = rep(1, length(a.vals))")
    phat.vals <- rep(1, length(a.vals))
  }
  
  # grid search bandwidth
  if (is.null(bw)) {
    risk.est <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                       psi = dat$psi, a = dat$a)
    bw <- c(bw.seq[which.min(risk.est)])
  }

  # KWLS Regression
  out <- sapply(a.vals, kern_est_eco, psi = dat$psi, a = dat$a, bw = bw, se.fit = TRUE)
  
  estimate <- out[1,]
  variance <- out[2,]
  
  return(list(estimate = estimate, variance = variance))

}

## kernel estimation
kern_est <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta)),
               weights * (a.std * k.std * (psi - eta)))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(c(mu = mu))
  
}

## kernel estimation
kern_est_eco <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- sqrt(weights)*(a - a.new) / bw
  k.std <- sqrt(weights)*dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std)$coefficients
  mu <- b[1]
  
  if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, k.std*g.std))
    V <- cbind((k.std * (psi - eta)),
               (a.std * k.std * (psi - eta)))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(c(mu = mu))
  
}

## Leave-one-out cross-validated bandwidth
w.fn <- function(h, a, a.vals) {
  
  w.avals <- sapply(a.vals, function(a.tmp, ...) {
    a.std <- (a - a.tmp) / h
    k.std <- dnorm(a.std) / h
    return(mean(a.std^2 * k.std) * (dnorm(0) / h) /
             (mean(k.std) * mean(a.std^2 * k.std) - mean(a.std * k.std)^2))
  })
  
  return(w.avals / length(a))
  
}

hatvals <- function(h, a, a.vals) {
  approx(a.vals, w.fn(h = h, a = a, a.vals = a.vals), xout = a)$y
}

cts.eff.fn <- function(psi, a, h) {
  approx(locpoly(a, psi, bandwidth = h), xout = a)$y
}

risk.fn <- function(h, psi, a, a.vals) {
  hats <- hatvals(h = h, a = a, a.vals = a.vals)
  sqrt(mean(((psi - cts.eff.fn(psi = psi, a = a, h = h)) / (1 - hats))^2, na.rm = TRUE))
}
