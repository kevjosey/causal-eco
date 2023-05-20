
## kernel estimation
kern_est_simple <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- c(g.std %*% b)
    
    g.mat <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.mat * kern.mat * int.mat
    int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(c(mu = mu))
  
}

## kernel estimation
kern_est_complex <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- sqrt(weights)*(a - a.new) / bw
  k.std <- sqrt(weights)*dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std)$coefficients
  mu <- b[1]
  n.std <- sum(k.std) 
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- c(g.std %*% b)
    
    g.mat <- sapply(1:length(a.vals), function(i) sqrt(weights)*c(a.vals[i] - a.new) / bw)
    kern.mat <- sqrt(weights)*dnorm(g.mat) / bw
    
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.mat * kern.mat * int.mat
    int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    U <- solve(crossprod(g.std, k.std*g.std))
    V <- cbind((k.std * (psi - eta) + int1),
               (a.std * k.std * (psi - eta) + int2))
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
