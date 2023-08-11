## kernel estimation - simple
kern_est <- function(a.new, a, psi, bw = 1, weights = NULL, se.fit = FALSE,
                     x = NULL, astar = NULL, astar2 = NULL, cmat = NULL, ipw = NULL) {
  
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
    
    m <- ncol(cmat)
    U <- matrix(0, ncol = m, nrow = m)
    V <- matrix(0, ncol = m + 2, nrow = 2)
    meat <- matrix(0, ncol = m + 2, nrow = m + 2)
    eta <- c(g.std%*%b)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - ipw[i] * tcrossprod(cmat[i,])
      V[,1:m] <- V[,1:m] - weights[i]*k.std[i]*psi[i]*tcrossprod(g.std[i,],cmat[i,])
      V[,(m + 1):(m + 2)] <- V[,(m + 1):(m + 2)] - weights[i]*k.std[i]*tcrossprod(g.std[i,])
      
      meat <- meat + 
        tcrossprod(esteq(p = ipw[i], x = x[i,], psi = psi[i],
                         g.std = g.std[i,], k.std = k.std[i],
                         astar = astar[i], astar2 = astar2[i], 
                         tm = colMeans(x), weights = weights[i], eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = m + 2, ncol = m + 2)
    invbread[1:m,1:m] <- U
    invbread[(m + 1):(m + 2), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[m + 2, m + 2]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
  } else
    return(mu)
  
}

## kernel estimation for ecological exposures
kern_est_eco <- function(a.new, a, psi, bw = 1, weights = NULL, se.fit = FALSE,
                         x = NULL, astar = NULL, astar2 = NULL, 
                         cmat = NULL, ipw = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))

  # Gaussian Kernel
  a.std <- sqrt(weights)*(a - a.new) / bw
  k.std <- sqrt(weights)*dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ a.std, weights = k.std)$coefficients
  mu <- unname(b[1])
  
  if (se.fit) {
    
    m <- ncol(cmat)
    U <- matrix(0, ncol = m, nrow = m)
    V <- matrix(0, ncol = m + 2, nrow = 2)
    meat <- matrix(0, ncol = m + 2, nrow = m + 2)
    eta <- c(g.std%*%b)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - ipw[i] * tcrossprod(cmat[i,])
      V[,1:m] <- V[,1:m] - k.std[i]*psi[i]*tcrossprod(g.std[i,],cmat[i,])
      V[,(m + 1):(m + 2)] <- V[,(m + 1):(m + 2)] - k.std[i]*tcrossprod(g.std[i,])
      
      meat <- meat + 
        tcrossprod(esteq(p = ipw[i], x = x[i,], psi = psi[i],
                         g.std = g.std[i,], k.std = k.std[i],
                         astar = astar[i], astar2 = astar2[i], 
                         tm = colMeans(x), weights = 1, eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = m + 2, ncol = m + 2)
    invbread[1:m,1:m] <- U
    invbread[(m + 1):(m + 2), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[m + 2, m + 2]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
  } else
    return(mu)
  
}

esteq <- function(p, x, psi, tm, g.std, k.std, astar, astar2, weights = 1, eta) {
  
  eq1 <- p*x*astar
  eq2 <- p*astar2
  eq3 <- p*x - tm
  eq4 <- weights*k.std*(psi - eta)*g.std
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
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
