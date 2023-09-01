
## kernel estimation for ecological exposures
kwls_est <- function(a.new, a, psi, bw = 1, weights = NULL,
                     se.fit = FALSE, eco = FALSE, sandwich = FALSE,
                     astar = NULL, astar2 = NULL, x = NULL, cmat = NULL, ipw = NULL) {
  
  n <- length(a)
  
  if (eco) {
    
    if (is.null(weights))
      weights <- rep(1, times = length(a))
    
    # Gaussian Kernel
    a.std <- sqrt(weights)*(a - a.new) / bw
    k.std <- sqrt(weights)*dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  } else {
    
    a.std <- (a - a.new) / bw
    k.std <- dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  }
  
  b <- lm(psi ~ a.std, weights = k.std)$coefficients
  mu <- unname(b[1])
  
  if (se.fit & sandwich) {
    
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
        tcrossprod(esteq_kwls(ipw = ipw[i], x = x[i,], psi = psi[i],
                              g.std = g.std[i,], k.std = k.std[i],
                              astar = astar[i], astar2 = astar2[i], eta = eta[i]))
      
      
    }
    
    invbread <- matrix(0, nrow = m + 2, ncol = m + 2)
    invbread[1:m,1:m] <- U
    invbread[(m + 1):(m + 2), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      Sig <- NA
      variance <- NA
      
    } else {
      
      Sig <- bread %*% meat %*% t(bread)
      variance <- Sig[m + 1, m + 1]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
  } else if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, k.std*g.std))
    V <- cbind((weights * (psi - eta)),
               (a.std * k.std * (psi - eta)))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
      
  } else
    return(mu)
  
}


## kernel estimation for ecological exposures
kwls_est2 <- function(a.new, a, y, bw = 1, weights = NULL, family = gaussian(),
                     se.fit = FALSE, eco = FALSE, sandwich = FALSE,
                     astar = NULL, astar2 = NULL, x = NULL, w = NULL,
                     muhat = NULL, cmat = NULL, ipw = NULL) {
  
  n <- length(a)
  psi <- (y - muhat)*ipw
  
  if (eco) {
    
    if (is.null(weights))
      weights <- rep(1, times = length(a))
    
    # Gaussian Kernel
    a.std <- sqrt(weights)*(a - a.new) / bw
    k.std <- sqrt(weights)*dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  } else {
    
    a.std <- (a - a.new) / bw
    k.std <- dnorm(a.std) / bw
    g.std <- cbind(1, a.std)
    
  }
  
  b <- lm(psi ~ a.std, weights = k.std)$coefficients
  mu <- unname(b[1])
  
  if (se.fit & sandwich) {
    
    m <- ncol(cmat)
    l <- ncol(w)
    U <- matrix(0, ncol = m + l, nrow = m + l)
    V <- matrix(0, ncol = m + l + 2, nrow = 2)
    meat <- matrix(0, ncol = m + l + 2, nrow = m + l + 2)
    eta <- c(g.std%*%b)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - ipw[i] * tcrossprod(cmat[i,])
      U[(m + 1):(m + l),(m + 1):(m + l)] <- U[(m + 1):(m + l),(m + 1):(m + l)] - 
        weights[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(w[i,])
      
      V[,1:m] <- V[,1:m] - k.std[i]*psi[i]*tcrossprod(g.std[i,],cmat[i,])
      V[,(m + 1):(m + l)] <- V[,(m + 1):(m + l)] + k.std[i]*ipw[i]*
        family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(g.std[i,],w[i,])
      V[,(m + l + 1):(m + l + 2)] <- V[,(m + l + 1):(m + l + 2)] - k.std[i]*tcrossprod(g.std[i,])
      
      meat <- meat + 
        tcrossprod(esteq_kwls2(ipw = ipw[i], x = x[i,], w = w[i,], y = y[i], muhat = muhat[i],
                               g.std = g.std[i,], k.std = k.std[i], weights = weights[i],
                               astar = astar[i], astar2 = astar2[i], eta = eta[i]))
      
      
    }
    
    invbread <- matrix(0, nrow = m + l + 2, ncol = m + l + 2)
    invbread[1:(m + l),1:(m + l)] <- U
    invbread[(m + l + 1):(m + l + 2), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      Sig <- NA
      variance <- NA
      
    } else {
      
      Sig <- bread %*% meat %*% t(bread)
      variance <- Sig[m + l + 1, m + l + 1]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
  } else if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, k.std*g.std))
    V <- cbind((weights * (psi - eta)),
               (a.std * k.std * (psi - eta)))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(mu)
  
}

## Estimating equation for meat of sandwich estiamtor
esteq_kwls <- function(ipw, x, psi, g.std, k.std, astar, astar2, eta) {
  
  eq1 <- ipw*x*astar
  eq2 <- ipw*astar2
  eq3 <- ipw*x - x
  eq4 <- k.std*(psi - eta)*g.std
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}

## Estimating equation for meat of sandwich estiamtor
esteq_kwls2 <- function(ipw, x, w, y, g.std, k.std, 
                        astar, astar2, eta, muhat, weights) {
  
  psi <- ipw*(y - muhat)
  eq1 <- ipw*x*astar
  eq2 <- ipw*astar2
  eq3 <- ipw*x - x
  eq4 <- weights*(y - muhat)*w
  eq5 <- k.std*(psi - eta)*g.std
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}

## Leave-one-out cross-validated bandwidth
w.fn <- function(h, a, a.vals, n) {
  
  w.avals <- sapply(a.vals, function(a.tmp, ...) {
    a.std <- sqrt(n)*(a - a.tmp) / h
    k.std <- sqrt(n)*dnorm(a.std) / h
    return(mean(a.std^2 * k.std) * mean(sqrt(n)*dnorm(0) / h) /
             (mean(k.std) * mean(a.std^2 * k.std) - mean(a.std * k.std)^2))
  })
  
  return(w.avals / length(a))
  
}

hatvals <- function(h, a, a.vals, n) {
  approx(a.vals, w.fn(h = h, a = a, a.vals = a.vals, n = n), xout = a)$y
}

cts.eff.fn <- function(psi, a, h, a.vals, n) {
  approx(x = a.vals, 
         y = sapply(a.vals, kwls_est, a = a, weights = n, 
                    psi = psi, eco = TRUE, se.fit = FALSE, bw = h), 
         xout = a)$y
}

risk.fn <- function(h, psi, a, a.vals, n) {
  hats <- hatvals(h = h, a = a, a.vals = a.vals, n = n)
  sqrt(mean(((psi - cts.eff.fn(psi = psi, a = a, h = h, a.vals = a.vals, n = n)) / (1 - hats))^2, na.rm = TRUE))
}
