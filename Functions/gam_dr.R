## GAM Estimation of the ERCs with weights for ecological regression
gam_dr <- function(a, y, family = gaussian(), weights = NULL,
                    a.vals = seq(min(a), max(a), length.out = 100), 
                    se.fit = FALSE, k = -1, ipw, muhat, eco = TRUE,
                    x = NULL, w = NULL, astar = NULL, astar2 = NULL) {
  
  # regression objects
  n <- length(a)
  psi <- ipw*(y - muhat)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  if (!eco) {
    q <- weights
  } else {
    q <- rep(1, times = length(a))
  }
  
  # GAMs
  mod <- scam(psi ~ s(a, bs = "tp", k = k), data = data.frame(a = a, psi = psi),
              weights = weights, family = gaussian()) # needs to be gaussian because of negative values

  # predictions
  mu.vals <- predict(mod, newdata = data.frame(a = a.vals), type = "response")
  
  # Robust Variance
  if (se.fit) {
    
    # more predictions
    g <- predict(mod, type = "lpmatrix")
    mu <- predict(mod, type = "response")
    g.vals <- predict(mod, newdata = data.frame(a = a.vals), type = "lpmatrix")
    
    # ipw constraint matrix
    cmat <- cbind(x*astar, astar2, x)
    
    # dimensions
    m <- ncol(cmat)
    l <- ncol(w)
    o <- ncol(g)
    
    # initialize matrices
    U <- matrix(0, ncol = m + l, nrow = m + l)
    V <- matrix(0, ncol = m + l + o, nrow = o)
    meat <- matrix(0, ncol = m + l + o, nrow = m + l + o)
    
    # sandwich estimator mess
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - q[i]*ipw[i]*tcrossprod(cmat[i,])
      U[(m + 1):(m + l),(m + 1):(m + l)] <- U[(m + 1):(m + l),(m + 1):(m + l)] - 
        weights[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(w[i,])
      
      V[,1:m] <- V[,1:m] - weights[i]*psi[i]*tcrossprod(g[i,],cmat[i,])
      V[,(m + 1):(m + l)] <- V[,(m + 1):(m + l)] - weights[i]*ipw[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(g[i,],w[i,])
      V[,(m + l + 1):(m + l + o)] <- V[,(m + l + 1):(m + l + o)] - weights[i]*tcrossprod(g[i,])
      
      meat <- meat + 
        tcrossprod(esteq_gam_dr(y = y[i], x = x[i,], w = w[i,], g = g[i,],
                                ipw = ipw[i], muhat = muhat[i], p = weights[i], q = q[i],
                                astar = astar[i], astar2 = astar2[i], mu = mu[i]))
      
    }
    
    # bread matrix
    invbread <- matrix(0, nrow = m + l + o, ncol = m + l + o)
    invbread[1:(m + l),1:(m + l)] <- U
    invbread[(m + l + 1):(m + l + o), ] <- V
    bread <- try(solve(invbread), silent = TRUE)
    
    # sandwich variance
    if (inherits(bread, "try-error")) {
      Sig <- NULL
    } else {
      Sigma <- bread %*% meat %*% t(bread)
      Sig <- Sigma[(m + 1):(m + l + o),(m + 1):(m + l + o)]
    }
    
    return(list(mu.vals = mu.vals, Sig = Sig, g.vals = g.vals))
    
  } else
    return(mu.vals)
  
}

## Estimating Equations for Robust Variance
esteq_gam_dr <- function(y, x, w, g, p, q,
                         ipw, muhat, astar, astar2, mu) {
  
  psi <- ipw*(y - muhat)
  
  eq1 <- q*(ipw*x*astar)
  eq2 <- q*ipw*astar2
  eq3 <- q*(ipw*x - x)
  eq4 <- p*(y - muhat)*w
  eq5 <- p*(psi - mu)*g
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}