## GAM Estimation of the ERCs
gam_est <- function(a, y, family = gaussian(), weights = NULL, se.fit = FALSE, 
                    a.vals = seq(min(a), max(a), length.out = 100),
                    ipw = NULL, muhat = NULL, x = NULL, w = NULL,
                    astar = NULL, astar2 = NULL, cmat = NULL) {
  
  n <- length(a)
  psi <- ipw*(y - muhat)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # GAM Models
  mod <- gam(psi ~ s(a), weights = weights, family = gaussian())
  
  g <- predict(mod, type = "lpmatrix")
  mu <- predict(mod, type = "response")
  g.vals <- predict(mod, type = "lpmatrix",
                    newdata = data.frame(a = a.vals),
                    newdata.guaranteed = TRUE)
  mu.vals <- predict(mod, type = "response", 
                     newdata = data.frame(a = a.vals),
                     newdata.guaranteed = TRUE)
  
  if (se.fit) {
    
    m <- ncol(cmat)
    l <- ncol(w)
    o <- ncol(g)
    U <- matrix(0, ncol = m + l, nrow = m + l)
    V <- matrix(0, ncol = m + l + o, nrow = o)
    meat <- matrix(0, ncol = m + l + o, nrow = m + l + o)
    eta <- mod$fitted.values
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - ipw[i] * tcrossprod(cmat[i,])
      U[(m + 1):(m + l),(m + 1):(m + l)] <- U[(m + 1):(m + l),(m + 1):(m + l)] - 
        weights[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(w[i,])
      
      V[,1:m] <- V[,1:m] - weights[i]*psi[i]*tcrossprod(g[i,],cmat[i,])
      V[,(m + 1):(m + l)] <- V[,(m + 1):(m + l)] - weights[i]*ipw[i]*
        family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(g[i,],w[i,])
      V[,(m + l + 1):(m + l + o)] <- V[,(m + l + 1):(m + l + o)] - weights[i]*tcrossprod(g[i,])
      
      meat <- meat + 
        tcrossprod(esteq_gam(y = y[i], x = x[i,], w = w[i,], g = g[i,],
                             ipw = ipw[i], muhat = muhat[i], weights = weights[i],
                             astar = astar[i], astar2 = astar2[i], eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = m + l + o, ncol = m + l + o)
    invbread[1:(m + l),1:(m + l)] <- U
    invbread[(m + l + 1):(m + l + o), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      Sig <- NA
      variance(rep(NA, length(a.vals)))
      
    } else {
      
      Sig <- bread %*% meat %*% t(bread)
      BV <- Sig[(m + 1):(m + l + o),(m + 1):(m + l + o)]
      
    }
    
    return(list(mu = mu.vals, Sig = BV, g.vals = g.vals))
    
  } else
    return(mu.vals)
  
}

esteq_gam <- function(y, x, w, g,
                      ipw, muhat, weights,
                      astar, astar2, eta) {
  
  psi <- ipw*(y - muhat)
  eq1 <- ipw*x*astar
  eq2 <- ipw*astar2
  eq3 <- ipw*x - x
  eq4 <- weights*(y - muhat)*w
  eq5 <- weights*(psi - eta)*g
  
  eq <- c(eq1, eq2, eq3, eq4, eq5) 
  return(eq)
  
}