## GAM Estimation of the ERCs
gam_est <- function(a, psi, family = gaussian(), weights = NULL, 
                    a.vals = seq(min(a), max(a), length.out = 100), se.fit = FALSE,
                    astar = NULL, astar2 = NULL, x = NULL, cmat = NULL, ipw = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # GAM Models
  mod <- gam(psi ~ s(a), weights = weights, family = family)
  
  g <- predict(mod, type = "lpmatrix")
  mu <- predict(mod, type = "response")
  g.vals <- predict(mod, type = "lpmatrix", newdata = data.frame(a = a.vals))
  mu.vals <- predict(mod, type = "response", newdata = data.frame(a = a.vals))
  
  if (se.fit) {
    
    m <- ncol(cmat)
    l <- ncol(g)
    U <- matrix(0, ncol = m, nrow = m)
    V <- matrix(0, ncol = m + l, nrow = l)
    meat <- matrix(0, ncol = m + l, nrow = m + l)
    eta <- mod$fitted.values
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - ipw[i]*tcrossprod(cmat[i,])
      V[,1:m] <- V[,1:m] - weights[i]*psi[i]*tcrossprod(g.std[i,],cmat[i,])
      V[,(m + 1):(m + l)] <- V[,(m + 1):(m + l)] - weights[i]*family$mu.eta(family$linkfun(mu))*tcrossprod(g.std[i,])
      
      meat <- meat + 
        tcrossprod(esteq_gam(p = ipw[i], x = x[i,], psi = psi[i],
                             g.std = g.std[i,], weights = weights[i],
                             astar = astar[i], astar2 = astar2[i], eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = m + l, ncol = m + l)
    invbread[1:m,1:m] <- U
    invbread[(m + 1):(m + l), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      Sig <- NA
      variance(rep(NA, length(a.vals)))
      
    } else {
      
      Sig <- bread %*% meat %*% t(bread)
      BV <- Sig[(m + 1):(m + l),(m + 1):(m + l)]
      variance <- diag((family$mu.eta(family$linkfun(mu.vals))*g.vals) %*% BV %*% t(family$mu.eta(family$linkfun(mu.vals))*g.vals))
      
    }
    
    return(rbind(mu = mu, sig2 = variance))
    
  } else
    return(mu)
  
}

esteq_gam <- function(p, x, psi, g.std, weights, astar, astar2, eta) {
  
  eq1 <- p*x*astar
  eq2 <- p*astar2
  eq3 <- p*x - x
  eq4 <- weights*(psi - eta)*g
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}