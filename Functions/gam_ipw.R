gam_ipw <- function(a, y, family = gaussian(), ipw, weights = NULL,
                    a.vals = seq(min(a), max(a), length.out = 100), se.fit = FALSE, 
                    x = NULL, astar = NULL, astar2 = NULL, cmat = NULL) {
  n <- length(a)
  psi <- ipw*y
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # GAM Models
  mod <- gam(psi ~ s(a), weights = weights, family = family)
  
  # Naive Variance
  # if (se.fit) {
  #   pred <- predict(mod, newdata = data.frame(a = a.vals), se.fit = TRUE, type = "response")
  #   return(list(mu = pred$fit, sig2 = (pred$se.fit)^2))
  # } else {
  #   return(predict(mod, newdata = data.frame(a = a.vals), se.fit = FALSE, type = "response"))
  # }
  
  # Robust Variance
  g <- predict(mod, type = "lpmatrix")
  mu <- family$linkinv(c(g %*% mod$coefficients))
  g.vals <- predict(mod, type = "lpmatrix", newdata = data.frame(a = a.vals))
  mu.vals <- family$linkinv(c(g.vals %*% mod$coefficients))
  
  if (se.fit) {
    
    m <- ncol(cmat)
    o <- ncol(g)
    U <- matrix(0, ncol = m, nrow = m)
    V <- matrix(0, ncol = m + o, nrow = o)
    meat <- matrix(0, ncol = m + o, nrow = m + o)
    
    for (i in 1:n) {
      
      U[1:m,1:m] <- U[1:m,1:m] - weights[i]*ipw[i]*tcrossprod(cmat[i,])
      V[,1:m] <- V[,1:m] - weights[i]*psi[i]*tcrossprod(g[i,],cmat[i,])
      V[,(m + 1):(m + o)] <- V[,(m + 1):(m + o)] - weights[i]*family$mu.eta(family$linkfun(mu[i]))*tcrossprod(g[i,])
      
      meat <- meat + 
        tcrossprod(esteq_gam_ipw(y = y[i], x = x[i,], g = g[i,],
                                 ipw = ipw[i], weights = weights[i],
                                 astar = astar[i], astar2 = astar2[i], mu = mu[i]))
      
      
    }
    
    invbread <- matrix(0, nrow = m + o, ncol = m + o)
    invbread[1:m,1:m] <- U
    invbread[(m + 1):(m + o), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sig2.vals <- rep(NA, length(a.vals))
      
    } else {
      
      Sigma <- bread %*% meat %*% t(bread)
      del.vals <- family$mu.eta(family$linkfun(mu.vals))
      sig2.vals <- diag((del.vals*g.vals) %*% Sigma[(m + 1):(m + o), (m + 1):(m + o)] %*% t(g.vals*del.vals))
      
    }
    
    return(rbind(mu.vals = mu.vals, sig2.vals = sig2.vals))
    
  } else
    return(mu)
  
}

## Estimating equation for meat of sandwich estiamtor
esteq_gam_ipw <- function(y, x, g, weights, ipw,
                         astar, astar2, mu) {
  
  psi <- ipw*y
  eq1 <- weights*ipw*x*astar
  eq2 <- weights*ipw*astar2
  eq3 <- weights*(ipw*x - x)
  eq4 <- weights*(psi - mu)*g
  
  eq <- c(eq1, eq2, eq3, eq4) 
  return(eq)
  
}
