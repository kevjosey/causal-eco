## GAM Estimation of the ERCs with weights for ecological regression
gam_om <- function(a, y, family = gaussian(), muhat, weights = NULL,
                   a.vals = seq(min(a), max(a), length.out = 100), 
                   se.fit = FALSE,  w = NULL) {
  n <- length(a)
  psi <- (y - muhat)
  
  if (is.null(weights))
    weights <- rep(1, times = length(a))
  
  # GAM Models
  mod <- gam(psi ~ s(a), weights = weights, family = gaussian()) # needs to be gaussian because of negative values
  
  # Naive Variance
  # if (se.fit) {
  #   pred <- predict(mod, newdata = data.frame(a = a.vals), se.fit = TRUE, type = "response")
  #   return(list(mu = pred$fit, sig2 = (pred$se.fit)^2))
  # } else {
  #   return(predict(mod, newdata = data.frame(a = a.vals), se.fit = FALSE, type = "response"))
  # }
  
  # Robust Variance
  g <- predict(mod, type = "lpmatrix")
  eta <- c(g %*% mod$coefficients)
  g.vals <- predict(mod, type = "lpmatrix", newdata = data.frame(a = a.vals))
  eta.vals <- c(g.vals %*% mod$coefficients)
  
  if (se.fit) {
    
    l <- ncol(w)
    o <- ncol(g)
    U <- matrix(0, ncol = l, nrow = l)
    V <- matrix(0, ncol = l + o, nrow = o)
    meat <- matrix(0, ncol = l + o, nrow = l + o)
    
    for (i in 1:n) {
      
      U[1:l,1:l] <- U[1:l,1:l] - weights[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(w[i,])
      V[,1:l] <- V[,1:l] - weights[i]*family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(g[i,],w[i,])
      V[,(l + 1):(l + o)] <- V[,(l + 1):(l + o)] - weights[i]*tcrossprod(g[i,])
      
      meat <- meat + 
        tcrossprod(esteq_gam_om(y = y[i], w = w[i,], g = g[i,],
                                muhat = muhat[i], weights = weights[i], eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = l + o, ncol = l + o)
    invbread[1:l,1:l] <- U
    invbread[(l + 1):(l + o), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      Sig <- NULL
      
    } else {
      
      Sig <- bread %*% meat %*% t(bread)
      
    }
    
    return(list(eta.vals = eta.vals, Sig = Sig, g.vals = g.vals))
    
  } else
    return(eta.vals)
  
}

# Estimating Equations for Robust Variance
esteq_gam_om <- function(y, w, g, weights, muhat, eta) {
  
  psi <- (y - muhat)
  eq1 <- weights*(y - muhat)*w
  eq2 <- weights*(psi - eta)*g
  
  eq <- c(eq1, eq2) 
  return(eq)
  
}