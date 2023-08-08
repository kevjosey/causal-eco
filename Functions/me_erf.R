
# Kernel weighted least squares
kern_naive <- function(a.new, a, x, psi, astar, astar2, cmat, ipw, bw = 1, se.fit = FALSE) {
  
  n <- length(a)
  m <- ncol(x)
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(a.std, 1)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std)$coefficients
  mu <- unname(b[2])
  
  if (se.fit) {
    
    U <- matrix(0, ncol = m + 1, nrow = m + 1)
    V <- matrix(0, ncol = m + 3, nrow = 2)
    meat <- matrix(0, ncol = m + 3, nrow = m + 3)
    eta <- c(g.std%*%b)
    
    for (i in 1:n) {
      
      U[1:(m + 1),1:(m + 1)] <- U[1:(m + 1),1:(m + 1)] - ipw[i] * tcrossprod(cmat[i,])
      
      V[,1:(m + 1)] <- V[,1:(m + 1)] - k.std[i]*psi[i]*tcrossprod(g.std[i,],cmat[i,])
      V[,(m + 2):(m + 3)] <- V[,(m + 2):(m + 3)] - k.std[i]*tcrossprod(g.std[i,])
      
      meat <- meat + tcrossprod(esteq_naive(p = ipw[i], x = x[i,], psi = psi[i],
                                            g.std = g.std[i,], k.std = k.std[i],
                                            astar = astar[i], astar2 = astar2[i], 
                                            eta = eta[i]))
      
    }
    
    invbread <- matrix(0, nrow = m + 3, ncol = m + 3)
    invbread[1:(m + 1),1:(m + 1)] <- U
    invbread[(m + 2):(m + 3), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[m + 3, m + 3]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
  } else
    return(mu)
  
}

esteq <- function(p, x, psi, g.std, k.std, astar, astar2, eta) {
  
  eq1 <- p*x*astar
  eq2 <- p*astar2
  
  eq3 <- k.std*(psi - eta)*g.std
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}