
multi_blp <- function(w, w.id) {
  
  wts <- c(unname(table(w.id)))
  
  # initialize exposures
  z_tmp <- aggregate(w, by = list(w.id), mean)
  id <- c(z_tmp[,1])
  z <- as.matrix(z_tmp[,-1])
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  
  ord <- order(w.id)
  z_tmp <- z[rep(1:nrow(z), wts),]
  z_w <- z_tmp[order(ord),]
  
  p <- ncol(w)

  mu_z <- colSums(wts*z)/m
  muMat_z <- matrix(rep(mu_z, n), nrow = n, byrow = TRUE)
  nu <- m - sum(wts^2)/m
    
  Omega <- crossprod(as.matrix(w - z_w))/(m - n)
  Sigma <- as.matrix(crossprod(wts*(z - muMat_z), (z - muMat_z)) - (n - 1)*Omega)/nu
    
  a <- t(sapply(1:n, function(i, ...) {
      
    V <- Sigma + Omega/wts[i]
    out <- c(mu_z + Sigma %*% solve(V) %*% c(t(z[i,]) - mu_z))
    
    return(out)
      
  }))
  
  colnames(a) <- colnames(w)
  return(data.frame(id = id, a))
  
}
