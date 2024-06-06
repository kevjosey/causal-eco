# generic calibration function
calibrate <- function(cmat, target, base_weights = NULL, coefs_init = NULL,
                      optim_ctrl = list(maxit = 500, reltol = 1e-6), ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun(lagrange_ent)
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, nrow(cmat))
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS", hessian = TRUE,
                      cmat = cmat, base_weights = base_weights,
                      target = target, control = optim_ctrl)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  weights <- c( base_weights*exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

# entropy objective function
lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

ebcf <- function(Y, X, A, A.new, id = NULL, family = gaussian(),
                 base_weights = rep(1, length(A)),
                 eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE, ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  X.mat <- model.matrix(~ ., data = X)
  n <- nrow(X.mat)
  
  if (is.null(id))
    id <- 1:nrow(X.mat)
  
  # Margins and matrix constraints
  Astar <- c(A - mean(A))/var(A)
  Astar2 <- c((A - mean(A))^2/var(A) - 1)
  W1 <- cbind(X.mat*Astar, Astar2, X.mat)
  Atilde <- c(A.new - mean(A))/var(A)
  Atilde2 <- c(A.new - mean(A))^2/var(A) - 1
  W0 <- cbind(X.mat*Atilde, Atilde2, X.mat)
  
  # set optimization settings
  settings <- list(reltol = eps_rel,
                   abstol = eps_abs)
  
  # construct linear term vector
  solution <- calibrate(cmat = W1, target = colSums(W0), 
                        base_weights = base_weights, 
                        optim_ctrl = settings)
  weights <- solution$weights
  
  # estimate nuisance outcome model with GAM
  fmla <- as.formula(paste0(c("y ~ s(a, bs = 'tp')", colnames(X), paste0("a:", colnames(X))), collapse = " + "))
  mumod <- scam(fmla, data = data.frame(a = A, y = Y, X), family = family)
  muhat <- predict(mumod, type = "response")
  mutilde <- predict(mumod, newdata = data.frame(a = A.new, X), type = "response")
  lphat <- predict(mumod, type = "lpmatrix")
  lptilde <- predict(mumod, newdata = data.frame(a = A.new, X), type = "lpmatrix")
  
  # Point Estimates
  psi <- c(Y - muhat)*weights + mutilde
  theta <- mean(psi - Y)
  
  # compute imbalances
  imbalance <- c(colMeans(W0)) - c(t(W1) %*% weights)/n
  names(imbalance) <- colnames(W1)
  
  # Robust Variance
  
  # dimensions
  m <- ncol(W1)
  l <- ncol(lphat)
  
  # initialize matrices
  U <- matrix(0, ncol = m + l, nrow = m + l)
  V <- vector(mode = "numeric", length = m + l + 1)
  meat <- matrix(0, ncol = m + l + 1, nrow = m + l + 1)
  
  # sandwich estimator mess
  for (i in 1:n) {
    
    U[1:m,1:m] <- U[1:m,1:m] - weights[i]*tcrossprod(W1[i,])
    U[(m + 1):(m + l),(m + 1):(m + l)] <- U[(m + 1):(m + l),(m + 1):(m + l)] - 
      family$mu.eta(family$linkfun(muhat[i]))*tcrossprod(lphat[i,])
    
    V[1:m] <- V[1:m] - c(weights[i]*(Y[i] - muhat[i])*W1[i,])
    V[(m + 1):(m + l)] <- V[(m + 1):(m + l)] - c(weights[i]*family$mu.eta(family$linkfun(muhat[i]))*lphat[i,]) +
      c(family$mu.eta(family$linkfun(mutilde[i]))*lptilde[i,])
    V[(m + l + 1)] <- V[(m + l + 1)] - 1
    
    meat <- meat + 
      tcrossprod(esteq(y = Y[i], w0 = W0[i,], w1 = W1[i,], lp = lphat[i,],
                       weights = weights[i], muhat = muhat[i], psi = psi[i], theta = theta))
    
  }
  
  # bread matrix
  invbread <- matrix(0, nrow = m + l + 1, ncol = m + l + 1)
  invbread[1:(m + l),1:(m + l)] <- U
  invbread[(m + l + 1), ] <- V
  bread <- try(solve(invbread), silent = TRUE)
  
  # sandwich variance
  if (inherits(bread, "try-error")) {
    Sig <- NULL
    sig2 <- NULL
  } else {
    Sigma <- bread %*% meat %*% t(bread)
    sig2 <- Sigma[(m + l + 1),(m + l + 1)]
  }
  
  # Efficient IF
  theta_eif <- psi - Y - theta
  
  # variances (needs work)
  m <- length(unique(id))
  omega2 <- var(vapply(split(theta_eif, id), function(z) mean(z), 1))/m
  
  # return output
  return(list(theta = theta, sig2 = sig2, omega2 = omega2,
              weights = weights, imbalance = imbalance, weights = weights))
  
}

esteq <- function(y, w0, w1, lp, ipw, weights, muhat, psi, theta) {
  
  eq1 <- weights*w1 - w0
  eq2 <- (y - muhat)*lp
  eq3 <- psi - y - theta
  
  eq <- c(eq1, eq2, eq3) 
  return(eq)
  
}
