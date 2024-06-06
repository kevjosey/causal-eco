gen_data <- function(n, sig = sqrt(2), nonlinear = FALSE,
                     scenario = c("base", "ps-mis", "out-mis", "mis")){
  
  # alternative
  x1 <- scale(runif(n, -1, 1))
  x2 <- scale(runif(n, 0, 1))
  # u1 <- scale(exp(1/2*x1 + x2))
  # u2 <- scale((x1 + x2)^2)
  u1 <- (2*x2 - 1)/(exp(1 + x1))
  u2 <- (exp(2*x1) - exp(-2)) / (exp(2) - exp(-2))
  
  X <- cbind(int = rep(1, n), x1, x2)
  colnames(X) <- c("(int)", "x1", "x2")
  U <- cbind(int = rep(1, n), u1, u2)
  
  if (scenario == "base"){
    
    a <- rnorm(n, 8 + 2*x1 - 2*x2, sig)
    mu <- plogis(-4 + 0.6*a - 0.2*a*x2 + 0.4*a*x1)
    mu_scale <- plogis(-4 + 0.6*(a*0.9) - 0.2*(a*0.9)*x2 + 0.4*(a*0.9)*x1)
    mu_diff <- plogis(-4 + 0.6*(a - 1) - 0.2*(a - 1)*x2 + 0.4*(a - 1)*x1)
    mu_thresh <- plogis(-4 + 0.6*pmin(a, 8) - 0.2*pmin(a, 8)*x2 + 0.4*pmin(a, 8)*x1)
    
  } else if (scenario == "ps-mis"){
    
    a <- rnorm(n, 8 + 2*u1 - 2*u2, sig)
    mu <- plogis(-4 + 0.6*a - 0.2*a*x2 + 0.4*a*x1)
    mu_scale <- plogis(-4 + 0.6*(a*0.9) - 0.2*(a*0.9)*x2 + 0.4*(a*0.9)*x1)
    mu_diff <- plogis(-4 + 0.6*(a - 1) - 0.2*(a - 1)*x2 + 0.4*(a - 1)*x1)
    mu_thresh <- plogis(-4 + 0.6*pmin(a, 8) - 0.2*pmin(a, 8)*x2 + 0.4*pmin(a, 8)*x1)
    
  } else if (scenario == "out-mis") {
    
    a <- rnorm(n, 8 + 2*x1 - 2*x2, sig)
    mu <- plogis(-4 + 0.6*a - 0.2*a*u2 + 0.4*a*u1)
    mu_scale <- plogis(-4 + 0.6*(a*0.9) - 0.2*(a*0.9)*u2 + 0.4*(a*0.9)*u1)
    mu_diff <- plogis(-4 + 0.6*(a - 1) - 0.2*(a - 1)*u2 + 0.4*(a - 1)*u1)
    mu_thresh <- plogis(-4 + 0.6*pmin(a, 8) - 0.2*pmin(a, 8)*u2 + 0.4*pmin(a, 8)*u1)
    
  } else if (scenario == "mis"){
    
    a <- rnorm(n, 8 + 2*u1 - 2*u2, sig)
    mu <- plogis(-4 + 0.6*a - 0.2*a*u2 + 0.4*a*u1)
    mu_scale <- plogis(-4 + 0.6*(a*0.9) - 0.2*(a*0.9)*u2 + 0.4*(a*0.9)*u1)
    mu_diff <- plogis(-4 + 0.6*(a - 1) - 0.2*(a - 1)*u2 + 0.4*(a - 1)*u1)
    mu_thresh <- plogis(-4 + 0.6*pmin(a, 8) - 0.2*pmin(a, 8)*u2 + 0.4*pmin(a, 8)*u1)
    
  }
  
  if (nonlinear) {
    
    mu <- mu + cos(pi*a/4)
    mu_scale <- mu_scale + cos(pi*(a*0.9)/2)
    mu_diff <- mu_diff + cos(pi*(a - 1)/2)
    mu_thresh <- mu_thresh + cos(pi*pmin(a, 8)/2)
    
  }
  
  theta_scale <- mean(mu_scale - mu)
  theta_diff <- mean(mu_diff - mu)
  theta_thresh <- mean(mu_thresh - mu)
  
  y <- rbinom(n, 1, mu)
  
  # create simulation dataset
  sim <- list(y = y, a = a, X = X, U = U, 
              theta_scale = theta_scale, 
              theta_diff = theta_diff, 
              theta_thresh = theta_thresh)
  
  return(sim)
  
}

# Fits the RKHS balancing weights using a variety of kernels
sim_fit <- function(idx = 1, simDat, ...) {
  
  print(idx)
  
  dat <- simDat[,idx]
  
  # data components
  Y <- dat$y
  A <- dat$a
  X <- data.frame(dat$X[,-1])
  
  # Counterfactural Exposures
  A.scale <- 0.9*A
  A.diff <- A - 1
  A.thresh <- pmin(A, 8)
  
  n <- nrow(X)
  
  theta_scale <- c(dat$theta_scale)
  theta_diff <- c(dat$theta_diff)
  theta_thresh <- c(dat$theta_thresh)
  
  # entropy balancing
  linear_scale <- ebcf(Y = Y, X = X, A = A, A.new = A.scale, family = binomial(),
                   eps_abs = 1e-6, eps_rel = 1e-6, verbose = FALSE)
  linear_diff <- ebcf(Y = Y, X = X, A = A, A.new = A.diff, family = binomial(),
                  eps_abs = 1e-6, eps_rel = 1e-6, verbose = FALSE)
  linear_thresh <- ebcf(Y = Y, X = X, A = A, A.new = A.thresh, family = binomial(),
                    eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE)

  
  # LMTP
  lmtp_none <- lmtp_tmle( data = data.frame(A = A, Y = Y, X),
                          trt = "A", outcome = "Y", 
                          baseline = c("x1", "x2"),
                          outcome_type = "binomial",
                          shifted = data.frame(A = A, Y = Y, X),
                          mtp = TRUE, folds = 1,
                          learners_outcome = c("SL.mean", "SL.glm", "SL.glm.interaction"),
                          learners_trt = c("SL.mean","SL.glm"))
  
  lmtp_scale <- lmtp_tmle( data = data.frame(A = A, Y = Y, X),
                           trt = "A", outcome = "Y", 
                           baseline = c("x1", "x2"),
                           outcome_type = "binomial",
                           shifted = data.frame(A = A.scale, Y = Y, X),
                           mtp = TRUE, folds = 1,
                           learners_outcome = c("SL.mean", "SL.glm", "SL.glm.interaction"),
                           learners_trt = c("SL.mean","SL.glm"))
  
  lmtp_diff <- lmtp_tmle( data = data.frame(A = A, Y = Y, X),
                          trt = "A", outcome = "Y", 
                          baseline = c("x1", "x2"),
                          outcome_type = "binomial",
                          shifted = data.frame(A = A.diff, Y = Y, X),
                          mtp = TRUE, folds = 1,
                          learners_outcome = c("SL.mean", "SL.glm", "SL.glm.interaction"),
                          learners_trt = c("SL.mean","SL.glm"))
  
  lmtp_thresh <- lmtp_tmle( data = data.frame(A = A, Y = Y, X),
                           trt = "A", outcome = "Y",
                           baseline = c("x1", "x2"),
                           outcome_type = "binomial",
                           shifted = data.frame(A = A.thresh, Y = Y, X),
                           mtp = TRUE, fodls = 1,
                           learners_outcome = c("SL.glm.interaction"),
                           learners_trt = c("SL.glm","SL.mean"))
  
  lmtp_contrast_scale = lmtp_contrast(lmtp_scale, ref = lmtp_none)
  lmtp_contrast_diff = lmtp_contrast(lmtp_diff, ref = lmtp_none)
  lmtp_contrast_thresh = lmtp_contrast(lmtp_thresh, ref = lmtp_none)
  
  # combine results
  est_scale <- c(theta = theta_scale,
                 linear = linear_scale$theta,
                 lmtp = lmtp_contrast_scale$vals$theta)
  est_diff <- c(theta = theta_diff,
                linear = linear_diff$theta,
                lmtp = lmtp_contrast_diff$vals$theta)
  est_thresh <- c(theta = theta_thresh,
                  linear = linear_thresh$theta,
                  lmtp = lmtp_contrast_thresh$vals$theta)
  
  se_scale <- c(linear = sqrt(linear_scale$omega2),
                mlinear = sqrt(linear_scale$sig2),
                lmtp = lmtp_contrast_scale$vals$std.error)
  se_diff <- c(linear = sqrt(linear_diff$omega2),
               mlinear = sqrt(linear_diff$sig2),
               lmtp = lmtp_contrast_diff$vals$std.error)
  se_thresh <- c(linear = sqrt(linear_thresh$omega2),
                 mlinear = sqrt(linear_thresh$sig2),
                 lmtp = lmtp_contrast_thresh$vals$std.error)
  
  est = data.frame(rbind(est_scale, est_diff, est_thresh), 
                   estimate = c("scalar", "difference", "threshold"))
  se = data.frame(rbind(se_scale, se_diff, se_thresh), 
                  estimate = c("scalar", "difference", "threshold"))
  rownames(est) <- rownames(se) <- NULL
  
  return(list(est = est, se = se))
  
}

