library(Matrix)
library(kernlab)
library(osqp)
library(scam)
library(SuperLearner)
library(earth)
library(glmnet)
library(lmtp)
library(dplyr)
library(parallel)
library(ggplot2)
library(tidyr)
library(stringr)

# Extra code for SRFs
source("~/Github/causal-eco/Functions/calibrate.R")

## Function for Generating Data
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

# Fits LMTP + Balancing Weights
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

## Run the models

scenarios <- expand.grid(n = c(500, 1000, 2000),
                         mis = c("base", "ps-mis", "out-mis"),
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)

n.iter <- 1000 # simulation iterations

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  # run simulation study
  simDat <- replicate(n.iter, gen_data(n = scenario$n, scenario = scenario$mis, sig = sqrt(2)))
  simFit <- mclapply(1:n.iter, sim_fit, simDat = simDat, mc.cores = 25)
  
  # summarize
  est_mat <- do.call(rbind, lapply(simFit, function(lst, ...) lst$est))
  se_mat <- do.call(rbind, lapply(simFit, function(lst, ...) lst$se))
  ci_mat <- data.frame(linear_lower = est_mat[,2] - 1.96*se_mat[,1],
                       linear_upper = est_mat[,2] + 1.96*se_mat[,1],
                       mlinear_lower = est_mat[,2] - 1.96*se_mat[,2],
                       mlinear_upper = est_mat[,2] + 1.96*se_mat[,2],
                       lmtp_lower = est_mat[,3] - 1.96*se_mat[,3],
                       lmtp_upper = est_mat[,3] + 1.96*se_mat[,3],
                       theta = est_mat[,1], estimate = est_mat[,4])
  
  ci_diff <- data.frame(linear = ci_mat$linear_upper - ci_mat$linear_lower,
                        mlinear = ci_mat$mlinear_upper - ci_mat$mlinear_lower,
                        lmtp = ci_mat$lmtp_upper - ci_mat$lmtp_lower,
                        estimate = est_mat[,4])
  
  # results
  mu <- est_mat %>% group_by(estimate) %>%
    summarize(linear = mean(linear), lmtp = mean(lmtp))
  
  bias <- est_mat %>% group_by(estimate) %>% 
    summarize(linear = mean(linear - theta)/mean(theta),
              lmtp = mean(lmtp - theta)/mean(theta))
  
  rmse <- est_mat %>% group_by(estimate) %>% 
    summarize(linear = sqrt(mean((linear - theta)^2)),
              lmtp = sqrt(mean((lmtp - theta)^2)))
  
  cp <- ci_mat %>% group_by(estimate) %>% 
    summarize(cp_linear = mean(linear_lower < mean(theta) & linear_upper > mean(theta)),
              cp_mlinear = mean(mlinear_lower < mean(theta) & mlinear_upper > mean(theta)),
              cp_lmtp = mean(lmtp_lower < mean(theta) & lmtp_upper > mean(theta)))
  
  ci_length <- ci_diff %>% group_by(estimate) %>% 
    summarize(linear = mean(linear), mlinear = mean(mlinear), lmtp = mean(lmtp))
  
  results <- list(mu = mu, bias = bias, rmse = rmse,
                  cp = cp, ci_length = ci_length, 
                  est_mat = est_mat, se_mat = se_mat)
  
  save(results, file = paste0("~/Github/causal-eco/Output/si_", scenario$n, "_", scenario$mis, ".RData"))
  
}

# create a single table and plot
dat <- data.frame()
box <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  load(file = paste0("~/Github/causal-eco/Output/si_", scenario$n, "_", scenario$mis, ".RData"))
  
  dat.tmp <- data.frame(n = rep(scenario$n, 9),
                        mis = rep(scenario$mis, 9),
                        method = rep(c("TMLE", "Calibration", "Calibration"), each = 3),
                        se_method = rep(c("EIF", "EIF", "M-Estimation"), each = 3),
                        estimand = rep(results$bias$estimate, times = 3),
                        bias = c(results$bias$lmtp, results$bias$linear, rep(NA, 3)),
                        rmse = c(results$rmse$lmtp, results$rmse$linear, rep(NA, 3)),
                        cp = c(results$cp$cp_lmtp, results$cp$cp_linear, results$cp$cp_mlinear),
                        ci_length = c(results$ci_length$lmtp, results$ci_length$linear, results$ci_length$mlinear))
  
  box.tmp <- data.frame(n = rep(scenario$n, 2*nrow(results$est_mat)),
                        mis = rep(scenario$mis, 2*nrow(results$est_mat)),
                        method = rep(c("TMLE", "Calibration"), each = nrow(results$est_mat)),
                        estimand = rep(results$est_mat$estimate, times = 2),
                        bias = c(results$est_mat$lmtp - results$est_mat$theta,
                                 results$est_mat$linear - results$est_mat$theta))
  
  dat <- rbind(dat, dat.tmp)
  box <- rbind(box, box.tmp)
  
}

# Cleanup Box Plot
box_long <- box %>%
  mutate(estimand=replace(estimand,estimand=="threshold", "Threshold"),
         estimand=replace(estimand,estimand=="scalar", "Scaled Shift"),
         estimand=replace(estimand,estimand=="difference", "Additive Shift")) %>%
  mutate(mis=replace(mis,mis=="ps-mis", "Propensity Score"),
         mis=replace(mis,mis=="out-mis", "Outcome Model"),
         mis=replace(mis,mis=="base", "Both Models Correct"))

# Box Plot
grid_plot <- box_long %>% filter(n == 1000) %>%
  ggplot(aes(x = mis, y = bias, color = method)) +
  geom_boxplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid( ~ as.factor(estimand), scales='free') +
  theme_bw() + 
  labs(x='Misspecification',
       y='Bias',
       color='Method') +
  theme(legend.position='bottom',
        legend.key.height = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold")) ; 

pdf("~/Github/causal-eco/Output/srf-simulation-plot.pdf", width = 12, height = 6)
grid_plot
dev.off()

