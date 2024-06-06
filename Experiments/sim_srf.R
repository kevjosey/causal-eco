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

rm(list = ls())
source("~/Github/causal-eco/Functions/simulate.R")
source("~/Github/causal-eco/Functions/calibrate.R")
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
    summarize(linear = mean(linear - theta)/theta,
              lmtp = mean(lmtp - theta)/theta)
  
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
