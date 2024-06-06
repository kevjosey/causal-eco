library(Matrix)
library(kernlab)
library(osqp)
library(scam)
library(SuperLearner)
library(lmtp)
library(dplyr)
library(parallel)

rm(list = ls())
source("~/Github/causal-eco/Functions/simulate.R")
source("~/Github/causal-eco/Functions/calibrate.R")

n.iter <- 100 # simulation iterations

# run simulation study
simDat <- replicate(100, gen_data(n = 1000, scenario = "base", sig = sqrt(2)))
simFit <- lapply(1:n.iter, sim_fit, simDat = simDat)

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
  summarize(linear = mean(linear - theta),
            lmtp = mean(lmtp - theta))

rmse <- est_mat %>% group_by(estimate) %>% 
  summarize(linear = sqrt(mean((linear - theta)^2)),
            lmtp = sqrt(mean((lmtp - theta)^2)))

cp <- ci_mat %>% group_by(estimate) %>% 
  summarize(cp_linear = mean(linear_lower < mean(theta) & linear_upper > mean(theta)),
            cp_mlinear = mean(mlinear_lower < mean(theta) & mlinear_upper > mean(theta)),
            cp_lmtp = mean(lmtp_lower < mean(theta) & lmtp_upper > mean(theta)))

ci_length <- ci_diff %>% group_by(estimate) %>% 
  summarize(linear = mean(linear), mlinear = mean(mlinear), lmtp = mean(lmtp))
