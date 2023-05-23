### Dependencies

library(dplyr)
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(KernSmooth)
library(ggplot2)

source('~/Github/erc-strata/Functions/gam_models.R')
source('~/Github/erc-strata/Functions/erf_models.R')
source('~/Github/erc-strata/Functions/calibrate.R')

### Simulation Function

fit_sim <- function(n, m, sig_gps = 2, gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a", "b")) {
  
  a.vals <- seq(4, 12, length.out = 81)
  bw.seq <- seq(0.8, 2, length.out = 16)
  
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  x4 <- rnorm(m)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp(x1/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  u4 <- as.numeric(scale((x2 + x4 + 20)^2))
  
  zip_data <- data.frame(zip = 1:m, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
                         u1 = u1, u2 = u2, u3 = u3, u4 = u4)
  
  if (gps_scen == "b") {
    mu_gps <- 8 - 0.25*u1 + 0.75*u2 - 0.75*u3 + 0.25*u4
  } else {
    mu_gps <- 8 - 0.25*x1 + 0.75*x2 - 0.75*x3 + 0.25*x4
  }
  
  zip_data$a <- a <- rnorm(m, mu_gps, sig_gps)
  
  w1 <- rbinom(n, 1, 0.3)
  w2 <- rbinom(n, 1, 0.7)
  
  if (ss_scen == "b"){
    
    mu_ss <- plogis(-0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4)
    
  } else {
    
    mu_ss <- plogis(-0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4)
    
  }
  
  prob <- mu_ss/sum(mu_ss)
  zip <- sample(1:m, n, replace = T, prob = prob)
  ind_data <- data.frame(zip, w1 = w1, w2 = w2)
  
  data <- merge(ind_data, zip_data, by = "zip")
  
  if (out_scen == "b") {
    mu_out <- with(data, plogis(-3 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 
                                  0.25*(a - 8)*u1 + 0.25*(a - 8)*w1 - 
                                  0.5*w1 + 0.5*w2))
    lambda <- sapply(a.vals, function(a.new) mean(plogis(-3 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                                                           0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                                                           0.25*(a.new - 8)*u1 + 0.25*(a.new - 8)*w1 -
                                                           0.5*w1 + 0.5*w2)))
  } else { # y_scen == "a"
    mu_out <- with(data, plogis(-3 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 
                                  0.25*(a - 8)*x1 + 0.25*(a - 8)*w1 - 
                                  0.5*w1 + 0.5*w2))
    lambda <- sapply(a.vals, function(a.new) mean(plogis(-3 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                                                           0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) - 
                                                           0.25*(a.new - 8)*x1 + 0.25*(a.new - 8)*w1 -
                                                           0.5*w1 + 0.5*w2)))
  }
  
  data$y <- rbinom(n, 1, mu_out)
  
  strata_data <- data %>% group_by(zip, w1, w2) %>% 
    summarise(x1 = mean(x1), x2 = mean(x2), 
              x3 = mean(x3), x4 = mean(x4),
              a = mean(a), y = sum(y), n = n())
  
  strata_data$ybar <- with(strata_data, y/n)
  
  data$y <- rbinom(n, 1, plogis(mu_out))
  
  ## LM GPS
  
  pimod <- lm(a ~ x1 + x2 + x3 + x4, data = zip_data)
  pimod.vals <- c(pimod$fitted.values)
  pimod.sd <- sigma(pimod)
  
  # nonparametric density
  pihat <- dnorm(zip_data$a, pimod.vals, pimod.sd) 
  # dens <- density(a.std)
  # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # ipw numerator
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    dnorm(a.tmp, pimod.vals, pimod.sd)
    # approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
  })
  
  phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  phat <- predict(smooth.spline(a.vals, phat.vals), x = zip_data$a)$y
  phat[phat < 0] <- .Machine$double.eps
  
  zip_data$ipw <- phat/pihat # LM GPS
  
  x.mat <- model.matrix(~ x1 + x2 + x3 + x4, data = data.frame(zip_data))
  astar <- c(zip_data$a - mean(zip_data$a))/var(zip_data$a)
  astar2 <- c((zip_data$a - mean(zip_data$a))^2/var(zip_data$a) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(nrow(x.mat), rep(0, ncol(x.mat) + 1)))
  
  zip_data$cal <- mod$weights # CALIBRATION
  
  strata_data <- merge(data.frame(zip = zip_data$zip, ipw = zip_data$ipw, cal = zip_data$cal),
                       strata_data, by = "zip")
  
  # fit gam outcome model
  w <- model.frame(~ x1 + x2 + x3 + x4 + w1 + w2, data = strata_data)[,-1]
  model_data <- gam_models(y = strata_data$y, a = strata_data$a, w = w,
                           log.pop = log(strata_data$n), id = strata_data$zip, 
                           weights = strata_data$ipw, a.vals = a.vals)
  
  # Separate Data into List
  mat.list <- with(model_data, split(cbind(exp(log.pop), resid, muhat.mat), id))
  wts <- do.call(c, lapply(split(exp(model_data$log.pop), model_data$id), sum))
  
  # Aggregate by ZIP-code-year
  mat <- do.call(rbind, lapply(mat.list, function(vec) {
    mat <- matrix(vec, ncol = length(a.vals) + 2)
    colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
  } ))
  
  mat.pool <- data.frame(id = as.numeric(names(mat.list)), mat)
  mhat.vals <- apply(mat.pool[,-(1:2)], 2, mean, na.rm = TRUE)
  resid.dat <- inner_join(mat.pool[,(1:2)], data.frame(a = zip_data$a, id = zip_data$zip), by = "id")
  resid.dat$mhat <- predict(smooth.spline(a.vals, mhat.vals), x = resid.dat$a)$y
  
  # Pseudo-Outcomes
  resid.dat$psi <- with(resid.dat, X1 + mhat)
  
  # grid search bandwidth
  risk.est <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                     psi = resid.dat$psi, a = resid.dat$a)
  bw <- c(bw.seq[which.min(risk.est)])
  
  mhat.mat <- matrix(rep(mhat.vals, nrow(mat.pool)), byrow = TRUE, nrow = nrow(mat.pool))
  phat.mat <- matrix(rep(phat.vals, nrow(mat.pool)), byrow = TRUE, nrow = nrow(mat.pool))
  int.mat <- (mat.pool[,-(1:2)] - mhat.mat)*phat.mat
  
  rm(mat, mat.list, mat.pool); gc()
  
  none <- sapply(a.vals, kern_est_simple, psi = resid.dat$psi, a = resid.dat$a, 
                   bw = bw[1], a.vals = a.vals, se.fit = TRUE, int.mat = int.mat)
  
  simple <- sapply(a.vals, kern_est_simple, psi = resid.dat$psi, a = resid.dat$a, weights = wts,
                       bw = bw[1], a.vals = a.vals, se.fit = TRUE, int.mat = int.mat)
  
  complex <- sapply(a.vals, kern_est_complex, psi = resid.dat$psi, a = resid.dat$a, weights = wts,
                       bw = bw[1], a.vals = a.vals, se.fit = TRUE, int.mat = int.mat)
  
  return(list(est.none = none[1,], se.none = sqrt(none[2,]),
              est.simple = simple[1,], se.simple = sqrt(simple[2,]),
              est.complex = complex[1,], se.complex = sqrt(complex[2,]), 
              lower.none = none[1,] - 1.96*sqrt(none[2,]), upper.none = none[1,] + 1.96*sqrt(none[2,]),
              lower.simple = simple[1,] - 1.96*sqrt(simple[2,]), upper.simple = simple[1,] + 1.96*sqrt(simple[2,]),
              lower.complex = complex[1,] - 1.96*sqrt(complex[2,]), upper.complex = complex[1,] + 1.96*sqrt(complex[2,]),
              lambda = lambda, a.vals = a.vals))
  
}

cl <- makePSOCKcluster(25)
clusterExport(cl, "fit_sim")

parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE)
  parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
            substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)

clusterEvalQ(cl, {
  
  ## set up each worker.  Could also use clusterExport()
  library(dplyr)
  library(parallel)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(splines)
  library(gam)
  library(KernSmooth)
  
  source('~/Github/erc-strata/Functions/gam_models.R')
  source('~/Github/erc-strata/Functions/erf_models.R')
  source('~/Github/erc-strata/Functions/calibrate.R')
  
})

### Run Simulation

scenarios = expand.grid(n = c(100000), m = c(10000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a"))
a.vals <- seq(4, 12, length.out = 81)
n.iter <- 500
dat <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scen <- scenarios[i,]
  m <- as.numeric(scen$m)
  n <- as.numeric(scen$n)
  gps_scen <- as.character(scen$gps_scen)
  out_scen <- as.character(scen$out_scen)
  ss_scen <- as.character(scen$ss_scen)
  
  clusterExport(cl, c("n","m","gps_scen","out_scen","ss_scen"))
  
  out <- parReplicate(cl, n.iter, fit_sim(n = n, m = m, gps_scen = gps_scen, 
                                          ss_scen = ss_scen, out_scen = out_scen))
  
  lambda <- rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda))
  
  # estimate
  est.none <- rowMeans(sapply(1:n.iter, function(i) out[,i]$est.none))
  est.simple <- rowMeans(sapply(1:n.iter, function(i) out[,i]$est.simple))
  est.complex <- rowMeans(sapply(1:n.iter, function(i) out[,i]$est.complex))
  
  # lower bound
  lower.none <- rowMeans(sapply(1:n.iter, function(i) out[,i]$lower.none))
  lower.simple <- rowMeans(sapply(1:n.iter, function(i) out[,i]$lower.simple))
  lower.complex <- rowMeans(sapply(1:n.iter, function(i) out[,i]$lower.complex))
  
  # upper bound
  upper.none <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.none))
  upper.simple <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.simple))
  upper.complex <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.complex))
  
  # absolute bias
  bias.none <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.none)))
  bias.simple <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.simple)))
  bias.complex <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.none)))
  
  # root mean squared error
  rmse.none <- sqrt(rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.none)^2))
  rmse.simple <- sqrt(rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.simple)^2))
  rmse.complex <- sqrt(rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) - sapply(1:n.iter, function(i) out[,i]$est.none)^2))

  # confidence length
  cl.none <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.none) - sapply(1:n.iter, function(i) out[,i]$lower.none))
  cl.simple <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.simple) - sapply(1:n.iter, function(i) out[,i]$lower.simple))
  cl.complex <- rowMeans(sapply(1:n.iter, function(i) out[,i]$upper.complex) - sapply(1:n.iter, function(i) out[,i]$lower.complex))
  
  # coverage probability
  cp.none <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) < sapply(1:n.iter, function(i) out[,i]$upper.none) &
                    rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) > sapply(1:n.iter, function(i) out[,i]$lower.none))
  
  cp.simple <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) < sapply(1:n.iter, function(i) out[,i]$upper.simple) &
                      rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) > sapply(1:n.iter, function(i) out[,i]$lower.simple))
  
  cp.complex <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) < sapply(1:n.iter, function(i) out[,i]$upper.complex) &
                       rowMeans(sapply(1:n.iter, function(i) out[,i]$lambda)) > sapply(1:n.iter, function(i) out[,i]$lower.complex))
  
  dat <- rbind(dat, data.frame(a.vals = rep(a.vals, times = 4),
                               est = c(lambda, est.none, est.simple, est.complex), 
                               lower = c(rep(NA, length(a.vals)), lower.none, lower.simple, lower.complex),
                               upper = c(rep(NA, length(a.vals)), upper.none, upper.simple, upper.complex),
                               bias = c(rep(NA, length(a.vals)), bias.none, bias.simple, bias.complex),
                               rmse = c(rep(NA, length(a.vals)), rmse.none, rmse.simple, rmse.complex),
                               cp = c(rep(NA, length(a.vals)), cp.none, cp.simple, cp.complex),
                               cl = c(rep(NA, length(a.vals)), cl.none, cl.simple, cl.complex),
                               adjust = rep(c("true", "none", "simple", "complex"), each = length(a.vals)),
                               gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))
  
}

save(dat, file = "~/Github/erc-strata/Output/simulation_results.RData")

### Plots

dat$label <- ifelse(dat$adjust == "true", "True ERF",
                    ifelse(dat$adjust == "none", "No Weighting", 
                           ifelse(dat$adjust == "simple", "Weighted Likelihood", "Scaled Kernel Weight")))
dat$scenario <- ifelse(dat$gps_scen == "a" & dat$out_scen == "a", "Correct Specification",
                       ifelse(dat$gps_scen == "a" & dat$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(dat$gps_scen == "b" & dat$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

dat$label <- factor(dat$label, levels = c("True ERF", "No Weighting", "Weighted Likelihood", "Scaled Kernel Weight"))
dat$scenario <- factor(dat$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

dat_tmp <- subset(dat, adjust %in% c("true", "none", "simple") & scenario == "Correct Specification")

first <- dat_tmp %>%
  ggplot(aes(x = a.vals, y = est, color = factor(label))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Absolute Risk of of Mortality",
       color = "Weighting Approach") +
  theme_bw() +
  coord_cartesian(xlim = c(4,12), 
                  ylim = c(0, 0.2)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#000000", "#D81B60", "#F57328")) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/erc-strata/Output/initial.pdf", width = 8, height = 8)
first
dev.off()

dat_tmp <- subset(dat, !(gps_scen == "b" & out_scen == "b"))

second <- dat_tmp %>%
  ggplot(aes(x = a.vals, y = est, color = factor(label))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ scenario) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Absolute Risk of Mortality",
       color = "Weighting Approach") +
  theme_bw() +
  coord_cartesian(xlim = c(4,12), 
                  ylim = c(0, 0.2)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#000000", "#D81B60", "#F57328", "#004D40")) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/erc-strata/Output/dr_property.pdf", width = 16, height = 8)
second
dev.off()

output <- dat %>% group_by(adjust, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse), cp = mean(cp), cl = mean(cl))

save(output, file = "~/Github/erc-strata/Output/simulation_summary.RData")

stopCluster(cl)

