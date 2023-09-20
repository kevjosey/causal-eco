## set up each worker.  Could also use clusterExport()
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(KernSmooth)
library(splines)
library(mgcv)
library(sandwich)

source('~/Github/erc-strata/Functions/gam_dr.R')
source('~/Github/erc-strata/Functions/gam_ipw.R')
source('~/Github/erc-strata/Functions/calibrate.R')

### Simulation Function

fit_sim <- function(i, n, m, sig_gps = 2, gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a", "b"),
                    a.vals = seq(4, 12, length.out = 81), bw.seq = seq(0.1, 2, length.out = 16)) {
  
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
    mu_ss <- exp(-0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4)
  } else {
    mu_ss <- exp(-0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4)
  }
  
  prob <- mu_ss/sum(mu_ss)
  zip <- sample(1:m, n, replace = TRUE, prob = prob)
  ind_data <- merge(data.frame(zip = zip, w1 = w1, w2 = w2), zip_data, by = "zip")
  
  if (out_scen == "b") {
    mu_out <- with(ind_data, plogis(-3 + u1 - u2 - u3 + u4 +
                                      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                      0.25*(a - 8)*u1 - 0.25*(a - 8)*u2))
    lambda <- with(ind_data, sapply(a.vals, function(a.new, ...) 
      mean(plogis(-3 + u1 - u2 - u3 + u4 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) - 
                    0.25*(a.new - 8)*u1 - 0.25*(a.new - 8)*u2))))
  } else { # y_scen == "a"
    mu_out <- with(ind_data, plogis(-3 + x1 - x2 - x3 + x4 +
                                      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                      0.25*(a - 8)*x1 - 0.25*(a - 8)*x2))
    lambda <- with(ind_data, sapply(a.vals, function(a.new, ...)
      mean(plogis(-3 + x1 - x2 - x3 + x4 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                    0.25*(a.new - 8)*x1 - 0.25*(a.new - 8)*x2))))
  }
  
  ind_data$y <- rbinom(n, 1, mu_out)
  
  strata_data <- ind_data %>% group_by(zip, w1, w2) %>% 
    summarise(x1 = mean(x1), x2 = mean(x2), 
              x3 = mean(x3), x4 = mean(x4),
              a = mean(a), y = sum(y), n = n())
  
  zip_data <- subset(zip_data, zip %in% unique(strata_data$zip))
  data <- merge(zip_data, strata_data %>% group_by(zip) %>% summarise(y = sum(y), n = sum(n)), by = "zip")
  
  ## Calibration weights
  x.mat <- model.matrix(~ x1 + x2 + x3 + x4, data = data.frame(data))
  astar <- c(data$a - mean(data$a))/var(data$a)
  astar2 <- c((data$a - mean(data$a))^2/var(data$a) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)

  # fit calibration model
  tm <- c(rep(0, ncol(x.mat) + 1), c(t(x.mat) %*% data$n))
  ipwmod <- calibrate(cmat = cmat, target = tm, base_weights = data$n)
  data$cal <- ipwmod$weights/ipwmod$base_weights
  
  ## GAM Outcome Model
  # data$ybar <- data$y/data$n
  # inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  # fmla <- as.formula(paste0("ybar ~ s(aa) + ", inner, " + aa:(", inner, ")"))
  # mumod <- bam(fmla, data = data.frame(aa = data$a, data),
  #              weights = data$n, family = quasipoisson())
  # w.mat <- predict(mumod, type = "lpmatrix")
  # Omega <- mumod$Vp
  
  ## Spline Outcome Model
  data$ybar <- data$y/data$n
  inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  nsa <- ns(data$a, df = 6)
  w.mat <- cbind(nsa, model.matrix(formula(paste0("~ ", inner, " + aa:(", inner, ")")),
                                   data = data.frame(aa = data$a, data)))
  mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = data$ybar, w.mat),
               weights = data$n, family = gaussian())
  Omega <- vcovHC(mumod)
  
  # fit GAM DR regression
  ipw <- gam_ipw(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, 
                 ipw = data$cal, a.vals = a.vals, se.fit = TRUE, 
                 x = x.mat, astar = astar, astar2 = astar2, cmat = cmat)
  
  dr <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, 
                ipw = data$cal, muhat = mumod$fitted.values, a.vals = a.vals,
                se.fit = TRUE,  x = x.mat, w = w.mat,
                astar = astar, astar2 = astar2, cmat = cmat)
  
  vals <- sapply(a.vals, function(a.tmp, ...) {
    
    ## preliminaries
    
    # GAM
    # w.tmp <- predict(mumod, type = "lpmatrix", newdata = data.frame(aa = a.tmp, data),
    #                  newdata.guaranteed = TRUE, block.size = nrow(data))
    
    # Splines
    nsa.tmp <- predict(nsa, newx = rep(a.tmp, nrow(data)))
    w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ ", inner, " + aa:(", inner, ")")),
                                         data = data.frame(aa = a.tmp, data)))
    
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    delta <- c(data$n*mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    Sig <- as.matrix(dr$Sig)
    
    # Outcome Model
    om.mu <- weighted.mean(mhat, w = data$n)
    om.sig2 <- c(t(delta) %*% w.tmp %*% Omega %*% t(w.tmp) %*% (delta))/(sum(data$n)^2 - sum(data$n^2)) 
    
    # Doubly-Robust Variance
    l <- ncol(w.tmp)
    o <- ncol(dr$g.vals)
    idx <- which(a.vals == a.tmp)
    g.val <- c(dr$g.vals[idx,])
    
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% (delta))/(sum(data$n)^2) + 
      2*c(t(delta) %*% w.tmp %*% Sig[1:l, (l + 1):(l + o)] %*% g.val)/sum(data$n)
    dr.sig2 <- first + c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    
    dr.mu <- om.mu + dr$eta.vals[idx]
    
    return(c(om.mu = om.mu, dr.mu = dr.mu, om.sig2 = om.sig2, dr.sig2 = dr.sig2))
    
  })
  
  return(list(est.ipw = ipw[1,], est.om = vals[1,], est.dr = vals[2,],
              se.ipw = sqrt(ipw[2,]), se.om = sqrt(vals[3,]), se.dr = sqrt(vals[4,]),
              lower.ipw = ipw[1,] - 1.96*sqrt(ipw[2,]), upper.ipw = ipw[1,] + 1.96*sqrt(ipw[2,]),
              lower.om = vals[1,] - 1.96*sqrt(vals[3,]), upper.om = vals[1,] + 1.96*sqrt(vals[3,]),
              lower.dr = vals[2,] - 1.96*sqrt(vals[4,]), upper.dr = vals[2,] + 1.96*sqrt(vals[4,]),
              lambda = lambda, a.vals = a.vals))
  
}

### Run Simulation

scenarios = expand.grid(n = c(10000), m = c(1000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a", "b"))
a.vals <- seq(4, 12, length.out = 81)
n.iter <- 200
df <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scen <- scenarios[i,]
  m <- as.numeric(scen$m)
  n <- as.numeric(scen$n)
  gps_scen <- as.character(scen$gps_scen)
  out_scen <- as.character(scen$out_scen)
  ss_scen <- as.character(scen$ss_scen)
  
  a.vals <- seq(4, 12, length.out = 81)
  bw.seq <- seq(0.1, 2, length.out = 16)
  sig_gps <- 2
  
  out <- mclapply(1:n.iter, fit_sim, n = n, m = m, gps_scen = gps_scen, 
                   ss_scen = ss_scen, out_scen = out_scen, sig_gps = sig_gps,
                   a.vals = a.vals, bw.seq = bw.seq, mc.cores = 25)
  
  lambda <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda))
  
  # estimate
  est.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.ipw))
  est.om <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.om))
  est.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr))
  
  # lower bound
  lower.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.ipw))
  lower.om <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.om))
  lower.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr))
  
  # upper bound
  upper.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw))
  upper.om <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.om))
  upper.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr))
  
  # absolute bias
  bias.ipw <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw)), na.rm = TRUE)
  bias.om <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.om)), na.rm = TRUE)
  bias.dr <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr)), na.rm = TRUE)
  
  # root mean squared error
  rmse.ipw <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw))^2, na.rm = TRUE))
  rmse.om <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.om))^2, na.rm = TRUE))
  rmse.dr <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr))^2, na.rm = TRUE))
  
  # confidence length
  cl.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw) - sapply(1:n.iter, function(i) out[[i]]$lower.ipw), na.rm = TRUE)
  cl.om <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.om) - sapply(1:n.iter, function(i) out[[i]]$lower.om), na.rm = TRUE)
  cl.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr) - sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  
  # coverage probability
  cp.ipw <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.ipw) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.ipw), na.rm = TRUE)
  cp.om <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.om) &
                      rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.om), na.rm = TRUE)
  cp.dr <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  
  df <- rbind(df, data.frame(a.vals = rep(a.vals, times = 4),
                               est = c(lambda, est.ipw, est.om, est.dr), 
                               lower = c(rep(NA, length(a.vals)), lower.ipw, lower.om, lower.dr),
                               upper = c(rep(NA, length(a.vals)), upper.ipw, upper.om, upper.dr),
                               bias = c(rep(NA, length(a.vals)), bias.ipw, bias.om, bias.dr),
                               rmse = c(rep(NA, length(a.vals)), rmse.ipw, rmse.om, rmse.dr),
                               cp = c(rep(NA, length(a.vals)), cp.ipw, cp.om, cp.dr),
                               cl = c(rep(NA, length(a.vals)), cl.ipw, cl.om, cl.dr),
                               adjust = rep(c("true", "ipw", "om", "dr"), each = length(a.vals)),
                               gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))

}

save(df, file = "~/Github/erc-strata/Output/simulation_results.RData")

### Plots

df$label <- ifelse(df$adjust == "true", "True ERF",
                    ifelse(df$adjust == "dr", "DR Estimator",
                           ifelse(df$adjust == "om", "OM Estimator", "IPW Estimator")))
df$scenario <- ifelse(df$gps_scen == "a" & df$out_scen == "a", "Correct Specification",
                       ifelse(df$gps_scen == "a" & df$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(df$gps_scen == "b" & df$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

df$label <- factor(df$label, levels = c("True ERF", "IPW Estimator", "OM Estimator", "DR Estimator"))
df$scenario <- factor(df$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

df_tmp <- subset(df, !(gps_scen == "b" & out_scen == "b") & ss_scen == "a")

plot <- df_tmp %>%
  ggplot(aes(x = a.vals, y = est, color = factor(label))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ scenario) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Absolute Risk of Mortality",
       color = "Method") +
  theme_bw() +
  coord_cartesian(xlim = c(4,12), 
                  ylim = c(0, 0.4)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#F57328", "#004D40", "#ffdb58", "#1E88E5")) +
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/erc-strata/Output/simulation_plot.pdf", width = 16, height = 8)
plot
dev.off()

output <- df %>% group_by(adjust, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse),  cp = mean(cp), cl = mean(cl))

save(output, file = "~/Github/erc-strata/Output/simulation_summary.RData")

