## set up each worker.  Could also use clusterExport()
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(KernSmooth)
library(splines)
library(mgcv)

source('~/Github/erc-strata/Functions/kwls.R')
source('~/Github/erc-strata/Functions/gam_test.R')
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
    mu_ss <- plogis(-0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4)
  } else {
    mu_ss <- plogis(-0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4)
  }
  
  prob <- mu_ss/sum(mu_ss)
  zip <- sample(1:m, n, replace = TRUE, prob = prob)
  ind_data <- merge(data.frame(zip = zip, w1 = w1, w2 = w2), zip_data, by = "zip")
  
  if (out_scen == "b") {
    mu_out <- with(ind_data, plogis(-2 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                  0.5*w1 + 0.5*w2 -  0.25*(a - 8)*u1 + 0.25*(a - 8)*w1))
    lambda <- with(ind_data, sapply(a.vals, function(a.new, ...) 
      mean(plogis(-2 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                    0.5*w1 + 0.5*w2 - 0.25*(a.new - 8)*u1 + 0.25*(a.new - 8)*w1))))
  } else { # y_scen == "a"
    mu_out <- with(ind_data, plogis(-2 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                  0.5*w1 + 0.5*w2 - 0.25*(a - 8)*x1 + 0.25*(a - 8)*w1))
    lambda <- with(ind_data, sapply(a.vals, function(a.new, ...)
      mean(plogis(-2 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                    0.5*w1 + 0.5*w2 - 0.25*(a.new - 8)*x1 + 0.25*(a.new - 8)*w1))))
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
  tm0 <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  ipwmod0 <- calibrate(cmat = cmat, target = tm0)
  data$cal0 <- ipwmod0$weights
  
  tm1 <- c(rep(0, ncol(x.mat) + 1), c(t(x.mat) %*% data$n))
  ipwmod1 <- calibrate(cmat = cmat, target = tm1, base_weights = data$n)
  data$cal1 <- ipwmod1$weights/data$n
  
  ## GAM Model
  data$ybar <- data$y/data$n
  inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  fmla <- as.formula(paste0("ybar ~ s(aa) + ", inner, " + aa:(", inner, ")"))
  mumod <- gam(fmla, data = data.frame(aa = data$a, data),
               weights = data$n, family = quasipoisson())
  w.mat <- predict(mumod, type = "lpmatrix")
  
  ## Spline Outcome Model
  # data$ybar <- data$y/data$n
  # inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  # nsa <- ns(data$a, df = 6)
  # w.mat <- cbind(nsa, model.matrix(formula(paste0("~ ", inner, " + a:(", inner, ")")),
  #                                  data = data.frame(aa = data$a, data)))
  # mumod <- glm(ybar ~ 0 + ., data = data.frame(ybar = data$ybar, w.mat),
  #              weights = data$n, family = quasipoisson())
  
  # fit GAM DR regression
  dr0 <- gam_est0(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, 
                    se.fit = TRUE, a.vals = a.vals, x = x.mat, w = w.mat,
                    ipw = data$cal0, muhat = mumod$fitted.values, 
                    astar = astar, astar2 = astar2, cmat = cmat)
  
  dr1 <- gam_est1(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, 
                    se.fit = TRUE, a.vals = a.vals, x = x.mat, w = w.mat,
                    ipw = data$cal1, muhat = mumod$fitted.values, 
                    astar = astar, astar2 = astar2, cmat = cmat)
  
  # linear algebra for variance
  vals0 <- sapply(a.vals, function(a.tmp, ...) {
    
    w.tmp <- predict(mumod, type = "lpmatrix", newdata = data.frame(aa = a.tmp, data))
    
    # nsa.tmp <- predict(nsa, newx = rep(a.tmp, nrow(data)))
    # w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ ", inner, " + aa:(", inner, ")")),
    #                                      data = data.frame(aa = a.tmp, data)))
    
    l <- ncol(w.tmp)
    o <- ncol(dr0$g.vals)
    idx <- which(a.vals == a.tmp)
    g.val <- c(dr0$g.vals[idx,])
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    one <- rep(1, times = nrow(w.mat))
    
    delta <- c(mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    first <- c(t(data$n*delta) %*% w.tmp %*% dr0$Sig[1:l,1:l] %*% t(w.tmp) %*% (data$n*delta))/(sum(data$n)^2) + 
                2*c(t(data$n*delta) %*% w.tmp %*% dr0$Sig[1:l, (l + 1):(l + o)] %*% g.val)/sum(data$n)
    sig2 <- first + c(t(g.val) %*% dr0$Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    
    mu <- weighted.mean(mhat, w = data$n) + dr0$mu[idx]
    
    return(c(mu = mu, sig2 = sig2))
    
  })
  
  vals1 <- sapply(a.vals, function(a.tmp, ...) {
    
    w.tmp <- predict(mumod, type = "lpmatrix", newdata = data.frame(aa = a.tmp, data),
                     newdata.guaranteed = TRUE, block.size = nrow(data))
    
    # nsa.tmp <- predict(nsa, newx = rep(a.tmp, nrow(data)))
    # w.tmp <- cbind(nsa.tmp, model.matrix(formula(paste0("~ ", inner, " + aa:(", inner, ")")),
    #                                      data = data.frame(aa = a.tmp, data)))
    
    l <- ncol(w.tmp)
    o <- ncol(dr1$g.vals)
    idx <- which(a.vals == a.tmp)
    g.val <- c(dr1$g.vals[idx,])
    mhat <- mumod$family$linkinv(c(w.tmp%*%mumod$coefficients))
    one <- rep(1, times = nrow(w.mat))
    
    delta <- c(mumod$family$mu.eta(mumod$family$linkfun(mhat)))
    first <- c(t(data$n*delta) %*% w.tmp %*% dr1$Sig[1:l,1:l] %*% t(w.tmp) %*% (data$n*delta))/(sum(data$n)^2) + 
      2*c(t(data$n*delta) %*% w.tmp %*% dr1$Sig[1:l, (l + 1):(l + o)] %*% g.val)/sum(data$n)
    sig2 <- first + c(t(g.val) %*% dr1$Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    
    mu <- weighted.mean(mhat, w = data$n) + dr1$mu[idx]
    
    return(c(mu = mu, sig2 = sig2))
    
  })
  
  return(list(est.dr0 = vals0[1,], se.dr0 = sqrt(vals0[2,]),
              est.dr1 = vals1[1,], se.dr0 = sqrt(vals1[2,]),
              lower.dr0 = vals0[1,] - 1.96*sqrt(vals0[2,]), upper.dr0 = vals0[1,] + 1.96*sqrt(vals0[2,]),
              lower.dr1 = vals1[1,] - 1.96*sqrt(vals1[2,]), upper.dr1 = vals1[1,] + 1.96*sqrt(vals1[2,]),
              lambda = lambda, a.vals = a.vals))
  
}

### Run Simulation

scenarios = expand.grid(n = c(100000), m = c(5000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a"))
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
  
  out <- mclapply(1:n.iter, fit_sim, n = n, m = m, gps_scen = gps_scen, 
                   ss_scen = ss_scen, out_scen = out_scen,
                   a.vals = a.vals, bw.seq = bw.seq, mc.cores = 16)
  
  lambda <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda))
  
  # estimate
  est.dr0 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr0))
  est.dr1 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr1))
  
  # lower bound
  lower.dr0 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr0))
  lower.dr1 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr1))
  
  # upper bound
  upper.dr0 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr0))
  upper.dr1 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr1))
  
  # absolute bias
  bias.dr0 <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr0)), na.rm = TRUE)
  bias.dr1 <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr1)), na.rm = TRUE)
  
  # root mean squared error
  rmse.dr0 <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr0))^2, na.rm = TRUE))
  rmse.dr1 <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr1))^2, na.rm = TRUE))
  
  # confidence length
  cl.dr0 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr0) - sapply(1:n.iter, function(i) out[[i]]$lower.dr0), na.rm = TRUE)
  cl.dr1 <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr1) - sapply(1:n.iter, function(i) out[[i]]$lower.dr1), na.rm = TRUE)
  
  # coverage probability
  cp.dr0 <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr0) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr0), na.rm = TRUE)
  cp.dr1 <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr1) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr1), na.rm = TRUE)
  
  df <- rbind(df, data.frame(a.vals = rep(a.vals, times = 3),
                               est = c(lambda, est.dr0, est.dr1), 
                               lower = c(rep(NA, length(a.vals)), lower.dr0, lower.dr1),
                               upper = c(rep(NA, length(a.vals)), upper.dr0, upper.dr1),
                               bias = c(rep(NA, length(a.vals)), bias.dr0, bias.dr1),
                               rmse = c(rep(NA, length(a.vals)), rmse.dr0, rmse.dr1),
                               cp = c(rep(NA, length(a.vals)), cp.dr0, cp.dr1),
                               cl = c(rep(NA, length(a.vals)), cl.dr0, cl.dr1),
                               adjust = rep(c("true", "dr0", "dr1"), each = length(a.vals)),
                               gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))

}

save(df, file = "~/Github/erc-strata/Output/simulation_results.RData")

### Plots

df$label <- ifelse(df$adjust == "true", "True ERF",
                    ifelse(df$adjust == "dr0", "Unscaled IPW", "Scaled IPW"))
df$scenario <- ifelse(df$gps_scen == "a" & df$out_scen == "a", "Correct Specification",
                       ifelse(df$gps_scen == "a" & df$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(df$gps_scen == "b" & df$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

df$label <- factor(df$label, levels = c("True ERF", "Unscaled IPW", "Scaled IPW"))
df$scenario <- factor(df$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

df_tmp <- subset(df, !(gps_scen == "b" & out_scen == "b" | m == 10000))

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
  scale_color_manual(values = c("#F57328", "#004D40", "#ffdb58")) +
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/erc-strata/Output/simulation_plot.pdf", width = 16, height = 8)
plot
dev.off()

output <- df %>% group_by(adjust, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse),  cp = mean(cp), cl = mean(cl))

save(output, file = "~/Github/erc-strata/Output/simulation_summary.RData")

