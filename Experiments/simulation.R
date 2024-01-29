library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(mgcv)
library(sandwich)
library(ggplot2)

source('~/Github/causal-eco/Functions/gam_dr.R')
source('~/Github/causal-eco/Functions/gam_ipw.R')
source('~/Github/causal-eco/Functions/calibrate.R')

### Simulation Function

fit_sim <- function(i, n, m, sig_gps = 2, gps_scen = c("a", "b"), out_scen = c("a", "b"), 
                    ss_scen = c("a", "b"), a.vals = seq(4, 12, length.out = 81)) {
  
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  x4 <- rnorm(m)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp((x1 + x4)/4)))
  u2 <- as.numeric(scale(x2/(1 + exp(x3))))
  u3 <- as.numeric(scale((x1 + x3)^2))
  u4 <- as.numeric(scale((x2 + x4)^2))
  
  zip_data <- data.frame(zip = 1:m, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
                         u1 = u1, u2 = u2, u3 = u3, u4 = u4)
  
  if (gps_scen == "b") {
    mu_gps <- 8 - 0.25*u1 + 0.75*u2 - 0.75*u3 + 0.25*u4
  } else {
    mu_gps <- 8 - 0.25*x1 + 0.75*x2 - 0.75*x3 + 0.25*x4
  }
  
  # combine zip data
  zip_data$a <- a <- rnorm(m, mu_gps, sig_gps)
  
  w1 <- rbinom(n, 1, 0.3)
  w2 <- rbinom(n, 1, 0.7)
  
  if (ss_scen == "b"){
    mu_ss <- exp(-0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4)
  } else {
    mu_ss <- exp(-0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4)
  }
  
  # simulate offset with mu_ss
  prob <- mu_ss/sum(mu_ss)
  zip <- sample(1:m, n, replace = TRUE, prob = prob)
  ind_data <- merge(data.frame(zip = zip, w1 = w1, w2 = w2), zip_data, by = "zip")
  
  # simulate outcomes
  if (out_scen == "b") {
    mu_out <- with(ind_data, plogis(-1.5 + 0.5*u1 - 0.5*u2 - 0.5*u3 +
                                      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                      0.25*(a - 8)*u3 + 0.5*(a - 8)*u4))
    theta <- with(ind_data, sapply(a.vals, function(a.new, ...) 
      mean(plogis(-1.5 + 0.5*u1 - 0.5*u2 - 0.5*u3 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) - 
                    0.25*(a.new - 8)*u3 + 0.5*(a.new - 8)*u4))))
  } else { # y_scen == "a"
    mu_out <- with(ind_data, plogis(-1.5 + 0.5*x1 - 0.5*x2 - 0.5*x3 +
                                      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) -
                                      0.25*(a - 8)*x3 + 0.5*(a - 8)*x4))
    theta <- with(ind_data, sapply(a.vals, function(a.new, ...)
      mean(plogis(-1.5 + 0.5*x1 - 0.5*x2 - 0.5*x3 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                    0.25*(a.new - 8)*x3 + 0.5*(a.new - 8)*x4))))
  }
  
  ind_data$y <- rbinom(n, 1, mu_out)
  
  # combine strata data
  strata_data <- ind_data %>% group_by(zip, w1, w2) %>% 
    summarise(x1 = mean(x1), x2 = mean(x2), 
              x3 = mean(x3), x4 = mean(x4),
              a = mean(a), y = sum(y), n = n())
  
  zip_data <- subset(zip_data, zip %in% unique(strata_data$zip))
  data <- merge(zip_data, strata_data %>% group_by(zip) %>% summarise(y = sum(y), n = sum(n)), by = "zip")
  
  # Ecological Regression Results
  if (out_scen == "b") {
    theta.eco <- with(data, sapply(a.vals, function(a.new, ...) 
      mean(plogis(-1.5 + 0.5*u1 - 0.5*u2 - 0.5*u3 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) - 
                    0.25*(a.new - 8)*u3 + 0.5*(a.new - 8)*u4))))
  } else { # y_scen == "a"
    theta.eco <- with(data, sapply(a.vals, function(a.new, ...)
      mean(plogis(-1.5 + 0.5*x1 - 0.5*x2 - 0.5*x3 +
                    0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                    0.25*(a.new - 8)*x3 + 0.5*(a.new - 8)*x4))))
  }
  
  ## Calibration weights
  x.mat <- model.matrix(~ x1 + x2 + x3 + x4, data = data.frame(data))
  astar <- c(data$a - mean(data$a))/sd(data$a)
  astar2 <- c((data$a - mean(data$a))^2/var(data$a) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)

  # fit individual calibration model
  tm.ind <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat*data$n))
  ipwmod.ind <- calibrate(cmat = cmat, target = tm.ind, base_weights = data$n)
  data$cal.ind <- ipwmod.ind$weights/ipwmod.ind$base_weights
  
  # fit ecological calibration model
  tm.eco <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  ipwmod.eco <- calibrate(cmat = cmat, target = tm.eco)
  data$cal.eco <- ipwmod.eco$weights/ipwmod.eco$base_weights
  
  ## GAM Outcome Model
  data$ybar <- data$y/data$n
  inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  fmla <- formula(paste0("ybar ~ s(aa, bs = 'cr', k = 6) + ", inner))
  mumod <- mgcv::gam(fmla, weights = data$n, family = quasipoisson(), 
                data = data.frame(aa = data$a, data))
  muhat <- predict(mumod, type = "response")
  w.mat <- predict(mumod, type = "lpmatrix")
  
  # fit GAM DR regression
  dr.ind <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, eco = FALSE,
                   ipw = data$cal.ind, muhat = mumod$fitted.values, a.vals = a.vals, se.fit = TRUE, 
                   x = x.mat, w = w.mat, astar = astar, astar2 = astar2)
  
  dr.eco <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, eco = TRUE,
                   ipw = data$cal.eco, muhat = mumod$fitted.values, a.vals = a.vals, se.fit = TRUE, 
                   x = x.mat, w = w.mat, astar = astar, astar2 = astar2)
  
  dr.ind.eco <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, eco = TRUE,
                       ipw = data$cal.eco, muhat = mumod$fitted.values, a.vals = a.vals, se.fit = TRUE, 
                       x = x.mat, w = w.mat, astar = astar, astar2 = astar2)
  
  dr.vals <- sapply(a.vals, function(a.tmp, ...) {
    
    ## preliminaries
    
    # GAM Predictions
    w.tmp <- predict(mumod, newdata = data.frame(aa = a.tmp, data), type = "lpmatrix")  
    muhat.tmp <- predict(mumod, newdata = data.frame(aa = a.tmp, data), type = "response")
    delta.eco <- mumod$family$mu.eta(mumod$family$linkfun(muhat.tmp))
    delta.ind <- c(data$n*delta.eco)
    
    # mean values
    mhat.val.ind <- weighted.mean(muhat.tmp, w = data$n)
    mhat.val.eco <- mean(muhat.tmp)

    idx <- which(a.vals == a.tmp)
    
    # Robust Variance
    o <- ncol(dr.ind$g.vals)
    l <- ncol(w.tmp)
    
    Sig <- as.matrix(dr.ind$Sig)
    g.val <- c(dr.ind$g.vals[idx,])
    
    first <- c(t(delta.ind) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta.ind)/(sum(data$n)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.ind <- first + second
    
    # Ecological Regression Variance
    o <- ncol(dr.eco$g.vals)
    l <- ncol(w.tmp)
    
    Sig <- as.matrix(dr.eco$Sig)
    g.val <- c(dr.eco$g.vals[idx,])
    
    first <- c(t(delta.eco) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta.eco)/(nrow(data)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.eco <- first + second
    
    # Eco IPW; Ind Outcome
    o <- ncol(dr.ind.eco$g.vals)
    l <- ncol(w.tmp)
    
    Sig <- as.matrix(dr.ind$Sig)
    g.val <- c(dr.ind$g.vals[idx,])
    
    first <- c(t(delta.ind) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta.ind)/(sum(data$n)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.ind.eco <- first + second
    
    # prdictions
    mu.ind <- mhat.val.ind + dr.ind$mu.vals[idx]
    mu.eco <- mhat.val.eco + dr.eco$mu.vals[idx]
    mu.ind.eco <- mhat.val.ind + dr.ind.eco$mu.vals[idx]
    
    return(c(mu.ind = mu.ind, sig2.ind = sig2.ind,
             mu.eco = mu.eco, sig2.eco = sig2.eco,
             mu.ind.eco = mu.ind.eco, sig2.ind.eco = sig2.ind.eco))
    
  })
  
  return(list(est.ind = dr.vals[1,], est.eco = dr.vals[3,], est.ind.eco = dr.vals[5,],
              se.ind = sqrt(dr.vals[2,]), se.eco = sqrt(dr.vals[4,]), se.ind.eco = sqrt(dr.vals[6,]),
              lower.ind = dr.vals[1,] - 1.96*sqrt(dr.vals[2,]), upper.ind = dr.vals[1,] + 1.96*sqrt(dr.vals[2,]),
              lower.eco = dr.vals[3,] - 1.96*sqrt(dr.vals[4,]), upper.eco = dr.vals[3,] + 1.96*sqrt(dr.vals[4,]),
              lower.ind.eco = dr.vals[5,] - 1.96*sqrt(dr.vals[6,]), upper.ind.eco = dr.vals[5,] + 1.96*sqrt(dr.vals[6,]),
              theta = theta, theta.eco = theta.eco, a.vals = a.vals))
  
}

### Run Simulation

scenarios = expand.grid(n = c(20000), m = c(1000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a"))
a.vals <- seq(4, 12, length.out = 81)
n.iter <- 1000
df <- perform <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  print(i)
  
  scen <- scenarios[i,]
  m <- as.numeric(scen$m)
  n <- as.numeric(scen$n)
  gps_scen <- as.character(scen$gps_scen)
  out_scen <- as.character(scen$out_scen)
  ss_scen <- as.character(scen$ss_scen)
  
  a.vals <- seq(4, 12, length.out = 81)
  sig_gps <- 2
  
  out <- mclapply(1:n.iter, fit_sim, n = n, m = m, gps_scen = gps_scen, 
                   ss_scen = ss_scen, out_scen = out_scen, sig_gps = sig_gps,
                   a.vals = a.vals, mc.cores = 15)
  
  theta <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta))
  theta.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco))
  
  # estimate
  est.ind <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.ind))
  est.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.eco))
  est.ind.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.ind.eco))
  
  # lower bound
  lower.ind <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.ind))
  lower.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.eco))
  lower.ind.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.ind.eco))
  

  # upper bound
  upper.ind <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ind))
  upper.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.eco))
  upper.ind.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ind.eco))
  
  # absolute bias
  bias.ind <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ind)), na.rm = TRUE)
  bias.eco <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.eco)), na.rm = TRUE)
  bias.ind.eco <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ind.eco)), na.rm = TRUE)
  bias.ind0 <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ind)), na.rm = TRUE)
  bias.eco0 <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.eco)), na.rm = TRUE)
  bias.ind.eco0 <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ind.eco)), na.rm = TRUE)
  
  
  # root mean squared error
  rmse.ind <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ind))^2, na.rm = TRUE))
  rmse.eco <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.eco))^2, na.rm = TRUE))
  rmse.ind.eco <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ind.eco))^2, na.rm = TRUE))
  rmse.ind0 <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ind))^2, na.rm = TRUE))
  rmse.eco0 <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.eco))^2, na.rm = TRUE))
  rmse.ind.eco0 <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ind.eco))^2, na.rm = TRUE))
  
  # coverage probability
  cp.ind <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) < sapply(1:n.iter, function(i) out[[i]]$upper.ind) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) > sapply(1:n.iter, function(i) out[[i]]$lower.ind), na.rm = TRUE)
  cp.eco <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) < sapply(1:n.iter, function(i) out[[i]]$upper.eco) &
                      rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) > sapply(1:n.iter, function(i) out[[i]]$lower.eco), na.rm = TRUE)
  cp.ind.eco <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) < sapply(1:n.iter, function(i) out[[i]]$upper.ind.eco) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) > sapply(1:n.iter, function(i) out[[i]]$lower.ind.eco), na.rm = TRUE)
  cp.ind0 <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) < sapply(1:n.iter, function(i) out[[i]]$upper.ind) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) > sapply(1:n.iter, function(i) out[[i]]$lower.ind), na.rm = TRUE)
  cp.eco0 <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) < sapply(1:n.iter, function(i) out[[i]]$upper.eco) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) > sapply(1:n.iter, function(i) out[[i]]$lower.eco), na.rm = TRUE)
  cp.ind.eco0 <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) < sapply(1:n.iter, function(i) out[[i]]$upper.ind.eco) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) > sapply(1:n.iter, function(i) out[[i]]$lower.ind.eco), na.rm = TRUE)
  
  
  df <- rbind(df, data.frame(a.vals = rep(a.vals, times = 5),
                   est = c(theta, theta.eco, est.ind, est.eco, est.ind.eco), 
                   lower = c(rep(NA, 2*length(a.vals)), lower.ind, lower.eco, lower.ind.eco),
                   upper = c(rep(NA, 2*length(a.vals)), upper.ind, upper.eco, upper.ind.eco),
                   ipw = rep(c(NA, NA, "Individual", "Ecological", "Ecological"), each = length(a.vals)),
                   outcome = rep(c(NA, NA, "Individual", "Ecological", "Individual"), each = length(a.vals)),
                   actual = rep(c("Individual", "Ecological", NA, NA, NA), each = length(a.vals)),
                   gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))
  
  perform <- rbind(perform, data.frame(a.vals = rep(a.vals, times = 3),
                                  bias = c(bias.ind, bias.eco, bias.ind.eco),
                                  rmse = c(rmse.ind, rmse.eco, rmse.ind.eco),
                                  cp = c(cp.ind, cp.eco, cp.ind.eco),
                                  bias0 = c(bias.ind0, bias.eco0, bias.ind.eco0),
                                  rmse0 = c(rmse.ind0, rmse.eco0, rmse.ind.eco0),
                                  cp0 = c(cp.ind0, cp.eco0, cp.ind.eco0),
                                  ipw = rep(c("Individual", "Ecological", "Ecological"), each = length(a.vals)),
                                  outcome = rep(c("Individual", "Ecological", "Individual"), each = length(a.vals)),
                                  gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))

}

save(df, file = "~/Github/causal-eco/Output/simulation_results.RData")
save(perform, file = "~/Github/causal-eco/Output/simulation_performance.RData")

### Plots

df$label <- with(df, ifelse(is.na(actual), paste("GPS: ", ipw, ", OM: ", outcome, sep = ""), NA))
df$truth <- with(df, ifelse(is.na(actual), "Estimated", 
                            ifelse(actual == "True Ecological", "Ecological")))

df$scenario <- ifelse(df$gps_scen == "a" & df$out_scen == "a", "Correct Specification",
                       ifelse(df$gps_scen == "a" & df$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(df$gps_scen == "b" & df$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

df$scenario <- factor(df$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

df_tmp <- subset(df, !(gps_scen == "b" & out_scen == "b") & label != "IPW Estimator" & ss_scen == "a")

plot <- df_tmp %>%
  ggplot(aes(x = a.vals, y = est, linetype = factor(label))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ scenario) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Absolute Risk of Mortality",
       color = "Estimand", linetype = "Method") +
  theme_bw() +
  coord_cartesian(xlim = c(4,12), 
                  ylim = c(0, 0.45)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#F57328", "#004D40", "#ffdb58", "#1E88E5")) +
  scale_y_continuous(breaks = seq(0, 0.45, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/causal-eco/Output/simulation_plot.pdf", width = 16, height = 8)
plot
dev.off()

output <- perform %>% group_by(ipw, outcome, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse), cp = mean(cp), 
                                                                                              bias0 = mean(bias0), rmse0 = mean(rmse0), cp0 = mean(cp0), )

save(output, file = "~/Github/causal-eco/Output/simulation_summary.RData")
