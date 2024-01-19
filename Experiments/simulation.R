library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(scam)
library(sandwich)
library(ggplot2)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_ipw.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')

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
  fmla <- formula(paste0("ybar ~ s(aa, bs = 'cr') + ", inner))
  mumod <- scam(fmla, weights = data$n, family = quasipoisson(), 
                data = data.frame(aa = data$a, data))
  muhat <- predict(mumod, type = "response")
  w.mat <- predict(mumod, type = "lpmatrix")
  
  # fit GAM IPW regression
  ipw.ind <- gam_ipw(a = data$a, y = data$ybar, family = gaussian(), weights = data$n, 
                 ipw = data$cal.eco, a.vals = a.vals, se.fit = TRUE, 
                 x = x.mat, astar = astar, astar2 = astar2)
  
  ipw.eco <- gam_ipw(a = data$a, y = data$ybar, family = gaussian(),
                     ipw = data$cal.eco, a.vals = a.vals, se.fit = TRUE, 
                     x = x.mat, astar = astar, astar2 = astar2)
  
  # fit GAM DR regression
  dr.ind <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, weights = data$n, eco = FALSE,
                   ipw = data$cal.ind, muhat = mumod$fitted.values, a.vals = a.vals, se.fit = TRUE, 
                   x = x.mat, w = w.mat, astar = astar, astar2 = astar2)
  
  dr.eco <- gam_dr(a = data$a, y = data$ybar, family = mumod$family, a.vals = a.vals,
                   ipw = data$cal.eco, muhat = mumod$fitted.values, se.fit = TRUE, 
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
    
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(data$n)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.ind <- first + second
    
    # Ecological Regression Variance
    o <- ncol(dr.eco$g.vals)
    l <- ncol(w.tmp)
    
    Sig <- as.matrix(dr.eco$Sig)
    g.val <- c(dr.eco$g.vals[idx,])
    
    first <- c(t(delta0) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta0)/(nrow(data)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.eco <- first + second
    
    # Mixture of the two
    o <- ncol(dr.ind.eco$g.vals)
    l <- ncol(w.tmp)
    
    Sig <- as.matrix(dr.ind$Sig)
    g.val <- c(dr.ind$g.vals[idx,])
    
    first <- c(t(delta) %*% w.tmp %*% Sig[1:l,1:l] %*% t(w.tmp) %*% delta)/(sum(data$n)^2)
    second <- c(t(g.val) %*% Sig[(l + 1):(l + o), (l + 1):(l + o)] %*% g.val)
    sig2.ind.eco <- first + second
    
    # prdictions
    mu.ind <- mhat.val.ind + dr.ind$eta.vals[idx]
    mu.eco <- mhat.val.eco + dr.eco$eta.vals[idx]
    mu.ind.eco <- mhat.val.ind + dr.ind.eco$eta.vals[idx]
    
    return(c(mu.ind = mu.ind, sig2.ind = sig2.ind,
             mu.eco = mu.eco, sig2.eco = sig2.eco,
             mu.eco = mu.eco, sig2.eco = sig2.eco))
    
  })
  
  return(list(est.ipw = ipw[1,], est.dr = dr.vals[1,], est.ipw.eco = ipw.eco[1,], est.dr.eco = dr.vals[3,],
              se.ipw = sqrt(ipw[2,]), se.dr = sqrt(dr.vals[2,]), se.ipw.eco = sqrt(ipw.eco[2,]), se.dr.eco = sqrt(dr.vals[4,]),
              lower.ipw = ipw[1,] - 1.96*sqrt(ipw[2,]), upper.ipw = ipw[1,] + 1.96*sqrt(ipw[2,]),
              lower.dr = dr.vals[1,] - 1.96*sqrt(dr.vals[2,]), upper.dr = dr.vals[1,] + 1.96*sqrt(dr.vals[2,]),
              lower.ipw.eco = ipw.eco[1,] - 1.96*sqrt(ipw.eco[2,]), upper.ipw.eco = ipw.eco[1,] + 1.96*sqrt(ipw.eco[2,]),
              lower.dr.eco = dr.vals[3,] - 1.96*sqrt(dr.vals[4,]), upper.dr.eco = dr.vals[3,] + 1.96*sqrt(dr.vals[4,]),
              theta = theta, theta.eco = theta.eco, a.vals = a.vals))
  
}

### Run Simulation

scenarios = expand.grid(n = c(10000), m = c(1000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a"))
a.vals <- seq(4, 12, length.out = 81)
n.iter <- 200
df <- data.frame()

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
                   a.vals = a.vals, mc.cores = 16)
  
  theta <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta))
  theta.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco))
  
  # estimate
  est.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.ipw))
  est.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr))
  est.ipw.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.ipw.eco))
  est.dr.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr.eco))
  
  # lower bound
  lower.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.ipw))
  lower.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr))
  lower.ipw.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.ipw.eco))
  lower.dr.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr.eco))
  
  # upper bound
  upper.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw))
  upper.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr))
  upper.ipw.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw.eco))
  upper.dr.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr.eco))
  
  # absolute bias
  bias.ipw <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw)), na.rm = TRUE)
  bias.dr <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.dr)), na.rm = TRUE)
  bias.ipw.eco <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw.eco)), na.rm = TRUE)
  bias.dr.eco <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.dr.eco)), na.rm = TRUE)
  
  # root mean squared error
  rmse.ipw <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw))^2, na.rm = TRUE))
  rmse.dr <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) - sapply(1:n.iter, function(i) out[[i]]$est.dr))^2, na.rm = TRUE))
  rmse.ipw.eco <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.ipw.eco))^2, na.rm = TRUE))
  rmse.dr.eco <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) - sapply(1:n.iter, function(i) out[[i]]$est.dr.eco))^2, na.rm = TRUE))
  
  # confidence length
  cl.ipw <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw) - sapply(1:n.iter, function(i) out[[i]]$lower.ipw), na.rm = TRUE)
  cl.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr) - sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  cl.ipw.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.ipw.eco) - sapply(1:n.iter, function(i) out[[i]]$lower.ipw.eco), na.rm = TRUE)
  cl.dr.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr.eco) - sapply(1:n.iter, function(i) out[[i]]$lower.dr.eco), na.rm = TRUE)
  
  # coverage probability
  cp.ipw <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) < sapply(1:n.iter, function(i) out[[i]]$upper.ipw) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) > sapply(1:n.iter, function(i) out[[i]]$lower.ipw), na.rm = TRUE)
  cp.dr <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr) &
                       rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  cp.ipw.eco <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) < sapply(1:n.iter, function(i) out[[i]]$upper.ipw.eco) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) > sapply(1:n.iter, function(i) out[[i]]$lower.ipw.eco), na.rm = TRUE)
  cp.dr.eco <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr.eco) &
                      rowMeans(sapply(1:n.iter, function(i) out[[i]]$theta.eco)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr.eco), na.rm = TRUE)
  
  df <- rbind(df, data.frame(a.vals = rep(a.vals, times = 6),
                   est = c(theta, theta.eco, est.ipw, est.dr, est.ipw.eco, est.dr.eco), 
                   lower = c(rep(NA, 2*length(a.vals)), lower.ipw, lower.dr, lower.ipw.eco, lower.dr.eco),
                   upper = c(rep(NA, 2*length(a.vals)), upper.ipw, upper.dr, upper.ipw.eco, upper.dr.eco),
                   bias = c(rep(NA, 2*length(a.vals)), bias.ipw, bias.dr, bias.ipw.eco, bias.dr.eco),
                   rmse = c(rep(NA, 2*length(a.vals)), rmse.ipw, rmse.dr, rmse.ipw.eco, rmse.dr.eco),
                             cp = c(rep(NA, 2*length(a.vals)), cp.ipw, cp.dr, cp.ipw.eco, cp.dr.eco),
                             cl = c(rep(NA, 2*length(a.vals)), cl.ipw, cl.dr, cl.ipw.eco, cl.dr.eco),
                             adjust = rep(c("true", "true", "ipw", "dr", "ipw", "dr"), each = length(a.vals)),
                             estimand = rep(c("Individual", "Ecological", "Individual", 
                                              "Individual", "Ecological", "Ecological"), each = length(a.vals)),
                             gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))

}

save(df, file = "~/Github/causal-eco/Output/simulation_results.RData")

### Plots

df$label <- ifelse(df$adjust == "true", "True ERF",
                    ifelse(df$adjust == "dr", "DR Estimator", "IPW Estimator"))

df$scenario <- ifelse(df$gps_scen == "a" & df$out_scen == "a", "Correct Specification",
                       ifelse(df$gps_scen == "a" & df$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(df$gps_scen == "b" & df$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

df$label <- factor(df$label, levels = c("True ERF", "IPW Estimator", "DR Estimator"))
df$estimand <- factor(df$estimand, levels = c("Individual", "Ecological"))
df$scenario <- factor(df$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

df_tmp <- subset(df, !(gps_scen == "b" & out_scen == "b") & label != "IPW Estimator" & ss_scen == "a")

plot <- df_tmp %>%
  ggplot(aes(x = a.vals, y = est, color = factor(estimand), linetype = factor(label))) +
  geom_line(size = 1) +
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

output <- df %>% group_by(adjust, estimand, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse), cp = mean(cp), cl = mean(cl))

save(output, file = "~/Github/causal-eco/Output/simulation_summary.RData")
