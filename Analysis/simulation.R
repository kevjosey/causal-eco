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
source('~/Github/erc-strata/Functions/gam.R')
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
  ind_data <- data.frame(zip, w1 = w1, w2 = w2)
  data <- merge(ind_data, zip_data, by = "zip")
  
  if (out_scen == "b") {
    mu_out <- with(data, plogis(-2 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 
                                  0.25*(a - 8)*u1 + 0.25*(a - 8)*w1 - 
                                  0.5*w1 + 0.5*w2))
    lambda <- sapply(a.vals, function(a.new) mean(plogis(-2 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4 +
                                                           0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) -
                                                           0.25*(a.new - 8)*u1 + 0.25*(a.new - 8)*w1 -
                                                           0.5*w1 + 0.5*w2)))
  } else { # y_scen == "a"
    mu_out <- with(data, plogis(-2 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                                  0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 
                                  0.25*(a - 8)*x1 + 0.25*(a - 8)*w1 - 
                                  0.5*w1 + 0.5*w2))
    lambda <- sapply(a.vals, function(a.new) mean(plogis(-2 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4 +
                                                           0.25*(a.new - 10) - 0.75*cos(pi*(a.new - 6)/4) - 
                                                           0.25*(a.new - 8)*x1 + 0.25*(a.new - 8)*w1 -
                                                           0.5*w1 + 0.5*w2)))
  }
  
  data$y <- rbinom(n, 1, mu_out)
  
  strata_data <- data %>% group_by(zip, w1, w2) %>% 
    summarise(x1 = mean(x1), x2 = mean(x2), 
              x3 = mean(x3), x4 = mean(x4),
              a = mean(a), y = sum(y), n = n())
  
  zip_data <- subset(zip_data, zip %in% unique(strata_data$zip))
  dat <- merge(zip_data, strata_data %>% group_by(zip) %>% summarise(y = sum(y), n = sum(n)), by = "zip")
    
  ## LM GPS
  # pimod <- lm(a ~ x1 + x2 + x3 + x4, data = dat)
  # pimod.vals <- c(pimod$fitted.values)
  # pimod.sd <- sigma(pimod)
  
  # nonparametric density
  # a.std <- (dat$a - pimod.vals)/pimod.sd
  # dens <- density(a.std)
  # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd
  
  # ipw numerator
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   a.std.tmp <- (a.tmp - pimod.vals)/pimod.sd
  #   approx(x = dens$x, y = dens$y, xout = a.std.tmp)$y / pimod.sd
  # })
  # 
  # phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
  # phat <- predict(smooth.spline(a.vals, phat.vals), x = dat$a)$y
  # phat[phat < 0] <- .Machine$double.eps
  # 
  # dat$ipw <- phat/pihat # LM GPS
  

  ## Calibration weights
  x.mat <- model.matrix(~ x1 + x2 + x3 + x4, data = data.frame(dat))
  astar <- c(dat$a - mean(dat$a))/var(dat$a)
  astar2 <- c((dat$a - mean(dat$a))^2/var(dat$a) - 1)
  cmat <- cbind(x.mat*astar, astar2, x.mat)
  tm <- c(rep(0, ncol(x.mat) + 1), colSums(x.mat))
  
  # fit calibration model
  mod <- calibrate(cmat = cmat, target = tm)
  dat$cal <- mod$weights
  
  dat$ybar <- dat$y/dat$n
  dat$ybar[dat$y > dat$n] <- 1 - .Machine$double.ep
  dat$psi <- dat$ybar*dat$cal
  
  # estimate nuisancße outcome model with gam
  inner <- paste(c("x1", "x2", "x3", "x4"), collapse = " + ")
  fmla <- as.formula(paste0("ybar ~ ns(aa, df = 6) + ", inner, " + aa:(", inner, ")"))
  mumod <- glm(fmla, data = data.frame(aa = dat$a, dat),
               weights = dat$n, family = quasipoisson())
  muhat <- mumod$fitted.values
  
  muhat.list <- lapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(aa = a.tmp, dat)
    return(predict(mumod, newdata = xa.tmp, type = "response", se.fit = TRUE))
  })
  
  muhat.mat <- do.call(cbind, lapply(muhat.list, function(z) z$fit))
  muhat.se <- do.call(cbind, lapply(muhat.list, function(z) z$se.fit))
  
  dat$psi2 <- (dat$ybar - muhat)*dat$cal
  
  # grid search bandwidth
  risk.est <- sapply(bw.seq, risk.fn, a.vals = a.vals, psi = dat$psi2, a = dat$a, n = dat$n)
  bw <- c(bw.seq[which.min(risk.est)])
  
  # fit kwls regressions
  erf <- sapply(a.vals, kwls_est, psi = dat$psi, a = dat$a, bw = bw[1], 
                se.fit = TRUE, sandwich = TRUE,
                x = x.mat, astar = astar, astar2 = astar2, cmat = cmat, ipw = dat$cal)
  erf.eco <- sapply(a.vals, kwls_est, psi = dat$psi, a = dat$a, bw = bw[1], 
                    se.fit = TRUE, sandwich = TRUE, eco = TRUE,
                    x = x.mat, astar = astar, astar2 = astar2, cmat = cmat, ipw = dat$cal)
  dr.eco <- sapply(a.vals, kwls_est, psi = dat$psi, a = dat$a, bw = bw[1], 
                   se.fit = TRUE, sandwich = TRUE, eco = TRUE,
                   x = x.mat, astar = astar, astar2 = astar2, cmat = cmat, ipw = dat$cal)
  sig2 <- c(dr.eco[2,] + colMeans(muhat.se^2)/nrow(x.mat))
  mu <- colMeans(muhat.mat) + dr.eco[1,]
  
  return(list(est.erf = erf[1,], se.erf = sqrt(erf[2,]),
              est.erf.eco = erf.eco[1,], se.erf.eco = sqrt(erf.eco[2,]), 
              est.dr = mu, se.dr = sqrt(sig2),
              lower.erf = erf[1,] - 1.96*sqrt(erf[2,]), upper.erf = erf[1,] + 1.96*sqrt(erf[2,]),
              lower.erf.eco = erf.eco[1,] - 1.96*sqrt(erf.eco[2,]), upper.erf.eco = erf.eco[1,] + 1.96*sqrt(erf.eco[2,]),
              lower.dr = mu - 1.96*sqrt(sig2), upper.dr = mu + 1.96*sqrt(sig2),
              lambda = lambda, a.vals = a.vals))
  
}

# cl <- makePSOCKcluster(25)
# clusterExport(cl, "fit_sim")
# 
# parReplicate <- function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE)
#   parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
#             substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)
# 
# clusterEvalQ(cl, {
#   
#   ## set up each worker.  Could also use clusterExport()
#   library(dplyr)
#   library(parallel)
#   library(data.table)
#   library(tidyr)
#   library(dplyr)
#   library(magrittr)
#   library(splines)
#   library(gam)
#   library(KernSmooth)
#   
#   source('~/Github/erc-strata/Functions/gam_models.R')
#   source('~/Github/erc-strata/Functions/erf_models.R')
#   source('~/Github/erc-strata/Functions/calibrate.R')
#   
# })

### Run Simulation

scenarios = expand.grid(n = c(100000), m = c(5000), gps_scen = c("a", "b"), out_scen = c("a", "b"), ss_scen = c("a"))
a.vals <- seq(4, 12, length.out = 81)
n.iter <- 200
dat <- data.frame()

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
                   a.vals = a.vals, bw.seq = bw.seq, mc.cores = 25)
  
  lambda <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda))
  
  # estimate
  est.erf <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.erf))
  est.erf.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.erf.eco))
  est.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$est.dr))
  
  # lower bound
  lower.erf <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.erf))
  lower.erf.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.erf.eco))
  lower.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$lower.dr))
  
  # upper bound
  upper.erf <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.erf))
  upper.erf.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.erf.eco))
  upper.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr))
  
  # absolute bias
  bias.erf <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.erf)), na.rm = TRUE)
  bias.erf.eco <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.erf.eco)), na.rm = TRUE)
  bias.dr <- rowMeans(abs(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr)), na.rm = TRUE)
  
  # root mean squared error
  rmse.erf <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.erf))^2, na.rm = TRUE))
  rmse.erf.eco <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.erf.eco))^2, na.rm = TRUE))
  rmse.dr <- sqrt(rowMeans((rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) - sapply(1:n.iter, function(i) out[[i]]$est.dr))^2, na.rm = TRUE))
  
  # confidence length
  cl.erf <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.erf) - sapply(1:n.iter, function(i) out[[i]]$lower.erf), na.rm = TRUE)
  cl.erf.eco <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.erf.eco) - sapply(1:n.iter, function(i) out[[i]]$lower.erf.eco), na.rm = TRUE)
  cl.dr <- rowMeans(sapply(1:n.iter, function(i) out[[i]]$upper.dr) - sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  
  # coverage probability
  cp.erf <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.erf) &
                        rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.erf), na.rm = TRUE)
  
  cp.erf.eco <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.erf.eco) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.erf.eco), na.rm = TRUE)
  
  cp.dr <- rowMeans(rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) < sapply(1:n.iter, function(i) out[[i]]$upper.dr) &
                           rowMeans(sapply(1:n.iter, function(i) out[[i]]$lambda)) > sapply(1:n.iter, function(i) out[[i]]$lower.dr), na.rm = TRUE)
  
  dat <- rbind(dat, data.frame(a.vals = rep(a.vals, times = 4),
                               est = c(lambda, est.erf, est.erf.eco, est.dr), 
                               lower = c(rep(NA, length(a.vals)), lower.erf, lower.erf.eco, lower.dr),
                               upper = c(rep(NA, length(a.vals)), upper.erf, upper.erf.eco, upper.dr),
                               bias = c(rep(NA, length(a.vals)), bias.erf, bias.erf.eco, bias.dr),
                               rmse = c(rep(NA, length(a.vals)), rmse.erf, rmse.erf.eco, rmse.dr),
                               cp = c(rep(NA, length(a.vals)), cp.erf, cp.erf.eco, cp.dr),
                               cl = c(rep(NA, length(a.vals)), cl.erf, cl.erf.eco, cl.dr),
                               adjust = rep(c("true", "erf", "erf.eco", "gam"), each = length(a.vals)),
                               gps_scen = gps_scen, out_scen = out_scen, ss_scen = ss_scen, m = m, n = n))

}

save(dat, file = "~/Github/erc-strata/Output/simulation_results.RData")

### Plots

dat$label <- ifelse(dat$adjust == "true", "True ERF",
                    ifelse(dat$adjust == "erf", "Unweighted",
                           ifelse(dat$adjust == "erf.eco", "Scaled Kernel Weight", "Doubly-Robust")))
dat$scenario <- ifelse(dat$gps_scen == "a" & dat$out_scen == "a", "Correct Specification",
                       ifelse(dat$gps_scen == "a" & dat$out_scen == "b", "Outcome Model Misspecification",
                              ifelse(dat$gps_scen == "b" & dat$out_scen == "a", "GPS Misspecification", "Incorrect Specification")))

dat$label <- factor(dat$label, levels = c("True ERF", "Unweighted", "Scaled Kernel Weight", "Doubly-Robust"))
dat$scenario <- factor(dat$scenario, levels = c("Correct Specification", "GPS Misspecification", "Outcome Model Misspecification", "Incorrect Specification"))

dat_tmp <- subset(dat, !(gps_scen == "b" & out_scen == "b" | m == 10000))

plot <- dat_tmp %>%
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
  scale_color_manual(values = c("#D81B60", "#F57328", "#004D40", "#ffdb58")) +
  scale_y_continuous(breaks = seq(0, 0.4, by = 0.05)) +
  scale_x_continuous(breaks = c(4,5,6,7,8,9,10,11,12))

pdf(file = "~/Github/erc-strata/Output/simulation_plot.pdf", width = 16, height = 8)
plot
dev.off()

output <- dat %>% group_by(adjust, gps_scen, out_scen, ss_scen, n, m) %>% summarise(bias = mean(bias), rmse = mean(rmse),  cp = mean(cp), cl = mean(cl))

save(output, file = "~/Github/erc-strata/Output/simulation_summary.RData")

