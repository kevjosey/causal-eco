library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(ggplot2)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
set.seed(42)

# exposure points
a.vals <- seq(2, 31, length.out = 146)
bw.seq <- seq(0.1, 3, length.out = 15)

bootstrap_data <- function(data, index, u.zip) {
  
  n.zip <- length(u.zip)
  boot <- data.frame()
  
  aa <- u.zip[index]
  aa <- aa[which(aa %in% data$zip)]
  bb <- table(aa)
  
  for (j in 1:max(bb)) {
    
    cc <- data[data$zip %in% names(bb[bb == j]),]
    
    for (k in 1:j) {
      cc$boot.id <- paste(cc$id, k, sep = "-")
      boot <- rbind(boot, cc)
    }
    
  }
  
  return(boot)
  
}

### Fit Balancing Weights

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_Data/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_All/'
load(paste0(dir_data,"aggregate_data.RData"))

# Outcome and Person-Years At-Risk
w <- data.table(zip = aggregate_data$zip, year = aggregate_data$year, race = aggregate_data$race,
                female = aggregate_data$female, dual = aggregate_data$dual, entry_age_break = aggregate_data$entry_age_break,
                followup_year = aggregate_data$followup_year, dead = aggregate_data$dead, time_count = aggregate_data$time_count)[
                  ,lapply(.SD, sum), by = c("zip", "year", "race", "female", "dual", "entry_age_break", "followup_year")]

# ZIP Code Covariates
zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
          "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region")

x <- data.table(zip = aggregate_data$zip, year = aggregate_data$year,
                model.matrix(~ ., data = aggregate_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]

x$id <- paste(x$zip, x$year, sep = "-")
u.zip <- unique(w$zip)
m <- ceiling(1/10*length(u.zip)) # for m out of n bootstrap

index <- sample(1:length(u.zip), m, replace = FALSE)  # initialize bootstrap  index

## GPS Model

# zip-code data
x <- bootstrap_data(data = x, index = index, u.zip = u.zip)

x.tmp <- subset(x, select = -c(zip, pm25, id, boot.id))
x.tmp$year <- factor(x.tmp$year)
x.tmp <- x.tmp %>% mutate_if(is.numeric, scale)

## LM GPS

pimod <- lm(a ~ ., data = data.frame(a = x$pm25, x.tmp))
pimod.vals <- c(pimod$fitted.values)
pimod.sd <- sigma(pimod)

# nonparametric density
a.std <- c(x$pm25 - pimod.vals) / pimod.sd
dens <- density(a.std)
pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd

# ipw numerator
pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  std <- c(a.tmp - pimod.vals) / pimod.sd
  approx(x = dens$x, y = dens$y, xout = std)$y / pimod.sd
})

phat.vals <- colMeans(pihat.mat, na.rm = TRUE)
phat <- predict(smooth.spline(a.vals, phat.vals), x = x$pm25)$y
phat[phat < 0] <- .Machine$double.eps

x$ipw <- phat/pihat # LM GPS

# truncation
trunc0 <- quantile(x$ipw, 0.005)
trunc1 <- quantile(x$ipw, 0.995)
x$ipw[x$ipw < trunc0] <- trunc0
x$ipw[x$ipw > trunc1] <- trunc1

# format variables
w$zip <- factor(w$zip)
w$year <- factor(w$year)
w$female <- as.numeric(w$female)
w$race <- factor(w$race)
w$dual <- as.numeric(w$dual)
w$entry_age_break <- factor(w$entry_age_break)
w$followup_year <- factor(w$followup_year)

x$zip <- factor(x$zip)
x$year <- factor(x$year)

### Fit Outcome Model

# merge in ZIP-level covariates
wx <- inner_join(w, x, by = c("zip", "year"))
wx <- bootstrap_data(data = wx, index = index, u.zip = u.zip)

w.tmp <- subset(wx, select = -c(zip, pm25, dead, time_count,
                                id, boot.id, ipw))

# fit gam outcome model
model_data <- gam_models(y = wx$dead, a = wx$pm25, w = w.tmp, log.pop = log(wx$time_count), 
                         id = wx$id, weights = wx$ipw, a.vals = a.vals)

save(model_data, w, x, phat.vals,
     file = paste0(dir_mod, "causal-eco.RData"))

### Fit Exposure Response Function

load(paste0(dir_mod, "causal-eco.RData"))

# Separate Data into List
mat.list <- with(model_data, split(cbind(exp(log.pop), resid, muhat.mat), id))
wts <- do.call(c, lapply(split(exp(model_data$log.pop), model_data$id), sum))

# Aggregate by ZIP-code-year
mat <- do.call(rbind, lapply(mat.list, function(vec) {
  mat <- matrix(vec, ncol = length(a.vals) + 2)
  colSums(mat[,1]*mat[,-1,drop = FALSE])/sum(mat[,1])
} ))

mat.pool <- data.frame(id = names(mat.list), mat)
mhat.vals <- apply(mat.pool[,-(1:2)], 2, mean, na.rm = TRUE)
resid.dat <- inner_join(mat.pool[,(1:2)], data.frame(a = x$pm25, id = x$id), by = "id")
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
                  bw = 30, a.vals = a.vals, se.fit = TRUE, int.mat = int.mat)

est_data <- data.frame(a.vals = rep(a.vals, times = 3),
                       adjust = rep(c("none", "simple", "complex"), each = length(a.vals)),
                       est = c(none[1,], simple[1,], complex[1,]), 
                       lower = c(none[1,] - 1.96*sqrt(none[2,]), simple[1,] - 1.96*sqrt(simple[2,]), complex[1,] - 1.96*sqrt(complex[2,])),
                       upper = c(none[1,] + 1.96*sqrt(none[2,]), simple[1,] + 1.96*sqrt(simple[2,]), complex[1,] + 1.96*sqrt(complex[2,])))

save(est_data, w, x, file = paste0(dir_out, "causal-eco.RData"))

## Plot Mortality Rates

est_data$label <- ifelse(est_data$adjust == "none", "No Weighting", 
                         ifelse(est_data$adjust == "simple", "Weighted Likelihood", "Scaled Kernel Weight"))
est_data$label <- factor(est_data$label, levels = c("No Weighting", "Weighted Likelihood", "Scaled Kernel Weight"))
est_data_tmp <- subset(est_data, a.vals <= 15 & a.vals >= 5)

ar_plot <- est_data_tmp %>%
  ggplot(aes(x = a.vals, y = est, color = factor(label))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", 
       y = "Absolute Rate of Mortality",
       color = "Weighting Approach") +
  theme_bw() +
  coord_cartesian(xlim = c(5.45,14.55), 
                  ylim = c(0.045, 0.051)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#D81B60", "#F57328", "#004D40")) +
  scale_y_continuous(breaks = seq(0.04, 0.06, by = 0.001)) +
  scale_x_continuous(breaks = c(5:15))

pdf(file = "~/Figures/causal-eco.pdf", width = 10, height = 8)
ar_plot
dev.off()