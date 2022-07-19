
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(mgcv)
library(ggplot2)
library(cobalt)

set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all", "white", "black", "hispanic", "asian"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(3, 17, length.out = 71)
n.boot <- 1000

# Load/Save models
dir_data = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/GAM_qd_fix/'

## Run Models QD

load(paste0(dir_data, "2_all_qd.RData"))
w.tmp <- setDF(new_data$w)
w.tmp$race_num <- w.tmp$race
w.tmp$race <- with(w.tmp, ifelse(race_num == 1, "white", 
                                 ifelse(race_num == 2, "black", 
                                        ifelse(race_num == 5, "hispanic",
                                               ifelse(race_num == 4, "asian", "other")))))

x.tmp <- setDF(new_data$x)
w.tmp$dual <- w.tmp$dual
w.tmp$race_dual <- factor(paste0(w.tmp$race, w.tmp$dual))
wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
u.zip <- unique(x.tmp$zip)
n.zip <- length(u.zip)

a_x <- x.tmp$pm25
a_w <- wx.tmp$pm25
y <- wx.tmp$dead
log.pop <- log(wx.tmp$time_count)
x <- subset(x.tmp, select = -c(zip, pm25))

w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
fmla <- formula(y ~ s(a, k = 5, by = race_dual) + year + sex + age_break + mean_bmi +
                  smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                  poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                  summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))

# estimate nuisance outcome model with gam
mumod <- bam(fmla, data = data.frame(y = y, a = a_w, lp = log.pop, w),
             family = poisson(link = "log"), samfrac = 0.1)

pe <- data.frame()

# point-estimate
for (i in 1:nrow(scenarios)) {
  
  if (scenarios[i,]$dual == 2 &  scenarios[i,]$race == "all") {
    
    w.sub <- w
    log.pop.sub <- log.pop
    
  } else if (scenarios[i,]$race == "all") {
    
    w.sub <- subset(w, dual == scenarios[i,]$dual)
    log.pop.sub <- subset(log.pop, w$dual == scenarios[i,]$dual)
    
  } else if (scenarios[i,]$dual == 2) {
    
    w.sub <- subset(w, race == scenarios[i,]$race)
    log.pop.sub <- subset(log.pop, w$race == scenarios[i,]$race)
    
  } else {
    
    w.sub <- subset(w, race == scenarios[i,]$race & dual == scenarios[i,]$dual)
    log.pop.sub <- subset(log.pop, w$race == scenarios[i,]$race & w$dual == scenarios[i,]$dual)
    
  }
  
  # predict potential outcomes and aggregate by person years
  target <- sapply(a.vals, function(a.tmp, mumod, w, log.pop, ...) {
    
    sum(predict(mumod, newdata = data.frame(a = a.tmp, lp = log.pop, w),
                type = "response"))/sum(exp(log.pop))
    
  }, mumod = mumod, w = w.sub, log.pop = log.pop.sub)

  print(paste0("Initial Fit Complete: Scenario ", i, " QD"))
  pe <- rbind(pe, data.frame(race = scenarios[i,]$race, dual = scenarios[i,]$dual,
                             a.vals = a.vals, estimate = target))
   
}

# bootstrap
boot_list <- mclapply(1:n.boot, mc.cores = 1, FUN = function(j, ...) {
  
  print(j)
  
  idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE)
  aa <- u.zip[idx]
  bb <- table(aa)
  wx.boot.tmp <- NULL
  
  for (k in 1:max(bb)) {
    cc <- wx.tmp[wx.tmp$zip %in% names(bb[which(bb == k)]),]
    dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
    for (l in 1:k) {
      wx.boot.tmp <- rbind(wx.boot.tmp, cc)
    }
  }
  
  a_w.boot <- wx.boot.tmp$pm25
  y.boot <- wx.boot.tmp$dead
  log.pop.boot <- log(wx.boot.tmp$time_count)
  
  w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
  
  # estimate nuisance outcome model with gam
  mumod.boot <- gam(fmla, data = data.frame(y = y.boot, a = a_w.boot, lp = log.pop.boot, w.boot),
                    family = poisson(link = "log"))
  
  out <- vector(mode = "numeric", length = nrow(scenarios)*length(a.vals))
  
  for (i in 1:nrow(scenarios)) {
    
    if (scenarios[i,]$dual == 2 &  scenarios[i,]$race == "all") {
      
      w.boot.sub <- w.boot
      log.pop.boot.sub <- log.pop.boot
      
    } else if (scenarios[i,]$race == "all") {
      
      w.boot.sub <- subset(w.boot, dual == scenarios[i,]$dual)
      log.pop.boot.sub <- subset(log.pop.boot, w$dual == scenarios[i,]$dual)
      
    } else if (scenarios[i,]$dual == 2) {
      
      w.boot.sub <- subset(w.boot, race == scenarios[i,]$race)
      log.pop.boot.sub <- subset(log.pop.boot, w$race == scenarios[i,]$race)
      
    } else {
      
      w.boot.sub <- subset(w.boot, race == scenarios[i,]$race & dual == scenarios[i,]$dual)
      log.pop.boot.sub <- subset(log.pop.boot, w$race == scenarios[i,]$race & w$dual == scenarios[i,]$dual)
      
    }
    
    # predict potential outcomes and aggregate by person years
    boot_target <- sapply(a.vals, function(a.tmp, mumod.boot, w.boot, log.pop.boot, ...) {
      
      sum(predict(mumod.boot, newdata = data.frame(a = a.tmp, lp = log.pop.boot, w.boot),
                  type = "response"))/sum(exp(log.pop.boot))
      
    }, mumod.boot = mumod.boot, w.boot = w.boot.sub, log.pop.boot = log.pop.boot.sub)
  
    out <- c(out, boot_target)
    
  }
  
  return(out)
  
})

# save data
individual_data <- wx.tmp
zip_data <- x.tmp
boot_data <- cbind(pe, Reduce(cbind, boot_list))

colnames(boot_data) <- c("race", "dual", "a.vals", "estimate", paste0("boot", 1:n.boot))

print(paste0("Bootstrap Complete: Scenario ", i, " QD"))
save(individual_data, zip_data, boot_data, n.zip,
     file = paste0(dir_out, scenario$dual, "_", scenario$race, "_qd.RData"))
