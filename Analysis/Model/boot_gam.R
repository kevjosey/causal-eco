
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(gam)
library(ggplot2)
library(cobalt)

set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all", "white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(3, 17, length.out = 71)
n.boot <- 1000

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/GAM_qd_new/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/GAM_rm_new/'

## Run Models QD

for (i in c(1,6,8)) {

  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd_new.RData"))

  if (scenario$race == "all") {
    w.tmp <- setDF(subset(new_data$w, race %in% c(1, 2)))
    w.tmp$race <- factor(w.tmp$race)
  } else {
    w.tmp <- setDF(new_data$w)
  }

  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))

  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)

  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))

  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + dual + race + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + dual + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  } else if (scenario$race == "all"){
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + race + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  }

  # estimate nuisance outcome model with gam
  mumod <- gam(fmla, data = data.frame(y = y, a = a_w, lp = log.pop, w),
               family = poisson(link = "log"))

  # predict potential outcomes and aggregate by person years
  target <- sapply(a.vals, function(a.tmp, mumod, w, log.pop, ...) {

    sum(predict(mumod, newdata = data.frame(a = a.tmp, lp = log.pop, w),
                type = "response"))/sum(exp(log.pop))

  }, mumod = mumod, w = w, log.pop = log.pop)

  print(paste0("Initial Fit Complete: Scenario ", i, " QD"))

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

    if (scenario$dual == 2 & scenario$race == "all") {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
    } else if (scenario$dual == 2) {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dead, time_count))
    } else if (scenario$race == "all") {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dual, dead, time_count))
    } else {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    }

    # estimate nuisance outcome model with gam
    mumod.boot <- gam(fmla, data = data.frame(y = y.boot, a = a_w.boot, lp = log.pop.boot, w.boot),
                      family = poisson(link = "log"))

    # predict potential outcomes and aggregate by person years
    boot_target <- sapply(a.vals, function(a.tmp, mumod.boot, w.boot, log.pop.boot, ...) {

      sum(predict(mumod.boot, newdata = data.frame(a = a.tmp, lp = log.pop.boot, w.boot),
                  type = "response"))/sum(exp(log.pop.boot))

    }, mumod.boot = mumod.boot, w.boot = w.boot, log.pop.boot = log.pop.boot)

    return(boot_target)

  })

  individual_data <- wx.tmp
  zip_data <- x.tmp
  boot_data <- data.frame(a.vals = a.vals, estimate = target, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))

  print(paste0("Bootstrap Complete: Scenario ", i, " QD"))
  save(individual_data, zip_data, boot_data, n.zip,
       file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))

}

## Run Models RM

for (i in c(1,6,8)) {

  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm_new.RData"))

  if (scenario$race == "all") {
    w.tmp <- setDF(subset(new_data$w, race %in% c(1, 2)))
    w.tmp$race <- factor(w.tmp$race)
  } else {
    w.tmp <- setDF(new_data$w)
  }

  x.tmp <- setDF(new_data$x)
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))

  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)

  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  y <- wx.tmp$dead
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))

  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + dual + race + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
    w$dual <- factor(w$dual)
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + dual + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
    w$dual <- factor(w$dual)
  } else if (scenario$race == "all"){
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + race + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    fmla <- formula(y ~ s(a, 5) + year + sex + age_break + mean_bmi +
                      smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue +
                      poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx +
                      summer_rmax + winter_rmax + regionNORTHEAST +regionSOUTH + regionWEST + offset(lp))
  }

  # estimate nuisance outcome model with gam
  mumod <- gam(fmla, data = data.frame(y = y, a = a_w, lp = log.pop, w),
               family = poisson(link = "log"))

  # predict potential outcomes and aggregate by person years
  target <- sapply(a.vals, function(a.tmp, mumod, w, log.pop, ...) {

    sum(predict(mumod, newdata = data.frame(a = a.tmp, lp = log.pop, w),
                type = "response"))/sum(exp(log.pop))

  }, mumod = mumod, w = w, log.pop = log.pop)

  print(paste0("Initial Fit Complete: Scenario ", i, " RM"))

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

    if (scenario$dual == 2 & scenario$race == "all") {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
    } else if (scenario$dual == 2) {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dead, time_count))
    } else if (scenario$race == "all") {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dual, dead, time_count))
    } else {
      w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
    }

    # estimate nuisance outcome model with gam
    mumod.boot <- gam(fmla, data = data.frame(y = y.boot, a = a_w.boot, lp = log.pop.boot, w.boot),
                      family = poisson(link = "log"))

    # predict potential outcomes and aggregate by person years
    boot_target <- sapply(a.vals, function(a.tmp, mumod.boot, w.boot, log.pop.boot, ...) {

      sum(predict(mumod.boot, newdata = data.frame(a = a.tmp, lp = log.pop.boot, w.boot),
                  type = "response"))/sum(exp(log.pop.boot))

    }, mumod.boot = mumod.boot, w.boot = w.boot, log.pop.boot = log.pop.boot)

    return(boot_target)

  })

  individual_data <- wx.tmp
  zip_data <- x.tmp
  boot_data <- data.frame(a.vals = a.vals, estimate = target, Reduce(cbind, boot_list))
  colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))

  print(paste0("Bootstrap Complete: Scenario ", i, " RM"))
  save(individual_data, zip_data, boot_data, n.zip,
       file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))

}
