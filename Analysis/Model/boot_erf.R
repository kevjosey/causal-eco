
library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(gam)
library(SuperLearner)
library(xgboost)
library(ggplot2)
library(cobalt)

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/erf.R')
source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white","black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(4, 16, length.out = 121)
n.boot <- 1000

# Load/Save models
dir_data_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_data_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/rm/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

## Run Models QD

for (i in 1:9) {

  scenario <- scenarios[i,]
  load(paste0(dir_data_qd, scenario$dual, "_", scenario$race, "_qd.RData"))

  x.tmp <- setDF(new_data$x)
  w.tmp <- new_data$w
  # w.tmp <- setDF(subset(new_data$w, age_break == "[85,95)"))
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))

  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")

  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)

  y <- wx.tmp$dead
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))

  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count, age_break))
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count, age_break))
  } else if (scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count, age_break))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count, age_break))
  }

  target <- count_erf(a_w = a_w, y = y, w = w, log.pop = log.pop,
                      a_x = a_x, x = x, w.id = w.id, x.id = x.id,
                      a.vals = a.vals, bw = 1, loess = FALSE,
                      sl.lib = c("SL.mean", "SL.glm", "SL.xgboost"))

  print(paste0("Initial Fit Complete: Scenario ", i, " QD"))

  # boot_list <- mclapply(1:n.boot, mc.cores = 1, FUN = function(j, ...) {
  #
  #   print(j)
  #
  #   idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE)
  #   aa <- u.zip[idx]
  #   bb <- table(aa)
  #   wx.boot.tmp <- x.boot.tmp <- NULL
  #   w.id.boot <- x.id.boot <- NULL
  #
  #   for (k in 1:max(bb)) {
  #     cc <- wx.tmp[wx.tmp$zip %in% names(bb[which(bb == k)]),]
  #     dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
  #     for (l in 1:k) {
  #       wx.boot.tmp <- rbind(wx.boot.tmp, cc)
  #       x.boot.tmp <- rbind(x.boot.tmp, dd)
  #       x.id.boot <- c(x.id.boot, paste(dd$zip, dd$year, l, sep = "-"))
  #       w.id.boot <- c(w.id.boot, paste(cc$zip, cc$year, l, sep = "-"))
  #     }
  #   }
  #
  #   y.boot <- wx.boot.tmp$dead
  #   a_x.boot <- x.boot.tmp$pm25
  #   a_w.boot <- wx.boot.tmp$pm25
  #   log.pop.boot <- log(wx.boot.tmp$time_count)
  #   x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
  #
  #   if (scenario$dual == 2 & scenario$race == "all") {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
  #   } else if (scenario$dual == 2) {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dead, time_count))
  #   } else if (scenario$race == "all") {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dual, dead, time_count))
  #   } else {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  #   }
  #
  #   boot_target <- count_erf(a_w = a_w.boot,  y = y.boot, w = w.boot, log.pop = log.pop.boot,
  #                            a_x = a_x.boot, x = x.boot, w.id = w.id.boot, x.id = x.id.boot,
  #                            a.vals = a.vals, bw = 1, se.fit = FALSE)
  #
  #   return(boot_target$estimate.cal)
  #
  # })

  individual_data <- data.frame(wx.tmp,
                                weights.lm = target$weights.lm_w,
                                weights.sl = target$weights.sl_w,
                                weights.cal = target$weights.cal_w)

  zip_data <- data.frame(x.tmp,
                         weights.lm = target$weights.lm_x,
                         weights.sl = target$weights.sl_x,
                         weights.cal = target$weights.cal_x)

  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.sl = target$estimate.sl, se.sl = sqrt(target$variance.sl),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal),
                         linear.lm = predict(target$fit.lm, newdata = data.frame(a = a.vals)),
                         linear.sl = predict(target$fit.sl, newdata = data.frame(a = a.vals)),
                         linear.cal = predict(target$fit.cal, newdata = data.frame(a = a.vals)))

  # boot_data <- data.frame(a.vals = a.vals,
  #                         estimate = target$estimate.cal,
  #                         Reduce(cbind, boot_list))
  # colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))

  extra <- list(n.zip = n.zip, lm.coef = target$fit.lm$coefficients,
               sl.coef = target$fit.sl$coefficients,
               cal.coef = target$fit.cal$coefficients,
               lm.vcov = vcov(target$fit.lm),
               sl.vcov = vcov(target$fit.sl),
               cal.vcov = vcov(target$fit.cal))

  print(paste0("Fit Complete: Scenario ", i, " QD"))
  save(individual_data, zip_data, est_data, extra,
       file = paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))

}

## Run Models RM

for (i in 1:9) {

  scenario <- scenarios[i,]
  load(paste0(dir_data_rm, scenario$dual, "_", scenario$race, "_rm.RData"))

  x.tmp <- setDF(new_data$x)
  w.tmp <- new_data$w
  wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))

  x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
  w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")

  u.zip <- unique(x.tmp$zip)
  n.zip <- length(u.zip)

  y <- wx.tmp$dead
  a_x <- x.tmp$pm25
  a_w <- wx.tmp$pm25
  log.pop <- log(wx.tmp$time_count)
  x <- subset(x.tmp, select = -c(zip, pm25))

  if (scenario$dual == 2 & scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))
  } else if (scenario$dual == 2) {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dead, time_count))
  } else if (scenario$race == "all") {
    w <- subset(wx.tmp, select = -c(zip, pm25, dual, dead, time_count))
  } else {
    w <- subset(wx.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  }

  target <- count_erf(a_w = a_w, y = y, w = w, log.pop = log.pop,
                      a_x = a_x, x = x, w.id = w.id, x.id = x.id,
                      a.vals = a.vals, bw = 1, loess = FALSE, 
                      sl.lib = c("SL.mean", "SL.glm", "SL.xgboost"))

  print(paste0("Initial Fit Complete: Scenario ", i, " RM"))

  # boot_list <- mclapply(1:n.boot, mc.cores = 1, FUN = function(j, ...) {
  #
  #   print(j)
  #
  #   idx <- sample(1:n.zip, n.zip/log(n.zip), replace = TRUE)
  #   aa <- u.zip[idx]
  #   bb <- table(aa)
  #   wx.boot.tmp <- x.boot.tmp <- NULL
  #   w.id.boot <- x.id.boot <- NULL
  #
  #   for (k in 1:max(bb)) {
  #     cc <- wx.tmp[wx.tmp$zip %in% names(bb[which(bb == k)]),]
  #     dd <- x.tmp[x.tmp$zip %in% names(bb[which(bb == k)]),]
  #     for (l in 1:k) {
  #       wx.boot.tmp <- rbind(wx.boot.tmp, cc)
  #       x.boot.tmp <- rbind(x.boot.tmp, dd)
  #       x.id.boot <- c(x.id.boot, paste(dd$zip, dd$year, l, sep = "-"))
  #       w.id.boot <- c(w.id.boot, paste(cc$zip, cc$year, l, sep = "-"))
  #     }
  #   }
  #
  #   y.boot <- wx.boot.tmp$dead
  #   a_x.boot <- x.boot.tmp$pm25
  #   a_w.boot <- wx.boot.tmp$pm25
  #   log.pop.boot <- log(wx.boot.tmp$time_count)
  #   x.boot <- subset(x.boot.tmp, select = -c(zip, pm25))
  #
  #   if (scenario$dual == 2 & scenario$race == "all") {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dead, time_count))
  #   } else if (scenario$dual == 2) {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dead, time_count))
  #   } else if (scenario$race == "all") {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, dual, dead, time_count))
  #   } else {
  #     w.boot <- subset(wx.boot.tmp, select = -c(zip, pm25, race, dual, dead, time_count))
  #   }
  #
  #   boot_target <- count_erf(a_w = a_w.boot,  y = y.boot, w = w.boot, log.pop = log.pop.boot,
  #                            a_x = a_x.boot, x = x.boot, w.id = w.id.boot, x.id = x.id.boot,
  #                            a.vals = a.vals, bw = 1, se.fit = FALSE)
  #
  #   return(boot_target$estimate.cal)
  #
  # })

  individual_data <- data.frame(wx.tmp,
                                weights.lm = target$weights.lm_w,
                                weights.sl = target$weights.sl_w,
                                weights.cal = target$weights.cal_w)

  zip_data <- data.frame(x.tmp,
                         weights.lm = target$weights.lm_x,
                         weights.sl = target$weights.sl_x,
                         weights.cal = target$weights.cal_x)

  est_data <- data.frame(a.vals = a.vals,
                         estimate.lm = target$estimate.lm, se.lm = sqrt(target$variance.lm),
                         estimate.sl = target$estimate.sl, se.sl = sqrt(target$variance.sl),
                         estimate.cal = target$estimate.cal, se.cal = sqrt(target$variance.cal),
                         linear.lm = predict(target$fit.lm, newdata = data.frame(a = a.vals)),
                         linear.sl = predict(target$fit.sl, newdata = data.frame(a = a.vals)),
                         linear.cal = predict(target$fit.cal, newdata = data.frame(a = a.vals)))

  # boot_data <- data.frame(a.vals = a.vals,
  #                         estimate = target$estimate.cal,
  #                         Reduce(cbind, boot_list))
  # colnames(boot_data) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))

  extra <- list(n.zip = n.zip, lm.coef = target$fit.lm$coefficients,
                sl.coef = target$fit.sl$coefficients,
                cal.coef = target$fit.cal$coefficients,
                lm.vcov = vcov(target$fit.lm),
                sl.vcov = vcov(target$fit.sl),
                cal.vcov = vcov(target$fit.cal))

  print(paste0("Fit Complete: Scenario ", i, " RM"))
  save(individual_data, zip_data, est_data, extra,
       file = paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))

}
