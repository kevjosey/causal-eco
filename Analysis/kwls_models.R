library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(splines)
library(KernSmooth)
library(ggplot2)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/gam_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/erf_models.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/erc-strata/Functions/calibrate.R')
set.seed(42)

## Setup

# scenarios
scenarios <- expand.grid(dual = c("high", "low", ""), race = c("white","black","hispanic","asian", ""),
                         age_break = c("\\[65,75)","\\[75,85)","\\[85,95)", ""))
scen_names <- expand.grid(dual = c("high", "low", "both"), race = c("white","black","hispanic","asian","all"),
                          age_break = c("[65,75)","[75,85)","[85,95)",""))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$age_break <- as.character(scenarios$age_break)
a.vals <- seq(2, 31, length.out = 146)
n.boot <- 1000

# Load/Save models
dir_mod = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/DR_Age/'

filenames <- list.files(dir_mod, full.names = TRUE)
fnames <- list.files(dir_mod, full.names = FALSE)

## Run Models QD

for (i in 1:length(scenarios)) {
  
  print(i)
  
  scenario <- scenarios[i,]
  sname <- scen_names[i,]
  
  grep1 <- grep(scenario[1], fnames)
  grep2 <- grep(scenario[2], fnames)
  grep3 <- grep(scenario[3], fnames)
  idx <- 1:length(filenames)
  
  fn <- filenames[(idx %in% grep1) & (idx %in% grep2) & (idx %in% grep3)]
  
  w.id <- log.pop <- resid <- NULL
  muhat.mat <- phat.tmp  <- NULL
  
  for (j in 1:length(fn)) {
    
    load(paste0(fn[j]))
    
    x.tmp <- setDF(new_data$x)
    w.tmp <- setDF(new_data$w)
    wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
  
    x.id <- paste(x.tmp$zip, x.tmp$year, sep = "-")
    w.id <- paste(wx.tmp$zip, wx.tmp$year, sep = "-")
    
    # data
    y <- wx.tmp$dead
    a_x <- x.tmp$pm25
    
    # fitted weights
    weights.cal0 <- wx.tmp$cal0
    weights.cal0_trunc <- wx.tmp$cal0_trunc
    weights.cal1 <- wx.tmp$cal1ß
    weights.cal1_trunc <- wx.tmp$cal1_trunc
    
    log.pop <- c(log.pop, log(wx.tmp$time_count))
    w.id <- c(w.id,  paste(wx.tmp$zip, wx.tmp$year, sep = "-"))
    psi0 <- c(psi0, cal0*ybar)
    psi1 <- c(psi1, cal1*ybar)
    
  }
  
    if (is.null(bw)) {
    risk.est <- sapply(bw.seq, risk.fn, a.vals = a.vals,
                       psi = dat$psi, a = dat$a)
    bw <- c(bw.seq[which.min(risk.est)])
  }
  
  # fit exposure response curves
  target <- count_erf(resid, log.pop = log.pop, w.id = w.id, 
                      a = zip_data$pm25, x.id = zip_data$id,
                      a.vals = a.vals, phat.vals = phat.tmp, se.fit = TRUE)
  
  # extract estimates
  est_data <- data.frame(a.vals = a.vals, estimate = target$estimate, se = sqrt(target$variance))
  
  print(paste0("Fit Complete: Scenario ", i))
  print(Sys.time())
  
  save(individual_data, zip_data, est_data,
       file = paste0(dir_out, sname$dual, "_", sname$race, "_", sname$age_break, ".RData"))
  
  rm(individual_data, zip_data, model_data, est_data, target); gc()
  
}