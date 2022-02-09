
library(survival)
library(parallel)
library(data.table)
library(CausalGPS)
library(dplyr)
library(mgcv)

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nasph_ci3/ci3_analysis/pdez_measurementerror/Output/GAM/Trimmed/'
scenarios <- expand.grid(Sex = c("Both", "Male", "Female"), Race = c("All", "White", "Black", "Hispanic", "Asian"))
a.vals <- seq(3, 18, length.out = 100)

wrapper <- function(data, a.vals, n.boot = 1000, sex = c("both", "male", "female"), race = c("all", "white", "black", "hispanic", "asian")) {
  
  if (sex == "male") {
    sex <- 1
  } else if (sex == "female") {
    sex <- 2
  } else {
    sex <- c(1,2)
  }
  
  if (race == "white") {
    race <- 1
  } else if (race == "black") {
    race <- 2
  } else if (race == "hispanic") {
    race <- 5
  } else if (race == "asian") {
    race <- 4
  } else {
    race <- c(1,2,3,4,5,6)
  }
  
  sub_data <- data %>% filter(data$race %in% race & data$sex %in% sex)   
  idx_mat <- data.frame(cbind(1:nrow(sub_data), replicate(n.boot, sample(1:nrow(sub_data), replace = TRUE))))
  idx_list <- as.list(idx_mat)
  
  out <- mclapply(idx_list, mc.cores = 64, function(idx, ...) {
    
    dat <- sub_data[idx,]
    
    if (sex == "Both" & race == "All"){
      
      x <- model.frame(~ as.factor(sex) + as.factor(race) + as.factor(dual) +
                         as.factor(entry_age_break)+as.factor(followup_year) +
                         mean_bmi + smoke_rate + hispanic + pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region), data = dat)
      
    } else if (sex == "Both") {
      
      x <- model.frame(~ as.factor(race) + as.factor(dual) +
                         as.factor(entry_age_break) + as.factor(followup_year) +
                         mean_bmi + smoke_rate + hispanic + pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region), data = dat)
      
    } else if (race == "All") {
      
      x <- model.frame(~ as.factor(sex) + as.factor(dual) +
                         as.factor(entry_age_break) + as.factor(followup_year) +
                         mean_bmi + smoke_rate + hispanic + pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region), data = dat)
      
    } else {
      
      x <- model.frame(~ as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year) +
                         mean_bmi + smoke_rate + hispanic + pct_blk +
                         medhouseholdincome + medianhousevalue +
                         poverty + education + popdensity + pct_owner_occ +
                         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                         as.factor(year) + as.factor(region), data = dat)
      
    }
    
    a <- dat$pm25
    y <- dat$dead
    offset <- log(dat$time_count)
    
    mods <- fit_models(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
    
    ## g-formula
    gform_out <- sapply(a.vals, function(a.star, ...) {
      mean(predict(mods[[1]], newdata = data.frame(a = rep(a.star, nrow(subdata)), subdata)))
    })
    
    match_out <- predict(mods[[2]], newdata = data.frame(w = a.vals))
    dr_out <- predict(mods[[2]], newdata = data.frame(w = a.vals))
    curves <- list(gform = gform_out, match = match_out, dr = dr_out)
    return(curves)
    
  })
  
  gform <- out[1,1]
  match <- out[1,2]
  dr <- out[1,3]
  
  # TEST THESE INDEXES
  gform_boot <- out[2:(n.boot + 1),1]
  match_boot <- out[2:(n.boot + 1),2]
  dr_boot <- out[2:(n.boot + 1),3]
  
  list(estimates = data.frame(a.vals = a.vals, gform = gform, match = match, dr = dr), 
       bootstrap = list(a.vals = a.vals, gform = gform_boot, match = match_boot, dr = dr_boot))
  
}

## Randall Martin

load(paste0(dir_data,"covariates_rm.RData"))
load(paste0(dir_data,"aggregate_data_rm.RData"))
aggregate_data_rm<-merge(aggregate_data_rm,covariates_rm,by=c("zip","year"),all.x=T)

rm(covariates_rm, GPS_mod_rm); gc()

rm_result <- apply(scenarios, 1, function(scenario, ...){
  
  output <- wrapper(data = aggregate_data_rm,
                    a.vals = a.vals, n.boot = 1000, 
                    sex = scenario[1], race = scenario[2])
  
  estimates <- output$estimates
  bootstrap <- output$bootstrap
  fname_est <- paste0("/nfs/nsaph_ci3/ci3_analysis/kjosey_erc_strata/output/estimates", scenario[1],"_", scenario[2],"_RM.RData")
  fname_boot <- paste0("/nfs/nsaph_ci3/ci3_analysis/kjosey_erc_strata/output/bootstrap", scenario[1],"_", scenario[2],"_RM.RData")
  save(estimates, file = fname_est)
  save(bootstrap, file = fname_boot)
  
})

rm(aggregate_data_rm); gc()

## Qian Di

load(paste0(dir_data,"covariates_qd.RData"))
load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd<-merge(aggregate_data_qd,covariates_qd, by=c("zip","year"),all.x=T)
aggregate_data_qd$pm25 <- aggregate_data_qd$pm25_ensemble

rm(covariates_qd, GPS_mod_qd); gc()

a.vals <- seq(quantile(aggregate_data_qd$pm25, 0.025), quantile(aggregate_data_qd$pm25, 0.975), length.out = 100)

qd_result <- apply(scenarios, 1, function(scenario, ...){
  
  output <- wrapper(data = aggregate_data_qd, a.vals = a.vals, n.boot = 1000, 
                    sex = scenario[1], race = scenario[2])
  
  estimates <- output$estimates
  bootstrap <- output$boostrap
  fname_est <- paste0("/nfs/nsaph_ci3/ci3_analysis/kjosey_erc_strata/output/estimates/", scenario[1],"_", scenario[2],"_QD.RData")
  fname_boot <- paste0("/nfs/nsaph_ci3/ci3_analysis/kjosey_erc_strata/output/bootstrap/", scenario[1],"_", scenario[2],"_QD.RData")
  save(estimates, file = fname_est)
  save(bootstrap, file = fname_boot)
  
})
