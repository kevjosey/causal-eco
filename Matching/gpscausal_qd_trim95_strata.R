library(devtools)
require(doParallel)
library(data.table)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)
library(devtools)
try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
install_github("fasrc/CausalGPS", ref="develop")
#install_github("fasrc/CausalGPS", ref="937810da2350f5c937c11eab16c60f2ee9a2783f")
library("CausalGPS")
library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
library(fst)
library("parallel")
require(ggplot2)
require(cowplot)
require(ggExtra)

set.seed(1)
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Matching/'

load(paste0(dir_data,"aggregate_data_qd.RData"))
# Matching on single exposure level a, a.vals selects the caliper
a.vals <- seq(min(aggregate_data_qd$pm25_ensemble), max(aggregate_data_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])
aggregate_data_qd$year<-as.factor(aggregate_data_qd$year)
aggregate_data_qd$region<-as.factor(aggregate_data_qd$region)

#All
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])


match_pop_all <- generate_pseudo_pop(Y=covariates_qd$zip,
                                     w=covariates_qd$pm25_ensemble,
                                     c=covariates_qd[, c(4:19)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = TRUE,
                                     transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list(xgb_nrounds=c(50)),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 5,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
match_pop_data <- match_pop_all$pseudo_pop
covariates_qd_trim <- subset(covariates_qd,
                             pm25_ensemble < quantile(covariates_qd$pm25_ensemble,0.95)&
                               pm25_ensemble > quantile(covariates_qd$pm25_ensemble,0.05))

match_pop_data <- cbind(match_pop_data, covariates_qd_trim[,1:2])
aggregate_data_qd2 <- merge(aggregate_data_qd, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_female_qd.RData")
covariates_white_female_qd$year<-as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region<-as.factor(covariates_white_female_qd$region)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
a.vals <- seq(min(white_female_qd$pm25_ensemble), max(white_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_white_female_qd1 <- generate_pseudo_pop(Y=covariates_white_female_qd$zip,
                                                  w=covariates_white_female_qd$pm25_ensemble,
                                                  c=covariates_white_female_qd[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transform = TRUE,
                                                  transformers = list("pow2", "pow3"),
                                                  sl_lib = c("m_xgboost"),
                                                  params = list(xgb_nrounds=c(50)),
                                                  nthread = 8, # number of cores, you can change,
                                                  covar_bl_method = "absolute",
                                                  covar_bl_trs = 0.1,
                                                  trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                  optimized_compile = TRUE, #created a column counter for how many times matched,
                                                  max_attempt = 5,
                                                  matching_fun = "matching_l1",
                                                  delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                  scale = 1.0)
match_pop_white_female_qd <- match_pop_white_female_qd1$pseudo_pop
covariates_white_female_qd_trim <- subset(covariates_white_female_qd,
                                          pm25_ensemble < quantile(covariates_white_female_qd$pm25_ensemble,0.95)&
                                            pm25_ensemble > quantile(covariates_white_female_qd$pm25_ensemble,0.05))
match_pop_white_female_qd <- cbind(match_pop_white_female_qd, covariates_white_female_qd_trim[,1:2])
match_pop_white_female_qd2 <- merge(white_female_qd, match_pop_white_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_female_qd)

#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_male_qd.RData")
covariates_white_male_qd$year<-as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region<-as.factor(covariates_white_male_qd$region)

white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
a.vals <- seq(min(white_male_qd$pm25_ensemble), max(white_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_white_male_qd1 <- generate_pseudo_pop(Y=covariates_white_male_qd$zip,
                                                w=covariates_white_male_qd$pm25_ensemble,
                                                c=covariates_white_male_qd[, c(4:19)],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = TRUE,
                                                transformers = list("pow2", "pow3"),
                                                sl_lib = c("m_xgboost"),
                                                params = list(xgb_nrounds=c(50)),
                                                nthread = 8, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.1,
                                                trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                optimized_compile = TRUE, #created a column counter for how many times matched,
                                                max_attempt = 5,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
match_pop_white_male_qd <- match_pop_white_male_qd1$pseudo_pop
covariates_white_male_qd_trim <- subset(covariates_white_male_qd,
                                        pm25_ensemble < quantile(covariates_white_male_qd$pm25_ensemble,0.95)&
                                          pm25_ensemble > quantile(covariates_white_male_qd$pm25_ensemble,0.05))
match_pop_white_male_qd <- cbind(match_pop_white_male_qd, covariates_white_male_qd_trim[,1:2])
match_pop_white_male_qd2 <- merge(white_male_qd, match_pop_white_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_male_qd)

#black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_female_qd.RData")
covariates_black_female_qd$year<-as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region<-as.factor(covariates_black_female_qd$region)

black_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
a.vals <- seq(min(black_female_qd$pm25_ensemble), max(black_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_black_female_qd1 <- generate_pseudo_pop(Y=covariates_black_female_qd$zip,
                                                  w=covariates_black_female_qd$pm25_ensemble,
                                                  c=covariates_black_female_qd[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transform = TRUE,
                                                  transformers = list("pow2", "pow3"),
                                                  sl_lib = c("m_xgboost"),
                                                  params = list(xgb_nrounds=c(50)),
                                                  nthread = 8, # number of cores, you can change,
                                                  covar_bl_method = "absolute",
                                                  covar_bl_trs = 0.1,
                                                  trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                  optimized_compile = TRUE, #created a column counter for how many times matched,
                                                  max_attempt = 5,
                                                  matching_fun = "matching_l1",
                                                  delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                  scale = 1.0)
match_pop_black_female_qd <- match_pop_black_female_qd1$pseudo_pop
covariates_black_female_qd_trim <- subset(covariates_black_female_qd,
                                          pm25_ensemble < quantile(covariates_black_female_qd$pm25_ensemble,0.95)&
                                            pm25_ensemble > quantile(covariates_black_female_qd$pm25_ensemble,0.05))
match_pop_black_female_qd <- cbind(match_pop_black_female_qd, covariates_black_female_qd_trim[,1:2])
match_pop_black_female_qd2 <- merge(black_female_qd, match_pop_black_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_female_qd)

#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_male_qd.RData")
covariates_black_male_qd$year<-as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region<-as.factor(covariates_black_male_qd$region)

black_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
a.vals <- seq(min(black_male_qd$pm25_ensemble), max(black_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])


match_pop_black_male_qd1 <- generate_pseudo_pop(Y=covariates_black_male_qd$zip,
                                                w=covariates_black_male_qd$pm25_ensemble,
                                                c=covariates_black_male_qd[, c(4:19)],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = TRUE,
                                                transformers = list("pow2", "pow3"),
                                                sl_lib = c("m_xgboost"),
                                                params = list(xgb_nrounds=c(50)),
                                                nthread = 8, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.1,
                                                trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                optimized_compile = TRUE, #created a column counter for how many times matched,
                                                max_attempt = 5,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
match_pop_black_male_qd <- match_pop_black_male_qd1$pseudo_pop
covariates_black_male_qd_trim <- subset(covariates_black_male_qd,
                                        pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble,0.95)&
                                          pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble,0.05))
match_pop_black_male_qd <- cbind(match_pop_black_male_qd, covariates_black_male_qd_trim[,1:2])
match_pop_black_male_qd2 <- merge(black_male_qd, match_pop_black_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_male_qd)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_female_qd.RData")
covariates_hispanic_female_qd$year<-as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region<-as.factor(covariates_hispanic_female_qd$region)

hispanic_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
a.vals <- seq(min(hispanic_female_qd$pm25_ensemble), max(hispanic_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_hispanic_female_qd1 <- generate_pseudo_pop(Y=covariates_hispanic_female_qd$zip,
                                                     w=covariates_hispanic_female_qd$pm25_ensemble,
                                                     c=covariates_hispanic_female_qd[, c(4:19)],
                                                     ci_appr = "matching",
                                                     pred_model = "sl",
                                                     gps_model = "parametric",
                                                     use_cov_transform = TRUE,
                                                     transformers = list("pow2", "pow3"),
                                                     sl_lib = c("m_xgboost"),
                                                     params = list(xgb_nrounds=c(50)),
                                                     nthread = 8, # number of cores, you can change,
                                                     covar_bl_method = "absolute",
                                                     covar_bl_trs = 0.1,
                                                     trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                                     max_attempt = 5,
                                                     matching_fun = "matching_l1",
                                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                     scale = 1.0)
match_pop_hispanic_female_qd <- match_pop_hispanic_female_qd1$pseudo_pop
covariates_hispanic_female_qd_trim <- subset(covariates_hispanic_female_qd,
                                             pm25_ensemble < quantile(covariates_hispanic_female_qd$pm25_ensemble,0.95)&
                                               pm25_ensemble > quantile(covariates_hispanic_female_qd$pm25_ensemble,0.05))
match_pop_hispanic_female_qd <- cbind(match_pop_hispanic_female_qd, covariates_hispanic_female_qd_trim[,1:2])
match_pop_hispanic_female_qd2 <- merge(hispanic_female_qd, match_pop_hispanic_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_female_qd)

#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_male_qd.RData")
covariates_hispanic_male_qd$year<-as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region<-as.factor(covariates_hispanic_male_qd$region)

hispanic_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
a.vals <- seq(min(hispanic_male_qd$pm25_ensemble), max(hispanic_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_hispanic_male_qd1 <- generate_pseudo_pop(Y=covariates_hispanic_male_qd$zip,
                                                   w=covariates_hispanic_male_qd$pm25_ensemble,
                                                   c=covariates_hispanic_male_qd[, c(4:19)],
                                                   ci_appr = "matching",
                                                   pred_model = "sl",
                                                   gps_model = "parametric",
                                                   use_cov_transform = TRUE,
                                                   transformers = list("pow2", "pow3"),
                                                   sl_lib = c("m_xgboost"),
                                                   params = list(xgb_nrounds=c(50)),
                                                   nthread = 8, # number of cores, you can change,
                                                   covar_bl_method = "absolute",
                                                   covar_bl_trs = 0.1,
                                                   trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                   optimized_compile = TRUE, #created a column counter for how many times matched,
                                                   max_attempt = 5,
                                                   matching_fun = "matching_l1",
                                                   delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                   scale = 1.0)
match_pop_hispanic_male_qd <- match_pop_hispanic_male_qd1$pseudo_pop
covariates_hispanic_male_qd_trim <- subset(covariates_hispanic_male_qd,
                                           pm25_ensemble < quantile(covariates_hispanic_male_qd$pm25_ensemble,0.95)&
                                             pm25_ensemble > quantile(covariates_hispanic_male_qd$pm25_ensemble,0.05))
match_pop_hispanic_male_qd <- cbind(match_pop_hispanic_male_qd, covariates_hispanic_male_qd_trim[,1:2])
match_pop_hispanic_male_qd2 <- merge(hispanic_male_qd, match_pop_hispanic_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_male_qd)

#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_female_qd.RData")
covariates_asian_female_qd$year<-as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region<-as.factor(covariates_asian_female_qd$region)

asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
a.vals <- seq(min(asian_female_qd$pm25_ensemble), max(asian_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_asian_female_qd1 <- generate_pseudo_pop(Y=covariates_asian_female_qd$zip,
                                                  w=covariates_asian_female_qd$pm25_ensemble,
                                                  c=covariates_asian_female_qd[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transform = TRUE,
                                                  transformers = list("pow2", "pow3"),
                                                  sl_lib = c("m_xgboost"),
                                                  params = list(xgb_nrounds=c(50)),
                                                  nthread = 8, # number of cores, you can change,
                                                  covar_bl_method = "absolute",
                                                  covar_bl_trs = 0.1,
                                                  trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                  optimized_compile = TRUE, #created a column counter for how many times matched,
                                                  max_attempt = 5,
                                                  matching_fun = "matching_l1",
                                                  delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                  scale = 1.0)
match_pop_asian_female_qd <- match_pop_asian_female_qd1$pseudo_pop
covariates_asian_female_qd_trim <- subset(covariates_asian_female_qd,
                                          pm25_ensemble < quantile(covariates_asian_female_qd$pm25_ensemble,0.95)&
                                            pm25_ensemble > quantile(covariates_asian_female_qd$pm25_ensemble,0.05))
match_pop_asian_female_qd <- cbind(match_pop_asian_female_qd, covariates_asian_female_qd_trim[,1:2])
match_pop_asian_female_qd2 <- merge(asian_female_qd, match_pop_asian_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_female_qd)

#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_male_qd.RData")
covariates_asian_male_qd$year<-as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region<-as.factor(covariates_asian_male_qd$region)

asian_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)
a.vals <- seq(min(asian_male_qd$pm25_ensemble), max(asian_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_asian_male_qd1 <- generate_pseudo_pop(Y=covariates_asian_male_qd$zip,
                                                w=covariates_asian_male_qd$pm25_ensemble,
                                                c=covariates_asian_male_qd[, c(4:19)],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = TRUE,
                                                transformers = list("pow2", "pow3"),
                                                sl_lib = c("m_xgboost"),
                                                params = list(xgb_nrounds=c(50)),
                                                nthread = 8, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.1,
                                                trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                                optimized_compile = TRUE, #created a column counter for how many times matched,
                                                max_attempt = 5,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
match_pop_asian_male_qd <- match_pop_asian_male_qd1$pseudo_pop
covariates_asian_male_qd_trim <- subset(covariates_asian_male_qd,
                                        pm25_ensemble < quantile(covariates_asian_male_qd$pm25_ensemble,0.95)&
                                          pm25_ensemble > quantile(covariates_asian_male_qd$pm25_ensemble,0.05))
match_pop_asian_male_qd <- cbind(match_pop_asian_male_qd, covariates_asian_male_qd_trim[,1:2])
match_pop_asian_male_qd2 <- merge(asian_male_qd, match_pop_asian_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_male_qd)

save(match_pop_all,match_pop_data, aggregate_data_qd2,
     match_pop_white_female_qd1, match_pop_white_female_qd, match_pop_white_female_qd2,
     match_pop_white_male_qd1, match_pop_white_male_qd, match_pop_white_male_qd2,
     match_pop_black_female_qd1, match_pop_black_female_qd, match_pop_black_female_qd2,
     match_pop_black_male_qd1, match_pop_black_male_qd, match_pop_black_male_qd2,
     match_pop_hispanic_female_qd1, match_pop_hispanic_female_qd, match_pop_hispanic_female_qd2,
     match_pop_hispanic_male_qd1, match_pop_hispanic_male_qd, match_pop_hispanic_male_qd2,
     match_pop_asian_female_qd1, match_pop_asian_female_qd, match_pop_asian_female_qd2,
     match_pop_asian_male_qd1, match_pop_asian_male_qd, match_pop_asian_male_qd2,
     file="/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/match_pseudo_pop_qd_strata.RData")

#Associations
#Calculating Associations
#All
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/match_pseudo_pop_qd_strata.RData")
main_gnm<-c()
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(sex)+as.factor(race)+as.factor(dual)+
                                 as.factor(entry_age_break)+as.factor(followup_year)),
                            data=aggregate_data_qd2, family=poisson(link="log"), weights=counter))
main_gnm[1]<-exp(10*matchingqd_gnm$coefficients[2])
main_gnm[1]

#White female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_white_female_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[2]<-exp(10*matchingqd_gnm$coefficients[2])


#White male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)) 
                            +(as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_white_male_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[3]<-exp(10*matchingqd_gnm$coefficients[2])

#Black female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_black_female_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[4]<-exp(10*matchingqd_gnm$coefficients[2])

#Black male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_black_male_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[5]<-exp(10*matchingqd_gnm$coefficients[2])

#Hispanic female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_hispanic_female_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[6]<-exp(10*matchingqd_gnm$coefficients[2])

#Hispanic male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_hispanic_male_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[7]<-exp(10*matchingqd_gnm$coefficients[2])

#Asian female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_asian_female_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[8]<-exp(10*matchingqd_gnm$coefficients[2])

#Asian male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_asian_male_qd2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[9]<-exp(10*matchingqd_gnm$coefficients[2])

#ERCS

matchingqd_gam <-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',k=3) +
                             as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                             offset(log(time_count))
                           , data=aggregate_data_qd2,family=poisson(link="log"), weights=counter)

white_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_white_female_qd2,family=poisson(link="log"), 
                                        weights=countrt)

white_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_white_male_qd2,family=poisson(link="log"), weights=counter)

black_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_black_female_qd2,family=poisson(link="log"), weights=counter)

black_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_black_male_qd2,family=poisson(link="log"), weights=counter)

hispanic_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_hispanic_female_qd2,family=poisson(link="log"), weights=counter)

hispanic_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_hispanic_male_qd2,family=poisson(link="log"), weights=counter)

asian_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_asian_female_qd2,family=poisson(link="log"))

asian_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_asian_male_qd2,family=poisson(link="log"), weights=counter)
#Plots
##QD
test.data.qd<-function(x, gm){
  test<-data.frame(pm25_ensemble = seq(min(x$pm25_ensemble), max(x$pm25_ensemble),length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   time_count=rep(1, 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50)

  )
  
  base<-data.frame(pm25_ensemble = rep(min(x$pm25),50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   time_count=rep(1, 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50)

  )
  
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
  base.vals<-predict(gm, newdata=base, type="link", se.fit = TRUE)
  dly<- pred.vals$fit - base.vals$fit
  return(gm$family$linkinv(dly))
  
}

#Plotting
#Plotting
testallqd<-test.data.qd(aggregate_data_qd2, matchingqd_gam)
test_white_female_qd<-test.data.qd(match_pop_white_female_qd2 , white_femaleqd_matching_gam)
test_white_male_qd<-test.data.qd(match_pop_white_male_qd2, white_maleqd_matching_gam)
test_black_female_qd<-test.data.qd(match_pop_black_female_qd2, black_femaleqd_matching_gam)
test_black_male_qd<-test.data.qd(match_pop_black_male_qd2, black_maleqd_matching_gam)
test_hispanic_female_qd<-test.data.qd(match_pop_hispanic_female_qd2, hispanic_femaleqd_matching_gam)
test_hispanic_male_qd<-test.data.qd(match_pop_hispanic_male_qd2, hispanic_maleqd_matching_gam)
test_asian_female_qd<-test.data.qd(match_pop_asian_female_qd2, asian_femaleqd_matching_gam)
test_asian_male_qd<-test.data.qd(match_pop_asian_male_qd2, asian_maleqd_matching_gam)

#Bootstrap
files<-list.files("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingqd/boots5strata/", 
                  pattern = "\\.RData",full.names=TRUE)
tallqd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_white_female_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_white_male_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_black_female_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_black_male_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_hispanic_female_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_hispanic_male_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_asian_female_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_asian_male_qd<-as.data.frame(matrix(ncol=50, nrow=length(files)))

t_gnm<-as.data.frame(matrix(nrow=length(files), ncol=9))
colnames(t_gnm)<-c("all", "white_female", "white_male", "black_female", "black_male",
                   "hispanic_female", "hispanic_male", "asian_female", "asian_male")

for(i in 1:length(files)){
  load(files[i])
       print(i)
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(sex)+as.factor(race)+as.factor(dual)+
                                   as.factor(entry_age_break)+as.factor(followup_year)),
                              data=aggregate_data_qd2, family=poisson(link="log"), weights=counter))
  t_gnm$all[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #White female
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_white_female_qd2, family=poisson(link="log"), 
                              weights=counter))
  t_gnm$white_female[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  
  #White male
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)) 
                              +(as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_white_male_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$white_male[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Black female
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_black_female_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$black_female[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Black male
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_black_male_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$black_male[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Hispanic female
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_hispanic_female_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$hispanic_female[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Hispanic male
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_hispanic_male_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$hispanic_male[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Asian female
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_asian_female_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$asian_female[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #Asian male
  matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_asian_male_qd2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingqd_gnm$coefficients[2])
  t_gnm$asian_male[i]=exp(10*matchingqd_gnm$coefficients[2])
  
  #ERC
  matchingqd_gam1 <-mgcv::bam(dead~ s(pm25_ensemble,bs='cr',k=3) +
                               as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                               offset(log(time_count))
                             , data=aggregate_data_qd2,family=poisson(link="log"), weights=counter)

  
  white_femaleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                            offset(log(time_count))
                                          , data=match_pop_white_female_qd2,family=poisson(link="log", weights=counter))
  
  white_maleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_white_male_qd2,family=poisson(link="log"), weights=counter)
  
  black_femaleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                            offset(log(time_count))
                                          , data=match_pop_black_female_qd2,family=poisson(link="log"), weights=counter)
  
  black_maleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_black_male_qd2,family=poisson(link="log"), weights=counter)
  
  hispanic_femaleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                               as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                               offset(log(time_count))
                                             , data=match_pop_hispanic_female_qd2,family=poisson(link="log"), weights=counter)
  
  hispanic_maleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_hispanic_male_qd2,family=poisson(link="log"), weights=counter)
  
  asian_femaleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr',k=3) +
                                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                            offset(log(time_count))
                                          , data=match_pop_asian_female_qd2,family=poisson(link="log", weights=counter))
  
  asian_maleqd_matching_gam1 <-mgcv::bam(dead~ s(pm25_ensemble, bs='cr', k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_asian_male_qd2,family=poisson(link="log", weights=counter))
  
  tallqd[i,]<-test.data.qd(aggregate_data_qd2, matchingqd_gam)
  t_white_female_qd[i,]<-test.data.qd(match_pop_white_female_qd2 , white_femaleqd_matching_gam)
  t_white_male_qd[i,]<-test.data.qd(match_pop_white_male_qd2, white_maleqd_matching_gam)
  t_black_female_qd[i,]<-test.data.qd(match_pop_black_female_qd2, black_femaleqd_matching_gam)
  t_black_male_qd[i,]<-test.data.qd(match_pop_black_male_qd2, black_maleqd_matching_gam)
  t_hispanic_female_qd[i,]<-test.data.qd(match_pop_hispanic_female_qd2, hispanic_femaleqd_matching_gam)
  t_hispanic_male_qd[i,]<-test.data.qd(match_pop_hispanic_male_qd2, hispanic_maleqd_matching_gam)
  t_asian_female_qd[i,]<-test.data.qd(match_pop_asian_female_qd2, asian_femaleqd_matching_gam)
  t_asian_male_qd[i,]<-test.data.qd(match_pop_asian_male_qd2, asian_maleqd_matching_gam)
}

#GNM
quantile(t_gnm$all, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$white_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$white_male, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$black_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$black_male, prob=c(0.025, 0.975, na.rm=TRUE))
quantile(t_gnm$hispanic_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_male, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_male, prob=c(0.025, 0.975), na.rm=TRUE)
                  
#ERC
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingqd/boots5/'

load(paste0(dir_data,"aggregate_data_qd.RData"))
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)
num_uniq_zip <- length(unique(covariates_qd$zip))
rm(covariates_qd)

 #All
require(matrixStats)
testallqd<-as.data.frame(testallqd)
#std<-apply(log(tallqd), 2, sd, na.rm=TRUE)
#testallqd$lower<- apply(tallqd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#testallqd$upper<- apply(tallqd,2, quantile, prob=c(0.975), na.rm=TRUE)
#testallqd$lower <-exp(log(testallqd$testallqd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
#testallqd$higher<-exp(log(testallqd$testallqd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
std<-apply((tallqd), 2, sd, na.rm=TRUE)
testallqd$median<-apply(tallqd, 2, mean, na.rm=TRUE)
testallqd$lower <-(testallqd$testallqd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
testallqd$higher<-(testallqd$testallqd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)


testallqd$pm25<-seq(min(aggregate_data_qd2$pm25_ensemble), 
                    max(aggregate_data_qd2$pm25_ensemble), length.out = 50)

pall_qd<-ggplot(data=testallqd, mapping=aes(x=pm25))+geom_line(aes(y=testallqd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="rm PM2.5", y="All HR")

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_female_qd.RData")
covariates_white_female_qd$year<-as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region<-as.factor(covariates_white_female_qd$region)
num_uniq_zip <- length(unique(covariates_white_female_qd$zip))
rm(covariates_white_female_qd)

test_white_female_qd<-as.data.frame(test_white_female_qd)
std<-apply((t_white_female_qd), 2, sd, na.rm=TRUE)
test_white_female_qd$median<-apply(t_white_female_qd, 2, mean, na.rm=TRUE)
test_white_female_qd$lower <-(test_white_female_qd$test_white_female_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_white_female_qd$higher<-(test_white_female_qd$test_white_female_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#test_white_female_qd$lower<- apply(t_white_female_qd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#test_white_female_qd$upper<- apply(t_white_female_qd,2, quantile, prob=c(0.975), na.rm=TRUE)

test_white_female_qd$pm25<-seq(min(match_pop_white_female_qd2$pm25_ensemble), max(match_pop_white_female_qd2$pm25_ensemble), length.out = 50)

#White male
#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_male_qd.RData")
covariates_white_male_qd$year<-as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region<-as.factor(covariates_white_male_qd$region)
num_uniq_zip <- length(unique(covariates_white_male_qd$zip))
rm(covariates_white_male_qd)

test_white_male_qd<-as.data.frame(test_white_male_qd)
std<-apply((t_white_male_qd), 2, sd, na.rm=TRUE)
test_white_male_qd$median<-apply(t_white_male_qd, 2, mean, na.rm=TRUE)
test_white_male_qd$lower <-(test_white_male_qd$test_white_male_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_white_male_qd$higher<-(test_white_male_qd$test_white_male_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#test_white_male_qd$lower<- apply(t_white_male_qd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#test_white_male_qd$upper<- apply(t_white_male_qd,2, quantile, prob=c(0.975), na.rm=TRUE)
test_white_male_qd$pm25<-seq(min(match_pop_white_male_qd2$pm25_ensemble), max(match_pop_white_male_qd2$pm25_ensemble), length.out = 50)

#Black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_female_qd.RData")
covariates_black_female_qd$year<-as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region<-as.factor(covariates_black_female_qd$region)
num_uniq_zip <- length(unique(covariates_black_female_qd$zip))
rm(covariates_black_female_qd)

test_black_female_qd<-as.data.frame(test_black_female_qd)
std<-apply((t_black_female_qd), 2, sd, na.rm=TRUE)
test_black_female_qd$median<-apply(t_black_female_qd, 2, mean, na.rm=TRUE)
test_black_female_qd$lower <-(test_black_female_qd$test_black_female_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_black_female_qd$higher<-(test_black_female_qd$test_black_female_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#test_black_female_qd$lower<- apply(t_black_female_qd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#test_black_female_qd$upper<- apply(t_black_female_qd,2, quantile, prob=c(0.975), na.rm=TRUE)
test_black_female_qd$pm25<-seq(min(match_pop_black_female_qd2$pm25_ensemble), max(match_pop_black_female_qd2$pm25_ensemble), length.out = 50)

#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_male_qd.RData")
covariates_black_male_qd$year<-as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region<-as.factor(covariates_black_male_qd$region)
num_uniq_zip <- length(unique(covariates_black_male_qd$zip))
rm(covariates_black_male_qd)

test_black_male_qd<-as.data.frame(test_black_male_qd)
std<-apply((t_black_male_qd), 2, sd, na.rm=TRUE)
test_black_male_qd$median<-apply(t_black_male_qd, 2, mean, na.rm=TRUE)
test_black_male_qd$lower <-(test_black_male_qd$test_black_male_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_black_male_qd$higher<-(test_black_male_qd$test_black_male_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_black_male_qd$pm25<-seq(min(match_pop_black_male_qd2$pm25_ensemble), max(match_pop_black_male_qd2$pm25_ensemble), length.out = 50)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_female_qd.RData")
covariates_hispanic_female_qd$year<-as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region<-as.factor(covariates_hispanic_female_qd$region)
num_uniq_zip <- length(unique(covariates_hispanic_female_qd$zip))
rm(covariates_hispanic_female_qd)

test_hispanic_female_qd<-as.data.frame(test_hispanic_female_qd)
std<-apply((t_hispanic_female_qd), 2, sd, na.rm=TRUE)
test_hispanic_female_qd$median<-apply(t_hispanic_female_qd, 2, mean, na.rm=TRUE)
test_hispanic_female_qd$lower <-(test_hispanic_female_qd$test_hispanic_female_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_hispanic_female_qd$higher<-(test_hispanic_female_qd$test_hispanic_female_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_hispanic_female_qd$pm25<-seq(min(match_pop_hispanic_female_qd2$pm25_ensemble), max(match_pop_hispanic_female_qd2$pm25_ensemble), length.out = 50)

#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_male_qd.RData")
covariates_hispanic_male_qd$year<-as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region<-as.factor(covariates_hispanic_male_qd$region)
num_uniq_zip <- length(unique(covariates_white_female_qd$zip))
rm(covariates_hispanic_male_qd)

test_hispanic_male_qd<-as.data.frame(test_hispanic_male_qd)
std<-apply((t_hispanic_male_qd), 2, sd, na.rm=TRUE)
test_hispanic_male_qd$median<-apply(t_hispanic_male_qd, 2, mean, na.rm=TRUE)
test_hispanic_male_qd$lower <-(test_hispanic_male_qd$test_hispanic_male_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_hispanic_male_qd$higher<-(test_hispanic_male_qd$test_hispanic_male_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_hispanic_male_qd$pm25<-seq(min(match_pop_hispanic_male_qd2$pm25_ensemble), max(match_pop_hispanic_male_qd2$pm25_ensemble), length.out = 50)

#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_female_qd.RData")
covariates_asian_female_qd$year<-as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region<-as.factor(covariates_asian_female_qd$region)
num_uniq_zip <- length(unique(covariates_asian_female_qd$zip))
rm(covariates_asian_female_qd)

test_asian_female_qd<-as.data.frame(test_asian_female_qd)
std<-apply((t_asian_female_qd), 2, sd, na.rm=TRUE)
test_asian_female_qd$median<-apply(t_asian_female_qd, 2, mean, na.rm=TRUE)
test_asian_female_qd$lower <-(test_asian_female_qd$test_asian_female_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_asian_female_qd$higher<-(test_asian_female_qd$test_asian_female_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#test_asian_female_qd$lower<- apply(t_asian_female_qd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#test_asian_female_qd$upper<- apply(t_asian_female_qd,2, quantile, prob=c(0.975), na.rm=TRUE)
test_asian_female_qd$pm25<-seq(min(match_pop_asian_female_qd2$pm25_ensemble), max(match_pop_asian_female_qd2$pm25_ensemble), length.out = 50)

#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_male_qd.RData")
covariates_asian_male_qd$year<-as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region<-as.factor(covariates_asian_male_qd$region)
num_uniq_zip <- length(unique(covariates_asian_male_qd$zip))
rm(covariates_asian_male_qd)

test_asian_male_qd<-as.data.frame(test_asian_male_qd)
std<-apply((t_asian_male_qd), 2, sd, na.rm=TRUE)
test_asian_male_qd$median<-apply(t_asian_male_qd, 2, mean, na.rm=TRUE)
test_asian_male_qd$lower <-(test_asian_male_qd$test_asian_male_qd) - 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)
test_asian_male_qd$higher<-(test_asian_male_qd$test_asian_male_qd)+ 1.96*std*sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)

#test_asian_male_qd$lower<- apply(t_asian_male_qd, 2, quantile,  prob=c(0.025), na.rm=TRUE)
#test_asian_male_qd$upper<- apply(t_asian_male_qd,2, quantile, prob=c(0.975), na.rm=TRUE)
test_asian_male_qd$pm25<-seq(min(match_pop_asian_male_qd2$pm25_ensemble), max(match_pop_asian_male_qd2$pm25_ensemble), length.out = 50)

pall_qd<-ggplot(data=testallqd, mapping=aes(x=pm25))+geom_line(aes(y=testallqd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="qd PM2.5", y="All HR")

p_white_female_qd<-ggplot(data=test_white_female_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_white_female_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="White female HR")

p_white_male_qd<-ggplot(data=test_white_male_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_white_male_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="White male HR")

p_black_female_qd<-ggplot(data=test_black_female_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_black_female_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Black female HR")

p_black_male_qd<-ggplot(data=test_black_male_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_black_male_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Black male HR")

p_hispanic_female_qd<-ggplot(data=test_hispanic_female_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_hispanic_female_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Hispanic female HR")

p_hispanic_male_qd<-ggplot(data=test_hispanic_male_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_hispanic_male_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Hispanic male HR")

p_asian_female_qd<-ggplot(data=test_asian_female_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_asian_female_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Asian female HR")

p_asian_male_qd<-ggplot(data=test_asian_male_qd, mapping=aes(x=pm25))+
  geom_line(aes(y=test_asian_male_qd))+
  geom_ribbon(aes(ymin=(lower), ymax=(higher)), alpha=0.1)+
  labs(x="PM2.5", y="Asian male HR")

require(cowplot)
plot_grid(p_white_female_qd+theme_bw(), p_white_male_qd + theme_bw(), p_black_female_qd + theme_bw(), p_black_male_qd + theme_bw(),
          p_hispanic_female_qd + theme_bw(), p_hispanic_male_qd + theme_bw(), p_asian_female_qd + theme_bw(), 
          p_asian_male_qd + theme_bw(),
          ncol=2, labels = "AUTO")
