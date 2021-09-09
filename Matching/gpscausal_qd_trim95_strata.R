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
aggregate_data_qd2 <- merge(aggregate_data_qd, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_white_female_qd.RData")
covariates_white_female_qd$year<-as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region<-as.factor(covariates_white_female_qd$region)

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
match_pop_white_female_qd2 <- merge(aggregate_data_qd, match_pop_white_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)

#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_white_male_qd.RData")
covariates_white_male_qd$year<-as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region<-as.factor(covariates_white_male_qd$region)


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
match_pop_white_male_qd2 <- merge(aggregate_data_qd, match_pop_white_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)

#black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_black_female_qd.RData")
covariates_black_female_qd$year<-as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region<-as.factor(covariates_black_female_qd$region)

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
match_pop_black_female_qd2 <- merge(aggregate_data_qd, match_pop_black_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)


#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_black_male_qd.RData")
covariates_black_male_qd$year<-as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region<-as.factor(covariates_black_male_qd$region)


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
match_pop_black_male_qd2 <- merge(aggregate_data_qd, match_pop_black_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_hispanic_female_qd.RData")
covariates_hispanic_female_qd$year<-as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region<-as.factor(covariates_hispanic_female_qd$region)

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
match_pop_hispanic_female_qd2 <- merge(aggregate_data_qd, match_pop_hispanic_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)


#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_hispanic_male_qd.RData")
covariates_hispanic_male_qd$year<-as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region<-as.factor(covariates_hispanic_male_qd$region)


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
match_pop_hispanic_male_qd2 <- merge(aggregate_data_qd, match_pop_hispanic_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)


#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_asian_female_qd.RData")
covariates_asian_female_qd$year<-as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region<-as.factor(covariates_asian_female_qd$region)

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
match_pop_asian_female_qd2 <- merge(aggregate_data_qd, match_pop_asian_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)


#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/balance_qd/covariates_asian_male_qd.RData")
covariates_asian_male_qd$year<-as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region<-as.factor(covariates_asian_male_qd$region)


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
match_pop_asian_male_qd2 <- merge(aggregate_data_qd, match_pop_asian_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.x =TRUE)

save(match_pop_all,match_pop_data, aggregate_data_qd2,
     match_pop_white_female_qd1, match_pop_white_female_qd, match_pop_white_female_qd2,
     match_pop_white_male_qd1, match_pop_white_male_qd, match_pop_white_male_qd2,
     match_pop_black_female_qd1, match_pop_black_female_qd, match_pop_black_female_qd2,
     match_pop_black_male, match_pop_black_male_qd, match_pop_black_male_qd2,
     match_pop_hispanic_female, match_pop_hispanic_female_qd, match_pop_hispanic_female_qd2,
     match_pop_hispanic_male, match_pop_hispanic_male_qd, match_pop_hispanic_male_qd2,
     match_pop_asian_female, match_pop_asian_female_qd, match_pop_asian_female_qd2,
     match_pop_asian_male, match_pop_asian_male_qd, match_pop_asian_male_qd2,
     file="/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerrror/Data/match_pseudo_pop_qd_strata.RData")

#Associations
#Calculating Associations
#All
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=aggregate_data_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#White female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_white_female_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#White male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_white_male_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Black female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_black_female_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Black male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_black_male_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_hispanic_female_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_hispanic_male_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Asian female
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_asian_female_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#Asian male
matchingqd_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)),
                            eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                            data=match_pop_asian_male_qd2, family=poisson(link="log"), weights=counter))
exp(10*matchingqd_gnm$coefficients[1])

#ERCS

matchingqd_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                             as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                             offset(log(time_count))
                           , data=aggregate_data_qd2,family=poisson(link="log"))

white_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_white_female_qd2,family=poisson(link="log"))

white_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_white_male_qd2,family=poisson(link="log"))

black_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_black_female_qd2,family=poisson(link="log"))

black_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_black_male_qd2,family=poisson(link="log"))

hispanic_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_hispanic_female_qd2,family=poisson(link="log"))

hispanic_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_hispanic_male_qd2,family=poisson(link="log"))
asian_femaleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_asian_female_qd2,family=poisson(link="log"))

asian_maleqd_matching_gam <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_asian_male_qd2,family=poisson(link="log"))

#Plots
##QD
test.data.qd<-function(x, gm){
  test<-data.frame(pm25_ensemble = seq(min(x$pm25_ensemble), max(x$pm25_ensemble),
                                       length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50),
                   time_count= rep(1, 50)) # This will give the appropriate RR estimate we desire
 
  pred.vals<-predict(gm, newdata=test, type="link", se.fit = TRUE)
 
  dly<- pred.vals$fit
  return(gm$family$linkinv(dly))
}



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

pall_qd<-ggplot(data=testallqd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="All HR")
palll_qd1<-ggMarginal(pall_qd, type="histogram")
p_white_female_qd<-ggplot(data=test_white_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White female HR")
p_white_female_qd1<-ggMarginal(p_white_female_qd, type="histogram")
p_white_male_qd<-ggplot(data=test_white_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="White male HR")
p_white_male_qd1<-ggMarginal(p_white_male_qd, type="histogram")
p_black_female_qd<-ggplot(data=test_black_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black female HR")
p_black_female_qd1<-ggMarginal(p_black_female_qd, type="histogram")
p_black_male_qd<-ggplot(data=test_black_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Black male HR")
p_black_male_qd1<-ggMarginal(p_black_male_qd, type="histogram")
p_hispanic_female_qd<-ggplot(data=test_hispanic_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic female HR")
p_hispanic_female_qd1<-ggMarginal(p_hispanic_female_qd, type="histogram")
p_hispanic_male_qd<-ggplot(data=test_hispanic_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Hispanic male HR")
p_hispanic_male_qd1<-ggMarginal(p_hispanic_male_qd, type="histogram")
p_asian_female_qd<-ggplot(data=test_asian_female_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian female HR")
p_asian_female_qd1<-ggMarginal(p_asian_female_qd, type="histogram")
p_asian_male_qd<-ggplot(data=test_asian_male_qd, mapping=aes(x=pm25_ensemble))+geom_point(aes(y=exp(y)))+
  geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), color="red")+
  labs(x="QD PM2.5", y="Asian male HR")
p_asian_male_qd1<-ggMarginal(p_asian_male_qd, type="histogram")

p_hist_all<-ggplot(data=matching_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="All QD PM2.5", y="Freq")
p_white_female_qd2<-ggplot(data=white_female_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="White female QD PM2.5", y="Freq")
p_white_male_qd2<-ggplot(data=white_male_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="White male QD PM2.5", y="Freq")
p_black_female_qd2<-ggplot(data=black_female_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Black female QD PM2.5", y="Freq")
p_black_male_qd2<-ggplot(data=black_male_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Black male QD PM2.5", y="Freq")
p_hispanic_female_qd2<-ggplot(data=hispanic_female_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Hispanic female QD PM2.5", y="Freq")
p_hispanic_male_qd2<-ggplot(data=hispanic_male_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Hispanic male QD PM2.5", y="Freq")
p_asian_female_qd2<-ggplot(data=asian_female_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Asian female QD PM2.5", y="Freq")
p_asian_male_qd2<-ggplot(data=asian_male_qd3, aes(x=pm25_ensemble))+geom_histogram()+
  labs(x="Asian male QD PM2.5", y="Freq")

plot_grid(p_white_female_qd, p_white_male_qd,
          p_black_female_qd, p_black_male_qd,
          p_hispanic_female_qd, p_hispanic_male_qd,
          p_asian_female_qd, p_asian_male_qd, ncol=2)
