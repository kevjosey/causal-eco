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

load(paste0(dir_data,"aggregate_data_rm.RData"))
# Matching on single exposure level a, a.vals selects the caliper
a.vals <- seq(min(aggregate_data_rm$pm25), max(aggregate_data_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])
aggregate_data_rm$year<-as.factor(aggregate_data_rm$year)
aggregate_data_rm$region<-as.factor(aggregate_data_rm$region)

#All
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_rm.RData")
covariates_rm$year<-as.factor(covariates_rm$year)
covariates_rm$region<-as.factor(covariates_rm$region)
a.vals <- seq(min(covariates_rm$pm25), max(covariates_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])


match_pop_all <- generate_pseudo_pop(Y=covariates_rm$zip,
                                     w=covariates_rm$pm25,
                                     c=covariates_rm[, c(4:19)],
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
covariates_rm_trim <- subset(covariates_rm,
                             pm25 < quantile(covariates_rm$pm25,0.95)&
                               pm25 > quantile(covariates_rm$pm25,0.05))

match_pop_data <- cbind(match_pop_data, covariates_rm_trim[,1:2])
aggregate_data_rm2 <- merge(aggregate_data_rm, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_white_female_rm.RData")
covariates_white_female_rm$year<-as.factor(covariates_white_female_rm$year)
covariates_white_female_rm$region<-as.factor(covariates_white_female_rm$region)

white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
a.vals <- seq(min(white_female_rm$pm25), max(white_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_white_female_rm1 <- generate_pseudo_pop(Y=covariates_white_female_rm$zip,
                                                  w=covariates_white_female_rm$pm25,
                                                  c=covariates_white_female_rm[, c(4:19)],
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
match_pop_white_female_rm <- match_pop_white_female_rm1$pseudo_pop
covariates_white_female_rm_trim <- subset(covariates_white_female_rm,
                                          pm25 < quantile(covariates_white_female_rm$pm25,0.95)&
                                            pm25 > quantile(covariates_white_female_rm$pm25,0.05))
match_pop_white_female_rm <- cbind(match_pop_white_female_rm, covariates_white_female_rm_trim[,1:2])
match_pop_white_female_rm2 <- merge(white_female_rm, match_pop_white_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_female_rm)

#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_white_male_rm.RData")
covariates_white_male_rm$year<-as.factor(covariates_white_male_rm$year)
covariates_white_male_rm$region<-as.factor(covariates_white_male_rm$region)

white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
a.vals <- seq(min(white_male_rm$pm25), max(white_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_white_male_rm1 <- generate_pseudo_pop(Y=covariates_white_male_rm$zip,
                                                w=covariates_white_male_rm$pm25,
                                                c=covariates_white_male_rm[, c(4:19)],
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
match_pop_white_male_rm <- match_pop_white_male_rm1$pseudo_pop
covariates_white_male_rm_trim <- subset(covariates_white_male_rm,
                                        pm25 < quantile(covariates_white_male_rm$pm25,0.95)&
                                          pm25 > quantile(covariates_white_male_rm$pm25,0.05))
match_pop_white_male_rm <- cbind(match_pop_white_male_rm, covariates_white_male_rm_trim[,1:2])
match_pop_white_male_rm2 <- merge(white_male_rm, match_pop_white_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_male_rm)

#black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_black_female_rm.RData")
covariates_black_female_rm$year<-as.factor(covariates_black_female_rm$year)
covariates_black_female_rm$region<-as.factor(covariates_black_female_rm$region)

black_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
a.vals <- seq(min(black_female_rm$pm25), max(black_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_black_female_rm1 <- generate_pseudo_pop(Y=covariates_black_female_rm$zip,
                                                  w=covariates_black_female_rm$pm25,
                                                  c=covariates_black_female_rm[, c(4:19)],
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
match_pop_black_female_rm <- match_pop_black_female_rm1$pseudo_pop
covariates_black_female_rm_trim <- subset(covariates_black_female_rm,
                                          pm25 < quantile(covariates_black_female_rm$pm25,0.95)&
                                            pm25 > quantile(covariates_black_female_rm$pm25,0.05))
match_pop_black_female_rm <- cbind(match_pop_black_female_rm, covariates_black_female_rm_trim[,1:2])
match_pop_black_female_rm2 <- merge(black_female_rm, match_pop_black_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_female_rm)

#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_black_male_rm.RData")
covariates_black_male_rm$year<-as.factor(covariates_black_male_rm$year)
covariates_black_male_rm$region<-as.factor(covariates_black_male_rm$region)

black_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
a.vals <- seq(min(black_male_rm$pm25), max(black_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])


match_pop_black_male_rm1 <- generate_pseudo_pop(Y=covariates_black_male_rm$zip,
                                                w=covariates_black_male_rm$pm25,
                                                c=covariates_black_male_rm[, c(4:19)],
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
match_pop_black_male_rm <- match_pop_black_male_rm1$pseudo_pop
covariates_black_male_rm_trim <- subset(covariates_black_male_rm,
                                        pm25 < quantile(covariates_black_male_rm$pm25,0.95)&
                                          pm25 > quantile(covariates_black_male_rm$pm25,0.05))
match_pop_black_male_rm <- cbind(match_pop_black_male_rm, covariates_black_male_rm_trim[,1:2])
match_pop_black_male_rm2 <- merge(black_male_rm, match_pop_black_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_male_rm)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_hispanic_female_rm.RData")
covariates_hispanic_female_rm$year<-as.factor(covariates_hispanic_female_rm$year)
covariates_hispanic_female_rm$region<-as.factor(covariates_hispanic_female_rm$region)

hispanic_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
a.vals <- seq(min(hispanic_female_rm$pm25), max(hispanic_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_hispanic_female_rm1 <- generate_pseudo_pop(Y=covariates_hispanic_female_rm$zip,
                                                     w=covariates_hispanic_female_rm$pm25,
                                                     c=covariates_hispanic_female_rm[, c(4:19)],
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
match_pop_hispanic_female_rm <- match_pop_hispanic_female_rm1$pseudo_pop
covariates_hispanic_female_rm_trim <- subset(covariates_hispanic_female_rm,
                                             pm25 < quantile(covariates_hispanic_female_rm$pm25,0.95)&
                                               pm25 > quantile(covariates_hispanic_female_rm$pm25,0.05))
match_pop_hispanic_female_rm <- cbind(match_pop_hispanic_female_rm, covariates_hispanic_female_rm_trim[,1:2])
match_pop_hispanic_female_rm2 <- merge(hispanic_female_rm, match_pop_hispanic_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_female_rm)

#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_hispanic_male_rm.RData")
covariates_hispanic_male_rm$year<-as.factor(covariates_hispanic_male_rm$year)
covariates_hispanic_male_rm$region<-as.factor(covariates_hispanic_male_rm$region)

hispanic_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
a.vals <- seq(min(hispanic_male_rm$pm25), max(hispanic_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_hispanic_male_rm1 <- generate_pseudo_pop(Y=covariates_hispanic_male_rm$zip,
                                                   w=covariates_hispanic_male_rm$pm25,
                                                   c=covariates_hispanic_male_rm[, c(4:19)],
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
match_pop_hispanic_male_rm <- match_pop_hispanic_male_rm1$pseudo_pop
covariates_hispanic_male_rm_trim <- subset(covariates_hispanic_male_rm,
                                           pm25 < quantile(covariates_hispanic_male_rm$pm25,0.95)&
                                             pm25 > quantile(covariates_hispanic_male_rm$pm25,0.05))
match_pop_hispanic_male_rm <- cbind(match_pop_hispanic_male_rm, covariates_hispanic_male_rm_trim[,1:2])
match_pop_hispanic_male_rm2 <- merge(hispanic_male_rm, match_pop_hispanic_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_male_rm)

#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_asian_female_rm.RData")
covariates_asian_female_rm$year<-as.factor(covariates_asian_female_rm$year)
covariates_asian_female_rm$region<-as.factor(covariates_asian_female_rm$region)

asian_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
a.vals <- seq(min(asian_female_rm$pm25), max(asian_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_asian_female_rm1 <- generate_pseudo_pop(Y=covariates_asian_female_rm$zip,
                                                  w=covariates_asian_female_rm$pm25,
                                                  c=covariates_asian_female_rm[, c(4:19)],
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
match_pop_asian_female_rm <- match_pop_asian_female_rm1$pseudo_pop
covariates_asian_female_rm_trim <- subset(covariates_asian_female_rm,
                                          pm25 < quantile(covariates_asian_female_rm$pm25,0.95)&
                                            pm25 > quantile(covariates_asian_female_rm$pm25,0.05))
match_pop_asian_female_rm <- cbind(match_pop_asian_female_rm, covariates_asian_female_rm_trim[,1:2])
match_pop_asian_female_rm2 <- merge(asian_female_rm, match_pop_asian_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_female_rm)

#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_asian_male_rm.RData")
covariates_asian_male_rm$year<-as.factor(covariates_asian_male_rm$year)
covariates_asian_male_rm$region<-as.factor(covariates_asian_male_rm$region)

asian_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)
a.vals <- seq(min(asian_male_rm$pm25), max(asian_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_asian_male_rm1 <- generate_pseudo_pop(Y=covariates_asian_male_rm$zip,
                                                w=covariates_asian_male_rm$pm25,
                                                c=covariates_asian_male_rm[, c(4:19)],
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
match_pop_asian_male_rm <- match_pop_asian_male_rm1$pseudo_pop
covariates_asian_male_rm_trim <- subset(covariates_asian_male_rm,
                                        pm25 < quantile(covariates_asian_male_rm$pm25,0.95)&
                                          pm25 > quantile(covariates_asian_male_rm$pm25,0.05))
match_pop_asian_male_rm <- cbind(match_pop_asian_male_rm, covariates_asian_male_rm_trim[,1:2])
match_pop_asian_male_rm2 <- merge(asian_male_rm, match_pop_asian_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_male_rm)

save(match_pop_all,match_pop_data, aggregate_data_rm2,
     match_pop_white_female_rm1, match_pop_white_female_rm, match_pop_white_female_rm2,
     match_pop_white_male_rm1, match_pop_white_male_rm, match_pop_white_male_rm2,
     match_pop_black_female_rm1, match_pop_black_female_rm, match_pop_black_female_rm2,
     match_pop_black_male_rm1, match_pop_black_male_rm, match_pop_black_male_rm2,
     match_pop_hispanic_female_rm1, match_pop_hispanic_female_rm, match_pop_hispanic_female_rm2,
     match_pop_hispanic_male_rm1, match_pop_hispanic_male_rm, match_pop_hispanic_male_rm2,
     match_pop_asian_female_rm1, match_pop_asian_female_rm, match_pop_asian_female_rm2,
     match_pop_asian_male_rm1, match_pop_asian_male_rm, match_pop_asian_male_rm2,
     file="/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/match_pseudo_pop_rm_strata.RData")

#Associations
#Calculating Associations
#All
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/match_pseudo_pop_rm_strata.RData")
main_gnm<-c()
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(sex)+as.factor(race)+as.factor(dual)+
                                 as.factor(entry_age_break)+as.factor(followup_year)),
                            data=aggregate_data_rm2, family=poisson(link="log"), weights=counter))
main_gnm[1]<-exp(10*matchingrm_gnm$coefficients[2])


#White female
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_white_female_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[2]<-exp(10*matchingrm_gnm$coefficients[2])


#White male
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)) 
                            +(as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_white_male_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[3]<-exp(10*matchingrm_gnm$coefficients[2])

#Black female
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_black_female_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[4]<-exp(10*matchingrm_gnm$coefficients[2])

#Black male
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_black_male_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[5]<-exp(10*matchingrm_gnm$coefficients[2])

#Hispanic female
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_hispanic_female_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[6]<-exp(10*matchingrm_gnm$coefficients[2])

#Hispanic male
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_hispanic_male_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[7]<-exp(10*matchingrm_gnm$coefficients[2])

#Asian female
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_asian_female_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[8]<-exp(10*matchingrm_gnm$coefficients[2])

#Asian male
matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                              (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                            data=match_pop_asian_male_rm2, family=poisson(link="log"), 
                            weights=counter))
main_gnm[9]<-exp(10*matchingrm_gnm$coefficients[2])

#ERCS

matchingrm_gam <-mgcv::bam(dead~ s(pm25,bs='cr',k=3) +
                             as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                             offset(log(time_count))
                           , data=aggregate_data_rm2,family=poisson(link="log"))

white_femalerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_white_female_rm2,family=poisson(link="log"))

white_malerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_white_male_rm2,family=poisson(link="log"))

black_femalerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_black_female_rm2,family=poisson(link="log"))

black_malerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_black_male_rm2,family=poisson(link="log"))

hispanic_femalerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_hispanic_female_rm2,family=poisson(link="log"))

hispanic_malerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_hispanic_male_rm2,family=poisson(link="log"))

asian_femalerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=match_pop_asian_female_rm2,family=poisson(link="log"))

asian_malerm_matching_gam <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                        as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                        offset(log(time_count))
                                      , data=match_pop_asian_male_rm2,family=poisson(link="log"))
#Plots
##rm
test.data.rm<-function(x, gm){
  test<-data.frame(pm25 = seq(min(x$pm25), max(x$pm25),length.out=50) ,
                   entry_age_break= rep(levels(as.factor(x$entry_age_break))[1], 50),
                   dual = rep(levels(as.factor(x$dual))[1],50),
                   sex = rep(levels(as.factor(x$sex))[1], 50),
                   race = rep(levels(as.factor(x$race))[1], 50),
                   time_count=rep(1, 50),
                   followup_year= rep(levels(as.factor(x$followup_year))[1], 50)
                   
  )
  
  base<-data.frame(pm25 = rep(min(x$pm25),50) ,
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
testallrm<-test.data.rm(aggregate_data_rm2, matchingrm_gam)
test_white_female_rm<-test.data.rm(match_pop_white_female_rm2 , white_femalerm_matching_gam)
test_white_male_rm<-test.data.rm(match_pop_white_male_rm2, white_malerm_matching_gam)
test_black_female_rm<-test.data.rm(match_pop_black_female_rm2, black_femalerm_matching_gam)
test_black_male_rm<-test.data.rm(match_pop_black_male_rm2, black_malerm_matching_gam)
test_hispanic_female_rm<-test.data.rm(match_pop_hispanic_female_rm2, hispanic_femalerm_matching_gam)
test_hispanic_male_rm<-test.data.rm(match_pop_hispanic_male_rm2, hispanic_malerm_matching_gam)
test_asian_female_rm<-test.data.rm(match_pop_asian_female_rm2, asian_femalerm_matching_gam)
test_asian_male_rm<-test.data.rm(match_pop_asian_male_rm2, asian_malerm_matching_gam)

#Bootstrap
files<-list.files("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingrm/boots5strata/", 
                  pattern = "\\.RData",full.names=TRUE)
tallrm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_white_female_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_white_male_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_black_female_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_black_male_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_hispanic_female_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_hispanic_male_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_asian_female_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))
t_asian_male_rm<-as.data.frame(matrix(ncol=50, nrow=length(files)))

t_gnm<-as.data.frame(matrix(nrow=length(files), ncol=9))
colnames(t_gnm)<-c("all", "white_female", "white_male", "black_female", "black_male",
                   "hispanic_female", "hispanic_male", "asian_female", "asian_male")

for(i in 1:length(files)){
  load(files[i])
  print(i)
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(sex)+as.factor(race)+as.factor(dual)+
                                   as.factor(entry_age_break)+as.factor(followup_year)),
                              data=aggregate_data_rm2, family=poisson(link="log"), weights=counter))
  t_gnm$all[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #White female
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_white_female_rm2, family=poisson(link="log"), 
                              weights=counter))
  t_gnm$white_female[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  
  #White male
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)) 
                              +(as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_white_male_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$white_male[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Black female
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_black_female_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$black_female[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Black male
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_black_male_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$black_male[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Hispanic female
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_hispanic_female_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$hispanic_female[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Hispanic male
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_hispanic_male_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$hispanic_male[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Asian female
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_asian_female_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$asian_female[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #Asian male
  matchingrm_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count))+
                                (as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)),
                              data=match_pop_asian_male_rm2, family=poisson(link="log"), 
                              weights=counter))
  exp(10*matchingrm_gnm$coefficients[2])
  t_gnm$asian_male[i]=exp(10*matchingrm_gnm$coefficients[2])
  
  #ERC
  matchingrm_gam1 <-mgcv::bam(dead~ s(pm25,bs='cr',k=3) +
                                as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                offset(log(time_count))
                              , data=aggregate_data_rm2,family=poisson(link="log"))
  
  
  white_femalerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_white_female_rm2,family=poisson(link="log"))
  
  white_malerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_white_male_rm2,family=poisson(link="log"))
  
  black_femalerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_black_female_rm2,family=poisson(link="log"))
  
  black_malerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_black_male_rm2,family=poisson(link="log"))
  
  hispanic_femalerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                                as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                                offset(log(time_count))
                                              , data=match_pop_hispanic_female_rm2,family=poisson(link="log"))
  
  hispanic_malerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                              as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                              offset(log(time_count))
                                            , data=match_pop_hispanic_male_rm2,family=poisson(link="log"))
  
  asian_femalerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr',k=3) +
                                             as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                             offset(log(time_count))
                                           , data=match_pop_asian_female_rm2,family=poisson(link="log"))
  
  asian_malerm_matching_gam1 <-mgcv::bam(dead~ s(pm25, bs='cr', k=3) +
                                           as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                           offset(log(time_count))
                                         , data=match_pop_asian_male_rm2,family=poisson(link="log"))
  
  tallrm[i,]<-test.data.rm(aggregate_data_rm2, matchingrm_gam)
  t_white_female_rm[i,]<-test.data.rm(match_pop_white_female_rm2 , white_femalerm_matching_gam)
  t_white_male_rm[i,]<-test.data.rm(match_pop_white_male_rm2, white_malerm_matching_gam)
  t_black_female_rm[i,]<-test.data.rm(match_pop_black_female_rm2, black_femalerm_matching_gam)
  t_black_male_rm[i,]<-test.data.rm(match_pop_black_male_rm2, black_malerm_matching_gam)
  t_hispanic_female_rm[i,]<-test.data.rm(match_pop_hispanic_female_rm2, hispanic_femalerm_matching_gam)
  t_hispanic_male_rm[i,]<-test.data.rm(match_pop_hispanic_male_rm2, hispanic_malerm_matching_gam)
  t_asian_female_rm[i,]<-test.data.rm(match_pop_asian_female_rm2, asian_femalerm_matching_gam)
  t_asian_male_rm[i,]<-test.data.rm(match_pop_asian_male_rm2, asian_malerm_matching_gam)
}

#GNM
quantile(t_gnm$all, prob=c(0.025, 0.975), na.rm==TRUE)
quantile(t_gnm$white_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$white_male, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$black_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$black_male, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_male, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_female, prob=c(0.025, 0.975), na.rm=TRUE)
quantile(t_gnm$hispanic_male, prob=c(0.025, 0.975), na.rm=TRUE)

#ERC
require(matrixStats)
testallrm<-as.data.frame(testallrm)
testallrm$lower<- apply(tallrm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
testallrm$upper<- apply(tallrm,2, quantile, prob=c(0.975), na.rm=TRUE)
testallrm$pm25<-seq(min(aggregate_data_rm2$pm25), max(aggregate_data_rm2$pm25), length.out = 50)

test_white_female_rm<-as.data.frame(test_white_female_rm)
test_white_female_rm$lower<- apply(t_white_female_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_white_female_rm$upper<- apply(t_white_female_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_white_female_rm$pm25<-seq(min(match_pop_white_female_rm2$pm25), max(match_pop_white_female_rm2$pm25), length.out = 50)

test_white_male_rm<-as.data.frame(test_white_male_rm)
test_white_male_rm$lower<- apply(t_white_male_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_white_male_rm$upper<- apply(t_white_male_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_white_male_rm$pm25<-seq(min(match_pop_white_male_rm2$pm25), max(match_pop_white_male_rm2$pm25), length.out = 50)

test_black_female_rm<-as.data.frame(test_black_female_rm)
test_black_female_rm$lower<- apply(t_black_female_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_black_female_rm$upper<- apply(t_black_female_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_black_female_rm$pm25<-seq(min(match_pop_black_female_rm2$pm25), max(match_pop_black_female_rm2$pm25), length.out = 50)

test_black_male_rm<-as.data.frame(test_black_male_rm)
test_black_male_rm$lower<- apply(t_black_male_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_black_male_rm$upper<- apply(t_black_male_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_black_male_rm$pm25<-seq(min(match_pop_black_male_rm2$pm25), max(match_pop_black_male_rm2$pm25), length.out = 50)

test_hispanic_female_rm<-as.data.frame(test_hispanic_female_rm)
test_hispanic_female_rm$lower<- apply(t_hispanic_female_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_hispanic_female_rm$upper<- apply(t_hispanic_female_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_hispanic_female_rm$pm25<-seq(min(match_pop_hispanic_female_rm2$pm25), max(match_pop_hispanic_female_rm2$pm25), length.out = 50)

test_hispanic_male_rm<-as.data.frame(test_hispanic_male_rm)
test_hispanic_male_rm$lower<- apply(t_hispanic_male_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_hispanic_male_rm$upper<- apply(t_hispanic_male_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_hispanic_male_rm$pm25<-seq(min(match_pop_hispanic_male_rm2$pm25), max(match_pop_hispanic_male_rm2$pm25), length.out = 50)

test_asian_female_rm<-as.data.frame(test_asian_female_rm)
test_asian_female_rm$lower<- apply(t_asian_female_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_asian_female_rm$upper<- apply(t_asian_female_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_asian_female_rm$pm25<-seq(min(match_pop_asian_female_rm2$pm25), max(match_pop_asian_female_rm2$pm25), length.out = 50)

test_asian_male_rm<-as.data.frame(test_asian_male_rm)
test_asian_male_rm$lower<- apply(t_asian_male_rm, 2, quantile,  prob=c(0.025), na.rm=TRUE)
test_asian_male_rm$upper<- apply(t_asian_male_rm,2, quantile, prob=c(0.975), na.rm=TRUE)
test_asian_male_rm$pm25<-seq(min(match_pop_asian_male_rm2$pm25), max(match_pop_asian_male_rm2$pm25), length.out = 50)

pall_rm<-ggplot(data=testallrm, mapping=aes(x=pm25))+geom_line(aes(y=testallrm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="All HR")

p_white_female_rm<-ggplot(data=test_white_female_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_white_female_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="White female HR")

p_white_male_rm<-ggplot(data=test_white_male_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_white_male_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="White male HR")

p_black_female_rm<-ggplot(data=test_black_female_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_black_female_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Black female HR")

p_black_male_rm<-ggplot(data=test_black_male_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_black_male_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Black male HR")

p_hispanic_female_rm<-ggplot(data=test_hispanic_female_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_hispanic_female_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Hispanic female HR")

p_hispanic_male_rm<-ggplot(data=test_hispanic_male_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_hispanic_male_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Hispanic male HR")

p_asian_female_rm<-ggplot(data=test_asian_female_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_asian_female_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Asian female HR")

p_asian_male_rm<-ggplot(data=test_asian_male_rm, mapping=aes(x=pm25))+
  geom_line(aes(y=test_asian_male_rm))+
  geom_ribbon(aes(ymin=(lower), ymax=(upper)), alpha=0.1)+
  labs(x="rm PM2.5", y="Asian male HR")

require(cowplot)
plot_grid(p_white_female_rm, p_white_male_rm, p_black_female_rm, p_black_male_rm,
          p_hispanic_female_rm, p_hispanic_male_rm, p_asian_female_rm, p_asian_male_rm,
          ncol=2, labels="AUTO")
