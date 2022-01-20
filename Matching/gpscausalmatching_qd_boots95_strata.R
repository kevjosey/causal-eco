library(devtools)
require(doParallel)
library(data.table)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)
#try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
#install_github("fasrc/CausalGPS", ref="master", force = TRUE)
require(CausalGPS)
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

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingqd/boots5/'

load(paste0(dir_data,"aggregate_data_qd.RData"))
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)

# Matching on single exposure level a, a.vals selects the caliper
aggregate_data_qd$year<-as.factor(aggregate_data_qd$year)
aggregate_data_qd$region<-as.factor(aggregate_data_qd$region)

white_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==2)
white_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==1 & aggregate_data_qd$sex==1)
black_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==2)
black_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==2 & aggregate_data_qd$sex==1)
hispanic_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==2)
hispanic_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==5 & aggregate_data_qd$sex==1)
asian_female_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==2)
asian_male_qd<-aggregate_data_qd %>% filter(aggregate_data_qd$race==4 & aggregate_data_qd$sex==1)

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_female_qd.RData")
covariates_white_female_qd$year<-as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region<-as.factor(covariates_white_female_qd$region)

#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_white_male_qd.RData")
covariates_white_male_qd$year<-as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region<-as.factor(covariates_white_male_qd$region)

#black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_female_qd.RData")
covariates_black_female_qd$year<-as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region<-as.factor(covariates_black_female_qd$region)

#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_black_male_qd.RData")
covariates_black_male_qd$year<-as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region<-as.factor(covariates_black_male_qd$region)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_female_qd.RData")
covariates_hispanic_female_qd$year<-as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region<-as.factor(covariates_hispanic_female_qd$region)

#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_hispanic_male_qd.RData")
covariates_hispanic_male_qd$year<-as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region<-as.factor(covariates_hispanic_male_qd$region)

#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_female_qd.RData")
covariates_asian_female_qd$year<-as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region<-as.factor(covariates_asian_female_qd$region)

#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_asian_male_qd.RData")
covariates_asian_male_qd$year<-as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region<-as.factor(covariates_asian_male_qd$region)

boots_id <- c(0:499)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]

#for(boots_id in 1:500){
set.seed(boots_id)
#All
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_qd, list(covariates_qd$zip))
num_uniq_zip <- length(unique(covariates_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))

match_pop_all <- generate_pseudo_pop(Y=cov$zip,
                                     w=cov$pm25_ensemble,
                                     c=cov[, c(4:19)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transfoqd = TRUE,
                                     transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                               pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))

match_pop_data <- cbind(match_pop_data, cov_trim[,1:2])
aggregate_data_qd2 <- merge(aggregate_data_qd, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)

#Strata
#White female
a.vals <- seq(min(white_female_qd$pm25_ensemble), max(white_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_white_female_qd, list(covariates_white_female_qd$zip))
num_uniq_zip <- length(unique(covariates_white_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))


match_pop_white_female_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25_ensemble,
                                                  c=cov[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transfoqd = TRUE,
                                                  transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                                          pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                                            pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_white_female_qd <- cbind(match_pop_white_female_qd, cov_trim[,1:2])
match_pop_white_female_qd2 <- merge(white_female_qd, match_pop_white_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_female_qd)

#White male
a.vals <- seq(min(white_male_qd$pm25_ensemble), max(white_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_white_male_qd, list(covariates_white_male_qd$zip))
num_uniq_zip <- length(unique(covariates_white_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))


match_pop_white_male_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25_ensemble,
                                                  c=cov[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transfoqd = TRUE,
                                                  transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_white_male_qd <- cbind(match_pop_white_male_qd, cov_trim[,1:2])
match_pop_white_male_qd2 <- merge(white_male_qd, match_pop_white_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_male_qd)

#Black female
a.vals <- seq(min(black_female_qd$pm25_ensemble), max(black_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_black_female_qd, list(covariates_black_female_qd$zip))
num_uniq_zip <- length(unique(covariates_black_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))

match_pop_black_female_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25_ensemble,
                                                  c=cov[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transfoqd = TRUE,
                                                  transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_black_female_qd <- cbind(match_pop_black_female_qd, cov_trim[,1:2])
match_pop_black_female_qd2 <- merge(black_female_qd, match_pop_black_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_female_qd)

#Black male
a.vals <- seq(min(black_male_qd$pm25_ensemble), max(black_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_black_male_qd, list(covariates_black_male_qd$zip))
num_uniq_zip <- length(unique(covariates_black_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))


match_pop_black_male_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                w=cov$pm25_ensemble,
                                                c=cov[, c(4:19)],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transfoqd = TRUE,
                                                transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_black_male_qd <- cbind(match_pop_black_male_qd, cov_trim[,1:2])
match_pop_black_male_qd2 <- merge(black_male_qd, match_pop_black_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_male_qd)

#Hispanic female
a.vals <- seq(min(hispanic_female_qd$pm25_ensemble), max(hispanic_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_hispanic_female_qd, list(covariates_hispanic_female_qd$zip))
num_uniq_zip <- length(unique(covariates_hispanic_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))

match_pop_hispanic_female_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25_ensemble,
                                                  c=cov[, c(4:19)],
                                                  ci_appr = "matching",
                                                  pred_model = "sl",
                                                  gps_model = "parametric",
                                                  use_cov_transfoqd = TRUE,
                                                  transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_hispanic_female_qd <- cbind(match_pop_hispanic_female_qd, cov_trim[,1:2])
match_pop_hispanic_female_qd2 <- merge(hispanic_female_qd, match_pop_hispanic_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_female_qd)

#Hispanic male
a.vals <- seq(min(hispanic_male_qd$pm25_ensemble), max(hispanic_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_hispanic_male_qd, list(covariates_hispanic_male_qd$zip))
num_uniq_zip <- length(unique(covariates_hispanic_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))


match_pop_hispanic_male_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                w=cov$pm25_ensemble,
                                                c=cov[, c(4:19)],
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transfoqd = TRUE,
                                                transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_hispanic_male_qd <- cbind(match_pop_hispanic_male_qd, cov_trim[,1:2])
match_pop_hispanic_male_qd2 <- merge(hispanic_male_qd, match_pop_hispanic_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_male_qd)

#Asian female
a.vals <- seq(min(asian_female_qd$pm25_ensemble), max(asian_female_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_asian_female_qd, list(covariates_asian_female_qd$zip))
num_uniq_zip <- length(unique(covariates_asian_female_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))

match_pop_asian_female_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                     w=cov$pm25_ensemble,
                                                     c=cov[, c(4:19)],
                                                     ci_appr = "matching",
                                                     pred_model = "sl",
                                                     gps_model = "parametric",
                                                     use_cov_transfoqd = TRUE,
                                                     transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_asian_female_qd <- cbind(match_pop_asian_female_qd, cov_trim[,1:2])
match_pop_asian_female_qd2 <- merge(asian_female_qd, match_pop_asian_female_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_female_qd)

#Asian male
a.vals <- seq(min(asian_male_qd$pm25_ensemble), max(asian_male_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_qd.list<-split(covariates_asian_male_qd, list(covariates_asian_male_qd$zip))
num_uniq_zip <- length(unique(covariates_asian_male_qd$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_qd.list[zip_sample]))


match_pop_asian_male_qd1 <- generate_pseudo_pop(Y=cov$zip,
                                                   w=cov$pm25_ensemble,
                                                   c=cov[, c(4:19)],
                                                   ci_appr = "matching",
                                                   pred_model = "sl",
                                                   gps_model = "parametric",
                                                   use_cov_transfoqd = TRUE,
                                                   transfoqders = list("pow2", "pow3"),
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
cov_trim <- subset(cov,
                   pm25_ensemble <= quantile(cov$pm25_ensemble,0.95)&
                     pm25_ensemble >= quantile(cov$pm25_ensemble,0.05))
match_pop_asian_male_qd <- cbind(match_pop_asian_male_qd, cov_trim[,1:2])
match_pop_asian_male_qd2 <- merge(asian_male_qd, match_pop_asian_male_qd[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_male_qd)

save(aggregate_data_qd2, match_pop_white_female_qd2, match_pop_white_male_qd2,
     match_pop_black_female_qd2, match_pop_black_male_qd2,
     match_pop_hispanic_female_qd2, match_pop_hispanic_male_qd2,
     match_pop_asian_female_qd2, match_pop_asian_male_qd2,
     file=paste0('/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingqd/boots5strata/gpscausal_qd_strata_trim5_', boots_id, '.RData'))


#}

