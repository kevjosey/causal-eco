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
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingrm/boots5strata/'

load(paste0(dir_data,"aggregate_data_rm.RData"))
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_rm.RData")
covariates_rm$year<-as.factor(covariates_rm$year)
covariates_rm$region<-as.factor(covariates_rm$region)

# Matching on single exposure level a, a.vals selects the caliper
aggregate_data_rm$year<-as.factor(aggregate_data_rm$year)
aggregate_data_rm$region<-as.factor(aggregate_data_rm$region)

white_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==2)
white_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==1 & aggregate_data_rm$sex==1)
black_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==2)
black_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==2 & aggregate_data_rm$sex==1)
hispanic_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==2)
hispanic_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==5 & aggregate_data_rm$sex==1)
asian_female_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==2)
asian_male_rm<-aggregate_data_rm %>% filter(aggregate_data_rm$race==4 & aggregate_data_rm$sex==1)

#White female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_white_female_rm.RData")
covariates_white_female_rm$year<-as.factor(covariates_white_female_rm$year)
covariates_white_female_rm$region<-as.factor(covariates_white_female_rm$region)

#White male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_white_male_rm.RData")
covariates_white_male_rm$year<-as.factor(covariates_white_male_rm$year)
covariates_white_male_rm$region<-as.factor(covariates_white_male_rm$region)

#black female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_black_female_rm.RData")
covariates_black_female_rm$year<-as.factor(covariates_black_female_rm$year)
covariates_black_female_rm$region<-as.factor(covariates_black_female_rm$region)

#black male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_black_male_rm.RData")
covariates_black_male_rm$year<-as.factor(covariates_black_male_rm$year)
covariates_black_male_rm$region<-as.factor(covariates_black_male_rm$region)

#hispanic female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_hispanic_female_rm.RData")
covariates_hispanic_female_rm$year<-as.factor(covariates_hispanic_female_rm$year)
covariates_hispanic_female_rm$region<-as.factor(covariates_hispanic_female_rm$region)

#hispanic male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_hispanic_male_rm.RData")
covariates_hispanic_male_rm$year<-as.factor(covariates_hispanic_male_rm$year)
covariates_hispanic_male_rm$region<-as.factor(covariates_hispanic_male_rm$region)

#asian female
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_asian_female_rm.RData")
covariates_asian_female_rm$year<-as.factor(covariates_asian_female_rm$year)
covariates_asian_female_rm$region<-as.factor(covariates_asian_female_rm$region)

#asian male
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_rm/covariates_asian_male_rm.RData")
covariates_asian_male_rm$year<-as.factor(covariates_asian_male_rm$year)
covariates_asian_male_rm$region<-as.factor(covariates_asian_male_rm$region)

boots_id <- c(0:499)[as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1]

#for(boots_id in 1:500){
set.seed(boots_id)
#All
a.vals <- seq(min(covariates_rm$pm25), max(covariates_rm$pm25), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_rm, list(covariates_rm$zip))
num_uniq_zip <- length(unique(covariates_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

match_pop_all <- generate_pseudo_pop(Y=cov$zip,
                                     w=cov$pm25,
                                     c=cov[, c(4:19)],
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
cov_trim <- subset(cov,pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))

match_pop_data <- cbind(match_pop_data, cov_trim[,1:2])
aggregate_data_rm2 <- merge(aggregate_data_rm, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)

#Strata
#White female
a.vals <- seq(min(white_female_rm$pm25), max(white_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_white_female_rm, list(covariates_white_female_rm$zip))
num_uniq_zip <- length(unique(covariates_white_female_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


match_pop_white_female_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25,
                                                  c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_white_female_rm <- cbind(match_pop_white_female_rm, cov_trim[,1:2])
match_pop_white_female_rm2 <- merge(white_female_rm, match_pop_white_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_female_rm)

#White male
a.vals <- seq(min(white_male_rm$pm25), max(white_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_white_male_rm, list(covariates_white_male_rm$zip))
num_uniq_zip <- length(unique(covariates_white_male_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


match_pop_white_male_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                w=cov$pm25,
                                                c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_white_male_rm <- cbind(match_pop_white_male_rm, cov_trim[,1:2])
match_pop_white_male_rm2 <- merge(white_male_rm, match_pop_white_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(white_male_rm)

#Black female
a.vals <- seq(min(black_female_rm$pm25), max(black_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_black_female_rm, list(covariates_black_female_rm$zip))
num_uniq_zip <- length(unique(covariates_black_female_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

match_pop_black_female_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25,
                                                  c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_black_female_rm <- cbind(match_pop_black_female_rm, cov_trim[,1:2])
match_pop_black_female_rm2 <- merge(black_female_rm, match_pop_black_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_female_rm)

#Black male
a.vals <- seq(min(black_male_rm$pm25), max(black_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_black_male_rm, list(covariates_black_male_rm$zip))
num_uniq_zip <- length(unique(covariates_black_male_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


match_pop_black_male_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                w=cov$pm25,
                                                c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_black_male_rm <- cbind(match_pop_black_male_rm, cov_trim[,1:2])
match_pop_black_male_rm2 <- merge(black_male_rm, match_pop_black_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(black_male_rm)

#Hispanic female
a.vals <- seq(min(hispanic_female_rm$pm25), max(hispanic_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_hispanic_female_rm, list(covariates_hispanic_female_rm$zip))
num_uniq_zip <- length(unique(covariates_hispanic_female_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

match_pop_hispanic_female_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                     w=cov$pm25,
                                                     c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_hispanic_female_rm <- cbind(match_pop_hispanic_female_rm, cov_trim[,1:2])
match_pop_hispanic_female_rm2 <- merge(hispanic_female_rm, match_pop_hispanic_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_female_rm)

#Hispanic male
a.vals <- seq(min(hispanic_male_rm$pm25), max(hispanic_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_hispanic_male_rm, list(covariates_hispanic_male_rm$zip))
num_uniq_zip <- length(unique(covariates_hispanic_male_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


match_pop_hispanic_male_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                   w=cov$pm25,
                                                   c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_hispanic_male_rm <- cbind(match_pop_hispanic_male_rm, cov_trim[,1:2])
match_pop_hispanic_male_rm2 <- merge(hispanic_male_rm, match_pop_hispanic_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(hispanic_male_rm)

#Asian female
a.vals <- seq(min(asian_female_rm$pm25), max(asian_female_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_asian_female_rm, list(covariates_asian_female_rm$zip))
num_uniq_zip <- length(unique(covariates_asian_female_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

match_pop_asian_female_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                  w=cov$pm25,
                                                  c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_asian_female_rm <- cbind(match_pop_asian_female_rm, cov_trim[,1:2])
match_pop_asian_female_rm2 <- merge(asian_female_rm, match_pop_asian_female_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_female_rm)

#Asian male
a.vals <- seq(min(asian_male_rm$pm25), max(asian_male_rm$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

all_rm.list<-split(covariates_asian_male_rm, list(covariates_asian_male_rm$zip))
num_uniq_zip <- length(unique(covariates_asian_male_rm$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
cov<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


match_pop_asian_male_rm1 <- generate_pseudo_pop(Y=cov$zip,
                                                w=cov$pm25,
                                                c=cov[, c(4:19)],
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
cov_trim <- subset(cov,
                   pm25 <= quantile(cov$pm25,0.95)&
                     pm25 >= quantile(cov$pm25,0.05))
match_pop_asian_male_rm <- cbind(match_pop_asian_male_rm, cov_trim[,1:2])
match_pop_asian_male_rm2 <- merge(asian_male_rm, match_pop_asian_male_rm[, c("year", "zip", "counter")], by = c("year", "zip"), all.y =TRUE)
rm(asian_male_rm)

save(aggregate_data_rm2, match_pop_white_female_rm2, match_pop_white_male_rm2,
     match_pop_black_female_rm2, match_pop_black_male_rm2,
     match_pop_hispanic_female_rm2, match_pop_hispanic_male_rm2,
     match_pop_asian_female_rm2, match_pop_asian_male_rm2,
     file=paste0('/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingrm/boots5strata/gpscausal_rm_strata_trim5_', boots_id, '.RData'))


#}

