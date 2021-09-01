library(devtools)
require(doParallel)
library(data.table)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)
try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
install_github("fasrc/CausalGPS", ref="develop")
#install_github("fasrc/CausalGPS", ref="937810da2350f5c937c11eab16c60f2ee9a2783f")
library("CausalGPS")
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

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

load(paste0(dir_data,"aggregate_data_rm.RData"))
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_rm.RData")
covariates_rm$year<-as.factor(covariates_rm$year)
covariates_rm$region<-as.factor(covariates_rm$region)
a.vals <- seq(min(covariates_rm$pm25), max(covariates_rm$pm25), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

# Matching on single exposure level a, a.vals selects the caliper
#a.vals <- seq(min(aggregate_data_rm$pm25), max(aggregate_data_rm$pm25), length.out = 50)
#delta_n <- (a.vals[2] - a.vals[1])
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
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_white_female_rm.RData")
covariates_white_female_rm$year<-as.factor(covariates_white_female_rm$year)
covariates_white_female_rm$region<-as.factor(covariates_white_female_rm$region)

#White male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_white_male_rm.RData")
covariates_white_male_rm$year<-as.factor(covariates_white_male_rm$year)
covariates_white_male_rm$region<-as.factor(covariates_white_male_rm$region)

#black female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_black_female_rm.RData")
covariates_black_female_rm$year<-as.factor(covariates_black_female_rm$year)
covariates_black_female_rm$region<-as.factor(covariates_black_female_rm$region)

#black male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_black_male_rm.RData")
covariates_black_male_rm$year<-as.factor(covariates_black_male_rm$year)
covariates_black_male_rm$region<-as.factor(covariates_black_male_rm$region)

#hispanic female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_hispanic_female_rm.RData")
covariates_hispanic_female_rm$year<-as.factor(covariates_hispanic_female_rm$year)
covariates_hispanic_female_rm$region<-as.factor(covariates_hispanic_female_rm$region)

#hispanic male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_hispanic_male_rm.RData")
covariates_hispanic_male_rm$year<-as.factor(covariates_hispanic_male_rm$year)
covariates_hispanic_male_rm$region<-as.factor(covariates_hispanic_male_rm$region)

#asian female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_asian_female_rm.RData")
covariates_asian_female_rm$year<-as.factor(covariates_asian_female_rm$year)
covariates_asian_female_rm$region<-as.factor(covariates_asian_female_rm$region)

#asian male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_rm/covariates_asian_male_rm.RData")
covariates_asian_male_rm$year<-as.factor(covariates_asian_male_rm$year)
covariates_asian_male_rm$region<-as.factor(covariates_asian_male_rm$region)

prematchdata <-function(x, y){
  dead_personyear<-aggregate(cbind(x$dead,
                                   x$time_count),
                             by=list(x$zip,
                                     x$year),
                             FUN=sum)
  colnames(dead_personyear)[1:4]<-c("zip","year","dead","time_count")
  dead_personyear[,"mortality"] <- dead_personyear[,"dead"]/dead_personyear[,"time_count"]
  
  prematch_data <- merge(dead_personyear, y,
                         by=c("zip", "year"))
  return(prematch_data)
}


for(boots_id in 1:500){
set.seed(boots_id)
#All
prematch_data1<-prematchdata(aggregate_data_rm, covariates_rm)

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]



c2=as.data.frame(data.matrix(c))

match_pop_all_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                       w=treat,
                                                       c=c2,
                                                       ci_appr = "matching",
                                                       pred_model = "sl",
                                                       gps_model = "parametric",
                                                       use_cov_transform = FALSE,
                                                       #transformers = list("pow2", "pow3"),
                                                       sl_lib = c("m_xgboost"),
                                                       params = list("xgb_nrounds"=50,
                                                                     "xgb_max_depth"=6,
                                                                     "xgb_eta"=0.3,
                                                                     "xgb_min_child_weight"=1),
                                                       nthread=16, # number of cores, you can change,
                                                       covar_bl_method = "absolute",
                                                       covar_bl_trs = 0.5,
                                                       trim_quantiles = c(0,1), # trimed, you can change
                                                       optimized_compile = FALSE, #created a column counter for how many times matched,
                                                       max_attempt = 1,
                                                       matching_fun = "matching_l1",
                                                       delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                       scale = 1.0)

match_pop_data_notrim<-match_pop_all_noncompile_notrim$pseudo_pop

erf_notrim_all<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                      matched_w=match_pop_data_notrim$w,
                                      bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                      w_vals=a.vals, nthread=16)


#trimed ERC
#plot(erf_notrim_all$erf[trimed_index])
#plot(a.vals[trimed_index], 
 #    erf_notrim_all$erf[trimed_index]/erf_notrim_all$erf[trimed_index][1])

#trimed matched dataset
match_pop_all_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                     w=treat,
                                                     c=c2,
                                                     ci_appr = "matching",
                                                     pred_model = "sl",
                                                     gps_model = "parametric",
                                                     use_cov_transform = FALSE,
                                                     #transformers = list("pow2", "pow3"),
                                                     sl_lib = c("m_xgboost"),
                                                     params = list("xgb_nrounds"=50,
                                                                   "xgb_max_depth"=6,
                                                                   "xgb_eta"=0.3,
                                                                   "xgb_min_child_weight"=1),
                                                     nthread=16, # number of cores, you can change,
                                                     covar_bl_method = "absolute",
                                                     covar_bl_trs = 0.1,
                                                     trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                     optimized_compile = FALSE, #created a column counter for how many times matched,
                                                     max_attempt = 1,
                                                     matching_fun = "matching_l1",
                                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                     scale = 1.0)

match_pop_data_trim<-match_pop_all_noncompile_trim$pseudo_pop

erf_trim_all<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                    matched_w = match_pop_data_trim$w,
                                    bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                    w_vals=a.vals,
                                    nthread=16)
#plot(erf_trim_all)
#plot(erf_trim_all$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_all$erf[trimed_index]/erf_trim_all$erf[trimed_index][1])

match_pop_all_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                       w = treat,
                                                       c = c,
                                                       ci_appr = "matching",
                                                       pred_model = "sl",
                                                       gps_model = "parametric",
                                                       use_cov_transform = FALSE,
                                                       sl_lib = c("m_xgboost"),
                                                       params = list("xgb_nrounds" = 50,
                                                                     "xgb_max_depth" = 6,
                                                                     "xgb_eta" = 0.3,
                                                                     "xgb_min_child_weight" = 1),
                                                       nthread=16, # number of cores, you can change,
                                                       covar_bl_method = "absolute",
                                                       covar_bl_trs = 0.1,
                                                       trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                       optimized_compile = FALSE, #created a column counter for how many times matched,
                                                       max_attempt = 1,
                                                       matching_fun = "matching_l1",
                                                       delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                       scale = 1.0)
match_pop_data_trim_onehot <- match_pop_all_noncompile_onehot$pseudo_pop

erf_trim_onehot_all <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                             matched_w = match_pop_data_trim_onehot$w,
                                             bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals = a.vals,
                                             nthread=16)



#Strata
#White female
prematch_data1<-prematchdata(white_female_rm, covariates_white_female_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_white_female_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                                w=treat,
                                                                c=c2,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                #transformers = list("pow2", "pow3"),
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds"=50,
                                                                              "xgb_max_depth"=6,
                                                                              "xgb_eta"=0.3,
                                                                              "xgb_min_child_weight"=1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.5,
                                                                trim_quantiles = c(0,1), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)

match_pop_data_notrim<-match_pop_white_female_noncompile_notrim$pseudo_pop

erf_notrim_white_female<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                               matched_w=match_pop_data_notrim$w,
                                               bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                               w_vals=a.vals, nthread=16)
#plot(erf_notrim_white_female)

#trimed ERC
#plot(erf_notrim_white_female$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_white_female$erf[trimed_index]/erf_notrim_white_female$erf[trimed_index][1])

#trimed matched dataset
match_pop_white_female_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_trim<-match_pop_white_female_noncompile_trim$pseudo_pop

erf_trim_white_female<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                             matched_w = match_pop_data_trim$w,
                                             bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals,
                                             nthread=16)
#plot(erf_trim_white_female)
#plot(erf_trim_white_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_white_female$erf[trimed_index]/erf_trim_white_female$erf[trimed_index][1])

match_pop_white_female_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                                w = treat,
                                                                c = c,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds" = 50,
                                                                              "xgb_max_depth" = 6,
                                                                              "xgb_eta" = 0.3,
                                                                              "xgb_min_child_weight" = 1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)
match_pop_data_trim_onehot <- match_pop_white_female_noncompile_onehot$pseudo_pop

erf_trim_onehot_white_female <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                      matched_w = match_pop_data_trim_onehot$w,
                                                      bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                      w_vals = a.vals,
                                                      nthread=16)
#plot(erf_trim_onehot_white_female)

#plot(erf_trim_onehot_white_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_white_female$erf[trimed_index]/erf_trim_onehot_white_female$erf[trimed_index][1])


#White male
prematch_data1<-prematchdata(white_male_rm, covariates_white_male_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_white_male_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.5,
                                                              trim_quantiles = c(0,1), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_notrim<-match_pop_white_male_noncompile_notrim$pseudo_pop

erf_notrim_white_male<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                             matched_w=match_pop_data_notrim$w,
                                             bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals, nthread=16)
#plot(erf_notrim_white_male)

#trimed ERC
#plot(erf_notrim_white_male$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_white_male$erf[trimed_index]/erf_notrim_white_male$erf[trimed_index][1])

#trimed matched dataset
match_pop_white_male_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                            w=treat,
                                                            c=c2,
                                                            ci_appr = "matching",
                                                            pred_model = "sl",
                                                            gps_model = "parametric",
                                                            use_cov_transform = FALSE,
                                                            #transformers = list("pow2", "pow3"),
                                                            sl_lib = c("m_xgboost"),
                                                            params = list("xgb_nrounds"=50,
                                                                          "xgb_max_depth"=6,
                                                                          "xgb_eta"=0.3,
                                                                          "xgb_min_child_weight"=1),
                                                            nthread=16, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                            optimized_compile = FALSE, #created a column counter for how many times matched,
                                                            max_attempt = 1,
                                                            matching_fun = "matching_l1",
                                                            delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                            scale = 1.0)

match_pop_data_trim<-match_pop_white_male_noncompile_trim$pseudo_pop

erf_trim_white_male<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                           matched_w = match_pop_data_trim$w,
                                           bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                           w_vals=a.vals,
                                           nthread=16)
#plot(erf_trim_white_male)
#plot(erf_trim_white_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_white_male$erf[trimed_index]/erf_trim_white_male$erf[trimed_index][1])

match_pop_white_male_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                              w = treat,
                                                              c = c,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds" = 50,
                                                                            "xgb_max_depth" = 6,
                                                                            "xgb_eta" = 0.3,
                                                                            "xgb_min_child_weight" = 1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)
match_pop_data_trim_onehot <- match_pop_white_male_noncompile_onehot$pseudo_pop

erf_trim_onehot_white_male <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                    matched_w = match_pop_data_trim_onehot$w,
                                                    bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                    w_vals = a.vals,
                                                    nthread=16)
#plot(erf_trim_onehot_white_male)

#plot(erf_trim_onehot_white_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_white_male$erf[trimed_index]/erf_trim_onehot_white_male$erf[trimed_index][1])

#Black female
prematch_data1<-prematchdata(black_female_rm, covariates_black_female_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]


#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_black_female_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                                w=treat,
                                                                c=c2,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                #transformers = list("pow2", "pow3"),
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds"=50,
                                                                              "xgb_max_depth"=6,
                                                                              "xgb_eta"=0.3,
                                                                              "xgb_min_child_weight"=1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.5,
                                                                trim_quantiles = c(0,1), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)

match_pop_data_notrim<-match_pop_black_female_noncompile_notrim$pseudo_pop

erf_notrim_black_female<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                               matched_w=match_pop_data_notrim$w,
                                               bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                               w_vals=a.vals, nthread=16)
#plot(erf_notrim_black_female)

#trimed ERC
#plot(erf_notrim_black_female$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_black_female$erf[trimed_index]/erf_notrim_black_female$erf[trimed_index][1])

#trimed matched dataset
match_pop_black_female_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_trim<-match_pop_black_female_noncompile_trim$pseudo_pop

erf_trim_black_female<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                             matched_w = match_pop_data_trim$w,
                                             bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals,
                                             nthread=16)
#plot(erf_trim_black_female)
#plot(erf_trim_black_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_black_female$erf[trimed_index]/erf_trim_black_female$erf[trimed_index][1])

match_pop_black_female_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                                w = treat,
                                                                c = c,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds" = 50,
                                                                              "xgb_max_depth" = 6,
                                                                              "xgb_eta" = 0.3,
                                                                              "xgb_min_child_weight" = 1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)
match_pop_data_trim_onehot <- match_pop_black_female_noncompile_onehot$pseudo_pop

erf_trim_onehot_black_female <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                      matched_w = match_pop_data_trim_onehot$w,
                                                      bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                      w_vals = a.vals,
                                                      nthread=16)
#plot(erf_trim_onehot_black_female)

#plot(erf_trim_onehot_black_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_black_female$erf[trimed_index]/erf_trim_onehot_black_female$erf[trimed_index][1])

#Black male
prematch_data1<-prematchdata(black_male_rm, covariates_black_male_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]


#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_black_male_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.5,
                                                              trim_quantiles = c(0,1), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_notrim<-match_pop_black_male_noncompile_notrim$pseudo_pop

erf_notrim_black_male<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                             matched_w=match_pop_data_notrim$w,
                                             bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals, nthread=16)
#plot(erf_notrim_black_male)

#trimed ERC
#plot(erf_notrim_black_male$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_black_male$erf[trimed_index]/erf_notrim_black_male$erf[trimed_index][1])

#trimed matched dataset
match_pop_black_male_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                            w=treat,
                                                            c=c2,
                                                            ci_appr = "matching",
                                                            pred_model = "sl",
                                                            gps_model = "parametric",
                                                            use_cov_transform = FALSE,
                                                            #transformers = list("pow2", "pow3"),
                                                            sl_lib = c("m_xgboost"),
                                                            params = list("xgb_nrounds"=50,
                                                                          "xgb_max_depth"=6,
                                                                          "xgb_eta"=0.3,
                                                                          "xgb_min_child_weight"=1),
                                                            nthread=16, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                            optimized_compile = FALSE, #created a column counter for how many times matched,
                                                            max_attempt = 1,
                                                            matching_fun = "matching_l1",
                                                            delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                            scale = 1.0)

match_pop_data_trim<-match_pop_black_male_noncompile_trim$pseudo_pop

erf_trim_black_male<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                           matched_w = match_pop_data_trim$w,
                                           bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                           w_vals=a.vals,
                                           nthread=16)
#plot(erf_trim_black_male)
#plot(erf_trim_black_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_black_male$erf[trimed_index]/erf_trim_black_male$erf[trimed_index][1])

match_pop_black_male_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                              w = treat,
                                                              c = c,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds" = 50,
                                                                            "xgb_max_depth" = 6,
                                                                            "xgb_eta" = 0.3,
                                                                            "xgb_min_child_weight" = 1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)
match_pop_data_trim_onehot <- match_pop_black_male_noncompile_onehot$pseudo_pop

erf_trim_onehot_black_male <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                    matched_w = match_pop_data_trim_onehot$w,
                                                    bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                    w_vals = a.vals,
                                                    nthread=16)
#plot(erf_trim_onehot_black_male)

#plot(erf_trim_onehot_black_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_black_male$erf[trimed_index]/erf_trim_onehot_black_male$erf[trimed_index][1])

#Hispanic female
prematch_data1 <-prematchdata(hispanic_female_rm, covariates_hispanic_female_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_hispanic_female_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                                   w=treat,
                                                                   c=c2,
                                                                   ci_appr = "matching",
                                                                   pred_model = "sl",
                                                                   gps_model = "parametric",
                                                                   use_cov_transform = FALSE,
                                                                   #transformers = list("pow2", "pow3"),
                                                                   sl_lib = c("m_xgboost"),
                                                                   params = list("xgb_nrounds"=50,
                                                                                 "xgb_max_depth"=6,
                                                                                 "xgb_eta"=0.3,
                                                                                 "xgb_min_child_weight"=1),
                                                                   nthread=16, # number of cores, you can change,
                                                                   covar_bl_method = "absolute",
                                                                   covar_bl_trs = 0.5,
                                                                   trim_quantiles = c(0,1), # trimed, you can change
                                                                   optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                   max_attempt = 1,
                                                                   matching_fun = "matching_l1",
                                                                   delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                   scale = 1.0)

match_pop_data_notrim<-match_pop_hispanic_female_noncompile_notrim$pseudo_pop

erf_notrim_hispanic_female<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                                  matched_w=match_pop_data_notrim$w,
                                                  bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                  w_vals=a.vals, nthread=16)
#plot(erf_notrim_hispanic_female)

#trimed ERC
#plot(erf_notrim_hispanic_female$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_hispanic_female$erf[trimed_index]/erf_notrim_hispanic_female$erf[trimed_index][1])

#trimed matched dataset
match_pop_hispanic_female_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                                 w=treat,
                                                                 c=c2,
                                                                 ci_appr = "matching",
                                                                 pred_model = "sl",
                                                                 gps_model = "parametric",
                                                                 use_cov_transform = FALSE,
                                                                 #transformers = list("pow2", "pow3"),
                                                                 sl_lib = c("m_xgboost"),
                                                                 params = list("xgb_nrounds"=50,
                                                                               "xgb_max_depth"=6,
                                                                               "xgb_eta"=0.3,
                                                                               "xgb_min_child_weight"=1),
                                                                 nthread=16, # number of cores, you can change,
                                                                 covar_bl_method = "absolute",
                                                                 covar_bl_trs = 0.1,
                                                                 trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                 optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                 max_attempt = 1,
                                                                 matching_fun = "matching_l1",
                                                                 delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                 scale = 1.0)

match_pop_data_trim<-match_pop_hispanic_female_noncompile_trim$pseudo_pop

erf_trim_hispanic_female<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                                matched_w = match_pop_data_trim$w,
                                                bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                w_vals=a.vals,
                                                nthread=16)
#plot(erf_trim_hispanic_female)
#plot(erf_trim_hispanic_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_hispanic_female$erf[trimed_index]/erf_trim_hispanic_female$erf[trimed_index][1])

match_pop_hispanic_female_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                                   w = treat,
                                                                   c = c,
                                                                   ci_appr = "matching",
                                                                   pred_model = "sl",
                                                                   gps_model = "parametric",
                                                                   use_cov_transform = FALSE,
                                                                   sl_lib = c("m_xgboost"),
                                                                   params = list("xgb_nrounds" = 50,
                                                                                 "xgb_max_depth" = 6,
                                                                                 "xgb_eta" = 0.3,
                                                                                 "xgb_min_child_weight" = 1),
                                                                   nthread=16, # number of cores, you can change,
                                                                   covar_bl_method = "absolute",
                                                                   covar_bl_trs = 0.1,
                                                                   trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                   optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                   max_attempt = 1,
                                                                   matching_fun = "matching_l1",
                                                                   delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                   scale = 1.0)
match_pop_data_trim_onehot <- match_pop_hispanic_female_noncompile_onehot$pseudo_pop

erf_trim_onehot_hispanic_female <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                         matched_w = match_pop_data_trim_onehot$w,
                                                         bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                         w_vals = a.vals,
                                                         nthread=16)
#plot(erf_trim_onehot_hispanic_female)

#plot(erf_trim_onehot_hispanic_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_hispanic_female$erf[trimed_index]/erf_trim_onehot_hispanic_female$erf[trimed_index][1])

#Hispanic male
prematch_data1<-prematchdata(hispanic_male_rm, covariates_hispanic_male_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_hispanic_male_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                                 w=treat,
                                                                 c=c2,
                                                                 ci_appr = "matching",
                                                                 pred_model = "sl",
                                                                 gps_model = "parametric",
                                                                 use_cov_transform = FALSE,
                                                                 #transformers = list("pow2", "pow3"),
                                                                 sl_lib = c("m_xgboost"),
                                                                 params = list("xgb_nrounds"=50,
                                                                               "xgb_max_depth"=6,
                                                                               "xgb_eta"=0.3,
                                                                               "xgb_min_child_weight"=1),
                                                                 nthread=16, # number of cores, you can change,
                                                                 covar_bl_method = "absolute",
                                                                 covar_bl_trs = 0.5,
                                                                 trim_quantiles = c(0,1), # trimed, you can change
                                                                 optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                 max_attempt = 1,
                                                                 matching_fun = "matching_l1",
                                                                 delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                 scale = 1.0)

match_pop_data_notrim<-match_pop_hispanic_male_noncompile_notrim$pseudo_pop

erf_notrim_hispanic_male<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                                matched_w=match_pop_data_notrim$w,
                                                bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                w_vals=a.vals, nthread=16)
#plot(erf_notrim_hispanic_male)

#trimed ERC
#plot(erf_notrim_hispanic_male$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_hispanic_male$erf[trimed_index]/erf_notrim_hispanic_male$erf[trimed_index][1])

#trimed matched dataset
match_pop_hispanic_male_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                               w=treat,
                                                               c=c2,
                                                               ci_appr = "matching",
                                                               pred_model = "sl",
                                                               gps_model = "parametric",
                                                               use_cov_transform = FALSE,
                                                               #transformers = list("pow2", "pow3"),
                                                               sl_lib = c("m_xgboost"),
                                                               params = list("xgb_nrounds"=50,
                                                                             "xgb_max_depth"=6,
                                                                             "xgb_eta"=0.3,
                                                                             "xgb_min_child_weight"=1),
                                                               nthread=16, # number of cores, you can change,
                                                               covar_bl_method = "absolute",
                                                               covar_bl_trs = 0.1,
                                                               trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                               optimized_compile = FALSE, #created a column counter for how many times matched,
                                                               max_attempt = 1,
                                                               matching_fun = "matching_l1",
                                                               delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                               scale = 1.0)

match_pop_data_trim<-match_pop_hispanic_male_noncompile_trim$pseudo_pop

erf_trim_hispanic_male<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                              matched_w = match_pop_data_trim$w,
                                              bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                              w_vals=a.vals,
                                              nthread=16)
#plot(erf_trim_hispanic_male)
#plot(erf_trim_hispanic_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_hispanic_male$erf[trimed_index]/erf_trim_hispanic_male$erf[trimed_index][1])

match_pop_hispanic_male_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                                 w = treat,
                                                                 c = c,
                                                                 ci_appr = "matching",
                                                                 pred_model = "sl",
                                                                 gps_model = "parametric",
                                                                 use_cov_transform = FALSE,
                                                                 sl_lib = c("m_xgboost"),
                                                                 params = list("xgb_nrounds" = 50,
                                                                               "xgb_max_depth" = 6,
                                                                               "xgb_eta" = 0.3,
                                                                               "xgb_min_child_weight" = 1),
                                                                 nthread=16, # number of cores, you can change,
                                                                 covar_bl_method = "absolute",
                                                                 covar_bl_trs = 0.1,
                                                                 trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                 optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                 max_attempt = 1,
                                                                 matching_fun = "matching_l1",
                                                                 delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                 scale = 1.0)
match_pop_data_trim_onehot <- match_pop_hispanic_male_noncompile_onehot$pseudo_pop

erf_trim_onehot_hispanic_male <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                       matched_w = match_pop_data_trim_onehot$w,
                                                       bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                       w_vals = a.vals,
                                                       nthread=16)
#plot(erf_trim_onehot_hispanic_male)

#plot(erf_trim_onehot_hispanic_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_hispanic_male$erf[trimed_index]/erf_trim_onehot_hispanic_male$erf[trimed_index][1])

#Asian female
prematch_data1<-prematchdata(asian_female_rm, covariates_asian_female_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))

treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_asian_female_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                                w=treat,
                                                                c=c2,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                #transformers = list("pow2", "pow3"),
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds"=50,
                                                                              "xgb_max_depth"=6,
                                                                              "xgb_eta"=0.3,
                                                                              "xgb_min_child_weight"=1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.5,
                                                                trim_quantiles = c(0,1), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)

match_pop_data_notrim<-match_pop_asian_female_noncompile_notrim$pseudo_pop

erf_notrim_asian_female<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                               matched_w=match_pop_data_notrim$w,
                                               bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                               w_vals=a.vals, nthread=16)
#plot(erf_notrim_asian_female)

#trimed ERC
#plot(erf_notrim_asian_female$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_asian_female$erf[trimed_index]/erf_notrim_asian_female$erf[trimed_index][1])

#trimed matched dataset
match_pop_asian_female_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_trim<-match_pop_asian_female_noncompile_trim$pseudo_pop

erf_trim_asian_female<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                             matched_w = match_pop_data_trim$w,
                                             bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals,
                                             nthread=16)
#plot(erf_trim_asian_female)
#plot(erf_trim_asian_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_asian_female$erf[trimed_index]/erf_trim_asian_female$erf[trimed_index][1])

match_pop_asian_female_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                                w = treat,
                                                                c = c,
                                                                ci_appr = "matching",
                                                                pred_model = "sl",
                                                                gps_model = "parametric",
                                                                use_cov_transform = FALSE,
                                                                sl_lib = c("m_xgboost"),
                                                                params = list("xgb_nrounds" = 50,
                                                                              "xgb_max_depth" = 6,
                                                                              "xgb_eta" = 0.3,
                                                                              "xgb_min_child_weight" = 1),
                                                                nthread=16, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                                max_attempt = 1,
                                                                matching_fun = "matching_l1",
                                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                                scale = 1.0)
match_pop_data_trim_onehot <- match_pop_asian_female_noncompile_onehot$pseudo_pop

erf_trim_onehot_asian_female <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                      matched_w = match_pop_data_trim_onehot$w,
                                                      bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                      w_vals = a.vals,
                                                      nthread=16)
#plot(erf_trim_onehot_asian_female)

#plot(erf_trim_onehot_asian_female$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_asian_female$erf[trimed_index]/erf_trim_onehot_asian_female$erf[trimed_index][1])

#Asian male
prematch_data1<-prematchdata(asian_male_rm, covariates_asian_male_rm)
all_rm.list<-split(prematch_data1, list(prematch_data1$zip))

q1<-quantile(prematch_data1$pm25, 0.05)
q2<-quantile(prematch_data1$pm25, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)

num_uniq_zip <- length(unique(prematch_data1$zip))
zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
prematch_data<-data.frame(Reduce(rbind, all_rm.list[zip_sample]))


treat=prematch_data["pm25"]$pm25
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

#Reproducing HEI report
c2=as.data.frame(data.matrix(c))

match_pop_asian_male_noncompile_notrim <- generate_pseudo_pop(Y=Y,
                                                              w=treat,
                                                              c=c2,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              #transformers = list("pow2", "pow3"),
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds"=50,
                                                                            "xgb_max_depth"=6,
                                                                            "xgb_eta"=0.3,
                                                                            "xgb_min_child_weight"=1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.5,
                                                              trim_quantiles = c(0,1), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)

match_pop_data_notrim<-match_pop_asian_male_noncompile_notrim$pseudo_pop

erf_notrim_asian_male<-estimate_npmetric_erf(matched_Y=match_pop_data_notrim$Y,
                                             matched_w=match_pop_data_notrim$w,
                                             bw_seq=seq(8*delta_n, 40*delta_n, 2*delta_n),
                                             w_vals=a.vals, nthread=16)
#plot(erf_notrim_asian_male)

#trimed ERC
#plot(erf_notrim_asian_male$erf[trimed_index])
#plot(a.vals[trimed_index], 
#     erf_notrim_asian_male$erf[trimed_index]/erf_notrim_asian_male$erf[trimed_index][1])

#trimed matched dataset
match_pop_asian_male_noncompile_trim <- generate_pseudo_pop(Y=Y,
                                                            w=treat,
                                                            c=c2,
                                                            ci_appr = "matching",
                                                            pred_model = "sl",
                                                            gps_model = "parametric",
                                                            use_cov_transform = FALSE,
                                                            #transformers = list("pow2", "pow3"),
                                                            sl_lib = c("m_xgboost"),
                                                            params = list("xgb_nrounds"=50,
                                                                          "xgb_max_depth"=6,
                                                                          "xgb_eta"=0.3,
                                                                          "xgb_min_child_weight"=1),
                                                            nthread=16, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                            optimized_compile = FALSE, #created a column counter for how many times matched,
                                                            max_attempt = 1,
                                                            matching_fun = "matching_l1",
                                                            delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                            scale = 1.0)

match_pop_data_trim<-match_pop_asian_male_noncompile_trim$pseudo_pop

erf_trim_asian_male<-estimate_npmetric_erf(matched_Y=match_pop_data_trim$Y,
                                           matched_w = match_pop_data_trim$w,
                                           bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                           w_vals=a.vals,
                                           nthread=16)
#plot(erf_trim_asian_male)
#plot(erf_trim_asian_male$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_asian_male$erf[trimed_index]/erf_trim_asian_male$erf[trimed_index][1])

match_pop_asian_male_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                              w = treat,
                                                              c = c,
                                                              ci_appr = "matching",
                                                              pred_model = "sl",
                                                              gps_model = "parametric",
                                                              use_cov_transform = FALSE,
                                                              sl_lib = c("m_xgboost"),
                                                              params = list("xgb_nrounds" = 50,
                                                                            "xgb_max_depth" = 6,
                                                                            "xgb_eta" = 0.3,
                                                                            "xgb_min_child_weight" = 1),
                                                              nthread=16, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.05,0.95), # trimed, you can change
                                                              optimized_compile = FALSE, #created a column counter for how many times matched,
                                                              max_attempt = 1,
                                                              matching_fun = "matching_l1",
                                                              delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                              scale = 1.0)
match_pop_data_trim_onehot <- match_pop_asian_male_noncompile_onehot$pseudo_pop

erf_trim_onehot_asian_male <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                                    matched_w = match_pop_data_trim_onehot$w,
                                                    bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                                    w_vals = a.vals,
                                                    nthread=16)
#plot(erf_trim_onehot_asian_male)

save(erf_notrim_all, erf_trim_all, erf_trim_onehot_all, 
     erf_notrim_white_female, erf_trim_white_female, erf_trim_onehot_white_female,
     erf_notrim_white_male, erf_trim_white_male, erf_trim_onehot_white_male,
     erf_notrim_black_female, erf_trim_black_female, erf_trim_onehot_black_female,
     erf_notrim_black_male, erf_trim_black_male, erf_trim_onehot_black_male,
     erf_notrim_hispanic_female, erf_trim_hispanic_female, erf_trim_onehot_hispanic_female,
     erf_notrim_hispanic_male, erf_trim_hispanic_male, erf_trim_onehot_hispanic_male,
     erf_notrim_asian_female, erf_notrim_asian_female, erf_trim_onehot_asian_female,
     erf_notrim_asian_male, erf_notrim_asian_male, erf_trim_onehot_asian_male,
     file=paste0('/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingrm/boots5/erf/erfgpscausal_rm_strata_trim5_', boots_id, '.RData'))


save(match_pop_all_noncompile_notrim, match_pop_all_noncompile_onehot, match_pop_all_noncompile_trim,
     match_pop_white_female_noncompile_notrim, match_pop_white_female_noncompile_onehot, match_pop_white_female_noncompile_trim,
     match_pop_white_male_noncompile_notrim, match_pop_white_male_noncompile_onehot, match_pop_white_male_noncompile_trim,
     match_pop_black_female_noncompile_notrim, match_pop_black_female_noncompile_onehot, match_pop_black_female_noncompile_trim,
     match_pop_black_male_noncompile_notrim, match_pop_black_male_noncompile_onehot, match_pop_black_male_noncompile_trim,
     match_pop_hispanic_female_noncompile_notrim, match_pop_hispanic_female_noncompile_onehot, match_pop_hispanic_female_noncompile_trim,
     match_pop_hispanic_male_noncompile_notrim, match_pop_hispanic_male_noncompile_onehot, match_pop_hispanic_male_noncompile_trim,
     match_pop_asian_female_noncompile_notrim, match_pop_asian_female_noncompile_onehot, match_pop_asian_female_noncompile_trim,
     match_pop_asian_male_noncompile_notrim, match_pop_asian_male_noncompile_onehot, match_pop_asian_male_noncompile_trim,
     file=paste0('/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/Bootstrap/matchingrm/boots5/gpscausal_rm_strata_trim5_', boots_id, '.RData'))
     

}
