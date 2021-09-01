library(devtools)
try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
install_github("fasrc/CausalGPS", ref="develop")
#install_github("fasrc/CausalGPS", ref="937810da2350f5c937c11eab16c60f2ee9a2783f")
library("CausalGPS")
library("dplyr")

set.seed(1)
dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

#All
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

load(paste0(dir_data,"aggregate_data_qd.RData"))
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
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_female_qd.RData")
covariates_white_female_qd$year<-as.factor(covariates_white_female_qd$year)
covariates_white_female_qd$region<-as.factor(covariates_white_female_qd$region)

#White male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_male_qd.RData")
covariates_white_male_qd$year<-as.factor(covariates_white_male_qd$year)
covariates_white_male_qd$region<-as.factor(covariates_white_male_qd$region)

#black female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_female_qd.RData")
covariates_black_female_qd$year<-as.factor(covariates_black_female_qd$year)
covariates_black_female_qd$region<-as.factor(covariates_black_female_qd$region)

#black male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_male_qd.RData")
covariates_black_male_qd$year<-as.factor(covariates_black_male_qd$year)
covariates_black_male_qd$region<-as.factor(covariates_black_male_qd$region)

#hispanic female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_female_qd.RData")
covariates_hispanic_female_qd$year<-as.factor(covariates_hispanic_female_qd$year)
covariates_hispanic_female_qd$region<-as.factor(covariates_hispanic_female_qd$region)

#hispanic male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_male_qd.RData")
covariates_hispanic_male_qd$year<-as.factor(covariates_hispanic_male_qd$year)
covariates_hispanic_male_qd$region<-as.factor(covariates_hispanic_male_qd$region)

#asian female
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_female_qd.RData")
covariates_asian_female_qd$year<-as.factor(covariates_asian_female_qd$year)
covariates_asian_female_qd$region<-as.factor(covariates_asian_female_qd$region)

#asian male
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_male_qd.RData")
covariates_asian_male_qd$year<-as.factor(covariates_asian_male_qd$year)
covariates_asian_male_qd$region<-as.factor(covariates_asian_male_qd$region)

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

#All
prematch_data<-prematchdata(aggregate_data_qd, covariates_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
#Reproducing HEI report
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
                                      nthread = 8, # number of cores, you can change,
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
                                  w_vals=a.vals, nthread=10)
plot(erf_notrim_all)

#trimed ERC
plot(erf_notrim_all$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_all$erf[trimed_index]/erf_notrim_all$erf[trimed_index][1])

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
                                                         nthread = 8, # number of cores, you can change,
                                                         covar_bl_method = "absolute",
                                                         covar_bl_trs = 0.1,
                                                         trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                nthread=10)
plot(erf_trim_all)
plot(erf_trim_all$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_all$erf[trimed_index]/erf_trim_all$erf[trimed_index][1])

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
                                                       nthread = 8, # number of cores, you can change,
                                                       covar_bl_method = "absolute",
                                                       covar_bl_trs = 0.1,
                                                       trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                         nthread = 10)
plot(erf_trim_onehot_all)

plot(erf_trim_onehot_all$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_all$erf[trimed_index]/erf_trim_onehot_all$erf[trimed_index][1])


#Strata
#White female
prematch_data<-prematchdata(white_female_qd, covariates_white_female_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                                nthread = 8, # number of cores, you can change,
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
                                               w_vals=a.vals, nthread=10)
plot(erf_notrim_white_female)

#trimed ERC
plot(erf_notrim_white_female$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_white_female$erf[trimed_index]/erf_notrim_white_female$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                             nthread=10)
plot(erf_trim_white_female)
plot(erf_trim_white_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_white_female$erf[trimed_index]/erf_trim_white_female$erf[trimed_index][1])

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
                                                                nthread = 8, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                      nthread = 10)
plot(erf_trim_onehot_white_female)

plot(erf_trim_onehot_white_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_white_female$erf[trimed_index]/erf_trim_onehot_white_female$erf[trimed_index][1])


#White male
prematch_data<-prematchdata(white_male_qd, covariates_white_male_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                              nthread = 8, # number of cores, you can change,
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
                                             w_vals=a.vals, nthread=10)
plot(erf_notrim_white_male)

#trimed ERC
plot(erf_notrim_white_male$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_white_male$erf[trimed_index]/erf_notrim_white_male$erf[trimed_index][1])

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
                                                            nthread = 8, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                           nthread=10)
plot(erf_trim_white_male)
plot(erf_trim_white_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_white_male$erf[trimed_index]/erf_trim_white_male$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                    nthread = 10)
plot(erf_trim_onehot_white_male)

plot(erf_trim_onehot_white_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_white_male$erf[trimed_index]/erf_trim_onehot_white_male$erf[trimed_index][1])

#Black female
prematch_data<-prematchdata(black_female_qd, covariates_black_female_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                                nthread = 8, # number of cores, you can change,
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
                                               w_vals=a.vals, nthread=10)
plot(erf_notrim_black_female)

#trimed ERC
plot(erf_notrim_black_female$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_black_female$erf[trimed_index]/erf_notrim_black_female$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                             nthread=10)
plot(erf_trim_black_female)
plot(erf_trim_black_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_black_female$erf[trimed_index]/erf_trim_black_female$erf[trimed_index][1])

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
                                                                nthread = 8, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                      nthread = 10)
plot(erf_trim_onehot_black_female)

plot(erf_trim_onehot_black_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_black_female$erf[trimed_index]/erf_trim_onehot_black_female$erf[trimed_index][1])

#Black male
prematch_data<-prematchdata(black_male_qd, covariates_black_male_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                              nthread = 8, # number of cores, you can change,
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
                                             w_vals=a.vals, nthread=10)
plot(erf_notrim_black_male)

#trimed ERC
plot(erf_notrim_black_male$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_black_male$erf[trimed_index]/erf_notrim_black_male$erf[trimed_index][1])

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
                                                            nthread = 8, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                           nthread=10)
plot(erf_trim_black_male)
plot(erf_trim_black_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_black_male$erf[trimed_index]/erf_trim_black_male$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                    nthread = 10)
plot(erf_trim_onehot_black_male)

plot(erf_trim_onehot_black_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_black_male$erf[trimed_index]/erf_trim_onehot_black_male$erf[trimed_index][1])

#Hispanic female
prematch_data<-prematchdata(hispanic_female_qd, covariates_hispanic_female_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                                   nthread = 8, # number of cores, you can change,
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
                                                  w_vals=a.vals, nthread=10)
plot(erf_notrim_hispanic_female)

#trimed ERC
plot(erf_notrim_hispanic_female$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_hispanic_female$erf[trimed_index]/erf_notrim_hispanic_female$erf[trimed_index][1])

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
                                                                 nthread = 8, # number of cores, you can change,
                                                                 covar_bl_method = "absolute",
                                                                 covar_bl_trs = 0.1,
                                                                 trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                nthread=10)
plot(erf_trim_hispanic_female)
plot(erf_trim_hispanic_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_hispanic_female$erf[trimed_index]/erf_trim_hispanic_female$erf[trimed_index][1])

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
                                                                   nthread = 8, # number of cores, you can change,
                                                                   covar_bl_method = "absolute",
                                                                   covar_bl_trs = 0.1,
                                                                   trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                         nthread = 10)
plot(erf_trim_onehot_hispanic_female)

plot(erf_trim_onehot_hispanic_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_hispanic_female$erf[trimed_index]/erf_trim_onehot_hispanic_female$erf[trimed_index][1])

#Hispanic male
prematch_data<-prematchdata(hispanic_male_qd, covariates_hispanic_male_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                                 nthread = 8, # number of cores, you can change,
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
                                                w_vals=a.vals, nthread=10)
plot(erf_notrim_hispanic_male)

#trimed ERC
plot(erf_notrim_hispanic_male$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_hispanic_male$erf[trimed_index]/erf_notrim_hispanic_male$erf[trimed_index][1])

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
                                                               nthread = 8, # number of cores, you can change,
                                                               covar_bl_method = "absolute",
                                                               covar_bl_trs = 0.1,
                                                               trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                              nthread=10)
plot(erf_trim_hispanic_male)
plot(erf_trim_hispanic_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_hispanic_male$erf[trimed_index]/erf_trim_hispanic_male$erf[trimed_index][1])

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
                                                                 nthread = 8, # number of cores, you can change,
                                                                 covar_bl_method = "absolute",
                                                                 covar_bl_trs = 0.1,
                                                                 trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                       nthread = 10)
plot(erf_trim_onehot_hispanic_male)

plot(erf_trim_onehot_hispanic_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_hispanic_male$erf[trimed_index]/erf_trim_onehot_hispanic_male$erf[trimed_index][1])

#Asian female
prematch_data<-prematchdata(asian_female_qd, covariates_asian_female_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                                nthread = 8, # number of cores, you can change,
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
                                               w_vals=a.vals, nthread=10)
plot(erf_notrim_asian_female)

#trimed ERC
plot(erf_notrim_asian_female$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_asian_female$erf[trimed_index]/erf_notrim_asian_female$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                             nthread=10)
plot(erf_trim_asian_female)
plot(erf_trim_asian_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_asian_female$erf[trimed_index]/erf_trim_asian_female$erf[trimed_index][1])

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
                                                                nthread = 8, # number of cores, you can change,
                                                                covar_bl_method = "absolute",
                                                                covar_bl_trs = 0.1,
                                                                trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                      nthread = 10)
plot(erf_trim_onehot_asian_female)

plot(erf_trim_onehot_asian_female$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_asian_female$erf[trimed_index]/erf_trim_onehot_asian_female$erf[trimed_index][1])

#Asian male
prematch_data<-prematchdata(asian_male_qd, covariates_asian_male_qd)

treat=prematch_data["pm25_ensemble"]$pm25_ensemble
c=prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "education", "popdensity", "pct_owner_occ",
                   "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", 
                   "region", "year")]
Y=prematch_data[,"mortality"]

q1<-quantile(prematch_data$pm25_ensemble, 0.01)
q2<-quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
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
                                                              nthread = 8, # number of cores, you can change,
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
                                             w_vals=a.vals, nthread=10)
plot(erf_notrim_asian_male)

#trimed ERC
plot(erf_notrim_asian_male$erf[trimed_index])
plot(a.vals[trimed_index], 
     erf_notrim_asian_male$erf[trimed_index]/erf_notrim_asian_male$erf[trimed_index][1])

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
                                                            nthread = 8, # number of cores, you can change,
                                                            covar_bl_method = "absolute",
                                                            covar_bl_trs = 0.1,
                                                            trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                           nthread=10)
plot(erf_trim_asian_male)
plot(erf_trim_asian_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_asian_male$erf[trimed_index]/erf_trim_asian_male$erf[trimed_index][1])

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
                                                              nthread = 8, # number of cores, you can change,
                                                              covar_bl_method = "absolute",
                                                              covar_bl_trs = 0.1,
                                                              trim_quantiles = c(0.01,0.99), # trimed, you can change
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
                                                    nthread = 10)
plot(erf_trim_onehot_asian_male)

plot(erf_trim_onehot_asian_male$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot_asian_male$erf[trimed_index]/erf_trim_onehot_asian_male$erf[trimed_index][1])

save(match_pop_all_noncompile_notrim, match_pop_all_noncompile_onehot, match_pop_all_noncompile_trim,
 match_pop_white_female_noncompile_notrim, match_pop_white_female_noncompile_onehot, match_pop_white_female_noncompile_trim,
 match_pop_white_male_noncompile_notrim, match_pop_white_male_noncompile_onehot, match_pop_white_male_noncompile_trim,
 match_pop_black_female_noncompile_notrim, match_pop_black_female_noncompile_onehot, match_pop_black_female_noncompile_trim,
 match_pop_black_male_noncompile_notrim, match_pop_black_male_noncompile_onehot, match_pop_black_male_noncompile_trim,
 match_pop_hispanic_female_noncompile_notrim, match_pop_hispanic_female_noncompile_onehot, match_pop_hispanic_female_noncompile_trim,
 match_pop_hispanic_male_noncompile_notrim, match_pop_hispanic_male_noncompile_onehot, match_pop_hispanic_male_noncompile_trim,
 match_pop_asian_female_noncompile_notrim, match_pop_asian_female_noncompile_onehot, match_pop_asian_female_noncompile_trim,
 match_pop_asian_male_noncompile_notrim, match_pop_asian_male_noncompile_onehot, match_pop_asian_male_noncompile_trim,
 file='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/gpscausal_qd_strata_trim1.RData'
)
#Plots
#all
require(ggplot2)
#plot(a.vals[trimed_index], erf_notrim_all$erf[trimed_index], ylim=c(0, 0.06))
#plot(a.vals[trimed_index], erf_trim_all$erf[trimed_index])
#plot(a.vals[trimed_index], erf_trim_onehot_all$erf[trimed_index])
#df_all$erf_notrim_all<-erf_notrim_all$erf[trimed_index]
#df_all$erf_trim_all<-erf_trim_all$erf[trimed_index]
#df_all$erf_trim_onehot_all<-erf_trim_onehot_all$erf[trimed_index]

#All
df_all<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_all$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_all$a.vals<-rep(a.vals[trimed_index], 3)
df_all$erf<-c(erf_notrim_all$erf[trimed_index],erf_trim_all$erf[trimed_index], erf_trim_onehot_all$erf[trimed_index] )
p_all<-ggplot(df_all, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF All")

#Strata
df_white_female<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_white_female$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_white_female$a.vals<-rep(a.vals[trimed_index], 3)
df_white_female$erf<-c(erf_notrim_white_female$erf[trimed_index],erf_trim_white_female$erf[trimed_index], erf_trim_onehot_white_female$erf[trimed_index] )
p_white_female<-ggplot(df_white_female, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF white_female")

df_white_male<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_white_male$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_white_male$a.vals<-rep(a.vals[trimed_index], 3)
df_white_male$erf<-c(erf_notrim_white_male$erf[trimed_index],erf_trim_white_male$erf[trimed_index], erf_trim_onehot_white_male$erf[trimed_index] )
p_white_male<-ggplot(df_white_male, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF white_male")

df_black_female<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_black_female$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_black_female$a.vals<-rep(a.vals[trimed_index], 3)
df_black_female$erf<-c(erf_notrim_black_female$erf[trimed_index],erf_trim_black_female$erf[trimed_index], erf_trim_onehot_black_female$erf[trimed_index] )
p_black_female<-ggplot(df_black_female, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF black_female")

df_black_male<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_black_male$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_black_male$a.vals<-rep(a.vals[trimed_index], 3)
df_black_male$erf<-c(erf_notrim_black_male$erf[trimed_index],erf_trim_black_male$erf[trimed_index], erf_trim_onehot_black_male$erf[trimed_index] )
p_black_male<-ggplot(df_black_male, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF black_male")

df_hispanic_female<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_hispanic_female$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_hispanic_female$a.vals<-rep(a.vals[trimed_index], 3)
df_hispanic_female$erf<-c(erf_notrim_hispanic_female$erf[trimed_index],erf_trim_hispanic_female$erf[trimed_index], erf_trim_onehot_hispanic_female$erf[trimed_index] )
p_hispanic_female<-ggplot(df_hispanic_female, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF hispanic_female")

df_hispanic_male<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_hispanic_male$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_hispanic_male$a.vals<-rep(a.vals[trimed_index], 3)
df_hispanic_male$erf<-c(erf_notrim_hispanic_male$erf[trimed_index],erf_trim_hispanic_male$erf[trimed_index], erf_trim_onehot_hispanic_male$erf[trimed_index] )
p_hispanic_male<-ggplot(df_hispanic_male, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF hispanic_male")

df_asian_female<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_asian_female$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_asian_female$a.vals<-rep(a.vals[trimed_index], 3)
df_asian_female$erf<-c(erf_notrim_asian_female$erf[trimed_index],erf_trim_asian_female$erf[trimed_index], erf_trim_onehot_asian_female$erf[trimed_index] )
p_asian_female<-ggplot(df_asian_female, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF asian_female")

df_asian_male<-as.data.frame(matrix(nrow=48*3, ncol=3))
df_asian_male$type<-c(rep("No Trim", 48), rep("Trim", 48), rep("Trim One Hot", 48))
df_asian_male$a.vals<-rep(a.vals[trimed_index], 3)
df_asian_male$erf<-c(erf_notrim_asian_male$erf[trimed_index],erf_trim_asian_male$erf[trimed_index], erf_trim_onehot_asian_male$erf[trimed_index] )
p_asian_male<-ggplot(df_asian_male, aes(x=a.vals, y=erf, color=type)) + geom_line()+ggtitle("ERF asian_male")

require(cowplot)
plot_grid(p_white_female, p_white_male, 
          p_black_female, p_black_male, 
          p_hispanic_female, p_hispanic_male,
          p_asian_female, p_asian_male)

save(match_pop_all_noncompile_notrim, match_pop_all_noncompile_onehot, match_pop_all_noncompile_trim,
     match_pop_white_female_noncompile_notrim, match_pop_white_female_noncompile_onehot, match_pop_white_female_noncompile_trim,
     match_pop_white_male_noncompile_notrim, match_pop_white_male_noncompile_onehot, match_pop_white_male_noncompile_trim,
     match_pop_black_female_noncompile_notrim, match_pop_black_female_noncompile_onehot, match_pop_black_female_noncompile_trim,
     match_pop_black_male_noncompile_notrim, match_pop_black_male_noncompile_onehot, match_pop_black_male_noncompile_trim,
     match_pop_hispanic_female_noncompile_notrim, match_pop_hispanic_female_noncompile_onehot, match_pop_hispanic_female_noncompile_trim,
     match_pop_hispanic_male_noncompile_notrim, match_pop_hispanic_male_noncompile_onehot, match_pop_hispanic_male_noncompile_trim,
     match_pop_asian_female_noncompile_notrim, match_pop_asian_female_noncompile_onehot, match_pop_asian_female_noncompile_trim,
     match_pop_asian_male_noncompile_notrim, match_pop_asian_male_noncompile_onehot, match_pop_asian_male_noncompile_trim,
     file='/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/gpscausal_qd_strata_trim1.RData'
)