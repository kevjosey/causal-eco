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
require(rlist)

#All
set.seed(1)
dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Matching/'

#All
load("/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Data/balance_qd/covariates_qd.RData")
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
q1<-quantile(prematch_data$pm25_ensemble, 0.05)
q2<-quantile(prematch_data$pm25_ensemble, 0.95)
trimed_index<-which(a.vals >=q1 & a.vals <=q2)
rm(aggregate_data_qd)

dir_data = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Matching/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/pdez_measurementerror/Output/Bootstrap/matchingrm/boots5/'
load(paste0(dir_data, "gpscausal_qd_strata_trim5.RData"))

match_pop_all_noncompile_trim$pseudo_pop$year<-as.factor(match_pop_all_noncompile_trim$pseudo_pop$year)
match_pop_all_noncompile_trim$pseudo_pop$region<-as.factor(match_pop_all_noncompile_trim$pseudo_pop$region)
temp<-cbind(match_pop_all_noncompile_trim$pseudo_pop, 
            prematch_data$dead[match_pop_all_noncompile_trim$pseudo_pop$row_index],
            prematch_data$time_count[match_pop_all_noncompile_trim$pseudo_pop$row_index])
#temp$mortality<-temp$V2/temp$V3
#Trimmed
#All
matchingqd_gnm<-summary(gnm(V2~w + +offset(log(V3)),
                            data=temp, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[2])

list1<-list.files(dir_out, full.names = TRUE)
res<-c()
for(i in 2:length(list)){
  f<-list1[i]
  load(f)
  match_pop_all_noncompile_trim$pseudo_pop$year<-as.factor(match_pop_all_noncompile_trim$pseudo_pop$year)
  match_pop_all_noncompile_trim$pseudo_pop$region<-as.factor(match_pop_all_noncompile_trim$pseudo_pop$region)
  temp<-cbind(match_pop_all_noncompile_trim$pseudo_pop, 
              prematch_data$dead[match_pop_all_noncompile_trim$pseudo_pop$row_index],
              prematch_data$time_count[match_pop_all_noncompile_trim$pseudo_pop$row_index])
  
  matchingqd_gnm<-summary(gnm(V2~w + +offset(log(V3)),
                              data=temp, 
                              family=poisson(link="log")))
  res<- list.append(res, exp(10*matchingqd_gnm$coefficients[2]))
  
}

#White female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_white_female_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#White male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_white_male_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Black female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_black_female_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Black male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_black_male_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_hispanic_female_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_hispanic_male_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Asian female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_asian_female_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Asian male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_asian_male_noncompile_trim$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])


#One hot
#All
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_all_noncompile_onehot$pseudo_pop[,], 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#White female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_white_female_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])exp(10*matchingqd_gnm$coefficients[1])

#White male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_white_male_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Black female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_black_female_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Black male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_black_male_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_hispanic_female_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Hispanic male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_hispanic_male_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])

#Asian female
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_asian_female_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])exp(10*matchingqd_gnm$coefficients[1])

#Asian male
matchingqd_gnm<-summary(gnm(Y~w,
                            data=match_pop_asian_male_noncompile_onehot$pseudo_pop, 
                            family=poisson(link="log")))
exp(10*matchingqd_gnm$coefficients[1])exp(10*matchingqd_gnm$coefficients[1])



#Plots
#all
require(ggplot2)

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