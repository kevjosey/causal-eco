library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
require(doParallel)
library(data.table)
library(fst)
library("parallel")
require(xgboost)

dir_data =  '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out =  '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
mc_cores=15

# Fit non-parametric kernel smoothing on matched set
matching_smooth<-function(matched_Y,
                          matched_w,
                          bw.seq = seq(0.2, 2, 0.2),
                          w.vals){
  ## The specified Gaussian kernel
  kern_fun <- function(t){ dnorm(t) }
  w_fun <- function(bw){
    w.avals <- NULL
    for (w.val in w.vals){
      w.std <- (matched_w-w.val) / bw
      kern.std <- kern_fun(w.std) / bw
      w.avals <- c(w.avals, mean(w.std^2*kern.std)*(kern_fun(0)/bw) /
                     (mean(kern.std) * mean(w.std^2 * kern.std) - mean(w.std * kern.std)^2))
    }
    return(w.avals / length(matched_w))
  }
  hatvals <- function(bw){approx(w.vals,
                                 w_fun(bw),
                                 xout = matched_w,
                                 rule = 2)$y}
  smooth_fun <- function(out,bw){
    approx(locpoly(matched_w,
                   out,bandwidth = bw,
                   gridsize = 1000),
           xout = matched_w, rule = 2)$y
  }
  ##
  risk_fun <- function(h){
    hats <- hatvals(h); mean(((matched_Y - smooth_fun(matched_Y, bw = h)) / (1 - hats))^2)
  }
  risk.val <- sapply(bw.seq, risk_fun)
  h.opt <- bw.seq[which.min(risk.val)]
  
  erf <- approx(locpoly(matched_w, matched_Y, bandwidth = h.opt), xout = w.vals)$y
  return(erf)
}

#National RM
f <- list.files(paste0(dir_data,"Matching_Output_rm/all/"),
                pattern = "\\.rds",
                full.names = TRUE)

matching_rm<-rbindlist(lapply(f, readRDS))
matching_rm2<-subset(matching_rm,time_count>0)
matching_rm2<-as.data.frame(lapply(matching_rm2, unlist))
matching_rm<- as.data.frame(lapply(matching_rm, unlist))

#Strata
#White female
f <- list.files(paste0(dir_data,"Matching_Output_rm/white_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

white_female_rm<-rbindlist(lapply(f, readRDS))
white_female_rm<-subset(white_female_rm,time_count>0)
white_female_rm<-as.data.frame(lapply(white_female_rm, unlist))
#White male
f <- list.files(paste0(dir_data,"Matching_Output_rm/white_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

white_male_rm<-rbindlist(lapply(f, readRDS))
white_male_rm<-subset(white_male_rm,time_count>0)
white_male_rm<-as.data.frame(lapply(white_male_rm, unlist))
#Black female
f <- list.files(paste0(dir_data,"Matching_Output_rm/black_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

black_female_rm<-rbindlist(lapply(f, readRDS))
black_female_rm<-subset(black_female_rm,time_count>0)
black_female_rm<-as.data.frame(lapply(black_female_rm, unlist))
#Black male
f <- list.files(paste0(dir_data,"Matching_Output_rm/black_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

black_male_rm<-rbindlist(lapply(f, readRDS))
black_male_rm<-subset(black_male_rm,time_count>0)
black_male_rm<-as.data.frame(lapply(black_male_rm, unlist))
#Hispanic female
f <- list.files(paste0(dir_data,"Matching_Output_rm/hispanic_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

hispanic_female_rm<-rbindlist(lapply(f, readRDS))
hispanic_female_rm<-subset(hispanic_female_rm,time_count>0)
hispanic_female_rm<-as.data.frame(lapply(hispanic_female_rm, unlist))
#Hispanic male
f <- list.files(paste0(dir_data,"Matching_Output_rm/hispanic_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

hispanic_male_rm<-rbindlist(lapply(f, readRDS))
hispanic_male_rm<-subset(hispanic_male_rm,time_count>0)
hispanic_male_rm<-as.data.frame(lapply(hispanic_male_rm, unlist))
#Asian female
f <- list.files(paste0(dir_data,"Matching_Output_rm/asian_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

asian_female_rm<-rbindlist(lapply(f, readRDS))
asian_female_rm<-subset(asian_female_rm,time_count>0)
asian_female_rm<-as.data.frame(lapply(asian_female_rm, unlist))
#Asian male
f <- list.files(paste0(dir_data,"Matching_Output_rm/asian_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

asian_male_rm<-rbindlist(lapply(f, readRDS))
asian_male_rm<-subset(asian_male_rm,time_count>0)
asian_male_rm<-as.data.frame(lapply(asian_male_rm, unlist))


#National QD
f <- list.files(paste0(dir_data,"Matching_Output_qd/all/"),
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd<-rbindlist(lapply(f, readRDS))
matching_qd2<-subset(matching_qd,time_count>0)
matching_qd2<-as.data.frame(lapply(matching_qd2, unlist))

#Strata
#White female
f <- list.files(paste0(dir_data,"Matching_Output_qd/white_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

white_female_qd<-rbindlist(lapply(f, readRDS))
white_female_qd<-subset(white_female_qd,time_count>0)
white_female_qd<-as.data.frame(lapply(white_female_qd, unlist))
#White male
f <- list.files(paste0(dir_data,"Matching_Output_qd/white_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

white_male_qd<-rbindlist(lapply(f, readRDS))
white_male_qd<-subset(white_male_qd,time_count>0)
white_male_qd<-as.data.frame(lapply(white_male_qd, unlist))
#Black female
f <- list.files(paste0(dir_data,"Matching_Output_qd/black_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

black_female_qd<-rbindlist(lapply(f, readRDS))
black_female_qd<-subset(black_female_qd,time_count>0)
black_female_qd<-as.data.frame(lapply(black_female_qd, unlist))
#Black male
f <- list.files(paste0(dir_data,"Matching_Output_qd/black_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

black_male_qd<-rbindlist(lapply(f, readRDS))
black_male_qd<-subset(black_male_qd,time_count>0)
black_male_qd<-as.data.frame(lapply(black_male_qd, unlist))
#Hispanic female
f <- list.files(paste0(dir_data,"Matching_Output_qd/hispanic_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

hispanic_female_qd<-rbindlist(lapply(f, readRDS))
hispanic_female_qd<-subset(hispanic_female_qd,time_count>0)
hispanic_female_qd<-as.data.frame(lapply(hispanic_female_qd, unlist))
#Hispanic male
f <- list.files(paste0(dir_data,"Matching_Output_qd/hispanic_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

hispanic_male_qd<-rbindlist(lapply(f, readRDS))
hispanic_male_qd<-subset(hispanic_male_qd,time_count>0)
hispanic_male_qd<-as.data.frame(lapply(hispanic_male_qd, unlist))
#Asian female
f <- list.files(paste0(dir_data,"Matching_Output_qd/asian_female/"),
                pattern = "\\.rds",
                full.names = TRUE)

asian_female_qd<-rbindlist(lapply(f, readRDS))
asian_female_qd<-subset(asian_female_qd,time_count>0)
asian_female_qd<-as.data.frame(lapply(asian_female_qd, unlist))
#Asian male
f <- list.files(paste0(dir_data,"Matching_Output_qd/asian_male/"),
                pattern = "\\.rds",
                full.names = TRUE)

asian_male_qd<-rbindlist(lapply(f, readRDS))
asian_male_qd<-subset(asian_male_qd,time_count>0)
asian_male_qd<-as.data.frame(lapply(asian_male_qd, unlist))


######RM
#ALl
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                  eliminate=(as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                  data=matching_rm2, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])


#White female
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=white_female_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#White Male
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=white_male_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])


#Black female
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=black_female_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Black Male
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=black_male_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])


#Hispanic female
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=hispanic_female_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Hispanic Male
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=hispanic_male_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Asian female
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=asian_female_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Asian Male
matching_gnm<-summary(gnm(dead~ pm25 + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=asian_male_rm, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#######QD
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=matching_qd2, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#White female
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=white_female_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#White Male
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=white_male_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])


#Black female
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=black_female_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Black Male
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=black_male_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])


#Hispanic female
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=hispanic_female_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Hispanic Male
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=hispanic_male_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Asian female
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=asian_female_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])

#Asian Male
matching_gnm<-summary(gnm(dead~ pm25_ensemble + offset(log(time_count)), 
                          eliminate=(as.factor(sex):as.factor(race):as.factor(entry_age_break):as.factor(dual):as.factor(followup_year)),
                          data=asian_male_qd, family=poisson(link="log")))
exp(10*matching_gnm$coefficients[1])              

#################GAM
#RM
matching_rm4<-matching_rm2[matching_rm2$pm25<quantile(matching_rm2$pm25, 0.95) & 
                      matching_rm2$pm25> quantile(matching_rm2$pm25, 0.05),]

matching_rm3<-subset(matching_rm2, matching_rm2$pm25<=15.41103 & matching_rm2$pm25>=4.341009)
matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                      as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                      offset(log(time_count))
                  , data=matching_rm3,family=poisson(link="log"))

white_female_rm3<-subset(white_female_rm, white_female_rm$pm25<=15.41103 & white_female_rm$pm25>=4.341009)
white_male_rm3<-subset(white_male_rm, white_male_rm$pm25<=15.41103 & white_male_rm$pm25>=4.341009)
black_female_rm3<-subset(black_female_rm, black_female_rm$pm25<=15.41103 & black_female_rm$pm25>=4.341009)
black_male_rm3<-subset(black_male_rm, black_male_rm$pm25<=15.41103 & black_male_rm$pm25>=4.341009)
hispanic_female_rm3<-subset(hispanic_female_rm, hispanic_female_rm$pm25<=15.41103 & hispanic_female_rm$pm25>=4.341009)
hispanic_male_rm3<-subset(hispanic_male_rm, hispanic_male_rm$pm25<=15.41103 & hispanic_male_rm$pm25>=4.341009)
asian_female_rm3<-subset(asian_female_rm, asian_female_rm$pm25<=15.41103 & asian_female_rm$pm25>=4.341009)
asian_male_rm3<-subset(asian_male_rm, asian_male_rm$pm25<=15.41103 & asian_male_rm$pm25>=4.341009)

white_female_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                            offset(log(time_count))
                          , data=white_female_rm3,family=poisson(link="log"))

white_male_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=white_male_rm3,family=poisson(link="log"))

black_female_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=black_female_rm3,family=poisson(link="log"))

black_male_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                       offset(log(time_count))
                                     , data=black_male_rm3,family=poisson(link="log"))

hispanic_female_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=hispanic_female_rm3,family=poisson(link="log"))

hispanic_male_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                       offset(log(time_count))
                                     , data=hispanic_male_rm3,family=poisson(link="log"))

asian_female_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                            offset(log(time_count))
                                          , data=asian_female_rm3,family=poisson(link="log"))

asian_male_matching_gam3 <-mgcv::bam(dead~ s(pm25,k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=asian_male_rm3,family=poisson(link="log"))


plot(matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")
par(mfrow=c(4,2))
plot(white_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="White female HR")
plot(white_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="White male HR")
plot(black_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black female HR")
plot(black_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black male HR")
plot(hispanic_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic female HR")
plot(hispanic_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic male HR")
plot(asian_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian female HR")
plot(asian_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian male HR")

#Kernel smoothinh
#rm2
w.vals=seq(min(matching_rm3$pm25),
           max(matching_rm3$pm25), length.out=50)
erf_rm2<-matching_smooth(matching_rm2$dead/matching_rm2$time_count, matching_rm2$pm25, 
                         bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df<-data.frame(matrix(nrow=50, ncol=2))
df[,1]<-w.vals
df[,2]<-erf_rm2
require(ggplot2)
p<-ggplot(df, aes(x=X1, y=X2))+geom_line()+labs(x="pm25", y="deaths/timecount")

#White
w.vals=seq(min(white_female_rm2$pm25),
           max(white_female_rm2$pm25), length.out=50)
erf_white_female_rm2<-matching_smooth(white_female_rm2$dead/white_female_rm2$time_count, white_female_rm2$pm25, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_white_female_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_white_female_rm2[,1]<-w.vals
df_white_female_rm2[,2]<-erf_white_female_rm2

w.vals=seq(min(white_male_rm2$pm25),
           max(white_male_rm2$pm25), length.out=50)
erf_white_male_rm2<-matching_smooth(white_male_rm2$dead/white_male_rm2$time_count, white_male_rm2$pm25, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_white_male_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_white_male_rm2[,1]<-w.vals
df_white_male_rm2[,2]<-erf_white_male_rm2


#Black
w.vals=seq(min(black_female_rm2$pm25),
           max(black_female_rm2$pm25), length.out=50)
erf_black_female_rm2<-matching_smooth(black_female_rm2$dead/black_female_rm2$time_count, black_female_rm2$pm25, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_black_female_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_black_female_rm2[,1]<-w.vals
df_black_female_rm2[,2]<-erf_black_female_rm2

w.vals=seq(min(black_male_rm2$pm25),
           max(black_male_rm2$pm25), length.out=50)
erf_black_male_rm2<-matching_smooth(black_male_rm2$dead/black_male_rm2$time_count, black_male_rm2$pm25, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_black_male_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_black_male_rm2[,1]<-w.vals
df_black_male_rm2[,2]<-erf_black_male_rm2

#Hispanic
w.vals=seq(min(hispanic_female_rm2$pm25),
           max(hispanic_female_rm2$pm25), length.out=50)
erf_hispanic_female_rm2<-matching_smooth(hispanic_female_rm2$dead/hispanic_female_rm2$time_count, hispanic_female_rm2$pm25, 
                                         bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_hispanic_female_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_hispanic_female_rm2[,1]<-w.vals
df_hispanic_female_rm2[,2]<-erf_hispanic_female_rm2

w.vals=seq(min(hispanic_male_rm2$pm25),
           max(hispanic_male_rm2$pm25), length.out=50)
erf_hispanic_male_rm2<-matching_smooth(hispanic_male_rm2$dead/hispanic_male_rm2$time_count, hispanic_male_rm2$pm25, 
                                       bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_hispanic_male_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_hispanic_male_rm2[,1]<-w.vals
df_hispanic_male_rm2[,2]<-erf_hispanic_male_rm2

#Asian
w.vals=seq(min(asian_female_rm2$pm25),
           max(asian_female_rm2$pm25), length.out=50)
erf_asian_female_rm2<-matching_smooth(asian_female_rm2$dead/asian_female_rm2$time_count, asian_female_rm2$pm25, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_asian_female_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_asian_female_rm2[,1]<-w.vals
df_asian_female_rm2[,2]<-erf_asian_female_rm2

w.vals=seq(min(asian_male_rm2$pm25),
           max(asian_male_rm2$pm25), length.out=50)
erf_asian_male_rm2<-matching_smooth(asian_male_rm2$dead/asian_male_rm2$time_count, asian_male_rm2$pm25, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_asian_male_rm2<-data.frame(matrix(nrow=50, ncol=2))
df_asian_male_rm2[,1]<-w.vals
df_asian_male_rm2[,2]<-erf_asian_male_rm2

par(mfrow=c(4,2))
plot(df_white_female_rm2$X1, df_white_female_rm2$X2,xlab="pm25", ylab="White female deaths", type="l")
plot(df_white_male_rm2$X1, df_white_male_rm2$X2,xlab="pm25", ylab="White male deaths", type="l")
plot(df_black_female_rm2$X1, df_black_female_rm2$X2,xlab="pm25", ylab="Black female deaths", type="l")
plot(df_black_male_rm2$X1, df_black_male_rm2$X2,xlab="pm25", ylab="Black male deaths", type="l")
plot(df_hispanic_female_rm2$X1, df_hispanic_female_rm2$X2,xlab="pm25", ylab="Hispanic female deaths", type="l")
plot(df_hispanic_male_rm2$X1, df_hispanic_male_rm2$X2,xlab="pm25", ylab="Hispanic male deaths", type="l")
plot(df_asian_female_rm2$X1, df_asian_female_rm2$X2,xlab="pm25", ylab="Asian female deaths", type="l")
plot(df_asian_male_rm2$X1, df_asian_male_rm2$X2,xlab="pm25", ylab="Asian male deaths", type="l")

#################GAM
#QD
#matching_qd3<-matching_qd2[matching_qd2$pm25<quantile(matching_qd2$pm25, 0.95) & 
#                             matching_qd2$pm25> quantile(matching_qd2$pm25, 0.05),]
matching_qd3<-subset(matching_qd2, matching_qd2$pm25_ensemble <=15.34556 & matching_qd2$pm25_ensemble>=4.581972)

white_female_qd3<-subset(white_female_qd, white_female_qd$pm25_ensemble<=15.34556 & white_female_qd$pm25_ensemble>=4.581972)
white_male_qd3<-subset(white_male_qd, white_male_qd$pm25_ensemble<=15.34556 & white_male_qd$pm25_ensemble>=4.581972)
black_female_qd3<-subset(black_female_qd, black_female_qd$pm25_ensemble<=15.34556 & black_female_qd$pm25_ensemble>=4.581972)
black_male_qd3<-subset(black_male_qd, black_male_qd$pm25_ensemble<=15.34556 & black_male_qd$pm25_ensemble>=4.581972)
hispanic_female_qd3<-subset(hispanic_female_qd, hispanic_female_qd$pm25_ensemble<=15.34556 & hispanic_female_qd$pm25_ensemble>=4.581972)
hispanic_male_qd3<-subset(hispanic_male_qd, hispanic_male_qd$pm25_ensemble<=15.34556 & hispanic_male_qd$pm25_ensemble>=4.581972)
asian_female_qd3<-subset(asian_female_qd, asian_female_qd$pm25_ensemble<=15.34556 & asian_female_qd$pm25_ensemble>=4.581972)
asian_male_qd3<-subset(asian_male_qd, asian_male_qd$pm25_ensemble<=15.34556 & asian_male_qd$pm25_ensemble>=4.581972)

matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                            as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                            offset(log(time_count))
                          , data=matching_qd3,family=poisson(link="log"))

white_female_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=white_female_qd3,family=poisson(link="log"))

white_male_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                       offset(log(time_count))
                                     , data=white_male_qd3,family=poisson(link="log"))

black_female_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=black_female_qd3,family=poisson(link="log"))

black_male_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                       offset(log(time_count))
                                     , data=black_male_qd3,family=poisson(link="log"))

hispanic_female_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                            as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                            offset(log(time_count))
                                          , data=hispanic_female_qd3,family=poisson(link="log"))

hispanic_male_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                          as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                          offset(log(time_count))
                                        , data=hispanic_male_qd3,family=poisson(link="log"))

asian_female_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                         as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                         offset(log(time_count))
                                       , data=asian_female_qd3,family=poisson(link="log"))

asian_male_matching_gam3 <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                                       as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                                       offset(log(time_count))
                                     , data=asian_male_qd3,family=poisson(link="log"))


plot(matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")
par(mfrow=c(4,2))
plot(white_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="White female HR")
plot(white_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="White male HR")
plot(black_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black female HR")
plot(black_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Black male HR")
plot(hispanic_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic female HR")
plot(hispanic_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Hispanic male HR")
plot(asian_female_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian female HR")
plot(asian_male_matching_gam3, all.terms = FALSE, trans=exp, shift=0.05, ylab="Asian male HR")


matching_gam3_qd <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                            as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                            offset(log(time_count))
                          , data=matching_qd3,family=poisson(link="log"))
matching_gam3_rm <-mgcv::bam(dead~ s(pm25,k=3) +
                            as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                            offset(log(time_count))
                          , data=matching_rm3,family=poisson(link="log"))

par(nfrow=c(1,2))
plot(matching_gam3_qd, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")
plot(matching_gam3_rm, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")

matching_qd4<-matching_qd2[matching_qd2$pm25_ensemble<quantile(matching_qd2$pm25_ensemble, 0.95) & 
                                                          matching_qd2$pm25_ensemble> quantile(matching_qd2$pm25_ensemble, 0.05),]

matching_rm4<-matching_rm2[matching_rm2$pm25<quantile(matching_rm2$pm25, 0.95) & 
                             matching_rm2$pm25> quantile(matching_rm2$pm25, 0.05),]


matching_gam4_qd <-mgcv::bam(dead~ s(pm25_ensemble,k=3) +
                               as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                               offset(log(time_count))
                             , data=matching_qd4,family=poisson(link="log"))
matching_gam4_rm <-mgcv::bam(dead~ s(pm25,k=3) +
                               as.factor(sex)+as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)+
                               offset(log(time_count))
                             , data=matching_rm4,family=poisson(link="log"))

par(mfrow=c(1,2))
plot(matching_gam4_qd, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")
plot(matching_gam4_rm, all.terms = FALSE, trans=exp, shift=0.05, ylab="HR")


test.data<-data.frame(pm25_ensemble = seq(min(matching_qd4$pm25_ensemble), max(matching_qd4$pm25_ensemble), length.out=50), 
                      entry_age_break= rep(levels(as.factor(matching_qd4$entry_age_break))[1], 50),
                      dual = rep(levels(as.factor(matching_qd4$dual))[1], 50),
                      sex = rep(levels(as.factor(matching_qd4$sex))[1], 50),
                      race = rep(levels(as.factor(matching_qd4$race))[1], 50),
                      followup_year= rep(levels(as.factor(matching_qd4$followup_year))[1], 50),
                      time_count= rep(mean(matching_qd4$time_count),50)
                    )

test.data1<-data.frame(pm25_ensemble = seq(min(matching_qd4$pm25_ensemble), max(matching_qd4$pm25_ensemble), length.out=50), 
                       entry_age_break= rep(levels(as.factor(matching_qd4$entry_age_break))[1], 50),
                       dual = rep(levels(as.factor(matching_qd4$dual))[1], 50),
                       followup_year= rep(levels(as.factor(matching_qd4$followup_year))[1], 50),
                       
                       time_count= rep(mean(matching_qd4$time_count),50))
                       
pred.fit<-predict(matching_gam4_qd, newdata=test.data, se=TRUE)
plot(y=exp(pred.fit$fit), x=test.data$pm25_ensemble)



test.data<-data.frame(pm25 = seq(min(matching_rm4$pm25), max(matching_rm4$pm25), length.out=50), 
                      entry_age_break= rep(levels(as.factor(matching_rm4$entry_age_break))[1], 50),
                      dual = rep(levels(as.factor(matching_rm4$dual))[1], 50),
                      sex = rep(levels(as.factor(matching_rm4$sex))[1], 50),
                      race = rep(levels(as.factor(matching_rm4$race))[1], 50),
                      followup_year= rep(levels(as.factor(matching_rm4$followup_year))[1], 50),
                      time_count= rep(mean(matching_rm4$time_count),50)
)
pred.fit<-predict(matching_gam4_rm, newdata=test.data, se=TRUE)
plot(y=exp(pred.fit$fit), x=test.data$pm25)


#QD kernel smoothing
#Kernel smoothinh
#qd2
w.vals=seq(min(matching_qd3$pm25_ensemble),
           max(matching_qd3$pm25_ensemble), length.out=50)
erf_qd2<-matching_smooth(matching_qd2$dead/matching_qd2$time_count, matching_qd2$pm25_ensemble, 
                         bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df<-data.frame(matrix(nrow=50, ncol=2))
df[,1]<-w.vals
df[,2]<-erf_qd2
require(ggplot2)
p<-ggplot(df, aes(x=X1, y=X2))+geom_line()+labs(x="pm25_ensemble", y="deaths/timecount")

#White
w.vals=seq(min(white_female_qd3$pm25_ensemble),
           max(white_female_qd3$pm25_ensemble), length.out=50)
erf_white_female_qd2<-matching_smooth(white_female_qd3$dead/white_female_qd3$time_count, white_female_qd3$pm25_ensemble, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_white_female_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_white_female_qd2[,1]<-w.vals
df_white_female_qd2[,2]<-erf_white_female_qd2

w.vals=seq(min(white_male_qd3$pm25_ensemble),
           max(white_male_qd3$pm25_ensemble), length.out=50)
erf_white_male_qd2<-matching_smooth(white_male_qd3$dead/white_male_qd3$time_count, white_male_qd3$pm25_ensemble, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_white_male_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_white_male_qd2[,1]<-w.vals
df_white_male_qd2[,2]<-erf_white_male_qd2


#Black
w.vals=seq(min(black_female_qd3$pm25_ensemble),
           max(black_female_qd3$pm25_ensemble), length.out=50)
erf_black_female_qd2<-matching_smooth(black_female_qd3$dead/black_female_qd3$time_count, black_female_qd3$pm25_ensemble, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_black_female_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_black_female_qd2[,1]<-w.vals
df_black_female_qd2[,2]<-erf_black_female_qd2

w.vals=seq(min(black_male_qd3$pm25_ensemble),
           max(black_male_qd3$pm25_ensemble), length.out=50)
erf_black_male_qd2<-matching_smooth(black_male_qd3$dead/black_male_qd3$time_count, black_male_qd3$pm25_ensemble, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_black_male_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_black_male_qd2[,1]<-w.vals
df_black_male_qd2[,2]<-erf_black_male_qd2

#Hispanic
w.vals=seq(min(hispanic_female_qd3$pm25_ensemble),
           max(hispanic_female_qd3$pm25_ensemble), length.out=50)
erf_hispanic_female_qd2<-matching_smooth(hispanic_female_qd3$dead/hispanic_female_qd3$time_count, hispanic_female_qd3$pm25_ensemble, 
                                         bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_hispanic_female_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_hispanic_female_qd2[,1]<-w.vals
df_hispanic_female_qd2[,2]<-erf_hispanic_female_qd2

w.vals=seq(min(hispanic_male_qd3$pm25_ensemble),
           max(hispanic_male_qd3$pm25_ensemble), length.out=50)
erf_hispanic_male_qd2<-matching_smooth(hispanic_male_qd3$dead/hispanic_male_qd3$time_count, hispanic_male_qd3$pm25_ensemble, 
                                       bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_hispanic_male_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_hispanic_male_qd2[,1]<-w.vals
df_hispanic_male_qd2[,2]<-erf_hispanic_male_qd2

#Asian
w.vals=seq(min(asian_female_qd3$pm25_ensemble),
           max(asian_female_qd3$pm25_ensemble), length.out=50)
erf_asian_female_qd2<-matching_smooth(asian_female_qd3$dead/asian_female_qd3$time_count, asian_female_qd3$pm25_ensemble, 
                                      bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_asian_female_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_asian_female_qd2[,1]<-w.vals
df_asian_female_qd2[,2]<-erf_asian_female_qd2

w.vals=seq(min(asian_male_qd3$pm25_ensemble),
           max(asian_male_qd3$pm25_ensemble), length.out=50)
erf_asian_male_qd2<-matching_smooth(asian_male_qd3$dead/asian_male_qd3$time_count, asian_male_qd3$pm25_ensemble, 
                                    bw.seq=seq(0.2,2,0.2), w.vals=w.vals)
df_asian_male_qd2<-data.frame(matrix(nrow=50, ncol=2))
df_asian_male_qd2[,1]<-w.vals
df_asian_male_qd2[,2]<-erf_asian_male_qd2

par(mfrow=c(4,2))
plot(df_white_female_qd2$X1, df_white_female_qd2$X2,xlab="pm25_ensemble", ylab="White female deaths", type="l")
plot(df_white_male_qd2$X1, df_white_male_qd2$X2,xlab="pm25_ensemble", ylab="White male deaths", type="l")
plot(df_black_female_qd2$X1, df_black_female_qd2$X2,xlab="pm25_ensemble", ylab="Black female deaths", type="l")
plot(df_black_male_qd2$X1, df_black_male_qd2$X2,xlab="pm25_ensemble", ylab="Black male deaths", type="l")
plot(df_hispanic_female_qd2$X1, df_hispanic_female_qd2$X2,xlab="pm25_ensemble", ylab="Hispanic female deaths", type="l")
plot(df_hispanic_male_qd2$X1, df_hispanic_male_qd2$X2,xlab="pm25_ensemble", ylab="Hispanic male deaths", type="l")
plot(df_asian_female_qd2$X1, df_asian_female_qd2$X2,xlab="pm25_ensemble", ylab="Asian female deaths", type="l")
plot(df_asian_male_qd2$X1, df_asian_male_qd2$X2,xlab="pm25_ensemble", ylab="Asian male deaths", type="l")






