library("wCorr")
library("parallel")
library(fst)
library(data.table)
library("xgboost")
require(polycor)
require(dplyr)
require(tidyr)

#White female
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_female",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_female_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_female/cor_matched.rds")


#White male
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_male",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_white_male_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/white_male/cor_matched.rds")


#Black female
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_female",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_female_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_female/cor_matched.rds")


#Black male
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_male",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_black_male_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/black_male/cor_matched.rds")



#Hispanic female
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/hispanic_female",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_hispanic_female_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/hispanic_female/cor_matched.rds")




#Asian female
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_female",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_female_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_female/cor_matched.rds")


#Asian male
f <- list.files("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_male",
                pattern = "\\.rds",
                full.names = TRUE)

matching_qd2<-rbindlist(lapply(f, readRDS), use.names = FALSE)
matching_qd2<-as.data.frame(matching_qd2)
matching_qd2<-matching_qd2[complete.cases(matching_qd2),]

load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_asian_male_qd.RData")
matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(covariates_black_male_qd$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(covariates_black_male_qd$pm25_ensemble, 0.01))

matching_qd2 <- subset(matching_qd2[complete.cases(matching_qd2), ],
                       pm25_ensemble < quantile(matching_qd2$pm25_ensemble, 0.99) &
                         pm25_ensemble > quantile(matching_qd2$pm25_ensemble, 0.01))

cor_matched2 <- c(abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$year_fac)),
                  abs(polyserial(matching_qd2$pm25_ensemble, matching_qd2$region)),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$mean_bmi, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$smoke_rate, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$hispanic, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_blk, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medhouseholdincome, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$medianhousevalue, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$poverty, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$education, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$popdensity, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$pct_owner_occ, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$summer_rmax, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_tmmx, method = c("spearman"))),
                  abs(cor(matching_qd2$pm25_ensemble, matching_qd2$winter_rmax, method = c("spearman"))))
cor_matched2
summary(cor_matched2)
saveRDS(cor_matched2, file="/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/CB_qd/asian_male/cor_matched.rds")
