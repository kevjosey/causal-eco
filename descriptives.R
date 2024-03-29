library(matrixStats)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cobalt)
library(fst)

## Person-Year Statistics

## write a function to add a row for continuous variables

add_cont <- function(x, w = rep(1, length(x)), nm_var, nm_level, ndig = 2) {
  
    n_x<-sum(!is.na(x))
    mean_x<-round(weightedMean(x = x, w = w, na.rm = T),ndig)
    sd_x<-round(weightedSd(x = x, w = w, na.rm=T), ndig)
    return(c(nm_var, nm_level, n_x, mean_x, sd_x))
  
}

## write a function to add a row for a level of a categorical variable
add_cat <- function(x, w = rep(1, length(x)), nm_var, nm_level, ndig = 2){

    n_x<-sum(!is.na(x))
    n_x1<-sum(w[x==1], na.rm=T)
    pct_x1<-round(100*n_x1/sum(w), ndig)
    return(c(nm_var, nm_level, n_x, n_x1, pct_x1))
  
}

## write a function to add a death-row
add_death <- function(x, w = rep(1, length(x)), nm_var, ndig = 2){
  
  n_x<-sum(!is.na(x))
  n_x1<-sum(x, na.rm=T)
  pct_x1<-round(100*n_x1/sum(w), ndig)
  return(c(nm_var, nm_level = '', n_x, n_x1, pct_x1))
  
}

# Save Location
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'
load(paste0(dir_data,"aggregate_data.RData"))

# scenarios
scenarios <- expand.grid(dual = c("both","high","low"), race = c("all","white","black","hispanic","asian"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

for (i in 1:nrow(scenarios)){
  
  scenario <- scenarios[i,]
  
  if (scenario$dual == "high") {
    dual0 <- 0
  } else if (scenario$dual == "low") {
    dual0 <- 1
  } else {
    dual0 <- c(0,1)
  }
  
  if (scenario$race == "white") {
    race0 <- 1
  } else if (scenario$race == "black") {
    race0 <- 2
  } else if (scenario$race == "asian") {
    race0 <- 4
  } else if (scenario$race == "hispanic") {
    race0 <- 5
  } else if (scenario$race == "other") {
    race0 <- 3
  } else {
    race0 <- c(0,1,2,3,4,5,6)
  }

  sub_data <- subset(aggregate_data, race %in% race0 & dual %in% dual0)

  zcov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
            "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
  
  x.tmp <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                  model.matrix(~ ., data = sub_data[,zcov])[,-1])[,lapply(.SD, min), by = c("zip", "year", "region")]
  x <- merge(x.tmp, data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                               y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region")],
             by = c("zip", "year", "region"))
  w <- data.table(zip = sub_data$zip, year = sub_data$year, region = sub_data$region,
                  dual = sub_data$dual, race = sub_data$race, sex = sub_data$sex, age_break = sub_data$age_break, 
                  y = sub_data$dead, n = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year", "region", "dual", "race", "sex", "age_break")]
  
  # Table 1 -----------------------------------------------------------------
  
  ## sex ##
  table1<-c(add_cat(x = as.numeric(w$sex==1), w = w$n, nm_var='Female', nm_level=''))
  
  ## region ##
  table1<-rbind(table1, add_cat(x = as.numeric(w$region == "MIDWEST"),
                                w = w$n, nm_var='Region', nm_level='Midwest'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$region == "NORTHEAST"),
                                w = w$n, nm_var='Region', nm_level='Northeast'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$region == "SOUTH"),
                                w = w$n, nm_var='Region', nm_level='South'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$region == "WEST"),
                                w = w$n, nm_var='Region', nm_level='West'))
  
  ## age ##
  table1<-rbind(table1, add_cat(x = as.numeric(w$age_break == "[65,75)"),
                                w = w$n, nm_var =' Age', nm_level = '65-74'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$age_break == "[75,85)"),
                                w = w$n, nm_var = 'Age', nm_level = '75-84'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$age_break == "[85,95)"),
                                w = w$n, nm_var = 'Age', nm_level = '85-94'))
  table1<-rbind(table1, add_cat(x = as.numeric(w$age_break == "[95,125)"),
                                w = w$n, nm_var = 'Age', nm_level = '95+'))
  
  ## race ##
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$race==1),
                                 w = w$n, nm_var='Race',nm_level='White')))
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$race==2),
                                 w = w$n, nm_var='Race',nm_level='Black')))
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$race==5),
                                 w = w$n, nm_var='Race',nm_level='Hispanic')))
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$race==4),
                                 w = w$n, nm_var='Race',nm_level='Asian')))
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$race==3),
                                 w = w$n, nm_var='Race',nm_level='Other or Unknown')))
  
  ## dual eligible ##
  table1<-rbind(table1,c(add_cat(x = as.numeric(w$dual==1), 
                                 w = w$n, nm_var='Medicaid Eligible',nm_level='')))
  
  ## deaths ##
  table1<-rbind(table1,c(add_death(x = as.numeric(w$y), w = w$n, nm_var = "Deaths")))
  table1<-rbind(table1,c(add_death(x = as.numeric(w$n), w = w$n, nm_var = "Person-years")))
  
  colnames(table1) <- c("Variable", "Level", "Rows", "N", "%")
  
  write_csv(data.frame(table1), path = paste0('~/Tables/table1_',scenario$dual, "_", scenario$race, ".csv"))
  
  # Table 2 -----------------------------------------------------------------
  
  ## demographics and health ##
  table2 <- c(add_cont(x = as.numeric(x$pm25), ndig=4,
                       w = x$n, nm_var='PM2.5', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$mean_bmi), ndig=4,
                                   w = x$n, nm_var='Average BMI', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$smoke_rate), ndig=4,
                                   w = x$n, nm_var='Smoking Rate', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$hispanic), ndig=4,
                                   w = x$n, nm_var='Percent Hispanic', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$pct_blk), ndig=4,
                                   w = x$n, nm_var='Percent Black', nm_level=''))
  
  ## socioeconomic variables ##
  table2 <- rbind(table2, add_cont(x = as.numeric(x$popdensity), ndig=4,
                                   w = x$n, nm_var='Population Density', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$medhouseholdincome), ndig=4,
                                   w = x$n, nm_var='Median Household Income', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$medianhousevalue), ndig=4,
                                   w = x$n, nm_var='Median House Value', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$poverty), ndig=4,
                                   w = x$n, nm_var='Poverty', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$education), ndig=4,
                                   w = x$n, nm_var='Less than High-School', nm_level=''))
  
  table2 <- rbind(table2, add_cont(x = as.numeric(x$pct_owner_occ), ndig=4,
                                   w = x$n, nm_var='Percenent Owner Occupied', nm_level=''))
  
  ## temperature/humidity ##
  table2 <- rbind(table2, add_cont(x = as.numeric(x$summer_tmmx), ndig=4,
                                   w = x$n, nm_var='Summer Temperature', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$winter_tmmx), ndig=4,
                                   w = x$n, nm_var='Winter Temperature', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$summer_rmax), ndig=4,
                                   w = x$n, nm_var='Summer Humidity', nm_level=''))
  table2 <- rbind(table2, add_cont(x = as.numeric(x$winter_rmax), ndig=4,
                                   w = x$n, nm_var='Winter Humidity', nm_level=''))
  
  colnames(table2) <- c("Variable", "Level", "Rows", "Mean", "SD")

  write_csv(data.frame(table2), path = paste0('~/Tables/table2_',scenario$dual, "_", scenario$race, ".csv"))
  
}

## Covariate Balance Plot

# scenarios
a.vals <- seq(4, 16, length.out = 121)
set.seed(42)

bal_dat <- function(a, x, weights){
  
  val <- bal.tab(x, treat = a, weights = weights, method = "weighting", continuous = "std", s.d.denom = "all")
  bal_df <- val$Balance
  labs <- rep(rownames(bal_df), 2)
  vals_tmp <- bal_df$Corr.Adj
  vals_year <- mean(abs(vals_tmp[1:17]))
  vals_region <- mean(abs(vals_tmp[18:21]))
  vals_tmp2 <- c(abs(vals_tmp[-c(1:21)]),
                 vals_year, vals_region)
  names(vals_tmp2) <- c("Calendar Year", "Region", "Mean BMI", "Smoking Rate", "% Hispanic", "% Black",
                        "Median Household Income", "Median House Value", "% Below Poverty Level",
                        "% Below High School Education", "Population Density", "% Owner-Occupied Housing",
                        "Summer Temperature","Winter Temperature", "Summer Humidity", "Winter Humidity")
  vals <- vals_tmp2[order(vals_tmp2, decreasing = TRUE)]
  adjust <- rep("Unadjusted", each = length(vals_tmp2))
  labs <- names(vals)
  df <- data.frame(labs = labs, vals = vals, adjust = adjust)
  df$labs <- factor(df$labs, levels = rev(names(vals)))
  
  return(df)
  
}

# Standalone Plots + ESS
load(paste0(dir_out, "all.RData"))
wx <- new_data$wx[,c(3:4,8:21)]

bdat_1 <- bal_dat(a = new_data$wx$pm25, x = wx, weights = new_data$wx$n)
bdat_2 <- bal_dat(a = new_data$wx$pm25, x = wx, weights = new_data$wx$ipw*new_data$wx$n)
bdat_3 <- bal_dat(a = new_data$wx$pm25, x = wx, weights = new_data$wx$trunc*new_data$wx$n)

bdat_2$adjust <- "Calibration"
bdat_3$adjust <- "Truncated Calibration"

df <- rbind(bdat_1, bdat_2, bdat_3)
df$adjust <- factor(df$adjust, levels = c("Unadjusted", "Calibration", "Truncated Calibration"))

# Effective Sample Size
ess_1 <- sum(new_data$wx$ipw*new_data$wx$n)^2/sum((new_data$wx$ipw*new_data$wx$n)^2)
ess_2 <- sum(new_data$wx$trunc*new_data$wx$n)^2/sum((new_data$wx$trunc*new_data$wx$n)^2)
ess <- data.frame(cal = ess_1, cal_trunc = ess_2)

balance_plot <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.35) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file = paste0("~/Figures/balance_plot.pdf"), width = 8, height = 8)
balance_plot
dev.off()

write_csv(ess, file = "~/Tables/ess.csv")

## Big Table

f <- list.files("/n/dominici_nsaph_l3/Lab/data/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst",
                full.names = TRUE)

myvars <- c("qid", "year","zip","sex","race","age","dual","entry_age_break","statecode",
            "followup_year","followup_year_plus_one","dead","pm25_ensemble",
            "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
            "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")

national_merged2016 <- rbindlist(lapply(f, read_fst, columns = myvars, as.data.table = TRUE))
national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)

NORTHEAST = c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")
SOUTH = c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST = c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST = c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

# creates region
national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST",
                                     ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                            ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                   ifelse(national_merged2016$state %in% WEST, "WEST", NA))))

national_merged2016 <- national_merged2016[complete.cases(national_merged2016[,c(1:28)]) ,]

national_merged2016$time_count <- rep(1, nrow(national_merged2016))
national_merged2016$age_break <- cut(national_merged2016$age, c(65,75,85,95,125), right = FALSE)
national_merged2016$sex <- national_merged2016$sex - 1

# label unknown and american/alaska natives as "other"
national_merged2016$race[national_merged2016$race == 6] <- 3
national_merged2016$race[national_merged2016$race == 0] <- 3

# unique persons
qid_data <- national_merged2016 %>% distinct(qid, .keep_all = TRUE)

# Individual-level data
national_merged2016_sub <- subset(national_merged2016, race == 4 & dual == 1) # need to change based on dual and race status
qid_data_sub <- subset(qid_data, race == 4 & dual == 1) # need to change based on dual and race status

nrow(qid_data_sub) # strata total at-risk
nrow(qid_data_sub)/nrow(qid_data)*100 # strata proportion of at-risk
nrow(national_merged2016_sub) # person-years at risk
nrow(national_merged2016_sub)/nrow(national_merged2016)*100 # proportion of person-years
sum(national_merged2016_sub$dead) # number of deaths
sum(national_merged2016_sub$dead)/sum(national_merged2016$dead)*100 # strata number of deaths
table(qid_data_sub$age_break)/nrow(qid_data_sub)*100 # age break proportions
table(qid_data_sub$sex)/nrow(qid_data_sub)*100 # female/male propotions
table(qid_data_sub$dual)/nrow(qid_data_sub)*100 # Medicaid eligibility

# calculate median follow-up year
followup_data_sub <- national_merged2016_sub %>% group_by(qid) %>% summarise(max_followup = max(followup_year))
median(followup_data_sub$max_followup)
