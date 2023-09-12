library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# scenarios
scenarios <- expand.grid(dual = c("high", "low","both"), race = c("white","black","hispanic","asian","all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(4, 16, length.out = 121)

# data directories
dir_data <- '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/qd/'
dir_out <- '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data/'

tea_8 <- tea_9 <- tea_10 <- tea_11 <- tea_12 <- data.frame()

# contrast indexes
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))

### Predict Total Events Avoided

## Under Absolute Risks
for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_both_all.RData"))
  est_data <- new_data$est_data
  wx <- new_data$wx
  
  tea_est_8 <- with(subset(wx, pm25 > 8), sum(y - n*as.numeric(est_data$estimate[idx8])))
  tea_est_9 <- with(subset(wx, pm25 > 9), sum(y - n*as.numeric(est_data$estimate[idx9])))
  tea_est_10 <- with(subset(wx, pm25 > 10), sum(y - n*as.numeric(est_data$estimate[idx10])))
  tea_est_11 <- with(subset(wx, pm25 > 11), sum(y - n*as.numeric(est_data$estimate[idx11])))
  tea_est_12 <- with(subset(wx, pm25 > 12), sum(y - n*as.numeric(est_data$estimate[idx12])))
  
  tea_var_8 <- with(subset(wx, pm25 > 8), sum(pm25 > 8)*(est_data$se[idx8]*n)^2)
  tea_var_9 <- with(subset(wx, pm25 > 9), sum(pm25 > 9)*(est_data$se[idx9]*n)^2)
  tea_var_10 <- with(subset(wx, pm25 > 10), sum(pm25 > 10)*(est_data$se[idx10]*n)^2)
  tea_var_11 <- with(subset(wx, pm25 > 11), sum(pm25 > 11)*(est_data$se[idx11]*n)^2)
  tea_var_12 <- with(subset(wx, pm25 > 12), sum(pm25 > 12)*(est_data$se[idx12]*n)^2)
  
  tea_lower_8 <- tea_est_8 - 1.96*sqrt(tea_var_8)
  tea_lower_9 <- tea_est_9 - 1.96*sqrt(tea_var_9)
  tea_lower_10 <- tea_est_10 - 1.96*sqrt(tea_var_10)
  tea_lower_11 <- tea_est_11 - 1.96*sqrt(tea_var_11)
  tea_lower_12 <- tea_est_12 - 1.96*sqrt(tea_var_12)
  
  tea_upper_8 <- tea_est_8 + 1.96*sqrt(tea_var_8)
  tea_upper_9 <- tea_est_9 + 1.96*sqrt(tea_var_9)
  tea_upper_10 <- tea_est_10 + 1.96*sqrt(tea_var_10)
  tea_upper_11 <- tea_est_11 + 1.96*sqrt(tea_var_11)
  tea_upper_12 <- tea_est_12 + 1.96*sqrt(tea_var_12)
  
  n_8 <- with(subset(wx, pm25 > 8), sum(y))
  n_9 <- with(subset(wx, pm25 > 9), sum(y))
  n_10 <- with(subset(wx, pm25 > 10), sum(y))
  n_11 <- with(subset(wx, pm25 > 11), sum(y))
  n_12 <- with(subset(wx, pm25 > 12), sum(y))
  
  tea_8 <- rbind(tea_8, data.frame(dual = scenario$dual, race = scenario$race,
                                   tea_8 = tea_est_8, n_8 = n_8, prop_8 = tea_est_8/n_8,
                                   tea_lower_8 = tea_lower_8, tea_upper_8 = tea_upper_8))
  tea_9 <- rbind(tea_9, data.frame(dual = scenario$dual, race = scenario$race,
                                   tea_9 = tea_est_9, n_9 = n_9, prop_9 = tea_est_9/n_9, 
                                   tea_lower_9 = tea_lower_9, tea_upper_9 = tea_upper_9))
  tea_10 <- rbind(tea_10, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_10 = tea_est_10, n_10 = n_10, prop_10 = tea_est_10/n_10, 
                                     tea_lower_10 = tea_lower_10, tea_upper_10 = tea_upper_10))
  tea_11 <- rbind(tea_11, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_11 = tea_est_11, n_11 = n_11, prop_11 = tea_est_11/n_11, 
                                     tea_lower_11 = tea_lower_11, tea_upper_11 = tea_upper_11))
  tea_12 <- rbind(tea_12, data.frame(dual = scenario$dual, race = scenario$race,
                                     tea_12 = tea_est_12, n_12 = n_12, prop_12 = tea_est_12/n_12,
                                     tea_lower_12 = tea_lower_12, tea_upper_12 = tea_upper_12))
  
}

save(tea_8, tea_9, tea_10, tea_11, tea_12, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/tea.RData')