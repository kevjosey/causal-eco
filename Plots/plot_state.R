library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scam)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_ipw.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/bootstrap.R')
set.seed(42)

dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
load(paste0(dir_data,"aggregate_data.RData"))

# scenarios
states <- as.character(unique(aggregate_data$statecode))
scenarios <- expand.grid(state = states, dual = c("both", "high", "low"))
scenarios$state <- as.character(scenarios$state)
scenarios$dual <- as.character(scenarios$dual)

# data directories
dir_erc = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF_RM/'
dir_srf = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_SRF_RM/'

erc_dat <- data.frame()
si_dat <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  ## ERC_output
  load(paste0(dir_erc, scenario$state, "_", scenario$dual, ".RData"))
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  wx <- new_data$wx
  
  a.vals <- new_data$est_data$a.vals
  erc_est <- est_data$estimate
  erc_se <- est_data$se
  
  erc_tmp <- data.frame(a.vals = c(est_data$a.vals),
                        estimate = erc_est,
                        lower = erc_est - 1.96*erc_se,
                        upper = erc_est + 1.96*erc_se,
                        dual = rep(scenario$dual, nrow(est_data)),
                        state = rep(scenario$state, nrow(est_data)),
                        deaths = rep(sum(new_data$wx$y), nrow(est_data)))
  
  ## SI Output
  load(paste0(dir_srf, scenario$state, "_", scenario$dual, ".RData"))
  
  delta <- est_data$delta
  si_est <- est_data$est
  si_se <- est_data$se
  
  si_tmp <- data.frame(delta = delta,
                       si_est = si_est,
                       si_excess = sum(wx$n)*si_est,
                       si_lower = si_est - 1.96*si_se,
                       si_upper = si_est + 1.96*si_se,
                       si_lower_excess = sum(wx$n)*(si_est - 1.96*si_se),
                       si_upper_excess = sum(wx$n)*(si_est + 1.96*si_se),
                       dual = rep(scenario$dual, nrow(est_data)),
                       state = rep(scenario$state, nrow(est_data)))
  
  erc_dat <- rbind(erc_dat, erc_tmp)
  si_dat <- rbind(si_dat, si_tmp)
  
}


### Contrast State Plot (delta = 9)

contr_state <- subset(si_dat, dual == "both" & delta == 6)
idx <- order(contr_state$si_est, decreasing = TRUE)
contr_state$state <- factor(contr_state$state, levels = contr_state$state[idx])

contrast_state_plot <- contr_state %>%
  ggplot(aes(x = state, y = 100*si_est)) +
  geom_pointrange(aes(ymin = 100*si_lower, ymax = 100*si_upper), 
                  position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  labs(x = ~ PM[2.5]*" Cutoffs", y = "% Reduction to Mortality", color = "Region") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "bottom",
        legend.key.height = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file = "~/Figures/contrast_state_plot_9.pdf", width = 10, height = 8)
contrast_state_plot
dev.off()
