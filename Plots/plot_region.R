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

# scenarios
scenarios <- expand.grid(dual = c("both", "high","low"), region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST", "US"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$region <- as.character(scenarios$region)

# data directories
dir_erc = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'
dir_srf = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_SRF/'

erc_dat <- data.frame()
si_dat <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  ## ERC_output
  load(paste0(dir_erc, scenario$region, "_", scenario$dual, ".RData"))
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
                        region = rep(scenario$region, nrow(est_data)),
                        deaths = rep(sum(new_data$wx$y), nrow(est_data)))
  
  ## SI Output
  load(paste0(dir_srf, scenario$region, "_", scenario$dual, ".RData"))
  
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
                       region = rep(scenario$region, nrow(est_data)))
  
  erc_dat <- rbind(erc_dat, erc_tmp)
  si_dat <- rbind(si_dat, si_tmp)
  
}

### Plot by Region
dat_ref <- subset(erc_dat, dual == "both" & region == "US")
dat_ref <- dat_ref[rep(seq_len(nrow(dat_ref)), 4), ]
dat_tmp <- subset(erc_dat, dual == "both" & region != "US")
dat_tmp$region2 <- dat_tmp$region
dat_ref$region2 <- rep(c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"), each = length(a.vals))
dat_tmp <- rbind(dat_tmp, dat_ref)

dat_tmp$region <- str_to_title(dat_tmp$region)
dat_tmp$region[dat_tmp$region == "Us"] <- "US"

# graph breaks
ylim <- c(0.03, 0.06)
ybreaks <- round(seq(ylim[1], ylim[2], length.out = 7), 4)

# dual ineligible + eligible
erf_region_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, color = region)) + 
  facet_wrap(~ region2, nrow = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.03,0.06)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
       color = "Region", title = "Mortality by Census Region") + 
  scale_y_continuous(breaks = ybreaks) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "~/Figures/erf_region_plot.pdf", width = 8, height = 6)
erf_region_plot
dev.off()

### Plot by Dual
dat_dual <- subset(erc_dat, region == "US")
dat_dual$dual <- ifelse(dat_dual$dual == "low", "Lower Income",
                        ifelse(dat_dual$dual == "high", "Higher Income", 
                               "All Medicare Recipients"))

# graph breaks
ylim <- c(0.03, 0.09)
ybreaks <- round(seq(ylim[1], ylim[2], length.out = 7), 4)


# dual ineligible + eligible
erf_dual_plot <- dat_dual %>% 
  ggplot(aes(x = a.vals, color = dual)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.03,0.09)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
       color = "Region", title = "Mortality by Census Region") + 
  scale_y_continuous(breaks = ybreaks) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "~/Figures/erf_dual_plot.pdf", width = 8, height = 6)
erf_dual_plot
dev.off()

### Contrast Region Plot
contr_region <- subset(si_dat, dual == "both")
contr_region$region <- str_to_title(contr_region$region)
contr_region$region[contr_region$region == "Us"] <- "US"
xbreaks = 6:12

contrast_region_plot <- contr_region %>%
  ggplot(aes(x = delta, y = 100*si_est, color = region)) +
  geom_pointrange(aes(ymin = 100*si_lower, ymax = 100*si_upper), 
                  position = position_dodge(width = 0.4)) +
  geom_line(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  labs(x = ~ PM[2.5]*" Cutoffs", y = "% Reduction to Mortality", color = "Region") +
  scale_x_continuous(breaks = xbreaks) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file = "~/Figures/contrast_region_plot.pdf", width = 10, height = 8)
contrast_region_plot
dev.off()

### Contrast Dual Plot

contr_dual <- subset(si_dat, region == "US")
contr_dual$dual <- ifelse(contr_dual$dual == "low", "Lower Income",
                          ifelse(contr_dual$dual == "high", "Higher Income", 
                                 "All Medicare Recipients"))

contrast_dual_plot <- contr_dual %>%
  ggplot(aes(x = delta, y = 100*si_est, color = dual)) +
  geom_pointrange(aes(ymin = 100*si_lower, ymax = 100*si_upper), 
                  position = position_dodge(width = 0.4)) +
  geom_line(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  labs(x = ~ PM[2.5]*" Cutoffs", y = "% Reduction to Mortality", color = "Dual Medicaid") +
  scale_x_continuous(breaks = xbreaks) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file = "~/Figures/contrast_dual_plot.pdf", width = 10, height = 8)
contrast_dual_plot
dev.off()