library(parallel)
library(data.table)
library(tidyr)
library(dplyr)
library(scam)
library(sandwich)

source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/calibrate.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_dr.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/gam_om.R')
source('/n/dominici_nsaph_l3/projects/kjosey-erc-strata/causal-eco/Functions/bootstrap.R')
set.seed(42)

# scenarios
scenarios <- expand.grid(race = c("all", "white","black","hispanic","asian"),
                         region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))
scenarios$race <- as.character(scenarios$race)
scenarios$region <- as.character(scenarios$region)
a.vals <- seq(4, 16, length.out = 121)

# data directories
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'

dat <- data.frame()
contr <- data.frame()
rr <- data.frame()

# contrast indexes
idx5 <- which.min(abs(a.vals - 5))
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))
idx15 <- which.min(abs(a.vals - 15))

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$race, "_", scenario$region, "_rti.RData"))
  
  u.zip <- unique(new_data$wx$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  n <- length(u.zip)
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  
  # asymptotics
  estimate <- est_data$estimate
  # se <- est_data$se
  excess_est <- excess_death$estimate
  # excess_se <- excess_death$se
  
  # bootstrapping
  boot_erc <- new_data$boot_erc
  boot_ed <- new_data$boot_ed
  # estimate <- colMeans(boot_erc, na.rm = T)
  se <- sqrt(m/n)*apply(boot_erc, 2, sd, na.rm = T)
  # excess_est <- colMeans(boot_ed, na.rm = T)
  excess_se <- sqrt(m/n)*apply(boot_ed, 2, sd, na.rm = T)
  
  ## Absoilute Risk
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals),
                        estimate = estimate,
                        excess = excess_est,
                        lower = estimate - 1.96*se,
                        upper = estimate + 1.96*se,
                        lower.ed = excess_est - 1.96*excess_se,
                        upper.ed = excess_est + 1.96*excess_se,
                        race = rep(scenario$race, nrow(est_data)),
                        region = rep(scenario$region, nrow(est_data)),
                        deaths = rep(sum(new_data$wx$y), nrow(est_data)))
  
  ## Relative Risks
  rr_tmp_8 <- c(as.numeric(estimate)/as.numeric(estimate[idx8]))
  rr_tmp_9 <- c(as.numeric(estimate)/as.numeric(estimate[idx9]))
  rr_tmp_10 <- c(as.numeric(estimate)/as.numeric(estimate[idx10]))
  rr_tmp_11 <- c(as.numeric(estimate)/as.numeric(estimate[idx11]))
  rr_tmp_12 <- c(as.numeric(estimate)/as.numeric(estimate[idx12]))
  
  # delta method standard errors
  log_rr_se_8 <- sqrt(c(se^2)/c(estimate^2) + c(se[idx8]^2)/c(estimate[idx8]^2))
  log_rr_se_9 <- sqrt(c(se^2)/c(estimate^2) + c(se[idx9]^2)/c(estimate[idx9]^2))
  log_rr_se_10 <- sqrt(c(se^2)/c(estimate^2) + c(se[idx10]^2)/c(estimate[idx10]^2))
  log_rr_se_11 <- sqrt(c(se^2)/c(estimate^2) + c(se[idx11]^2)/c(estimate[idx11]^2))
  log_rr_se_12 <- sqrt(c(se^2)/c(estimate^2) + c(se[idx12]^2)/c(estimate[idx12]^2))
  
  rr_tmp <- data.frame(a.vals = rep(est_data$a.vals, 5),
                       estimate = c(rr_tmp_8, rr_tmp_9, rr_tmp_10, rr_tmp_11, rr_tmp_12),
                       lower = c(exp(log(rr_tmp_8) - 1.96*log_rr_se_8), 
                                 exp(log(rr_tmp_9) - 1.96*log_rr_se_9), 
                                 exp(log(rr_tmp_10) - 1.96*log_rr_se_10),
                                 exp(log(rr_tmp_11) - 1.96*log_rr_se_11),
                                 exp(log(rr_tmp_12) - 1.96*log_rr_se_12)),
                       upper = c(exp(log(rr_tmp_8) + 1.96*log_rr_se_8), 
                                 exp(log(rr_tmp_9) + 1.96*log_rr_se_9), 
                                 exp(log(rr_tmp_10) + 1.96*log_rr_se_10),
                                 exp(log(rr_tmp_11) + 1.96*log_rr_se_11),
                                 exp(log(rr_tmp_12) + 1.96*log_rr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       race = rep(scenario$race, nrow(est_data)),
                       region = rep(scenario$region, nrow(est_data)))
  
  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

### All Subjects

# factor race
dat_tmp <- subset(rr, pm0 == 12 & race == "all")
dat_tmp$region <- str_to_title(dat_tmp$region)

# graph breaks
ylim <- c(0.85, 1.1)
breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)

# dual ineligible + eligible
erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, color = region)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  coord_cartesian(xlim = c(5,15), ylim = ylim) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
       color = "Region", title = "Mortality by Census Region") + 
  scale_y_continuous(breaks = breaks) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

erf_plot

### Plot by Region

plot_list_race <- list()
situations <- expand.grid(region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  main <- str_to_title(situation)
  
  # factor race
  dat_tmp <- subset(rr, region == situation & pm0 == 12 & race != "all")
  dat_tmp$race <- str_to_title(dat_tmp$race)
  
  # graph breaks
  ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
            max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  ylim <- c(0.8, 1.25)
  breaks <- round(seq(ylim[1], ylim[2], length.out = 7), 4)
  
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = race)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate)) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Race", title = main) + 
    scale_color_manual(values = c("#75bad3", "#489f8c","#ea3323","#ea8832")) +
    scale_y_continuous(breaks = breaks) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list_race[[i]] <- erf_strata_plot
  
}

strata_plot_race <- ggarrange(plotlist = plot_list_race[1:4], ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

pdf(file = "~/Figures/strata_plot.pdf", width = 12, height = 6)
strata_plot_race
dev.off()

#### Plot by Race

plot_list_region <- list()
situations <- expand.grid(race = c("asian","black", "white", "hispanic"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  main <- str_to_title(situation)
  
  # factor race
  dat_tmp <- subset(rr, race == situation & pm0 == 12)
  dat_tmp$region <- str_to_title(dat_tmp$region)
  
  # graph breaks
  # ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
  #           max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  ylim <- c(0.8, 1.25)
  breaks <- round(seq(ylim[1], ylim[2], length.out = 7), 4)
  
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = region)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate)) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Region", title = main) + 
    scale_y_continuous(breaks = breaks) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list_region[[i]] <- erf_strata_plot
  
}

strata_plot_region <- ggarrange(plotlist = plot_list_region[1:4], ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

pdf(file = "~/Figures/strata_plot_region.pdf", width = 12, height = 6)
strata_plot_region
dev.off()