library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white", "black", "asian", "hispanic"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals = seq(4, 16, length.out = 121)

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

### Create Data

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_rti.RData"))
  
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
                        dual = rep(scenario$dual, nrow(est_data)),
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
                       dual = rep(scenario$dual, nrow(est_data)))

  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

save(dat, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/estimate_eco.RData')
save(rr, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/relative_eco.RData')

### Main Plot

# histogram and ERF data
load(paste0(dir_out, "both_all_rti.RData"))
rr_tmp <- subset(rr, dual == "both" & race == "all" & pm0 == 12)
a_dat <- rep(new_data$wx$pm25, new_data$wx$n)

# exposure response curve
erf_plot <- rr_tmp %>% 
  ggplot(aes(x = a.vals)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  geom_hline(yintercept = 1, color = "dodgerblue3", size = 1.1, linetype = "dotted") +
  coord_cartesian(xlim = c(5,15), ylim = c(0.925,1.05)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio",
       title = "Exposure Response Curve for\n All Medicare Recipients (2009-2014)") + 
  scale_y_continuous(breaks = c(0.925,0.95,0.975,1.0,1.025, 1.05)) +
  scale_x_continuous(breaks = c(6,8,10,12,14)) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  grids(linetype = "dashed")

# histogram
a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a)) + 
  geom_density(fill = "darkorange3", alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(5,15), ylim = c(0.0,0.2)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15, 0.2)) +
  scale_x_continuous(breaks = c(6,8,10,12,14)) +
  guides(fill = "none") +
  theme_cowplot()

align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
main_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "~/Figures/erc_plot.pdf", width = 8, height = 8)
main_plot
dev.off()

### ERC by Race

plot_list <- list()
dual.vals <- c("both", "high", "low")

for (i in 1:length(dual.vals)){
  
  if (dual.vals[i] == "low") {
    main <- "Lower Income Only"
  } else if (dual.vals[i] == "high") {
    main <- "Higher Income Only"
  } else {
    main <- "All Participants"
  }
  
  rr_tmp <- subset(rr, dual == dual.vals[i] & race != "all" & pm0 == 12)
  rr_tmp$race <- str_to_title(rr_tmp$race)
  ylim <- c(0.8,1.1)
  
  rr_tmp$race <- factor(rr_tmp$race)
    
  # dual ineligible ERCs
  erf_strata_tmp <- rr_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = race)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate)) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Race", title = main) +
    scale_color_manual(values = c("#75bad3", "#489f8c","#ea3323","#ea8832")) +
    scale_y_continuous(breaks = seq(from = 0.8, to = 1.1, length.out = 7)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) 
  
  plot_list[[i]] <- erf_strata_tmp
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

pdf(file = "~/Figures/strata_plot.pdf", width = 12, height = 6)
strata_plot
dev.off()

### Contrast Plot

contr <- subset(dat, a.vals == 8)

contrast_plot <- contr %>% 
  ggplot(aes(x = str_to_upper(race), y = 100*excess/deaths, color = str_to_upper(dual))) + 
  geom_pointrange(aes(ymin = 100*lower.ed/deaths, ymax = 100*upper.ed/deaths), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Percent of Deaths Avoidable (%)", 
       title = ~"Excess Death Estimates from Reducing "*PM[2.5]*" to 8 "*mu*g*"/"*m^3,
       color = "Socioeconomic Position") +
  theme(legend.position = c(0.12, 0.85),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#008080", "#FF00FF","#FFD700")) +
  grids(linetype = "dashed")

pdf(file = "~/Figures/contrast_plot.pdf", width = 8, height = 8)
contrast_plot
dev.off()