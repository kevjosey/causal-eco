library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(cobalt)

# scenarios
scenarios <- expand.grid(dual = c("low", "high", "both"), race = c("all","white", "black", "asian", "hispanic"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(4, 16, length.out = 121)

# data directories
dir_out_qd = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data/'
dir_out_rm = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data_RM/'

### Exposure Assessment Sensitivity Plot

scenario <- scenarios[3,]

# QD
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_both_all.RData"))
dat_qd <- data.frame(a.vals = c(new_data$est_data$a.vals), 
                     estimate = c(new_data$est_data[,6]),
                     lower = c(new_data$est_data[,6] - 1.96*new_data$est_data[,7]),
                     upper = c(new_data$est_data[,6] + 1.96*new_data$est_data[,7]),
                     exposure = rep("Di et al. (2019)", nrow(new_data$est_data)),
                     race = rep(scenario$race, nrow(new_data$est_data)),
                     dual = rep(scenario$dual, nrow(new_data$est_data)))

a_dat_sens <- data.frame(a = rep(new_data$wx$pm25, new_data$wx$time_count), exposure = "Di et al. (2019)")

# RM
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_both_all.RData"))
dat_rm <- data.frame(a.vals = c(new_data$est_data$a.vals),
                     estimate = c(new_data$est_data[,6]),
                     lower = c(new_data$est_data[,6] - 1.96*new_data$est_data[,7]),
                     upper = c(new_data$est_data[,6] + 1.96*new_data$est_data[,7]),
                     exposure = rep("van Donkelaar et al. (2016)", nrow(new_data$est_data)),
                     race = rep(scenario$race, nrow(new_data$est_data)),
                     dual = rep(scenario$dual, nrow(new_data$est_data)))

a_dat_sens <- rbind(a_dat_sens, data.frame(a = rep(new_data$wx$pm25, new_data$wx$time_count), exposure = "van Donkelaar et al. (2016)"))

# combine
dat_qd_rm <- rbind(dat_qd, dat_rm)

erf_plot_sens <- dat_qd_rm %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       color = "Exposure Assessment", title = "Exposure Response Curve under\n Different Exposure Assessments") + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.044,0.045,0.046,0.047,0.048,0.049)) +
  grids(linetype = "dashed")

leg_plot <- dat_qd_rm %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       color = "Exposure Assessment") +
  theme_bw() +
  theme(legend.background = element_rect(colour = "black"))

leg <- gtable_filter(ggplot_gtable(ggplot_build(leg_plot)), "guide-box")

a_hist_sens <- ggplot(a_dat_sens, mapping = aes(x = a, fill = exposure)) + 
  geom_density(alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
  theme(panel.grid = element_blank()) +
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
  guides(fill="none") +
  theme_cowplot()

align <- align_plots(a_hist_sens, erf_plot_sens +
                       annotation_custom(leg, xmin = 6, xmax = 7.5, ymin = 0.048, ymax = 0.049), 
                     align = "hv", axis = "tblr")
sens_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "~/Figures/sens_plot.pdf", width = 8, height = 8)
sens_plot
dev.off()
