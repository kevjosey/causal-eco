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
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(4, 16, length.out = 121)

# data directories
dir_data_qd <- '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/'
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'

dat <- data.frame()
contr <- data.frame()

# contrast indexes
idx5 <- which.min(abs(a.vals - 5))
idx8 <- which.min(abs(a.vals - 8))
idx10 <- which.min(abs(a.vals - 10))
idx12 <- which.min(abs(a.vals - 12))
idx15 <- which.min(abs(a.vals - 15))

### Create Data

for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate.cal),
                        linear = c(est_data$linear.cal),
                        lower = c(est_data[,6] - 1.96*est_data[,7]),
                        upper = c(est_data[,6] + 1.96*est_data[,7]),
                        exposure = rep("Di et al. (2019)", nrow(est_data)),
                        race = rep(scenario$race, nrow(est_data)),
                        dual = rep(scenario$dual, nrow(est_data)))
  
  tmp_1 <- as.numeric(est_data[idx10,6]) - as.numeric(est_data[idx5,6])
  tmp_2 <- as.numeric(est_data[idx12,6]) - as.numeric(est_data[idx8,6])
  tmp_3 <- as.numeric(est_data[idx15,6]) - as.numeric(est_data[idx10,6])
  tmp_4 <- sqrt(as.numeric(est_data[idx10,7])^2 + as.numeric(est_data[idx5,7])^2)
  tmp_5 <- sqrt(as.numeric(est_data[idx12,7])^2 + as.numeric(est_data[idx8,7])^2)
  tmp_6 <- sqrt(as.numeric(est_data[idx15,7])^2 + as.numeric(est_data[idx10,7])^2)
  
  contr_tmp <- data.frame(estimate = c(tmp_1, tmp_2, tmp_3),
                          lower = c(tmp_1 - 1.96*tmp_4, tmp_2 - 1.96*tmp_5, tmp_3 - 1.96*tmp_6),
                          upper = c(tmp_1 + 1.96*tmp_4, tmp_2 + 1.96*tmp_5, tmp_3 + 1.96*tmp_6),
                          pm0 = c(5, 8, 10),
                          pm1 = c(10, 12, 15),
                          race = scenario$race,
                          dual = scenario$dual)
  
  contr_tmp$contrast <- paste0(contr_tmp$pm1, " vs. ", contr_tmp$pm0)
  
  dat <- rbind(dat, dat_tmp)
  contr <- rbind(contr, contr_tmp)
  
}

save(dat, file = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/estimate.RData')
save(contr, file = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/contr.RData')

### Main Plot

scenario <- scenarios[3,]

# histogram and ERF data
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
dat_tmp <- subset(dat, dual == scenario$dual & race == scenario$race)
a_dat <- rep(individual_data$pm25, individual_data$time_count)

# exposure response curve
erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, y = estimate)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate",
       title = "Exposure Response Curve for\n All Medicare Recipients") + 
  theme(legend.position = c(0.02, 0.8),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.044,0.045,0.046,0.047,0.048,0.049)) +
  grids(linetype = "dashed")

# histogram
a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a_dat)) + 
  geom_density(fill = "grey", alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
  labs(x = "Annual Average PM2.5", y = "Exposure Density") + 
  theme(panel.grid = element_blank()) +
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
  guides(fill = "none") +
  theme_cowplot() +
  grids(linetype = "dashed")

align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
main_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/erc_plot.pdf", width = 8, height = 8)
main_plot
dev.off()

### ERC by Race

plot_list <- list()
dual.vals <- c(2, 0, 1)

for (i in 1:length(dual.vals)){
  
  if (dual.vals[i] == 1) {
    main <- "Dual Eligible"
    d <- dual.vals[i]
  } else if (dual.vals[i] == 0) {
    main <- "Dual Ineligible"
    d <- dual.vals[i]
  } else {
    main <- "Dual Eligible + Ineligible"
    d <- c(0,1)
  }
  
  dat_tmp <- subset(dat, dual == as.numeric(dual.vals[i]) & race != "all")
  dat_tmp$race <- str_to_title(dat_tmp$race)
  ylim <- c(min(dat_tmp$lower), max(dat_tmp$upper))
  
  dat_tmp$race <- factor(dat_tmp$race)
  
  # black data
  load(paste0(dir_out_qd, dual.vals[i], "_black_qd.RData"))
  black_data <- subset(individual_data, dual %in% as.numeric(d) & race == 2)
  a_dat_tmp <- data.frame(a = rep(black_data$pm25, black_data$time_count), race = "Black")
  
  # white data
  load(paste0(dir_out_qd, dual.vals[i], "_white_qd.RData"))
  white_data <- subset(individual_data, dual %in% as.numeric(d) & race == 1)
  a_dat <- rbind(a_dat_tmp, data.frame(a = rep(white_data$pm25, white_data$time_count), race = "White"))
    
  if (dual.vals[i] == 2) {
    
    # dual eligible + dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
      geom_line(size = 1) +
      coord_cartesian(xlim = c(5.45,14.55), ylim = c(0.043,0.053)) +
      labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate", 
           color = "Race", title = main) +
      theme(legend.position = c(0.02, 0.9),
            legend.background = element_rect(colour = "black"),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_manual(values = c("#367E18", "#F57328")) +
      scale_y_continuous(breaks = c(0.043,0.045,0.047,0.049,0.051,0.053)) +
      grids(linetype = "dashed")
    
  } else if (dual.vals[i] == 0) {
  
    # dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
      geom_line(size = 1) +
      coord_cartesian(xlim = c(5.45,14.55), ylim = c(0.034, 0.044)) +
      labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate", title = main) +
      theme(legend.position = "none",
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_manual(values = c("#367E18", "#F57328")) +
      scale_y_continuous(breaks = c(0.034,0.036,0.038,0.04,0.042,0.044)) +
      grids(linetype = "dashed")
    
  } else {
    
    # dual eligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
      geom_line(size = 1) +
      coord_cartesian(xlim = c(5.45,14.55), ylim = c(0.065, 0.105)) +
      labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate", title = main) + 
      theme(legend.position = "none",
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_manual(values = c("#367E18", "#F57328")) +
      scale_y_continuous(breaks = c(0.065,0.075,0.085,0.095,0.105)) +
      grids(linetype = "dashed")
    
  }
  
  # histogram
  a_hist_tmp <- ggplot(a_dat, mapping = aes(x = a, fill = race)) + 
    geom_density(alpha = 0.3, adjust = 3)+
    coord_cartesian(xlim = c(5.45,14.55), ylim = c(0,0.15)) +
    labs(x = "Annual Average PM2.5", y = "Exposure Density") + 
    theme(panel.grid = element_blank()) +
    scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
    scale_fill_manual(values = c("#367E18", "#F57328")) +
    guides(fill = "none") +
    theme_cowplot()
  
  align_tmp <- align_plots(a_hist_tmp, erf_strata_tmp, align = "hv", axis = "tblr")
  erf_strata_plot <- ggdraw(align_tmp[[1]]) + draw_plot(align_tmp[[2]])
  
  plot_list[[i]] <- erf_strata_plot
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/strata_plot.pdf", width = 16, height = 8)
strata_plot
dev.off()

### Contrast Plot

contr <- subset(contr, !(dual %in% c(0,1) & race == "all"))
contr$race_dual <- paste(str_to_title(contr$race), ifelse(contr$dual == 0, "- Dual\n Ineligible", 
                                                          ifelse(contr$dual == 1, "- Dual\n Eligible", "- All")))

contr$race_dual <- ifelse(contr$race_dual == "All - All", "All", contr$race_dual)
contr$race_dual <- factor(contr$race_dual, levels = c("All", "White - All", "Black - All", 
                                                      "White - Dual\n Ineligible",  "Black - Dual\n Ineligible",
                                                      "White - Dual\n Eligible", "Black - Dual\n Eligible"))

contrast_plot <- contr %>% 
  ggplot(aes(x = race_dual, y = 100*estimate, color = contrast)) + 
  geom_pointrange(aes(ymin = 100*lower, ymax = 100*upper), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Risk Difference (%)") +
  ggtitle(expression("Risk Difference Estimates")) +
  guides(color = guide_legend(title = "Contrast (mcg/m^3)")) +
  theme_bw() +
  scale_color_manual(values = c("#008080", "#FF00FF","#FFD700")) +
  scale_y_continuous(breaks = round(seq(0, max(100*contr$upper), by = 0.1),1)) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  grids(linetype = "dashed")

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/contrast_plot.pdf", width = 8, height = 8)
contrast_plot
dev.off()

### Linear Plot

dat$dual.name <- ifelse(dat$dual == 1, "Dual Eligible", 
                        ifelse(dat$dual == 0, "Dual Ineligible", "Dual Eligible + Ineligible"))
dat$race <- ifelse(dat$race == "all", "All Races", 
                   ifelse(dat$race == "black", "Black", "White"))

dat$race <- factor(dat$race)
dat$dual.name <- factor(dat$dual.name, levels = c("Dual Eligible + Ineligible", "Dual Ineligible",
                                                  "Dual Eligible"))

linear_plot <- dat %>% 
  ggplot(aes(x = a.vals, y = linear, linetype = dual.name, color = race)) + 
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(min(dat$linear), max(dat$linear))) +
  labs(x = "Annual Average PM2.5", y = "Approximate All-cause Mortality Rate", 
       color = "Race", linetype = "Dual Medicaid/Medicare Eligibility",
       title = "Linear Approximations of the Absolute Risk") +
  theme(legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#000000", "#367E18", "#F57328")) +
  grids(linetype = "dashed")

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/linear_plot.pdf", width = 8, height = 8)
linear_plot
dev.off()
