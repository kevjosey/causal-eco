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
  load(paste0(dir_out, scenario$race, "_", scenario$dual, ".RData"))
  
  u.zip <- unique(new_data$wx$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  n <- nrow(new_data$wx)
  
  dat_tmp <- data.frame(a.vals = c(new_data$est_data$a.vals), 
                        estimate = c(new_data$est_data$estimate),
                        excess = c(new_data$excess_death$estimate),
                        lower = c(new_data$est_data$estimate) - 
                          1.96*c(new_data$est_data$se),
                        upper = c(new_data$est_data$estimate) + 
                          1.96*c(new_data$est_data$se),
                        lower.ed = c(new_data$excess_death$estimate) - 
                          1.96*c(new_data$excess_death$se),
                        upper.ed = c(new_data$excess_death$estimate) +
                          1.96*c(new_data$excess_death$se),
                        race = rep(scenario$race, nrow(new_data$est_data)),
                        dual = rep(scenario$dual, nrow(new_data$est_data)),
                        n = rep(sum(new_data$wx$y), nrow(new_data$est_data)))
  
  tmp_1 <- as.numeric(new_data$est_data[idx10,2]) - as.numeric(new_data$est_data[idx5,2])
  tmp_2 <- as.numeric(new_data$est_data[idx12,2]) - as.numeric(new_data$est_data[idx8,2])
  tmp_3 <- as.numeric(new_data$est_data[idx15,2]) - as.numeric(new_data$est_data[idx10,2])
  tmp_4 <- sqrt(as.numeric(new_data$est_data[idx10,3])^2 + as.numeric(new_data$est_data[idx5,3])^2)
  tmp_5 <- sqrt(as.numeric(new_data$est_data[idx12,3])^2 + as.numeric(new_data$est_data[idx8,3])^2)
  tmp_6 <- sqrt(as.numeric(new_data$est_data[idx15,3])^2 + as.numeric(new_data$est_data[idx10,3])^2)
  
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

save(dat, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/estimate.RData')
save(contr, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/contrast.RData')

### Main Plot

# histogram and ERF data
load(paste0(dir_out, "all_both.RData"))
dat_tmp <- subset(dat, dual == "both" & race == "all")
a_dat <- rep(new_data$wx$pm25, new_data$wx$n)

# exposure response curve
erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  coord_cartesian(xlim = c(6,14), ylim = c(0.042,0.047)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       title = "Exposure Response Curve for\n All Medicare Recipients") + 
  scale_y_continuous(breaks = c(0.042,0.043,0.044,0.045,0.046,0.047)) +
  scale_x_continuous(breaks = c(6,8,10,12,14)) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  grids(linetype = "dashed")

# histogram
a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a)) + 
  geom_density(fill = "grey", alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(6,14), ylim = c(0,0.15)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
  scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
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
  
  dat_tmp <- subset(dat, dual == dual.vals[i] & race != "all")
  dat_tmp$race <- str_to_title(dat_tmp$race)
  ylim <- c(min(dat_tmp$lower), max(dat_tmp$upper))
  
  dat_tmp$race <- factor(dat_tmp$race)
  
  # # black data
  # load(paste0(dir_out, dual.vals[i], "_black.RData"))
  # a_dat_tmp <- data.frame(a = rep(new_data$wx$pm25, new_data$wx$n), race = "Black")
  # 
  # # white data
  # load(paste0(dir_out, dual.vals[i], "_white.RData"))
  # a_dat_tmp <- rbind(a_dat_tmp, data.frame(a = rep(new_data$wx$pm25, new_data$wx$n), race = "White"))
  # 
  # # asian data
  # load(paste0(dir_out, dual.vals[i], "_asian.RData"))
  # a_dat_tmp <- rbind(a_dat_tmp, data.frame(a = rep(new_data$wx$pm25, new_data$wx$n), race = "Asian"))
  # 
  # # hispanic data
  # load(paste0(dir_out, dual.vals[i], "_hispanic.RData"))
  # a_dat_tmp <- rbind(a_dat_tmp, data.frame(a = rep(new_data$wx$pm25, new_data$wx$n), race = "Hispanic"))
  
  if (dual.vals[i] == "both") {
    
    # dual eligible + dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(size = 1, aes(y = estimate)) +
      coord_cartesian(xlim = c(5,15), ylim = c(0.025,0.055)) +
      labs(x = " ", y = "All-cause Mortality Rate", color = "Race", title = main) +
      scale_color_manual(values = c("#75bad3", "#ea8832","#ea3323","#489f8c")) +
      scale_y_continuous(breaks = c(0.025,0.03,0.035,0.04,0.045,0.05, 0.055)) +
      theme_cowplot() +
      grids(linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = c(0.02, 0.8),
            legend.background = element_rect(colour = "black", fill = "white")) 
    
  } else if (dual.vals[i] == "high") {
    
    # dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(size = 1, aes(y = estimate)) +
      coord_cartesian(xlim = c(5,15), ylim = c(0.02, 0.045)) +
      labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
           color = "Race", title = main) +
      scale_color_manual(values = c("#75bad3", "#ea8832","#ea3323","#489f8c")) +
      scale_y_continuous(breaks = c(0.02,0.025,0.03,0.035,0.04,0.045,0.05)) +
      theme_cowplot() +
      grids(linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none") 
    
  } else {
    
    # dual eligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(size = 1, aes(y = estimate)) +
      coord_cartesian(xlim = c(5,15), ylim = c(0.035, 0.105)) +
      labs(x = " ", y = "All-cause Mortality Rate", color = "Race", title = main) + 
      scale_color_manual(values = c("#75bad3", "#ea8832","#ea3323","#489f8c")) +
      scale_y_continuous(breaks = c(0.035,0.045,0.055,0.065,0.075,0.085,0.095,0.105)) +
      theme_cowplot() +
      grids(linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")
    
  }
  
  leg <- get_legend(erf_strata_tmp)
  
  if (dual.vals[i] == "both") {
    erf_strata_tmp <- erf_strata_tmp + theme(legend.position = "none") +
      annotation_custom(leg, xmin = 6, xmax = 9, ymin = 0.051, ymax = 0.053)
  }
  
  # histogram
  # a_hist_tmp <- ggplot(a_dat_tmp, mapping = aes(x = a, fill = race)) + 
  #   geom_density(alpha = 0.3, adjust = 3)+
  #   coord_cartesian(xlim = c(6,14), ylim = c(0,0.15)) +
  #   labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
  #   theme(panel.grid = element_blank()) +
  #   scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
  #   scale_x_continuous(breaks = c(6,8,10,12,14)) +
  #   scale_fill_manual(values = c("#75bad3", "#ea8832","#ea3323","#489f8c")) +
  #   guides(fill = "none") +
  #   theme_cowplot() +
  #   grids(linetype = "dashed") 
  
  # if (dual.vals[i] == "both") {
  #   align_tmp <- align_plots(a_hist_tmp, erf_strata_tmp + 
  #                              theme(legend.position = "none") +
  #                              annotation_custom(leg, xmin = 7, xmax = 9, ymin = 0.051, ymax = 0.053), 
  #                            align = "hv", axis = "tblr")
  # } else {
  #     align_tmp <- align_plots(a_hist_tmp, erf_strata_tmp, align = "hv", axis = "tblr")
  # }
  
  # erf_strata_plot <- ggdraw(align_tmp[[1]]) + draw_plot(align_tmp[[2]])
  
  plot_list[[i]] <- erf_strata_tmp
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1)

pdf(file = "~/Figures/strata_plot.pdf", width = 12, height = 6)
strata_plot
dev.off()

### Contrast Plot

contr <- subset(dat, a.vals == 8)

contrast_plot <- contr %>% 
  ggplot(aes(x = str_to_upper(race), y = 100*excess/n, color = str_to_upper(dual))) + 
  geom_pointrange(aes(ymin = 100*lower.ed/n, ymax = 100*upper.ed/n), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Percent of Deaths Avoidable (%)", title = "Excess Death Estimates",
       color = "Socioeconomic Position") +
  theme(legend.position = c(0.12, 0.9),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#008080", "#FF00FF","#FFD700")) +
  grids(linetype = "dashed")

pdf(file = "~/Figures/contrast_plot.pdf", width = 8, height = 8)
contrast_plot
dev.off()