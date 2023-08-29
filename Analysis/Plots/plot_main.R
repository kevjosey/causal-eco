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
scenarios <- expand.grid(dual = c("high", "low", "both"), race = c("white", "black", "asian", "hispanic", "all"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals = seq(2, 31, length.out = 146)

# data directories
dir_data = '/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/'
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Age_Strata_Data/'

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
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_both_all.RData"))

  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate),
                        lower = c(est_data[,2] - 1.96*est_data[,3]),
                        upper = c(est_data[,2] + 1.96*est_data[,3]),
                        exposure = rep("Di et al. (2019)", nrow(est_data)),
                        race = rep(scenario$race, nrow(est_data)),
                        dual = rep(scenario$dual, nrow(est_data)))
  
  tmp_1 <- as.numeric(est_data[idx10,2]) - as.numeric(est_data[idx5,2])
  tmp_2 <- as.numeric(est_data[idx12,2]) - as.numeric(est_data[idx8,2])
  tmp_3 <- as.numeric(est_data[idx15,2]) - as.numeric(est_data[idx10,2])
  tmp_4 <- sqrt(as.numeric(est_data[idx10,3])^2 + as.numeric(est_data[idx5,3])^2)
  tmp_5 <- sqrt(as.numeric(est_data[idx12,3])^2 + as.numeric(est_data[idx8,3])^2)
  tmp_6 <- sqrt(as.numeric(est_data[idx15,3])^2 + as.numeric(est_data[idx10,3])^2)
  
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

# histogram and ERF data
load(paste0(dir_out, "both_all_both_all.RData"))
dat_tmp <- subset(dat, dual == "both" & race == "all")
agg_data <- merge(data.frame(id = individual_data$id, n = individual_data$n), 
               data.frame(id = zip_data$id, pm25 = zip_data$pm25), by = "id")
a_dat <- rep(agg_data$pm25, agg_data$time_count)

# exposure response curve
erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, y = estimate)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15), ylim = c(0.044,0.049)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate",
       title = "Exposure Response Curve for\n All Medicare Recipients") + 
  theme(legend.position = c(0.02, 0.8),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.044,0.045,0.046,0.047,0.048,0.049)) +
  grids(linetype = "dashed")

# histogram
a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a)) + 
  geom_density(fill = "grey", alpha = 0.3, adjust = 3)+
  coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
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
dual.vals <- c("both", "high", "low")

for (i in 1:length(dual.vals)){
  
  if (dual.vals[i] == "low") {
    main <- "Low SEP"
  } else if (dual.vals[i] == "high") {
    main <- "High SEP"
  } else {
    main <- "High + Low SEP"
  }
  
  dat_tmp <- subset(dat, dual == dual.vals[i] & race != "all")
  dat_tmp$race <- str_to_title(dat_tmp$race)
  ylim <- c(min(dat_tmp$lower), max(dat_tmp$upper))
  
  dat_tmp$race <- factor(dat_tmp$race)
  
  # black data
  load(paste0(dir_out_qd, dual.vals[i], "_black_both_all.RData"))
  agg_data <- merge(data.frame(id = individual_data$id, n = individual_data$n), 
                    data.frame(id = zip_data$id, pm25 = zip_data$pm25), by = "id")
  a_dat_tmp <- data.frame(rep(agg_data$pm25, agg_data$time_count), race = "Black")
  
  # white data
  load(paste0(dir_out_qd, dual.vals[i], "_white_both_all.RData"))
  agg_data <- merge(data.frame(id = individual_data$id, n = individual_data$n), 
                    data.frame(id = zip_data$id, pm25 = zip_data$pm25), by = "id")
  a_dat_tmp <- rbind(a_dat_tmp, data.frame(rep(agg_data$pm25, agg_data$time_count), race = "White"))
  
  # asian data
  load(paste0(dir_out_qd, dual.vals[i], "_asian_both_all.RData"))
  agg_data <- merge(data.frame(id = individual_data$id, n = individual_data$n), 
                    data.frame(id = zip_data$id, pm25 = zip_data$pm25), by = "id")
  a_dat_tmp <- rbind(a_dat_tmp, data.frame(rep(agg_data$pm25, agg_data$time_count), race = "Asian"))
  
  # hispanic data
  load(paste0(dir_out_qd, dual.vals[i], "_asian_both_all.RData"))
  agg_data <- merge(data.frame(id = individual_data$id, n = individual_data$n), 
                    data.frame(id = zip_data$id, pm25 = zip_data$pm25), by = "id")
  a_dat_tmp <- rbind(a_dat_tmp, data.frame(rep(agg_data$pm25, agg_data$time_count), race = "Asian"))
    
  if (dual.vals[i] == "both") {
    
    # dual eligible + dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
      geom_line(size = 1) +
      coord_cartesian(xlim = c(5,15), ylim = c(0.043,0.053)) +
      labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
           color = "Race", title = main) +
      theme(legend.position = "none",
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_manual(values = c("#367E18", "#F57328")) +
      scale_y_continuous(breaks = c(0.043,0.045,0.047,0.049,0.051,0.053)) +
      grids(linetype = "dashed")
    
  } else if (dual.vals[i] == "high") {
  
    # dual ineligible ERCs
    erf_strata_tmp <- dat_tmp %>% 
      ggplot(aes(x = a.vals, y = estimate, color = race)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
      geom_line(size = 1) +
      coord_cartesian(xlim = c(5,15), ylim = c(0.034, 0.044)) +
      labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", title = main) +
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
      coord_cartesian(xlim = c(5,15), ylim = c(0.065, 0.105)) +
      labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", title = main) + 
      theme(legend.position = "none",
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_manual(values = c("#367E18", "#F57328")) +
      scale_y_continuous(breaks = c(0.065,0.075,0.085,0.095,0.105)) +
      grids(linetype = "dashed")
    
  }
  
  leg_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = race)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(5,15), ylim = c(0.043,0.053)) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
         color = "Race", title = main) +
    theme_bw() +
    theme(legend.background = element_rect(colour = "black")) +
    scale_color_manual(values = c("#367E18", "#F57328"))
  
  leg <- gtable_filter(ggplot_gtable(ggplot_build(leg_plot)), "guide-box")
  
  # histogram
  a_hist_tmp <- ggplot(a_dat, mapping = aes(x = a, fill = race)) + 
    geom_density(alpha = 0.3, adjust = 3)+
    coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
    theme(panel.grid = element_blank()) +
    scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
    scale_fill_manual(values = c("#367E18", "#F57328")) +
    guides(fill = "none") +
    theme_cowplot()
  
  align_tmp <- align_plots(a_hist_tmp, erf_strata_tmp +
                             annotation_custom(leg, xmin = 5, xmax = 7, ymin = 0.051, ymax = 0.053), 
                           align = "hv", axis = "tblr")
  erf_strata_plot <- ggdraw(align_tmp[[1]]) +  draw_plot(align_tmp[[2]])
  
  plot_list[[i]] <- erf_strata_plot
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/strata_plot.pdf", width = 16, height = 8)
strata_plot
dev.off()

### Contrast Plot

contr <- subset(contr, !(dual %in% c("low","high") & race == "all"))
contr$race_dual <- paste0(str_to_title(contr$race), ifelse(contr$dual == "high", " -\n High SEP", 
                                                          ifelse(contr$dual == "low", " -\n Low SEP",
                                                                 " -\n High + Low SEP")))

contr$race_dual <- ifelse(contr$race_dual == "All -\n High + Low SEP", "All - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "White -\n High + Low SEP", "White - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "Black -\n High + Low SEP", "Black - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "Asian -\n High + Low SEP", "Asian - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "Hispanic -\n High + Low SEP", "Hispanic - High\n + Low SEP", contr$race_dual)
contr$race_dual <- factor(contr$race_dual, levels = c("All - High\n + Low SEP",
                                                      "White - High\n + Low SEP", 
                                                      "Black - High\n + Low SEP", 
                                                      "Hispanic - High\n + Low SEP",
                                                      "Asian - High\n + Low SEP",
                                                      "White -\n High SEP", "Black -\n High SEP",
                                                      "Hispanic -\n High SEP", "Asian -\n High SEP",
                                                      "White -\n Low SEP",  "Black -\n Low SEP",
                                                      "Hispanic -\n Low SEP",  "Asian -\n Low SEP"))

contrast_plot <- contr %>% 
  ggplot(aes(x = race_dual, y = 100*estimate, color = contrast)) + 
  geom_pointrange(aes(ymin = 100*lower, ymax = 100*upper), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Risk Difference (%)", title = "Risk Difference Estimates") +
  guides(color = guide_legend(title = ~ PM[2.5]*" Contrasts ("*mu*g*"/"*m^3*")")) +
  theme(legend.position = c(0.17, 0.87),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#008080", "#FF00FF","#FFD700")) +
  scale_y_continuous(breaks = round(seq(0, max(100*contr$upper), by = 0.1),1)) +
  grids(linetype = "dashed")

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/contrast_plot.pdf", width = 8, height = 8)
contrast_plot
dev.off()
