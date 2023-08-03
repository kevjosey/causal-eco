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

load("/Users/kevin/Library/CloudStorage/Dropbox/Projects/ERC-Strata/Output/estimate.RData")

contr <- rr <- data.frame()
a.vals <- dat$a.vals
dat$se <- (dat$upper - dat$estimate)/1.96

# contrast indexes

### Create Data

for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  
  dat_tmp <- subset(dat, race == scenario$race & dual == scenario$dual)
  
  idx8 <- which.min(abs(dat_tmp$a.vals - 8))
  idx9 <- which.min(abs(dat_tmp$a.vals - 9))
  idx10 <- which.min(abs(dat_tmp$a.vals - 10))
  idx11 <- which.min(abs(dat_tmp$a.vals - 11))
  idx12 <- which.min(abs(dat_tmp$a.vals - 12))
  
  tmp_1 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx10])
  tmp_2 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx9])
  tmp_3 <- as.numeric(dat_tmp$estimate[idx12]) - as.numeric(dat_tmp$estimate[idx8])
  tmp_4 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx10])^2)
  tmp_5 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx9])^2)
  tmp_6 <- sqrt(as.numeric(dat_tmp$se[idx12])^2 + as.numeric(dat_tmp$se[idx8])^2)
  
  contr_tmp <- data.frame(estimate = c(tmp_1, tmp_2, tmp_3),
                          lower = c(tmp_1 - 1.96*tmp_4, tmp_2 - 1.96*tmp_5, tmp_3 - 1.96*tmp_6),
                          upper = c(tmp_1 + 1.96*tmp_4, tmp_2 + 1.96*tmp_5, tmp_3 + 1.96*tmp_6),
                          pm0 = c(10, 9, 8),
                          pm1 = c(12, 12, 12),
                          race = scenario$race,
                          dual = scenario$dual)
  
  contr_tmp$contrast <- paste0(contr_tmp$pm1, " vs. ", contr_tmp$pm0)
  contr <- rbind(contr, contr_tmp)
  
  # relative risks
  
  rr_tmp_1 <- c(as.numeric(dat_tmp$estimate))/as.numeric(dat_tmp$estimate[idx8])
  rr_tmp_2 <- c(as.numeric(dat_tmp$estimate))/as.numeric(dat_tmp$estimate[idx9])
  rr_tmp_3 <- c(as.numeric(dat_tmp$estimate))/as.numeric(dat_tmp$estimate[idx10])
  rr_tmp_4 <- c(as.numeric(dat_tmp$estimate))/as.numeric(dat_tmp$estimate[idx11])
  
  rr_tmp <- data.frame(estimate = c(rr_tmp_1, rr_tmp_2, rr_tmp_3, rr_tmp_4),
                       pm0 = c(rep(8, nrow(dat_tmp)),
                               rep(9, nrow(dat_tmp)),
                               rep(10, nrow(dat_tmp)),
                               rep(11, nrow(dat_tmp))),
                       pm1 = rep(dat_tmp$a.vals, 4),
                       race = scenario$race,
                       dual = scenario$dual)
  
  rr <- rbind(rr, rr_tmp)

  
}

### Main Plot

scenario <- scenarios[3,]

# histogram and ERF data
rr_tmp <- subset(rr, dual == scenario$dual & race == scenario$race)

# exposure response curve
erf_plot <- rr_tmp %>% 
  ggplot(aes(x = pm1, y = estimate, color = factor(pm0))) + 
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(x = 8, y = 0.85, xend = 8, yend = 1, linetype = "dashed", color = "#D81B60") +
  geom_segment(x = 9, y = 0.85, xend = 9, yend = 1, linetype = "dashed", color = "#1E88E5") +
  geom_segment(x = 10, y = 0.85, xend = 10, yend = 1, linetype = "dashed", color = "#FFC107") +
  geom_segment(x = 11, y = 0.85, xend = 11, yend = 1, linetype = "dashed", color = "#004D40") +
  coord_cartesian(xlim = c(5,15), ylim = c(0.9,1.05)) +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Relative Risks",
       title = "Relative Risks for\n All Medicare Recipients", color = ~ "Baseline "*PM[2.5]) + 
  theme_bw() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = c(0.9,0.925,0.95,0.975,1,1.025,1.05)) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15)) +
  scale_color_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
                     labels = c(~ "8 "*mu*g*"/"*m^3, ~ "9 "*mu*g*"/"*m^3, ~ "10 "*mu*g*"/"*m^3, ~ "11 "*mu*g*"/"*m^3)) +
  grids(linetype = "dashed")

# histogram
# a_hist <- ggplot(data.frame(a = a_dat), mapping = aes(x = a_dat)) + 
#   geom_density(fill = "grey", alpha = 0.3, adjust = 3)+
#   coord_cartesian(xlim = c(5,15), ylim = c(0,0.15)) +
#   labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Exposure Density") + 
#   theme(panel.grid = element_blank()) +
#   scale_y_continuous(position = "right", breaks = c(0, 0.05, 0.10, 0.15)) +
#   guides(fill = "none") +
#   theme_cowplot() +
#   grids(linetype = "dashed")
# 
# align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
# main_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "~/Documents/erc_plot.pdf", width = 8, height = 8)
erf_plot
dev.off()

### ERC by Race

plot_list <- list()
dual.vals <- c(2, 0, 1)

rr$dual_label <- ifelse(rr$dual == 0, "High SEP", ifelse(rr$dual == 1, "Low SEP", "High + Low SEP"))
rr$dual_label <- factor(rr$dual_label, levels = c("High + Low SEP", "High SEP", "Low SEP"))
rr$race <- str_to_title(rr$race)

rr_tmp <- subset(rr, race %in% c("All", "Black", "White") & pm0 == 8)

erf_strata <- rr_tmp %>% 
  ggplot(aes(x = a.vals, y = estimate, color = race)) + 
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  geom_segment(x = 8, y = 0.85, xend = 8, yend = 1, linetype = "dashed", color = "#000000") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ dual_label) +
  coord_cartesian(xlim = c(5,15)) +
  theme_bw() +
  labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Hazard", 
       color = "Race", title = ~ "Stratified Relative Risk Estimates Evaluated with a Baseline "*PM[2.5]*" Concentration of "*8*" "*mu*g*"/"*m^3) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("#000000", "#367E18", "#F57328")) +
  scale_x_continuous(breaks = c(5,6,7,8,9,10,11,12,13,14,15)) +
  grids(linetype = "dashed")

pdf(file = "~/Documents/erc_strata.pdf", width = 16, height = 8)
erf_strata
dev.off()

### N > PM

scenarios <- expand.grid(dual = c("High SEP", "Low SEP", "High + Low SEP"), 
                         race = c("All","White", "Black"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)

r <- c(rep(0, times = 100000*0.85), rep(1, times = 100000*0.15))
d <- rbinom(100000, 1, 0.1 + 0.1*r)
pm <- abs(rnorm(100000, 6 + 0.5*r + d, 2))

fake_dat <- data.frame(r = r, d = d, pm = pm)
fake_dat$race <- ifelse(fake_dat$r == 1, "Black", "White")
fake_dat$dual <- ifelse(fake_dat$d == 1, "Low SEP", "High SEP")
prop <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  
  if (scenario$dual != "High + Low SEP") 
    f_dat <- subset(fake_dat, dual == scenario$dual)
  
  if (scenario$race != "All")
    f_dat <- subset(f_dat, race == scenario$race)
  
  prop8 <- mean(f_dat$pm > 8)  
  prop9 <- mean(f_dat$pm > 9)
  prop10 <- mean(f_dat$pm > 10)
  prop11 <- mean(f_dat$pm > 11)
  
  prop_tmp <- data.frame(estimate = c(prop8, prop9, prop10, prop11),
                         pm = c(8,9,10,11),
                         race = scenario$race,
                         dual = scenario$dual)
  
  prop <- rbind(prop, prop_tmp)
  
}

prop$label <- paste(prop$race, prop$dual, sep = " - ")

prop_plot <- prop %>% 
  ggplot(aes(x = label, y = estimate, fill = factor(pm))) + 
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = ~ "Subpopulation", y = "Proportion at Risk",
       title = ~ "Proportion of Medicare Participants with Exceeding Levels of "*PM[2.5], 
       fill = ~ PM[2.5]*" Cutoffs ("*mu*g*"/"*m^3*")") + 
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_manual(values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40"),
                    labels = c(~ PM[2.5]*" > 8", ~ PM[2.5]*" > 9", ~ PM[2.5]*" > 10", ~ PM[2.5]*" > 11")) +
  grids(linetype = "dashed")

pdf(file = "~/Documents/pm_prop.pdf", width = 16, height = 8)
prop_plot
dev.off()

### Contrast Plot

contr <- subset(contr, !(dual %in% c(0,1) & race == "all"))
contr$contrast <- factor(contr$contrast, levels = c("12 vs. 10", "12 vs. 9", "12 vs. 8"))
contr$race_dual <- paste0(str_to_title(contr$race), ifelse(contr$dual == 0, " -\n High SEP", 
                                                           ifelse(contr$dual == 1, " -\n Low SEP",
                                                                  " -\n High + Low SEP")))

contr$race_dual <- ifelse(contr$race_dual == "All -\n High + Low SEP", "All - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "White -\n High + Low SEP", "White - High\n + Low SEP", contr$race_dual)
contr$race_dual <- ifelse(contr$race_dual == "Black -\n High + Low SEP", "Black - High\n + Low SEP", contr$race_dual)
contr$race_dual <- factor(contr$race_dual, levels = c("All - High\n + Low SEP",
                                                      "White - High\n + Low SEP", 
                                                      "Black - High\n + Low SEP", 
                                                      "White -\n High SEP", "Black -\n High SEP",
                                                      "White -\n Low SEP",  "Black -\n Low SEP"))

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

contrast_plot
