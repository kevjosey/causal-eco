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
scenarios <- c("Ecological", "Individual", "Individual_Ecological")
a.vals = seq(4, 16, length.out = 121)

# data directories
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Eco/'

dat <- data.frame()
contr <- data.frame()

# contrast indexes
idx5 <- which.min(abs(a.vals - 5))
idx8 <- which.min(abs(a.vals - 8))
idx10 <- which.min(abs(a.vals - 10))
idx12 <- which.min(abs(a.vals - 12))
idx15 <- which.min(abs(a.vals - 15))

### Create Data

for (i in 1:length(scenarios)) {
  
  # QD
  scenario <- scenarios[i]
  load(paste0(dir_out, scenario, ".RData"))
  
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
                        estimand = scenario,
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
                          estimand = scenario)
  
  contr_tmp$contrast <- paste0(contr_tmp$pm1, " vs. ", contr_tmp$pm0)
  
  dat <- rbind(dat, dat_tmp)
  contr <- rbind(contr, contr_tmp)
  
}

save(dat, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/estimate.RData')
save(contr, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/contrast.RData')

### ERC by Estimand

plot_list <- list()

ylim <- c(min(dat$lower), max(dat$upper))
dat$estimand <- factor(dat$estimand)

# dual eligible + dual ineligible ERCs
erf_compare <- dat %>% 
  ggplot(aes(x = a.vals, color = estimand)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1, aes(y = estimate)) +
  coord_cartesian(xlim = c(5,15), ylim = ylim) +
  labs(x = " ", y = "All-cause Mortality Rate", color = "Method", title = "Comparison of ERCs") +
  scale_color_manual(values = c("#75bad3", "#ea8832","#489f8c")) +
  scale_y_continuous(breaks = c(0.025,0.03,0.035,0.04,0.045,0.05, 0.055)) +
  theme_cowplot() +
  grids(linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.background = element_rect(colour = "black", fill = "white")) 

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