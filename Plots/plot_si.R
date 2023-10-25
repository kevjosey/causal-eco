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
a.vals = seq(4, 16, length.out = 121)

# data directories
dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_SRF/'

dat <- data.frame()
contr <- data.frame()

### Create Data

for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, ".RData"))
  
  dat_tmp <- data.frame(d.vals = c(new_data$est_data$d), 
                        estimate = c(new_data$est_data$estimate),
                        lower = c(new_data$est_data[,2] - 1.96*new_data$est_data[,3]),
                        upper = c(new_data$est_data[,2] + 1.96*new_data$est_data[,3]),
                        race = rep(scenario$race, nrow(new_data$est_data)),
                        dual = rep(scenario$dual, nrow(new_data$est_data)))
  
  tmp_8 <- as.numeric(new_data$est_data[2,2]) - as.numeric(new_data$est_data[1,2])
  tmp_9 <- as.numeric(new_data$est_data[3,2]) - as.numeric(new_data$est_data[1,2])
  tmp_10 <- as.numeric(new_data$est_data[4,2]) - as.numeric(new_data$est_data[1,2])
  tmp_11 <- as.numeric(new_data$est_data[5,2]) - as.numeric(new_data$est_data[1,2])
  tmp_12 <- as.numeric(new_data$est_data[6,2]) - as.numeric(new_data$est_data[1,2])
  se_8 <- sqrt(as.numeric(new_data$est_data[2,3])^2 + as.numeric(new_data$est_data[1,3])^2)
  se_9 <- sqrt(as.numeric(new_data$est_data[3,3])^2 + as.numeric(new_data$est_data[1,3])^2)
  se_10 <- sqrt(as.numeric(new_data$est_data[4,3])^2 + as.numeric(new_data$est_data[1,3])^2)
  se_11 <- sqrt(as.numeric(new_data$est_data[5,3])^2 + as.numeric(new_data$est_data[1,3])^2)
  se_12 <- sqrt(as.numeric(new_data$est_data[6,3])^2 + as.numeric(new_data$est_data[1,3])^2)
  
  contr_tmp <- data.frame(cutoff = c(8,9,10,11,12),
                          estimate = c(tmp_8, tmp_9, tmp_10, tmp_11, tmp_12),
                          se = c(se_8, se_9, se_10, se_11, se_12),
                          lower = c(tmp_8 - 1.96*se_8, tmp_9 - 1.96*se_9, tmp_10 - 1.96*tmp_10,
                                    tmp_11 - 1.96*se_11, tmp_12 - 1.96*se_12),
                          upper = c(tmp_8 + 1.96*se_8, tmp_9 + 1.96*se_9, tmp_10 + 1.96*tmp_10,
                                     tmp_11 + 1.96*se_11, tmp_12 + 1.96*se_12),
                          race = scenario$race,
                          dual = scenario$dual)
  
  dat <- rbind(dat, dat_tmp)
  contr <- rbind(contr, contr_tmp)
  
}

save(dat, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/estimate_srf.RData')
save(contr, file = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/contrast_srf.RData')

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
  ggplot(aes(x = race_dual, y = 100*estimate, color = factor(cutoff))) + 
  geom_pointrange(aes(ymin = 100*lower, ymax = 100*upper), position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Risk Difference (%)", title = "Risk Difference Estimates") +
  guides(color = guide_legend(title = ~ PM[2.5]*" Contrasts ("*mu*g*"/"*m^3*")")) +
  theme(legend.position = c(0.17, 0.87),
        legend.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(breaks = round(seq(-0.5, max(100*contr$upper), by = 0.1),1)) +
  grids(linetype = "dashed")

pdf(file = "~/Figures/contrast_plot.pdf", width = 8, height = 8)
contrast_plot
dev.off()
