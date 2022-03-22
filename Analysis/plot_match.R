library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(cobalt)

# scenarios
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios <- rbind(c(dual = 2, race = "all"), scenarios)
a.vals <- seq(5, 15, length.out = 51)
knot.list <- list(pm25 = c(8,12))
n.boot <- 1000

# Data Directories
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/Match_rm/'

dat_qd <- data.frame()
dat_rm <- data.frame()
contr <- data.frame()

# contrast indexes
idx5 <- which(a.vals == 5)
idx8 <- which(a.vals == 8)
idx10 <- which(a.vals == 10)
idx12 <- which(a.vals == 12)

# Race or dual Plot
for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  adj <- sqrt(2*sqrt(n.zip))/sqrt(n.zip)
  dat_qd_tmp <- data.frame(a.vals = boot_data$a.vals, estimate = boot_data$estimate,
                           lower = sapply(1:nrow(boot_data), function(j,...) exp(log(boot_data[j,2]) - 1.96*sd(log(boot_data[j,3:(n.boot + 2)]))*adj)),
                           upper = sapply(1:nrow(boot_data), function(j,...) exp(log(boot_data[j,2]) + 1.96*sd(log(boot_data[j,3:(n.boot + 2)]))*adj)),
                           exposure = rep("QD", nrow(boot_data)),
                           race = rep(scenario$race, nrow(boot_data)),
                           dual = rep(scenario$dual, nrow(boot_data)))
  
  dat_qd <- rbind(dat_qd, dat_qd_tmp)
  
  qd_tmp_1 <- as.numeric(boot_data[idx10,2:(n.boot + 2)])/as.numeric(boot_data[idx5,2:(n.boot + 2)])
  qd_tmp_2 <- as.numeric(boot_data[idx12,2:(n.boot + 2)])/as.numeric(boot_data[idx8,2:(n.boot + 2)])
  
  contr_qd_tmp <- data.frame(estimate = c(qd_tmp_1[1], qd_tmp_2[1]),
                             lower = c(exp(log(qd_tmp_1[1]) - 1.96*sd(log(qd_tmp_1[2:n.boot + 1]))*adj),
                                       exp(log(qd_tmp_2[1]) - 1.96*sd(log(qd_tmp_2[2:n.boot + 1]))*adj)),
                             upper = c(exp(log(qd_tmp_1[1]) + 1.96*sd(log(qd_tmp_1[2:n.boot + 1]))*adj), 
                                       exp(log(qd_tmp_2[1]) + 1.96*sd(log(qd_tmp_2[2:n.boot + 1]))*adj)),
                             pm0 = c(5, 8),
                             pm1 = c(10, 12),
                             exposure = c("QD"),
                             race = scenario$race,
                             dual = scenario$dual)
  
  # RM
  load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  adj <- sqrt(2*sqrt(n.zip))/sqrt(n.zip)
  dat_rm_tmp <- data.frame(a.vals = boot_data$a.vals, 
                           estimate = boot_data$estimate,
                           lower = sapply(1:nrow(boot_data), function(j,...) exp(log(boot_data[j,2]) - 1.96*sd(log(boot_data[j,3:(n.boot + 2)]))*adj)),
                           upper = sapply(1:nrow(boot_data), function(j,...) exp(log(boot_data[j,2]) + 1.96*sd(log(boot_data[j,3:(n.boot + 2)]))*adj)),
                           exposure = rep("RM", nrow(boot_data)),
                           race = rep(scenario$race, nrow(boot_data)),
                           dual = rep(scenario$dual, nrow(boot_data)))
  dat_rm <- rbind(dat_rm, dat_rm_tmp)
  
  rm_tmp_1 <- as.numeric(boot_data[idx10,2:(n.boot + 2)])/as.numeric(boot_data[idx5,2:(n.boot + 2)])
  rm_tmp_2 <- as.numeric(boot_data[idx12,2:(n.boot + 2)])/as.numeric(boot_data[idx8,2:(n.boot + 2)])
  
  contr_rm_tmp <- data.frame(estimate = c(rm_tmp_1[1], rm_tmp_2[1]),
                             lower = c(exp(log(rm_tmp_1[1]) - 1.96*sd(log(rm_tmp_1[2:n.boot + 1]))*adj), 
                                       exp(log(rm_tmp_2[1]) - 1.96*sd(log(rm_tmp_2[2:n.boot + 1]))*adj)),
                             upper = c(exp(log(rm_tmp_1[1]) + 1.96*sd(log(rm_tmp_1[2:n.boot + 1]))*adj), 
                                       exp(log(rm_tmp_2[1]) + 1.96*sd(log(rm_tmp_2[2:n.boot + 1]))*adj)),
                             pm0 = c(5, 8),
                             pm1 = c(10, 12),
                             exposure = c("RM"),
                             race = scenario$race,
                             dual = scenario$dual)
  
  contr <- rbind(contr, contr_qd_tmp, contr_rm_tmp)
  
}

### Main Plot

i <- 1
scenario <- scenarios[i,]

# QD
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
dat_qd_tmp <- subset(dat_qd, dual == scenario$dual & race == scenario$race)
a_dat <- data.frame(a = match_data[!duplicated(paste0(match_data$zip, match_data$year)),]$pm25, exposure = "QD")

# RM
load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
dat_rm_tmp <- subset(dat_rm, dual == scenario$dual & race == scenario$race)
a_dat <- rbind(a_dat, data.frame(a = match_data[!duplicated(paste0(match_data$zip, match_data$year)),]$pm25, exposure = "RM"))

# combine
dat_tmp <- rbind(dat_qd_tmp, dat_rm_tmp)

erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(5,15)) +
  labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate",
       color = "Exposure Model") + 
  theme(legend.position = c(0.02, 0.9),
        legend.background = element_rect(colour = "black"),
        panel.grid=element_blank())

a_dat <- subset(a_dat, a >= 5 & a <= 15)

a_hist <- ggplot(a_dat, mapping = aes(x = a, fill = exposure)) + 
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.25)+
  coord_cartesian(xlim = c(5,15)) +
  labs(x = "Annual Average PM2.5", y = "Exposure Density") + 
  theme(panel.grid=element_blank()) +
  scale_y_continuous(position = "right") +
  guides(fill="none") +
  theme_cowplot()

align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
big_bad_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/erc_plot.pdf", width = 8, height = 8)
big_bad_plot
dev.off()

### Plot by Race

plot_list <- list()
situations <- expand.grid(dual = c(2, 0, 1), exposure = c("QD", "RM"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  if (situation$dual == 1)
    dname <- "Dual Eligible"
  else if (situation$dual == 0)
    dname <- "Dual Ineligible"
  else
    dname <- "All"
  
  if (situation$exposure == "QD"){
    
    main <- paste(dname, "QD")
    dat_tmp <- subset(dat_qd, dual == as.numeric(situation$dual) & race != "all")
    
  } else {
    
    main <- paste(dname, "RM")
    dat_tmp <- subset(dat_rm, dual == as.numeric(situation$dual) & race != "all")
    
  }
  
  dat_tmp$race <- str_to_title(dat_tmp$race)
  
  if (situation$dual == 0)
    ylim <- c(0.03, 0.045)
  else if (situation$dual == 1)
    ylim <- c(0.06, 0.11)
  else
    ylim <- c(0.04, 0.0525)
  
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = race)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(3,18), ylim = ylim) +
    labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate", 
         color = "Race", title = main) + 
    theme_bw() +
    guides(color = guide_legend(title = "Race")) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9")) +
    grids(linetype = "dashed")
  
  plot_list[[i]] <- erf_strata_plot
  
}

strata_plot <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, legend = "bottom", common.legend = TRUE)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/strata_plot.pdf", width = 10, height = 10)
strata_plot
dev.off()

### Contrast Plot

contr$race_dual <- paste(str_to_title(contr$race), ifelse(contr$dual == 0, "- Dual\n Ineligible", 
                                                          ifelse(contr$dual == 1, "- Dual\n Eligible", "- All")))
contr$race_dual <- ifelse(contr$race_dual == "All - All", "All", contr$race_dual)
contr$race_dual <- factor(contr$race_dual, levels = c("All", "White - All", "White - Dual\n Ineligible", "White - Dual\n Eligible",
                                                      "Black - All", "Black - Dual\n Ineligible", "Black - Dual\n Eligible"))
contr_1 <- subset(contr, pm0 == 5) 
contr_2 <- subset(contr, pm0 == 8) 

contrast_plot_1 <- contr_1 %>% 
  ggplot(aes(x = race_dual, y = estimate, color = exposure)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25)) +
  labs(x = "", y = "Risk Ratio") +
  ggtitle(expression("Changing PM2.5 from 5 " ~ mu * "g/m3 to 10 " ~ mu * "g/m3")) +
  theme_bw() +
  guides(color = guide_legend(title = "Exposure Model")) + 
  grids(linetype = "dashed")

contrast_plot_2 <- contr_2 %>% 
  ggplot(aes(x = race_dual, y = estimate, color = exposure)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25)) +
  labs(x = "", y = "Risk Ratio") +
  ggtitle(expression("Changing PM2.5 from 8 " ~ mu * "g/m3 to 12 " ~ mu * "g/m3")) + 
  theme_bw() +
  guides(color = guide_legend(title = "Exposure Model")) + 
  grids(linetype = "dashed")

contrast_plot <- ggarrange(contrast_plot_1 + theme(legend.position="none"), 
                           contrast_plot_2 + theme(legend.position="none"),
                           ncol = 1, nrow = 2, 
                           legend = "bottom", common.legend = TRUE)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/contrast_plot.pdf", width = 12, height = 6)
contrast_plot
dev.off()

### Balance Plot

i <- 1
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
corr_data <- corr_data[order(corr_data[,1], decreasing = FALSE),]
df <- data.frame(covars = rep(rownames(corr_data), 2), vals = c(corr_data$original, corr_data$adjusted),
                 adjust = rep(c("Unadjusted", "Adjusted"), each = nrow(corr_data)))
df$covars <- factor(df$covars, levels = rownames(corr_data))

bplot_qd <- ggplot(data = df, aes(x = covars, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) + 
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.5) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() + # use a white background
  theme(axis.text.y = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  ggtitle("QD")

load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
corr_data <- corr_data[order(corr_data[,1], decreasing = FALSE),]
df <- data.frame(covars = rep(rownames(corr_data), 2), vals = c(corr_data$original, corr_data$adjusted),
                 adjust = rep(c("Unadjusted", "Adjusted"), each = nrow(corr_data)))
df$covars <- factor(df$covars, levels = rownames(corr_data))

bplot_rm <- ggplot(data = df, aes(x = covars, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) + 
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.4) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() + # use a white background
  theme(axis.text.y = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  ggtitle("RM")


balance_plot <- ggarrange(bplot_qd + theme(legend.position="none"), 
                          bplot_rm + theme(legend.position="none"),
                          ncol = 2, nrow = 1, 
                          legend = "bottom", common.legend = TRUE)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/balance_plot.pdf", width = 8, height = 8)
balance_plot
dev.off()