library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(cobalt)

# scenarioss
scenarios <- expand.grid(dual = c(0, 1, 2), race = c("all","white", "black"))
scenarios$dual <- as.numeric(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
a.vals <- seq(3, 17, length.out = 106)

# Data Directories
dir_out_qd = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/'
dir_out_rm = '/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/'

dat_qd <- data.frame()
dat_rm <- data.frame()
contr <- data.frame()

# contrast indexes
idx5 <- which.min(abs(a.vals - 5))
idx8 <- which.min(abs(a.vals - 8))
idx10 <- which.min(abs(a.vals - 10))
idx12 <- which.min(abs(a.vals - 12))

# Race or dual Plot
for (i in 1:nrow(scenarios)) {
  
  # QD
  scenario <- scenarios[i,]
  load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
  dat_qd_tmp <- data.frame(a.vals = rep(est_data$a.vals, 2), 
                           estimate = c(est_data$estimate.lm, est_data$estimate.sl),
                           lower = c(est_data[,2] - 1.96*est_data[,3], est_data[,4] - 1.96*est_data[,5]),
                           upper = c(est_data[,2] + 1.96*est_data[,3], est_data[,4] + 1.96*est_data[,5]),
                           exposure = rep("Di et al. (2019)", nrow(est_data)),
                           gps = rep(c("LM","SL"), each = nrow(est_data)),
                           race = rep(scenario$race, 2*nrow(est_data)),
                           dual = rep(scenario$dual, 2*nrow(est_data)))
  
  dat_qd <- rbind(dat_qd, dat_qd_tmp)
  
  # RM
  load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
  dat_rm_tmp <- data.frame(a.vals = rep(est_data$a.vals, 2), 
                           estimate = c(est_data$estimate.lm, est_data$estimate.sl),
                           lower = c(est_data[,2] - 1.96*est_data[,3], est_data[,4] - 1.96*est_data[,5]),
                           upper = c(est_data[,2] + 1.96*est_data[,3], est_data[,4] + 1.96*est_data[,5]),
                           exposure = rep("van Donkelaar et al. (2019)", nrow(est_data)),
                           gps = rep(c("LM","SL"), each = nrow(est_data)),
                           race = rep(scenario$race, nrow(est_data)),
                           dual = rep(scenario$dual, nrow(est_data)))
  dat_rm <- rbind(dat_rm, dat_rm_tmp)
  
}

### Main Plot

i <- 3
scenario <- scenarios[i,]

# QD
load(paste0(dir_out_qd, scenario$dual, "_", scenario$race, "_qd.RData"))
dat_qd_tmp <- subset(dat_qd, dual == scenario$dual & race == scenario$race)
a_dat <- data.frame(a = zip_data$pm25, exposure = "Di et al. (2019)")

# RM
load(paste0(dir_out_rm, scenario$dual, "_", scenario$race, "_rm.RData"))
dat_rm_tmp <- subset(dat_rm, dual == scenario$dual & race == scenario$race)
a_dat <- rbind(a_dat, data.frame(a = zip_data$pm25, exposure = "van Donkelaar et al. (2016)"))

# combine
dat_tmp <- rbind(dat_qd_tmp, dat_rm_tmp)

erf_plot <- dat_tmp %>% 
  ggplot(aes(x = a.vals, y = estimate, color = exposure, linetype = gps)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(3,17)) +
  labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate",
       color = "Exposure Assessment", linetype = "GPS Method") + 
  theme(legend.position = c(0.02, 0.8),
        legend.background = element_rect(colour = "black"),
        panel.grid=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

a_dat <- subset(a_dat, a >= 3 & a <= 17)

a_hist <- ggplot(a_dat, mapping = aes(x = a, fill = exposure)) + 
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.25)+
  coord_cartesian(xlim = c(3,17)) +
  labs(x = "Annual Average PM2.5", y = "Exposure Density") + 
  theme(panel.grid=element_blank()) +
  scale_y_continuous(position = "right") +
  guides(fill="none") +
  theme_cowplot()

align <- align_plots(a_hist, erf_plot, align = "hv", axis = "tblr")
main_plot <- ggdraw(align[[1]]) + draw_plot(align[[2]])

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/erc_plot.pdf", width = 8, height = 8)
main_plot
dev.off()

### Plot by Race

plot_list <- list()
situations <- expand.grid(dual = c(2, 0, 1), exposure = c("Di et al. (2019)", "van Donkelaar et al. (2016)"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  if (situation$dual == 1)
    main <- "Dual Eligible"
  else if (situation$dual == 0)
    main <- "Dual Ineligible"
  else
    main <- "All"
  
  if (situation$exposure == "Di et al. (2019)"){
    
    dat_tmp <- subset(dat_qd, dual == as.numeric(situation$dual) & race != "all" & gps == "SL")
    
  } else {
    
    dat_tmp <- subset(dat_rm, dual == as.numeric(situation$dual) & race != "all" & gps == "SL")
    
  }
  
  dat_tmp$race <- str_to_title(dat_tmp$race)
  
  if (situation$dual == 0)
    ylim <- c(0.033, 0.045)
  else if (situation$dual == 1)
    ylim <- c(0.06, 0.105)
  else
    ylim <- c(0.039, 0.053)
  
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, y = estimate, color = race)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted") +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(3,17), ylim = ylim) +
    labs(x = "Annual Average PM2.5", y = "All-cause Mortality Rate", 
         color = "Race", title = main) + 
    theme_bw() +
    guides(color = guide_legend(title = "Race")) + 
    scale_color_manual(values=c("#FFAF40", "#1CADBA")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    grids(linetype = "dashed")
  
  plot_list[[i]] <- erf_strata_plot
  
}

strata_plot_tmp1 <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp2 <- ggarrange(plotlist = plot_list[4:6], ncol = 3, nrow = 1, legend = "bottom", common.legend = TRUE)

strata_plot1 <- annotate_figure(strata_plot_tmp1, top = text_grob("Di et al. (2019)", face = "bold", size = 14))
strata_plot2 <- annotate_figure(strata_plot_tmp2, top = text_grob("van Donkelaar et al. (2016)", face = "bold", size = 14))

strata_plot <- ggarrange(strata_plot1, strata_plot2, nrow = 2, ncol = 1, legend = "bottom", common.legend = TRUE)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/strata_plot.pdf", width = 10, height = 10)
strata_plot
dev.off()

# Covariate Balance Plot

bal_dat <- function(a, x, weights){
  
  val <- bal.tab(x, treat = a, weights = weights, method = "weighting", continuous = "raw", s.d.denom = "pooled")
  bal_df <- val$Balance
  labs <- rep(rownames(bal_df), 2)
  vals_tmp <- cbind(bal_df$Corr.Un, bal_df$Corr.Adj)
  vals_year <- c(mean(abs(vals_tmp[1:17,1])), mean(abs(vals_tmp[1:17,2])))
  vals_region <- c(mean(abs(vals_tmp[32:34,1])), mean(abs(vals_tmp[32:34,2])))
  vals_tmp2 <- rbind(cbind(abs(vals_tmp[-c(1:17, 32:34),1]), 
                           abs(vals_tmp[-c(1:17, 32:34),2])),
                     vals_year, vals_region)
  rownames(vals_tmp2) <- c("Mean BMI", "Smoking Rate", "% Hispanic", "% Black", 
                           "Median Household Income", "Median House Value", "% Below Poverty Level", 
                           "% Below High School Education", "Population Density", "% Owner-Occupied Housing", 
                           "Summer Temperature","Winter Temperature", "Summer Humidity", "Winter Humidity",
                           "Calendar Year","Census Region")
  vals_tmp2 <- vals_tmp2[order(vals_tmp2[,1], decreasing = TRUE),]
  
  vals <- c(vals_tmp2[,1], vals_tmp2[,2])
  adjust <- rep(c("Unadjusted", "LM"), each = nrow(vals_tmp2))
  labs <- rep(rownames(vals_tmp2), times = 2)
  df <- data.frame(labs = labs, vals = vals, adjust = adjust)
  df$labs <- factor(df$labs, levels = rev(rownames(vals_tmp2)))
  
  return(df)
  
}

load("/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_qd/2_all_qd.RData")
bdat_1 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.lm)
bdat_2 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.sl)
bdat_tmp <- subset(bdat_2, adjust == "LM")
bdat_tmp$adjust <- "SuperLearner"

df1 <- rbind(bdat_1, bdat_tmp) 
df1$adjust <- factor(df1$adjust, levels = c("Unadjusted", "LM", "SuperLearner"))

bplot_1 <- ggplot(data = df1, aes(x = labs, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) + 
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.35) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() + # use a white background
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values=c("#F77452","#82F739", "#6867AA")) +
  ggtitle("Di et al. (2019)")

load("/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/DR_rm/0_all_rm.RData")
bdat_1 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.lm)
bdat_2 <- bal_dat(a = zip_data$pm25, x = zip_data[,c(2,4:20)], weights = zip_data$weights.sl)
bdat_tmp <- subset(bdat_2, adjust == "LM")
bdat_tmp$adjust <- "SuperLearner"

df2 <- rbind(bdat_1, bdat_tmp) 
df2$adjust <- factor(df2$adjust, levels = c("Unadjusted", "LM", "SuperLearner"))

bplot_2 <- ggplot(data = df2, aes(x = labs, y = vals, color = adjust)) +
  geom_point(pch = 21, size = 2) +
  geom_line(aes(group = adjust)) + 
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Covariates") + ylab("Absolute Correlation") +
  ylim(0, 0.35) +
  guides(color = guide_legend(title = "Implementation")) +
  theme_bw() + # use a white background
  theme(axis.text.y = element_text(angle = 30, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values=c("#F77452","#82F739", "#6867AA")) +
  ggtitle("van Donkelaar et al. (2016)")

balance_plot <- ggarrange(bplot_1 + theme(legend.position="none"), 
                          bplot_2 + theme(legend.position="none"),
                          ncol = 2, nrow = 1, 
                          legend = "bottom", common.legend = TRUE)

pdf(file = "/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Output/balance_plot.pdf", width = 10, height = 10)
balance_plot
dev.off()


