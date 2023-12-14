library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gtable)
library(cowplot)

# scenarios
scenarios <- expand.grid(dual = c("both", "high", "low"), race = c("all", "white","black","hispanic","asian"),
                         region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))
scenarios$dual <- as.character(scenarios$dual)
scenarios$race <- as.character(scenarios$race)
scenarios$region <- as.character(scenarios$region)
a.vals = seq(4, 16, length.out = 121)

dat <- data.frame()

dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'

# race by dual by region plots
for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario$dual, "_", scenario$race, "_", scenario$region, ".RData"))
  
  dat_tmp <- data.frame(a.vals = c(new_data$est_data$a.vals), 
                        estimate = c(new_data$est_data$estimate),
                        ls = c(new_data$est_data$ls),
                        spl = c(new_data$est_data$spl),
                        lower = c(new_data$est_data[,2] - 1.96*new_data$est_data[,3]),
                        upper = c(new_data$est_data[,2] + 1.96*new_data$est_data[,3]),
                        race = rep(scenario$race, nrow(new_data$est_data)),
                        sample_size = rep(paste0(str_to_title(scenario$race), " = ", 
                                                 formatC(sum(new_data$wx$n),
                                                         format = "d", big.mark = ",")), 
                                          nrow(new_data$est_data)),
                        dual = rep(scenario$dual, nrow(new_data$est_data)),
                        region = rep(scenario$region, nrow(new_data$est_data)))
  
  dat <- rbind(dat, dat_tmp)
  
}

plot_list <- list()
situations <- expand.grid(dual = c("both", "high", "low"), region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))
situations$dual <- as.character(situations$dual)
situations$region <- as.character(situations$region)

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  if (situation$dual == "low") {
    main <- "Low SEP"
  } else if (situation$dual == "high") {
    main <- "High SEP"
 } else {
    main <- "High + Low SEP"
 }  
  
  # factor race
  dat_tmp <- subset(dat, dual == situation$dual & race != "all" & region == situation$region)
  dat_tmp$race <- str_to_title(dat_tmp$race)
  dat_tmp$sample_size <- factor(as.character(dat_tmp$sample_size))
  
  # graph breaks
  ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
            max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)
    
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = sample_size)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate, linetype = "solid")) +
    geom_line(size = 1, aes(y = ls, linetype = "dashed")) +
    geom_line(size = 1, aes(y = spl, linetype = "dotted")) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Rate", 
         color = "Person-Years at Risk", title = main, linetype = "Model") + 
    guides(linetype = "none") +
    scale_y_continuous(breaks = breaks) +
    scale_color_manual(values = c("#75bad3", "#ea8832","#ea3323","#489f8c")) +
    scale_linetype_manual(labels = c("Linear", "Spline", "Shape-Constrained"), 
                          values = c("dashed", "dotted", "solid")) +   
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(colour = "black"))
  
  leg <- gtable_filter(ggplot_gtable(ggplot_build(erf_strata_plot)), "guide-box")
  
  if (i %in% c(1,2,4))
    plot_list[[i]] <- erf_strata_plot + theme(legend.position = "none") +
    annotation_custom(leg, xmin = 4.3, xmax = 8, ymin = breaks[5], ymax = breaks[6])
  else
    plot_list[[i]] <- erf_strata_plot + theme(legend.position = "none") +
    annotation_custom(leg, xmin = 12, xmax = 14.7, ymin = breaks[1], ymax = breaks[2])
  
    
}

strata_plot_tmp1 <- ggarrange(plotlist = plot_list[1:3], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp2 <- ggarrange(plotlist = plot_list[4:6], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp3 <- ggarrange(plotlist = plot_list[7:9], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)
strata_plot_tmp4 <- ggarrange(plotlist = plot_list[10:12], ncol = 3, nrow = 1, legend = "none", common.legend = TRUE)

strata_plot1 <- annotate_figure(strata_plot_tmp1, top = text_grob("Midwest", face = "bold", size = 14))
strata_plot2 <- annotate_figure(strata_plot_tmp2, top = text_grob("Northeast", face = "bold", size = 14))
strata_plot3 <- annotate_figure(strata_plot_tmp3, top = text_grob("South", face = "bold", size = 14))
strata_plot4 <- annotate_figure(strata_plot_tmp4, top = text_grob("West", face = "bold", size = 14))

pdf(file = "~/Figures/region_strata_plot.pdf", width = 16, height = 20)
ggarrange(strata_plot1, strata_plot2, strata_plot3, strata_plot4, nrow = 4, ncol = 1, common.legend = FALSE)
dev.off()
