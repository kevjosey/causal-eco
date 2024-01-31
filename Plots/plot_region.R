library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(splines)
library(ggplot2)
library(ggpubr)
library(gtable)
library(cowplot)

### Race Only ###

scenarios <- expand.grid(race = c("all", "white","black"))
scenarios$race <- as.character(scenarios$race)

# contrast indexes
idx5 <- which.min(abs(a.vals - 5))
idx8 <- which.min(abs(a.vals - 8))
idx9 <- which.min(abs(a.vals - 9))
idx10 <- which.min(abs(a.vals - 10))
idx11 <- which.min(abs(a.vals - 11))
idx12 <- which.min(abs(a.vals - 12))
idx15 <- which.min(abs(a.vals - 15))

dat <- data.frame()
rr <- data.frame()

for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i,]
  load(paste0(dir_out, scenario, ".RData"))
  
  u.zip <- unique(new_data$wx$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  n <- nrow(new_data$wx)
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate),
                        excess = c(excess_death$estimate),
                        lower = c(est_data$estimate) - 1.96*c(est_data$se),
                        upper = c(est_data$estimate) + 1.96*c(est_data$se),
                        lower.ed = c(excess_death$estimate) - 1.96*c(excess_death$se),
                        upper.ed = c(excess_death$estimate) + 1.96*c(excess_death$se),
                        race = rep(scenario, nrow(est_data)),
                        region = rep("All Regions", nrow(est_data)),
                        year = rep("All Years", nrow(est_data)),
                        n = rep(sum(new_data$wx$y), nrow(new_data$est_data)))
  
  # hazard ratios
  hr_tmp_8 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx8]))
  hr_tmp_9 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx9]))
  hr_tmp_10 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx10]))
  hr_tmp_11 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx11]))
  hr_tmp_12 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx12]))
  
  # delta method standard errors
  log_hr_se_8 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx8]^2)/c(est_data$estimate[idx8]^2))
  log_hr_se_9 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx9]^2)/c(est_data$estimate[idx9]^2))
  log_hr_se_10 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx10]^2)/c(est_data$estimate[idx10]^2))
  log_hr_se_11 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx11]^2)/c(est_data$estimate[idx11]^2))
  log_hr_se_12 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx12]^2)/c(est_data$estimate[idx12]^2))
  
  rr_tmp <- data.frame(a.vals = rep(est_data$a.vals, 5),
                       estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       lower = c(exp(log(hr_tmp_8) - 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) - 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) - 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) - 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) - 1.96*log_hr_se_12)),
                       upper = c(exp(log(hr_tmp_8) + 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) + 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) + 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) + 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) + 1.96*log_hr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       race = scenario,
                       region = "All Regions",
                       year = "All Years")
  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

### Race, Region and Years ###

scenarios_region_years <- expand.grid(race = c("all", "white","black"),
                                      region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"),
                                      years = c("2001-2004", "2005-2008", "2009-2012", "2013-2016"))
scenarios_region_years$years <- as.character(scenarios_region_years$years)
scenarios_region_years$race <- as.character(scenarios_region_years$race)
scenarios_region_years$region <- as.character(scenarios_region_years$region)
a.vals = seq(4, 16, length.out = 121)

dir_out = '/n/dominici_nsaph_l3/projects/kjosey-erc-strata/Output/Strata_ERF/'

# race by dual by region plots
for (i in 1:nrow(scenarios_region_years)) {
  
  scenario <- scenarios_region_years[i,]
  load(paste0(dir_out, scenario$race, "_", scenario$region, "_", scenario$years, ".RData"))
  
  u.zip <- unique(new_data$wx$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  n <- nrow(new_data$wx)
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate),
                        excess = c(excess_death$estimate),
                        lower = c(est_data$estimate) - 1.96*c(est_data$se),
                        upper = c(est_data$estimate) + 1.96*c(est_data$se),
                        lower.ed = c(excess_death$estimate) - 1.96*c(excess_death$se),
                        upper.ed = c(excess_death$estimate) + 1.96*c(excess_death$se),
                        race = rep(scenario$race, nrow(est_data)),
                        region = rep(scenario$region, nrow(est_data)),
                        year = rep(scenario$years, nrow(est_data)),
                        n = rep(sum(new_data$wx$y), nrow(new_data$est_data)))
  
  # hazard ratios
  hr_tmp_8 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx8]))
  hr_tmp_9 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx9]))
  hr_tmp_10 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx10]))
  hr_tmp_11 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx11]))
  hr_tmp_12 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx12]))
  
  # delta method standard errors
  log_hr_se_8 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx8]^2)/c(est_data$estimate[idx8]^2))
  log_hr_se_9 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx9]^2)/c(est_data$estimate[idx9]^2))
  log_hr_se_10 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx10]^2)/c(est_data$estimate[idx10]^2))
  log_hr_se_11 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx11]^2)/c(est_data$estimate[idx11]^2))
  log_hr_se_12 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx12]^2)/c(est_data$estimate[idx12]^2))
  
  rr_tmp <- data.frame(a.vals = rep(est_data$a.vals, 5),
                       estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       lower = c(exp(log(hr_tmp_8) - 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) - 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) - 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) - 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) - 1.96*log_hr_se_12)),
                       upper = c(exp(log(hr_tmp_8) + 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) + 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) + 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) + 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) + 1.96*log_hr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       race = scenario$race,
                       region = scenario$region,
                       year = scenario$years)
  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

### Race and Region ###

scenarios_region <- expand.grid(race = c("all", "white","black"),
                                region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))
scenarios_region$race <- as.character(scenarios_region$race)
scenarios_region$region <- as.character(scenarios_region$region)
a.vals = seq(4, 16, length.out = 121)

for (i in 1:nrow(scenarios_region)) {
  
  scenario <- scenarios_region[i,]
  load(paste0(dir_out, scenario$race, "_", scenario$region, ".RData"))
  
  n <- nrow(new_data$wx)
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate),
                        excess = c(excess_death$estimate),
                        lower = c(est_data$estimate) - 1.96*c(est_data$se),
                        upper = c(est_data$estimate) + 1.96*c(est_data$se),
                        lower.ed = c(excess_death$estimate) - 1.96*c(excess_death$se),
                        upper.ed = c(excess_death$estimate) + 1.96*c(excess_death$se),
                        race = rep(scenario$race, nrow(est_data)),
                        region = rep(scenario$region, nrow(est_data)),
                        year = rep("All Years", nrow(est_data)),
                        n = rep(sum(new_data$wx$y), nrow(new_data$est_data)))
  
  # hazard ratios
  hr_tmp_8 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx8]))
  hr_tmp_9 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx9]))
  hr_tmp_10 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx10]))
  hr_tmp_11 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx11]))
  hr_tmp_12 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx12]))
  
  # delta method standard errors
  log_hr_se_8 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx8]^2)/c(est_data$estimate[idx8]^2))
  log_hr_se_9 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx9]^2)/c(est_data$estimate[idx9]^2))
  log_hr_se_10 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx10]^2)/c(est_data$estimate[idx10]^2))
  log_hr_se_11 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx11]^2)/c(est_data$estimate[idx11]^2))
  log_hr_se_12 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx12]^2)/c(est_data$estimate[idx12]^2))
  
  rr_tmp <- data.frame(a.vals = rep(est_data$a.vals, 5),
                       estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       lower = c(exp(log(hr_tmp_8) - 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) - 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) - 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) - 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) - 1.96*log_hr_se_12)),
                       upper = c(exp(log(hr_tmp_8) + 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) + 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) + 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) + 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) + 1.96*log_hr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       race = scenario$race,
                       region = scenario$region,
                       year = "All Years")
  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

### Race and Years ###

scenarios_years <- expand.grid(race = c("all", "white","black"),
                               years = c("2001-2004", "2005-2008", "2009-2012", "2013-2016"))
scenarios_years$race <- as.character(scenarios_years$race)
scenarios_years$years <- as.character(scenarios_years$years)
a.vals = seq(4, 16, length.out = 121)

# race by dual by region plots
for (i in 1:nrow(scenarios_years)) {
  
  scenario <- scenarios_years[i,]
  load(paste0(dir_out, scenario$race, "_", scenario$years, ".RData"))
  
  n <- nrow(new_data$wx)
  est_data <- new_data$est_data
  excess_death <- new_data$excess_death
  
  dat_tmp <- data.frame(a.vals = c(est_data$a.vals), 
                        estimate = c(est_data$estimate),
                        excess = c(excess_death$estimate),
                        lower = c(est_data$estimate) - 1.96*c(est_data$se),
                        upper = c(est_data$estimate) + 1.96*c(est_data$se),
                        lower.ed = c(excess_death$estimate) - 1.96*c(excess_death$se),
                        upper.ed = c(excess_death$estimate) + 1.96*c(excess_death$se),
                        race = rep(scenario$race, nrow(est_data)),
                        region = rep("All Regions", nrow(est_data)),
                        year = rep(scenario$years, nrow(est_data)),
                        n = rep(sum(new_data$wx$y), nrow(new_data$est_data)))
  
  # hazard ratios
  hr_tmp_8 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx8]))
  hr_tmp_9 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx9]))
  hr_tmp_10 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx10]))
  hr_tmp_11 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx11]))
  hr_tmp_12 <- c(as.numeric(est_data$estimate)/as.numeric(est_data$estimate[idx12]))
  
  # delta method standard errors
  log_hr_se_8 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx8]^2)/c(est_data$estimate[idx8]^2))
  log_hr_se_9 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                        c(est_data$se[idx9]^2)/c(est_data$estimate[idx9]^2))
  log_hr_se_10 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx10]^2)/c(est_data$estimate[idx10]^2))
  log_hr_se_11 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx11]^2)/c(est_data$estimate[idx11]^2))
  log_hr_se_12 <- sqrt(c(est_data$se^2)/c(est_data$estimate^2) +
                         c(est_data$se[idx12]^2)/c(est_data$estimate[idx12]^2))
  
  rr_tmp <- data.frame(a.vals = rep(est_data$a.vals, 5),
                       estimate = c(hr_tmp_8, hr_tmp_9, hr_tmp_10, hr_tmp_11, hr_tmp_12),
                       lower = c(exp(log(hr_tmp_8) - 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) - 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) - 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) - 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) - 1.96*log_hr_se_12)),
                       upper = c(exp(log(hr_tmp_8) + 1.96*log_hr_se_8), 
                                 exp(log(hr_tmp_9) + 1.96*log_hr_se_9), 
                                 exp(log(hr_tmp_10) + 1.96*log_hr_se_10),
                                 exp(log(hr_tmp_11) + 1.96*log_hr_se_11),
                                 exp(log(hr_tmp_12) + 1.96*log_hr_se_12)),
                       pm0 = c(rep(8, nrow(est_data)),
                               rep(9, nrow(est_data)),
                               rep(10, nrow(est_data)), 
                               rep(11, nrow(est_data)),
                               rep(12, nrow(est_data))),
                       race = scenario$race,
                       region = "All Regions",
                       year = scenario$years)
  
  dat <- rbind(dat, dat_tmp)
  rr <- rbind(rr, rr_tmp)
  
}

### All Subjects

plot_list <- list()
situations <- expand.grid(region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  main <- str_to_title(situation)
  
  # factor race
  dat_tmp <- subset(dat, race == "all" & region == situation & year != "All Years")
  dat_tmp$year <- str_to_title(dat_tmp$year)
  
  # graph breaks
  ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
            max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)
    
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = year)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate, linetype = "solid")) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Years", title = main, linetype = "Model") + 
    guides(linetype = "none") +
    scale_y_continuous(breaks = breaks) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(colour = "black"))
  
  plot_list[[i]] <- erf_strata_plot
    
}

strata_plot <- ggarrange(plotlist = plot_list[1:4], ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE)

pdf(file = "~/Figures/region_strata_plot.pdf", width = 16, height = 16)
strata_plot
dev.off()

### Black Participants

plot_list <- list()
situations <- expand.grid(region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  main <- str_to_title(situation)
  
  # factor race
  dat_tmp <- subset(dat, race == "black" & region == situation & year != "All Years")
  dat_tmp$year <- str_to_title(dat_tmp$year)
  
  # graph breaks
  ylim <- c(min(dat_tmp$lower[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]),
            max(dat_tmp$upper[dat_tmp$a.vals >= 5 & dat_tmp$a.vals <= 15]))
  breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)
  
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = year)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate, linetype = "solid")) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Years", title = main, linetype = "Model") + 
    guides(linetype = "none") +
    scale_y_continuous(breaks = breaks) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(colour = "black"))
  
  plot_list[[i]] <- erf_strata_plot
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:4], ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE)

pdf(file = "~/Figures/region_strata_black_plot.pdf", width = 16, height = 16)
strata_plot
dev.off()

### White Participants

plot_list <- list()
situations <- expand.grid(region = c("MIDWEST", "NORTHEAST", "SOUTH", "WEST"))

for (i in 1:nrow(situations)){
  
  situation <- situations[i,]
  
  main <- str_to_title(situation)
  
  # factor race
  dat_tmp <- subset(dat, race == "all" & region == situation)
  dat_tmp$year <- str_to_title(dat_tmp$year)
  
  # graph breaks
  ylim <- c(0.039,0.06)
  breaks <- round(seq(ylim[1], ylim[2], length.out = 6), 4)
  
  # dual ineligible + eligible
  erf_strata_plot <- dat_tmp %>% 
    ggplot(aes(x = a.vals, color = year)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(size = 1, aes(y = estimate, linetype = "solid")) +
    coord_cartesian(xlim = c(5,15), ylim = ylim) +
    labs(x = ~ "Annual Average "*PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "All-cause Mortality Risk Ratio", 
         color = "Years", title = main, linetype = "Model") + 
    guides(linetype = "none") +
    scale_y_continuous(breaks = breaks) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.background = element_rect(colour = "black"))
  
  plot_list[[i]] <- erf_strata_plot
  
  
}

strata_plot <- ggarrange(plotlist = plot_list[1:4], ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE)

pdf(file = "~/Figures/region_strata_plot_white.pdf", width = 16, height = 16)
strata_plot
dev.off()

