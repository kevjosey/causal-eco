## Load libraries --------------------------------------------------------------

library(CausalGPS)
library(data.table)
library(dplyr)
library(splines)
library(mgcv)

## Load input data -------------------------------------------------------------

load("/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Data/qd/2_all_qd.RData")
print(paste0("Length of data:  ", nrow(new_data$x)))
a.vals <- seq(5, 15, length.out = 51)

## Pre-processing --------------------------------------------------------------

w <- setDF(new_data$w)
x.tmp <- setDF(new_data$x)
zip <- x.tmp$zip
a <- x.tmp$pm25
x <- subset(x.tmp, select = -c(zip, pm25))
delta_n <- a.vals[2] - a.vals[1]
trim <- 0.05

## Generate Pseudo Population --------------------------------------------------

set_logger(logger_level="DEBUG")

s_t <- proc.time()

set.seed(892)
match_pop <- generate_pseudo_pop(Y = zip,
                                 w = a,
                                 c = x,
                                 ci_appr = "matching",
                                 pred_model = "sl",
                                 gps_model = "parametric",
                                 use_cov_transform = TRUE,
                                 transformers = list("pow2", "pow3"),
                                 sl_lib = c("m_xgboost"),
                                 params = list(xgb_nrounds = seq(10,100)),
                                 nthread = 24,
                                 covar_bl_method = "absolute",
                                 covar_bl_trs = 0.1,
                                 covar_bl_trs_type = "mean",
                                 trim_quantiles = c(trim, 1 - trim),
                                 optimized_compile = TRUE,
                                 max_attempt = 10,
                                 matching_fun = "matching_l1",
                                 delta_n = delta_n,
                                 scale = 1.0)

e_t <- proc.time()

## Merge individual- and zip- level data ---------------------------------------

pseudo <- match_pop$pseudo_pop
match_data <- merge(w, data.frame(zip = pseudo$Y, year = pseudo$year,
                                  pm25 = pseudo$w, counter = pseudo$counter), 
                    by = c("zip", "year"), all = FALSE)
match_data <- subset(match_data, counter > 0)

## Fit cluster matched model conditional on individual covariates --------------

match_curve <- mgcv::bam(dead ~ s(pm25, bs = 'cr', k = 3) + factor(sex) + factor(race) + factor(dual) + factor(age_break), 
                         data = match_data, offset = log(time_count), 
                         family = poisson(link = "log"), weights = counter)

## Marginalize out individual level covariates ---------------------------------

# being overly-caustious about offsets
wts <- match_data$time_count
covar <- subset(match_data, select = c(sex, race, dual, age_break))

estimate <- sapply(a.vals, function(a.tmp, ...) {
  print(which(a.tmp == a.vals))
  match_estimate <- predict(match_curve, type = "response", newdata = data.frame(pm25 = a.tmp, covar, time_count = 1))
  return(weighted.mean(match_estimate, w = wts, na.rm = TRUE))
})

## Try with TMLE ----------------------------------------------------------------

source('/nfs/nsaph_ci3/ci3_analysis/josey_erc_strata/Code/R/tmle_glm.R')

w.tmp <- setDF(new_data$w)
x.tmp <- setDF(new_data$x)
wx.tmp <- merge(w.tmp, x.tmp, by = c("zip", "year"))
a_x <- x.tmp$pm25
a_w <- wx.tmp$pm25
y <- wx.tmp$dead
offset <- log(wx.tmp$time_count)
x <- subset(x.tmp, select = -c(zip, pm25))
w <- subset(wx.tmp, select = -c(zip, pm25, dead, time_count))

# note df = 4 because of intercept, trunc != trim
target <- tmle_glm(a_w = a_w, a_x = a_x, w = w, x = x,
                   y = y, offset = offset, df = 4,
                   family = poisson(link = "log"), 
                   a.vals = a.vals, trunc = 0.01)

## Compare Curves ----------------------------------------------------------------

plot(a.vals, estimate, type = "l", col = 2, ylim = c(0.04, 0.05), xlab = "PM2.5", ylab = "Risk", lwd = 2)
lines(a.vals, target$estimate, col = 3, lwd = 2)
legend(5, 0.05, legend = c("Matching", "TMLE"), col = c(2, 3), lty = c(1, 1), lwd = 2)