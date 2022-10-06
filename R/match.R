
match_estimate <- function(a, w, x, zip, a.vals, fmla, trim = 0.01, attempts = 5) {
  
  if (trim < 0 | trim > 0.5)
    stop("trim < 0 | trim > 0.5")
  
  if (length(a) != nrow(x))
    stop("length(a) != nrow(x)")
  
  if (length(a) != length(zip))
    stop("length(a) != length(zip)")
  
  # matching estimator
  match_pop <- generate_pseudo_pop(Y = zip, w = a, c = x,
                                   ci_appr = "matching",
                                   pred_model = "sl",
                                   gps_model = "non-parametric",
                                   use_cov_transform = TRUE,
                                   transformers = list("pow2", "pow3"),
                                   sl_lib = c("m_xgboost"),
                                   params = list(xgb_nrounds = c(200)),
                                   nthread = 12, # number of cores, you can change,
                                   covar_bl_method = "absolute",
                                   covar_bl_trs = 0.1,
                                   covar_bl_trs_type = "mean",
                                   trim_quantiles = c(trim, 1 - trim), # trimed, you can change,
                                   optimized_compile = TRUE, #created a column counter for how many times matched,
                                   max_attempt = attempts,
                                   matching_fun = "matching_l1",
                                   delta_n = 0.2, # you can change this to the one you used in previous analysis,
                                   scale = 0.5)
  
  # merge individual level data
  pseudo <- match_pop$pseudo_pop
  match_data <- merge(w, data.frame(zip = pseudo$Y,
                                    year = pseudo$year,
                                    pm25 = pseudo$w,
                                    counter = pseudo$counter), 
                      by = c("zip", "year"), all = FALSE)
  match_data <- subset(match_data, counter > 0)
  
  # fit model conditional on individual level covariates
  match_curve <- mgcv::bam(fmla, data = match_data, offset = log(time_count), 
                           family = poisson(link = "log"), weights = counter)
  
  # cautionary about offsets
  wts <- w$time_count
  covar.df <- subset(w, select = -c(zip, time_count, dead))
  
  # marginalize
  estimate <- sapply(a.vals, function(a.tmp, ...) {
    
    df.tmp <- data.frame(pm25 = a.tmp, covar.df,counter = 1)
    match_estimate <- predict(match_curve, type = "response", newdata = df.tmp)
    return(weighted.mean(match_estimate, w = wts, na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, match_data = match_data,
              adjusted_corr_results = match_pop$adjusted_corr_results, 
              original_corr_results = match_pop$original_corr_results))
  
}
