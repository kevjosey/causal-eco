wrapper <- function(data, a.vals, n.boot = 1000) {
  
  new_data.list <- split(new_data, list(new_data$zip))
  n_zip <- length(unique(new_data$zip))
  
  tmle_out <- mclapply(1:(n.boot + 1), mc.cores = 10, function(j, new_data, new_data.list, a.vals) {
    
    if (j == 1){
      boot_data <- new_data
    } else{  
      idx <- sample(1:n_zip, n_zip, replace = TRUE) 
      boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    }
    
    x <- setDF(subset(boot_data, select = -c(zip, dead, time_count, pm25)))
    a <- boot_data$pm25
    y <- boot_data$dead
    offset <- log(boot_data$time_count)
    
    estimate <- tmle_models(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
    
    return(estimate)
    
  }, new_data = new_data, new_data.list = new_data.list, a.vals = a.vals)
  
  estimate <- tmle_out[[1]]
  boot <- tmle_out[2:(n.boot + 1)]
  
  out <- data.frame(a.vals = a.vals, estimate = estimate, Reduce(cbind, boot))
  colnames(out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  return(out)
  
}

ipw_models <- function(a, x, y, offset, a.vals) {
  
  # IPW Estimator
  ipw_pop <- generate_pseudo_pop(Y = y, w = a, c = x,
                                 ci_appr = "weighting", pred_model = "sl",
                                 gps_model = "parametric", use_cov_transform = TRUE,
                                 transformers = list("pow2", "pow3"), 
                                 nthread = 8, sl_lib = c("m_xgboost"), 
                                 params = list(xgb_nrounds = c(50), xgb_eta = 0.25),
                                 covar_bl_method = "absolute", covar_bl_trs = 0.1,
                                 trim_quantiles = c(0.025,0.975), # trimmed, you can change,
                                 max_attempt = 5, scale = 1.0)
  
  ipw_data <- ipw_pop$pseudo_pop
  ipw_data$offset <- offset[ipw_data$row_index]
  ipw_curve <- mgcv::bam(Y ~ s(w, bs = 'tp'), data = ipw_data, offset = ipw_data$offset,
                         family = poisson(link = "log"), weights = ipw_data$ipw)
  ipw_estimate <- predict(ipw_curve, newdata = data.frame(w = a.vals), type = "response")
  return(ipw_estimate)
  
}

match_models <- function(a, x, y, offset, a.vals) {
  
  ## Matching Estimator
  match_pop <- generate_pseudo_pop(Y = y, w = a, c = x,
                                   ci_appr = "matching", pred_model = "sl",
                                   gps_model = "parametric", use_cov_transform = TRUE,
                                   transformers = list("pow2", "pow3"),
                                   nthread = 8, sl_lib = c("m_xgboost"), 
                                   params = list(xgb_nrounds = c(50), xgb_eta = 0.25),
                                   covar_bl_method = "absolute", covar_bl_trs = 0.1,
                                   trim_quantiles = c(0.025,0.975), # trimmed, you can change,
                                   optimized_compile = TRUE, #created a column counter for how many times matched,
                                   max_attempt = 5, matching_fun = "matching_l1",
                                   delta_n = c(a.vals[2] - a.vals[1]), scale = 1.0)

  match_data <- match_pop$pseudo_pop
  match_data$offset <- offset[match_data$row_index]
  match_data <- subset(match_data, counter > 0)
  match_curve <- mgcv::bam(Y ~ s(w, bs = 'tp'), data = match_data, offset = match_data$offset,
                           family = poisson(link = "log"), weights = match_data$counter)
  match_estimate <- predict(match_curve, newdata = data.frame(w = a.vals), type = "response")
  
  return(match_estimate)
  
}

tmle_models <- function(a, x, y, offset, a.vals){
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  xa <- data.frame(x, a = a)
  
  # estimate nuisance outcome model with splines
  fmla <- formula(paste0("y ~ ns(a, df = 4) +", paste0(colnames(x), collapse = "+")))
  mumod <- glm(fmla, data = xa, offset = offset, family =  poisson(link = "log"))
  muhat <- predict(mumod, newdata = xa, type = "response")
  
  # estimate nuisance GPS parameters with lm
  pimod <- lm(a ~ ., data = xa)
  pimod.vals <- c(pimod$fitted.values)
  pi2mod.vals <- sigma(pimod)^2
  
  # exposure models
  pihat <- dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) 
    dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals)))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  
  # outcomes models given a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    xa.tmp <- data.frame(x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa) 
    return(predict(mumod, newdata = xa.tmp, type = "response"))
    
  })
  
  # find marginal outcome estimates at a
  mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  
  # TMLE update
  nsa <- ns(a, df = 4, intercept = TRUE)
  base <- nsa*phat/pihat
  new_mod <- glm(y ~ 0 + base, offset = log(muhat)+ offset, 
                 family = poisson(link = "log"))
  param <- coef(new_mod)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    muhat.tmp <- muhat.mat[,k]
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    mat <- predict(nsa, newx = rep(a.tmp, n))*c(mean(pihat.tmp)/pihat.tmp)
    return(mean(exp(log(muhat.tmp) + c(mat%*%param))))
    
  })
  
  return(estimate)
  
}
