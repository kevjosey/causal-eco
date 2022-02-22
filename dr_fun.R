
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

tmle_glm <- function(a, x, y, offset, a.vals){
  
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
  a.std <- c(a - pimod.vals) / sqrt(pi2mod.vals)
  dens <- density(a.std)
  pihat <- approx(dens$x, dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / sqrt(pi2mod.vals)
    approx(dens$x, dens$y, xout = std)$y / sqrt(pi2mod.vals)
  } )
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat, na.rm = T)), x = a)$y
  phat[phat<0] <- 1e-4
  
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
  weights <- phat/pihat
  weights[weights > quantile(weights, 0.99)] <- quantile(weights, 0.99)
  base <- nsa*weights
  new_mod <- glm(y ~ 0 + base, offset = log(muhat)+ offset, 
                 family = poisson(link = "log"))
  param <- coef(new_mod)
  trim <- max(weights)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    muhat.tmp <- muhat.mat[,k]
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp, na.rm = TRUE)/pihat.tmp)
    wts[wts > trim] <- trim
    mat <- predict(nsa, newx = rep(a.tmp, n))*wts
    return(mean(exp(log(muhat.tmp) + c(mat%*%param)), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, weights = weights))
  
}

tmle_ranger <- function(a, x, y, offset, a.vals){
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  xa <- data.frame(x, a = a)
  
  # estimate nuisance outcome model with splines
  fmla <- formula(paste0("y ~ ns(a, df = 4) +", paste0(colnames(x), collapse = "+")))
  mumod <- glm(fmla, data = xa, offset = offset, family =  poisson(link = "log"))
  muhat <- predict(mumod, newdata = xa, type = "response")
  
  # estimate nuisance GPS parameters with lm
  pimod <- ranger(a ~ ., data = xa, num.trees = 500, min.node.size = 20)
  pimod.vals <- c(pimod$predictions)
  pi2mod.vals <- mean((a - pimod.vals)^2)
  
  # exposure models
  a.std <- c(a - pimod.vals) / sqrt(pi2mod.vals)
  dens <- density(a.std)
  pihat <- approx(dens$x, dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / sqrt(pi2mod.vals)
    approx(dens$x, dens$y, xout = std)$y / sqrt(pi2mod.vals)
  } )
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat, na.rm = T)), x = a)$y
  phat[phat<0] <- 1e-4
  
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
  weights <- phat/pihat
  weights[weights > quantile(weights, 0.99)] <- quantile(weights, 0.99)
  base <- nsa*weights
  new_mod <- glm(y ~ 0 + base, offset = log(muhat)+ offset, 
                 family = poisson(link = "log"))
  param <- coef(new_mod)
  trim <- max(weights)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    muhat.tmp <- muhat.mat[,k]
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp, na.rm = TRUE)/pihat.tmp)
    wts[wts > trim] <- trim
    mat <- predict(nsa, newx = rep(a.tmp, n))*wts
    return(mean(exp(log(muhat.tmp) + c(mat%*%param)), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, weights = weights))
  
}

bal_plot <- function(a, x, weights, main = "All QD"){
  
  val <- bal.tab(x, treat = a, weights = weights, method = "weighting")
  bal_df <- val$Balance[order(abs(val$Balance$Corr.Un), decreasing = TRUE),]
  labs <- rep(rownames(bal_df), 2)
  vals <- c(bal_df$Corr.Un, bal_df$Corr.Adj)
  adjust <- rep(c("Unadjusted", "Adjusted"), each = nrow(bal_df))
  df <- data.frame(labs = labs, vals = abs(vals), adjust = adjust)
  df$labs <- factor(df$labs, levels = rev(rownames(bal_df)))
  
  fp <- ggplot(data = df, aes(x = labs, y = vals, color = adjust)) +
    geom_point(pch = 21, size = 2) +
    geom_line(aes(group = adjust)) + 
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(yintercept = 0.1, lty = 3, colour = "black") +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Covariates") + ylab("Standardized Absolute Mean Difference") +
    ylim(0, 0.5) +
    guides(color = guide_legend(title = "GPS Adjusting")) +
    theme_bw() + # use a white background
    ggtitle(main)
  
  return(fp)
  
}
