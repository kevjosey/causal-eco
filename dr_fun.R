
match_models <- function(a, w, x, zip, a.vals, fmla, trim = 0.05) {
  
  if (trim < 0 | trim > 0.5)
    stop("trim < 0 | trim > 0.5")
  
  ## Matching Estimator
  match_pop <- generate_pseudo_pop(Y = zip, w = a, c = x,
                                   ci_appr = "matching",
                                   pred_model = "sl",
                                   gps_model = "parametric",
                                   use_cov_transform = TRUE,
                                   transformers = list("pow2", "pow3"),
                                   sl_lib = c("m_xgboost"),
                                   params = list(xgb_nrounds = c(50)),
                                   nthread = 12, # number of cores, you can change,
                                   covar_bl_method = "absolute",
                                   covar_bl_trs = 0.1,
                                   covar_bl_trs_type = "mean",
                                   trim_quantiles = c(trim, 1 - trim), # trimed, you can change,
                                   optimized_compile = TRUE, #created a column counter for how many times matched,
                                   max_attempt = 5,
                                   matching_fun = "matching_l1",
                                   delta_n = 0.2, # you can change this to the one you used in previous analysis,
                                   scale = 1.0)
  
  # merge individual level data
  pseudo <- match_pop$pseudo_pop
  match_data <- merge(w, data.frame(zip = pseudo$zip,
                                    year = pseudo$year,
                                    a = pseudo$w,
                                    counter = pseudo$counter), 
                      by = c("zip", "year"), all = FALSE)
  match_data <- subset(match_data, counter > 0)
  match_curve <- mgcv::bam(fmla, data = match_data, offset = log(time_count), family = poisson(link = "log"), weights = counter)
  
  estimate <- sapply(a.vals, function(a.tmp, ...) {
    
    match_estimate <- predict(match_curve, newdata = data.frame(a = a.tmp, w), type = "response")
    return(weighted.mean(match_estimate, w = exp(w$time_count), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, match_data = match_data,
              adjusted_corr_results = match_pop$adjusted_corr_results, 
              original_corr_results = match_pop$original_corr_results))
  
}

tmle_glm <- function(a, w, x, y, offset, a.vals, trim = 0.01){
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  pm_id <- which(colnames(w) == "pm25")
  
  # estimate nuisance outcome model with splines
  fmla <- formula(paste0( "y ~ ns(a, 4) ", paste0(colnames(w[,-pm_id]), collapse = "+")))
  mumod <- glm(fmla, data = data.frame(w, a = w[,pm_id]), offset = offset, family = poisson(link = "log"))
  muhat <- exp(log(mumod$fitted.values) - offset)
  
  # estimate nuisance GPS parameters with lm
  pimod <- lm(a ~ 0 + ., data = data.frame(x))
  pimod.vals <- c(pimod$fitted.values, predict(pimod, newdata = data.frame(w)))
  pi2mod.vals <- sigma(pimod)^2
  
  # parametric density
  # pihat <- dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  # })
  # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat, na.rm = T)), x = a)$y
  # phat[phat<0] <- 1e-4
  
  # nonparametric denisty
  a.std <- c(c(a, w[,pm_id]) - pimod.vals) / sqrt(pi2mod.vals)
  dens <- density(a.std[1:n])
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
  
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - pimod.vals) / sqrt(pi2mod.vals)
    approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
  })
  
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n,], na.rm = T)), x = c(a, w[,pm_id]))$y
  phat[phat<0] <- 1e-4
  
  # TMLE update
  nsa <- ns(a, df = 4, intercept = TRUE)
  weights <- phat/pihat
  trim0 <- quantile(weights[1:n], trim)
  trim1 <- quantile(weights[1:n], 1 - trim)
  weights[weights < trim0] <- trim0
  weights[weights > trim1] <- trim1
  base <- predict(nsa, newx = w[,pm_id])*weights[-(1:n)]
  new_mod <- glm(y ~ 0 + base, offset = log(muhat) + offset, 
                 family = poisson(link = "log"))
  param <- coef(new_mod)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    print(k)
    w$a <- a.vals[k]
    muhat.tmp <- predict(mumod, newdata = data.frame(w))
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp[1:n], na.rm = TRUE)/pihat.tmp[-(1:n)])
    wts[wts < trim0] <- trim0
    wts[wts > trim1] <- trim1
    mat <- predict(nsa, newx = rep(a.tmp, length(wts)))*wts
    return(weighted.mean(exp(log(muhat.tmp) + c(mat%*%param)), w = exp(offset), na.rm = TRUE))
  })
  
  return(list(estimate = estimate, weights = weights[-(1:n)]))
  
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
    xlab("Covariates") + ylab("Absolute Pearson Correlation") +
    ylim(0, 0.35) +
    guides(color = guide_legend(title = "GPS Adjusting")) +
    theme_bw() + # use a white background
    ggtitle(main)
  
  return(fp)
  
}
