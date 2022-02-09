
fit_models <- function(a, x, y, offset, a.vals) {
    
    gam_curve <- mgcv::bam(y ~ s(a, bs = 'cr', k = 3) - a + .,
                           data = data.frame(y = y, a = a, x), offset = offset,
                           family=poisson(link = "log"))
    
    match_pop <- generate_pseudo_pop(Y = y, w = a, c = x,
                                     ci_appr = "matching", pred_model = "sl",
                                     gps_model = "parametric", use_cov_transfoqd = TRUE,
                                     transfoqders = list("pow2", "pow3"), sl_lib = c("m_xgboost"),
                                     params = list(xgb_nrounds=c(50)), nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute", covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.025,0.975), # trimmed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 5, matching_fun = "matching_l1",
                                     delta.n = (max(a.vals) - min(a.vals)), scale = 1.0)
    
    match_pop$pseudo_pop$offset <- offset[match_pop$pseudo_pop$row_index]
    
    match_curve <- mgcv::bam(Y ~ s(w, bs = 'cr', k = 3), data = match_pop$pseudo_pop,
                             offset = offset, family = poisson(link = "log"), weights = counter)
    
    np <- np_est(y = y, a = a, x = x, offset = offset)
    
    muhat <- np$muhat
    mhat <- np$mhat
    pihat <- np$pihat
    phat <- np$phat
    
    pihat[pihat <= quantile(pihat, 0.025)] <- quantile(pihat, 0.025)
    pihat[pihat > quantile(pihat, 0.975)] <- quantile(pihat, 0.975)
    phat[phat <= quantile(phat, 0.025)] <- quantile(phat, 0.025)
    phat[phat > quantile(phat, 0.975)] <- quantile(phat, 0.975)

    pseudo <- exp(offset)*((y/exp(offset) - muhat)/(pihat/phat) + mhat)
    pseudo[pseudo < 0] <- 0

    dr_curve <- mgcv::bam(pseudo ~ s(a, bs = 'cr', k = 3), offset = offset, family = quasipoisson(link = "log"))
    
    return(list(gform = gam_curve, match = match_curve, dr_curve = dr_curve))
  
}

# Nonparametric estimation
np_est <- function(a, y, x, offset){
  
  a.vals <- seq(min(a.vals), max(a.vals), length.out = 100)
  
  # set up evaluation points & matrices for predictions
  xa.new <- rbind(cbind(x, a), 
                  cbind(x[rep(1:n, length(a.vals)), ],
                        a = rep(a.vals, rep(n, length(a.vals)))))
  x.new <- xa.new[, -dim(xa.new)[2]]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x) <- colnames(x.new)
  xa.new <- data.frame(xa.new)
  
  # estimate nuisance functions via super learner
  pimod <- SuperLearner(Y = a, X = data.frame(x), SL.library = "SL.xgboost", newX = x.new)
  pimod.vals <- c(pimod$SL.predict)
  pi2mod <- SuperLearner(Y = (a - pimod.vals[1:n])^2, X = x, SL.library = sl.lib, newX = x.new)
  pi2mod.vals <- c(pi2mod$SL.predict)
  mumod <- glm(y ~ ., data = data.frame(x, y = y, a = a), offset = offset, family = poisson(link = "log"))
  mumod.vals <- predict(mumod, newdata = xa.new, type = "response")
  
  # construct estimated pi/p and mu/m values
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
  pihat.vals <- approx(density(a.std)$x, density(a.std[1:n])$y, xout = a.std)$y / sqrt(pi2mod.vals)
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  mhat <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  
  out <- list(muhat = muhat, mhat = mhat, pihat = pihat, phat = phat)
  
  return(out)
  
}
