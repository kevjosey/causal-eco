# estimate glm and combine with weights
gam_models <- function(y, a, w, weights, id, a.vals, log.pop = NULL, trunc = 0.01, ...) {
  
  if (is.null(log.pop))
    log.pop <- rep(0, nrow(x))
  
  ybar <- y/exp(log.pop)
  ybar[y > exp(log.pop)] <- 1 - .Machine$double.eps
  
  # estimate nuisance outcome model with glm + splines
  mumod <- gam(ybar ~ s(a, df = 6) + . + (a:.) - a, weights = exp(log.pop), model = FALSE,
               data = data.frame(ybar = ybar, a = a, w), family = quasipoisson())
  muhat <- mumod$fitted.values
  
  # predictions along a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    wa.tmp <- data.frame(a = a.tmp, w)
    predict(mumod, newdata = wa.tmp, type = "response")
    
  })
  
  # pseudo outcome
  resid <- c(ybar - muhat)*weights
  
  out <- list(resid = resid, muhat.mat = muhat.mat, 
              id = id, log.pop = log.pop)
  
  return(out)
  
}
