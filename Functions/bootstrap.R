# function for creating cluster bootstrap samples
bootstrap_data <- function(data, index, u.zip) {
  
  n.zip <- length(u.zip)
  boot <- data.frame()
  
  aa <- u.zip[index]
  aa <- aa[which(aa %in% data$zip)]
  bb <- table(aa)
  
  for (j in 1:max(bb)) {
    
    cc <- data[data$zip %in% names(bb[bb == j]),]
    
    for (k in 1:j) {
      cc$boot.id <- paste(cc$id, k, sep = "-")
      boot <- rbind(boot, cc)
    }
    
  }
  
  return(boot)
  
}

# function that implements model_erc for a bootstrap sample
bootable <- function(w, x, z, a.vals, region) {
  
  u.zip <- unique(w$zip)
  m <- length(u.zip)/log(length(u.zip)) # for m out of n bootstrap
  index <- sample(1:length(u.zip), ceiling(m), replace = TRUE)  # initialize bootstrap  index
  x_boot <- try(bootstrap_data(data = x, index = index, u.zip = u.zip)) # bootstrap site-level data
  w_boot <- try(bootstrap_data(data = w, index = index, u.zip = u.zip)) # bootstrap individual-level data
  
  if (!inherits(x_boot, "try-error") & !inherits(w_boot, "try-error")) {
    
    boot <- try(gam_om(x = x_boot, w = w_boot, z = z, a.vals = a.vals,
                       se.fit = FALSE, boot = TRUE, region = region))
    
    if (!inherits(boot, "try-error"))
      return(list(a.vals = a.vals, erc = boot$est_data[,2], ed = boot$excess_death[,2]))
    else
      return(list(erc = rep(a.vals = a.vals, NA, length = length(a.vals)), ed = rep(NA, length = length(a.vals))))
    
  } else
    return(list(a.vals = a.vals, erc = rep(NA, length = length(a.vals)), ed = rep(NA, length = length(a.vals))))
  
}
