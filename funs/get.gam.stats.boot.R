
## This function calculates the estimates and standard errors of efficacy and potency via bootstrapping.

get.gam.stats.boot = function(fit, B = 10000, bootn = 1000)
{
  d <- fit$model
  d$ypred <- predict(fit)
  sig <- sqrt(fit$sig2) 
  x.seq = seq(min(d$logConc), max(d$logConc), length = B)
  
  boot_lst <- list()
  k <- 1
  for (i in 1:bootn) {
    ## resample data
    d$yboot <- rnorm(nrow(d), 
                     mean = d$ypred,
                     sd = sig)
    
    ## check if there is more than one time point
    if (length(unique(d$pert_time)) > 1) {
      boot_fit <- gam(yboot ~ pert_time + s(logConc, k = 4,
                                            by = pert_time),
                      data = d)
    } else {
      boot_fit <- gam(yboot ~ s(logConc, k = 4,
                                by = pert_time),
                      data = d)
    }
    
    for ( ft in levels(d$pert_time) ){
      ## obtain model fitted to simulated data
      pred = data.frame(fit = predict(boot_fit, 
                                      newdata = data.frame(logConc = x.seq, 
                                                           pert_time = ft)))
      
      ## Get difference between "best" (smallest/largest) predicted value and predicted value @ dose = 0
      index = which.max(abs( pred$fit - pred$fit[1] ))
      
      boot_lst[[k]] = data.frame(iter = i,
                                 pert_time = ft, 
                                 diff = pred$fit[index] - pred$fit[1])
      
      ## Get y50
      y50 = 0.5*( pred$fit[1] + pred$fit[index] )
      
      ## Get lx50
      index.50 = which.min(abs(pred$fit - y50))
      while ((!all(pred$fit[1:(index.50-1)] < pred$fit[index.50])) & 
             (!all(pred$fit[1:(index.50-1)] > pred$fit[index.50])) &
             (index.50 > 1)) {
        index.50 = which.min(
          abs(pred$fit[1:(index.50-1)] - y50)
        )
      }
      boot_lst[[k]]$lx50 = x.seq[index.50]    
      k <- k+1
    }
  }
  
  ## Get estimate and standard error of the best max log2 fold change and lx50
  boot_df <- do.call('rbind', boot_lst)
  st = boot_df %>% 
    group_by(pert_time) %>%
    dplyr::summarise(mean.diff = mean(diff),
                     se.diff = sd(diff),
                     mean.lx50 = mean(lx50),
                     se.lx50 = sd(lx50))
  
  return(st)
}
