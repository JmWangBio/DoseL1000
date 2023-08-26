
## This function calculates the estimates and standard errors of efficacy and potency.

get.gam.stats = function(fit, B = 10001)
{
  design.mat = function(lConc, fit, ft)
  {
    predict(fit, newdata=data.frame(logConc = lConc, pert_time = ft), 
            type="lpmatrix")
  }   
  
  d = fit$model
  x.seq = seq(min(d$logConc), max(d$logConc), length=B)
  
  st = list()
  for ( ft in levels(d$pert_time) ){
    pred = as.data.frame(predict(fit, 
                                 newdata = data.frame(logConc = x.seq, 
                                                      pert_time = ft), 
                                 se.fit = TRUE))
    
    ## Get difference between "best" (smallest/largest) predicted value and predicted value @ dose = 0
    index = which.max(abs( pred$fit - pred$fit[1] ))
    
    st[[ft]] = data.frame(pert_time = ft, 
                          Diff = pred$fit[index] - pred$fit[1])
    
    ## Get standard error of the best difference
    X = predict(fit, newdata=data.frame(logConc=x.seq[c(1, index)], 
                                        pert_time=ft), 
                type = "lpmatrix")
    X = X[2, ] - X[1, ]
    st[[ft]]$SE.diff = as.vector( sqrt( X%*%fit$Vp%*%X ) )
    
    ## Get t-stat of max log2 fold change
    st[[ft]]$tstat.diff <- st[[ft]]$Diff / st[[ft]]$SE.diff
    
    ## Get p-value of max log2 fold change
    st[[ft]]$pval.diff <- 2*pt(-abs(st[[ft]]$tstat.diff), 
                               fit$df.residual)
    
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
    st[[ft]]$lx50 = x.seq[index.50]    
    
    ## Get standard error of lx50    
    X = as.vector(predict(fit, 
                          newdata = data.frame(logConc = st[[ft]]$lx50, 
                                             pert_time = ft), 
                          type = "lpmatrix"))    
    
    h = 1e-8 * ifelse( abs(st[[ft]]$lx50) == 0, 1, abs(st[[ft]]$lx50) )
    H = as.vector( design.mat(lConc = st[[ft]]$lx50 + h, 
                              fit = fit, ft = ft) - 
                     design.mat(lConc = st[[ft]]$lx50 - h, 
                                fit = fit, ft = ft) )/(2*h)
    G = -X / sum( coef(fit) * H )      
    st[[ft]]$SE.lx50 = as.vector( sqrt( G%*%fit$Vp%*%G ) )
  }
  stats = do.call('rbind', st)
  
  return(st)
}
