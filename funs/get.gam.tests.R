
## This function conducts t-tests to compare the mean responses for generalized additive models.

get.gam.tests = function(all.fits, d.grid)
{
  stats = d.grid
  for ( v in c("Diff", "SE", "tStat", "pval", "degf") ) {
    stats[,v] = NA
  }
    
  stats = by(d.grid, d.grid$pert_time, function(d0){
    X = predict(all.fits, newdata=d0, type="lpmatrix")
    X = sweep(X, 2, X[1,], FUN="-")[-1, , drop=FALSE]
    out = d0[-1, , drop=FALSE]
    
    out$Diff = as.vector( X%*%coef(all.fits) )      
    out$SE = sqrt( diag(X%*%all.fits$Vp%*%t(X)) )
    out$tStat = out$Diff/out$SE
    out$pval = 2*pt(-abs(out$tStat), 
                    all.fits$df.residual)
    out$degf = all.fits$df.residual
    return(out)
  })
  
  stats = do.call("rbind", stats)
  return(stats)
}
