
## This function conducts t-tests to compare the mean responses for linear models.

get.lm.tests = function(all.fits, d.grid)
{
  stats = d.grid
  for ( v in c("Diff", "SE", "tStat", "pval", "degf") ) {
    stats[,v] = NA
  }

  stats = by( d.grid, d.grid$pert_time, function(d0){
    d0$abundance <- 0
    X = model.matrix(all.fits$terms, data=d0)
    X = sweep(X, 2, X[1, ], FUN = "-")[-1, , drop = FALSE]
    out = d0[-1, -which(colnames(d0) == "abundance"), drop = FALSE]
    
    out$Diff = as.vector( X%*%coef(all.fits) )
    out$SE = sqrt( diag(X%*%vcov(all.fits)%*%t(X)) )
    out$tStat = out$Diff/out$SE
    out$pval = 2*pt(-abs(out$tStat), 
                    summary(all.fits)$df[2])
    out$degf = summary(all.fits)$df[2]
    return(out)
  })
  
  stats = do.call("rbind", stats)
  return(stats)
}
