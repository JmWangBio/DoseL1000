
## This function estimates log2 fold change for generalized additive mixed models.

get.gamer.log2fc = function(all.fits, d.grid)
{
  stats = d.grid
  for ( v in c("Diff") ) {
    stats[,v] = NA
  }
  
  stats = by( d.grid, d.grid$pert_time, function(d0){
    X = predict(all.fits, newdata=d0, type="lpmatrix")
    X = sweep(X, 2, X[1,], FUN="-")[-1, , drop=FALSE]
    out = d0[-1, , drop=FALSE]
    
    out$Diff = as.vector( X%*%coef(all.fits) )
    return(out)
  })
  
  stats = do.call("rbind", stats)
  return(stats)
}
