
## This function estimates the log2 fold change for robust linear mixed models.

get.rlmer.log2fc = function(all.fits, d.grid)
{
  stats = d.grid
  for ( v in c("Diff") ) {
    stats[,v] = NA
  }

  stats = by( d.grid, d.grid$pert_time, function(d0){
    d0$abundance <- 0
    X = model.matrix(terms(all.fits), data=d0)
    X = sweep(X, 2, X[1, ], FUN = "-")[-1, , drop = FALSE]
    out = d0[-1, -which(colnames(d0) == "abundance"), drop = FALSE]

    out$Diff = as.vector( X%*%unlist(coef(all.fits)[[1]][1, ]) )
    return(out)
  })
  
  stats = do.call("rbind", stats)
  return(stats)
}
