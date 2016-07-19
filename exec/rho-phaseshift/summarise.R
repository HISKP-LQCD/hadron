compute.weights <- function(err, pvalues) {
    return(pvalues^2 * min(err, na.rm=TRUE)^2/err^2)
  }

compute.boots <- function(res, index=1) {
  pvalues <- as.vector((1-2*abs(res[1,,4]-0.5)))
  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  return(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE))
}


summarise.rho <- function(ens, frame, irrep, PC, threshold=-0.1, shift=pi) {
  filelist <- Sys.glob(paste("./", "rhoana.", PC, ".*", ens, ".", frame, ".", irrep, ".Rdata", sep=""))
  
  N <- length(filelist)
  cat("processing", N, "data files for ensemble", ens, "and", PC, "\n")
  load(filelist[1])
  R <- pc.matrixfit$boot.R
  rr <- c(2:(R+1))
  
  res <- array(0, dim=c(R+1, N, 4))
  
  for(i in c(1:N)) {
    load(filelist[i])
    
    ## ECM
    res[1,i,1] <- gs$Ecm
    res[rr,i,1] <- gs$Ecmboot
    
    ##delta
    res[1,i,2] <- gs$delta
    res[rr,i,2] <- gs$deltaboot

    ## shift where neccessary
    ii <- which(res[,i,2] < threshold)
    res[ii,i,2] <- res[ii,i,2]+shift

    ##tan(delta)
    res[1,i,3] <- gs$tandelta
    res[rr,i,3] <- gs$tandeltaboot

    ## Qvals
    res[1,i,4] <- pc.matrixfit$Qval
    res[2,i,4] <- pion.matrixfit$Qval
  }

  save(res, file=paste("res", ens, frame, irrep, "Rdata", sep="."))
  return(invisible(res))
}

estimate.error <- function(res, index=1, prob=c(0.1573, 0.8427), main) {
  pvalues <- as.vector((1-2*abs(res[1,,4]-0.5)))
  
  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  x <- weighted.quantile(as.vector(res[1,,index]), prob=c(0.5), w=w, na.rm=TRUE)
  
  if(interactive()) X11()
  weighted.hist(x=as.vector(res[1,,index]), w=w, main=main, na.rm=TRUE)
  ## statistical
  x[2] <- sd(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE), na.rm=TRUE)
  ## systematic error
  ## lower and upper
  x[c(3:4)] <- weighted.quantile(as.vector(res[1,,index]), w=w, prob=prob, na.rm=TRUE)-x[1]
  return(x)
}

compute.error.rho <- function(res, prob=c(0.1573, 0.8427), PC) {
  Ecm <- estimate.error(res, index=1, prob=prob, main="Ecm")
  delta <- estimate.error(res, index=2, prob=prob, main="delta")
  tandelta <- estimate.error(res, index=3, prob=prob, main="tandelta")
  
  return(list(Ecm=Ecm, delta=delta, tandelta=tandelta, PC=PC))
}

shift.delta <- function(res, threshold=-0.1, shift=pi) {
  ii <- which(res[,,2] < threshold)
  res[ii] <- res[ii] + shift
  return(res)
}
