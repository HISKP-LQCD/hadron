compute.weights <- function(err, pvalues) {
    return(pvalues^2 * min(err, na.rm=TRUE)^2/err^2)
  }

compute.boots <- function(res, index=1) {
  j <- 5
  if(index %% 2 == 0) j <- 6
  pvalues <- as.vector((1-2*abs(res[1,,j]-0.5)))
  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  return(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE))
}


summarise.rho <- function(ens, frame) {
  filelist <- Sys.glob(paste("./", "rhoana.*", ens, frame, ".Rdata", sep=""))
  
  N <- length(filelist)
  cat("processing", N, "data files for ensemble", ens, "\n")
  load(filelist[1])
  R <- pc.matrixfit$boot.R
  rr <- c(2:(R+1))
  
  res <- array(0, dim=c(R+1, N, 6))
  
  for(i in c(1:N)) {
    load(filelist[i])
    
    ## ECM1
    res[1,i,1] <- gs$Ecm
    res[rr,i,1] <- gs$Ecmboot
    ## ECM2
    res[1,i,2] <- fes$Ecm
    res[rr,i,2] <- fes$Ecmboot
    
    ##delta1
    res[1,i,3] <- gs$delta
    res[rr,i,3] <- gs$deltaboot
    ##delta1
    res[1,i,4] <- fes$delta
    res[rr,i,4] <- fes$deltaboot
    ## Qval1
    res[1,i,5] <- pc.matrixfit$Qval
    ## Qval2
    res[1,i,6] <- pc2.matrixfit$Qval
  }

  save(res, file=paste(path, "res.", ens, frame, ".Rdata", sep=""))
  return(invisible(res))
}

##source("../../pipi/analyse_pipi.R")

estimate.error <- function(res, index=1, prob=c(0.1573, 0.8427), main) {
  j <- 5
  if(index %% 2 == 0) j <- 6
  pvalues <- as.vector((1-2*abs(res[1,,j]-0.5)))
  
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

compute.error.rho <- function(res, prob=c(0.1573, 0.8427)) {
  Ecm1 <- estimate.error(res, index=1, prob=prob, main="Ecm1")
  Ecm2 <- estimate.error(res, index=2, prob=prob, main="Ecm2")
  delta1 <- estimate.error(res, index=3, prob=prob, main="delta1")
  delta2 <- estimate.error(res, index=4, prob=prob, main="delta2")
  
  return(list(Ecm1=Ecm1, Ecm2=Ecm2, delta1=delta1, delta2=delta2))
}
