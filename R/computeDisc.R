computeDisc <- function(cf, cf2,
                        real=TRUE, real2 = TRUE,
                        subtract.vev=TRUE, subtract.vev2=TRUE,
                        subtract.equal = TRUE,
                        use.samples) {
  T <- cf$Time
  ## the real part of the correlator matrix
  tcf <- cf$cf
  if(!real) tcf <- cf$icf
  
  nrSamples <- cf$nrSamples
  if(!missing(use.samples) && !(use.samples > nrSamples) && (use.samples > 0) ) {
    nrSamples <- use.samples
  }
  sindex <- c(1:nrSamples)

  
  ## the imaginary part of the correlator matrix

  ## number of gauges
  N <- dim(tcf)[3]
  ## index array for t
  i <- c(1:T)
  ## index array for t'
  i2 <- i
  ## space for the correlator
  Cf <- array(0., dim=c(N, T/2+1))
  vev <- 0.

  if(subtract.vev) {
    ## compute vev first
    ## mean over all gauges and times
    if(nrSamples == 1) vev <- mean(cf$cf)
    else vev <- mean(cf$cf[,sindex,])
  }

  if(missing(cf2)) {
    ## here we compute the actual correlation
    if(nrSamples == 1) {
      tcf <- tcf - vev
      for(dt in c(0:(T/2))) {
        ## shift the index array by 1 to the left
        Cf[,1+dt] <- apply(tcf[i,1,]*tcf[i2,1,], 2, mean)
        i2 <- (i2) %% T + 1
      }
    }
    else {
      ## re-order data
      mtcf <- tcf - vev
      ## average over samples, tcf has dim(T,N)
      tcf <- apply(mtcf[,sindex,], c(1,3), sum)
      for(dt in c(0:(T/2))) {
        Cf[,1+dt] <- apply(tcf[i,]*tcf[i2,], 2, mean)
        ## subtract product of equal samples
        if(subtract.equal) Cf[,1+dt] <- Cf[,1+dt] - apply(apply(mtcf[i,sindex,]*mtcf[i2,sindex,], c(2,3), mean), 2, sum)
        ## shift the index array by 1 to the left
        i2 <- (i2) %% T + 1
      }
      Cf <- Cf/nrSamples/(nrSamples-1)
      if(!subtract.equal) CF <- CF/nrSamples
    }
  }
  ## now the case of cross-correlators
  else {
    ## sanity checks
    if(cf2$nrSamples != nrSamples) {
      stop("sample numbers in two loops do not agree... Aborting...!\n")
    }
    if(cf2$Time != T) {
      stop("time extend in two loops does not agree... Aborting...!\n")
    }
    if(dim(cf2$cf)[3] != N) {
      stop("number of gauges for the two loops does not agree... Aborting...!\n")
    }
    if(!real2) tcf2 <- cf2$icf
    else tcf2 <- cf2$cf
    vev2 <- 0.
    if(subtract.vev2) {
      ## compute vev first
      ## mean over all gauges and times
      if(nrSamples == 1) vev2 <- mean(tcf2)
      else vev2 <- mean(cf2$cf[,sindex,])
    }

    if(nrSamples == 1) {
      tcf <- tcf - vev
      tcf2 <- tcf2 - vev2
      for(dt in c(0:(T/2))) {
        Cf[,1+dt] <- apply(tcf[i,1,]*tcf2[i2,1,], 2, mean)
        ## shift the index array by 1 to the left
        i2 <- (i2) %% T + 1
      }
    }
    else {
      ## re-order data
      mtcf <- tcf - vev
      mtcf2 <- tcf2 - vev2
      ## average over samples, tcf has dim(T,N)
      tcf <- apply(mtcf[,sindex,], c(1,3), sum)
      tcf2 <- apply(mtcf2[,sindex,], c(1,3), sum)
      for(dt in c(0:(T/2))) {
        Cf[,1+dt] <- apply(tcf[i,]*tcf2[i2,], 2, mean)
        ## subtract product of equal samples
        if(subtract.equal) Cf[,1+dt] <- Cf[,1+dt] - apply(apply(mtcf[i,sindex,]*mtcf2[i2,sindex,], c(2,3), mean), 2, sum)
        ## shift the index array by 1 to the left
        i2 <- (i2) %% T + 1
      }
      Cf <- Cf/nrSamples/(nrSamples-1)
      if(!subtract.equal) CF <- CF/nrSamples
    }


  }
  cf <- list(cf=Cf, icf=NULL, Time=T, nrStypes=1, nrObs=1, nrSamples=nrSamples)
  return(invisible(cf))
}
