computeDisc <- function(cf, real=TRUE, subtract.vev=TRUE) {
  T <- cf$Time
  ## the real part of the correlator matrix
  tcf <- cf$cf
  nrSamples <- cf$nrSamples
  ## the imaginary part of the correlator matrix
  if(!real) tcf <- cf$icf
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
    vev <- mean(cf$cf)
  }
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
    tcf <- apply(mtcf, c(1,3), sum)
    for(dt in c(0:(T/2))) {
      ## shift the index array by 1 to the left
      Cf[,1+dt] <- apply(tcf[i,]*tcf[i2,], 2, mean)
      Cf[,1+dt] <- Cf[,1+dt] - apply(apply(mtcf[i,,]*mtcf[i2,,], c(2,3), mean), 2, sum)
      i2 <- (i2) %% T + 1
    }
    Cf <- Cf/nrSamples/(nrSamples-1)
  }
  cf <- list(cf=Cf, icf=NULL, Time=T, nrStypes=1, nrObs=1, nrSamples=nrSamples)
  return(invisible(cf))
}
