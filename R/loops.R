loops <- function() {
  v <- 0.
  ngauge <- 10
  by <- 10
  L <- 48
  basename <- "cvc_2pt_disc_vv.0"
  vev <- FALSE
  if(vev) {
    for(i in seq(from=500, to=500+(ngauge-1)*by, by=by)) {
      filename <- paste(basename, i, sep="")
      cat(filename, "\n")
      data <- read.table(filename)
      pi0 <- data[data$V2 == 8,]
      v <- v + mean(pi0$V5)
    }
    v <- v/ngauge
    cat("v=", v, "\n")
  }
  S <- 24
  T <- 96
  res <- rep(0., times=T)
  ii1 <- c(0:(T-1))
  
  
  for(i in seq(from=500, to=500+(ngauge-1)*by, by=by)) {
    filename <- paste(basename, i, sep="")
    data <- read.table(filename)
    pi0 <- data[data$V2 == 8,]
    temp <- pi0$V6-v
    
    for(s1 in 0:(S-1)) {
      s2 <- s1+1
      while(s2 < S) {
        for(dt in 0:(T-1)) {
          ii2 <- (ii1+dt)%%T
          res[dt+1] <- res[dt+1] + mean((temp[ii1+s1*T+1])*(temp[ii2+s2*T+1]))
        }
        s2 <- s2 + 1
      }
    }
  }
  res <- res/((S-1)*S/2)/ngauge/L^3
  cat(res, "\n")
}

computeDisc <- function(cf, real=TRUE, subtract.vev=TRUE, avoid.equal=TRUE, L) {
  T <- cf$Time
  if(missing(L)) L <- T/2
  ## the real part of the correlator matrix
  tcf <- cf$cf
  nrSamples <- cf$nrSamples
  if(nrSamples == 1) avoid.equal <- FALSE
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
    Cf <- Cf/nrSamples/(nrSamples-1)/L^3
  }
  cf <- list(cf=Cf, icf=NULL, Time=T, nrStypes=1, nrObs=1, nrSamples=nrSamples, L=L)
  return(invisible(cf))
}
