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

computeDisc <- function(cf, real=TRUE, subtract.vev=TRUE) {
  T <- cf$Time
  tcf <- cf$cf
  if(!real) tcf <- cf$icf
  N <- dim(tcf)[1]

  i <- c(1:T)
  i2 <- i
  Cf <- array(0, dim=c(N, T/2+1))
  
  if(subtract.vev) {
    ## compute vev first
    ## mean over all gauges and times
    vev <- mean(tcf)
  }
  Cf[,1] <- apply((tcf - vev)*(tcf - vev), 1, mean)
  for(dt in c(1:(T/2))) {
    i2 <- (i2) %% 48 + 1
    Cf[,1+dt] <- apply((tcf[,i] - vev)*(tcf[,i2] - vev), 1, mean)
  }
  return(invisible(Cf))
}
