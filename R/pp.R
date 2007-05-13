pp <- function(filename, skip=0, from, to, S=1.5, A=0.01, m=0.01, plot=FALSE, debug=FALSE) {

  if(!missing(filename)) {
    psscar <- read.table(file=filename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option psfilename is mandatory!")
  }
  T2 <- (max(psscar$t)-min(psscar$t)+1)
  if(missing(from) || from < 2) {
    from <- 2
  }
  if(missing(to) || to > T2/2) {
    to <- T2/2
  }

  Z <- array(psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  result <- data.frame(t = array(0., dim=c((to-from+1))), mass = array(0., dim=c((to-from+1))), dmass = array(0., dim=c((to-from+1))),
                       amp = array(0., dim=c((to-from+1))), damp = array(0., dim=c((to-from+1))), ddmass = array(0., dim=c((to-from+1))),
                       masstauint = array(0., dim=c((to-from+1))), massdtauint = array(0., dim=c((to-from+1))), ddamp = array(0., dim=c((to-from+1))),
                       amptauint = array(0., dim=c((to-from+1))), ampdtauint = array(0., dim=c((to-from+1))))
  i=1
  cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")
  for(cutoff in from:to) {
    options(error = function() {cat("\nAn Error occured in uwerrderived function! Continuing with next timeslice...\nResults of this timeslice must not be used!!\n")})

    cat("fitting from timeslice", cutoff, "!\n")
    try(mass <- uwerrderived(getmass, data=Z[(cutoff):(T2+2-cutoff),skip:(length(psscar$ps)/T2)], S, pl=debug, T2, cutoff, A, m, debug=debug))
    try(amp <- uwerrderived(getamp, data=Z[(cutoff):(T2+2-cutoff),skip:(length(psscar$ps)/T2)], S, pl=debug, T2, cutoff, A, m, debug=debug))

    result$t[i] <- cutoff-1
    result$mass[i] <- mass$value[1]
    result$dmass[i] <- mass$dvalue[1]
    result$amp[i] <- amp$value[1]
    result$damp[i] <- amp$dvalue[1]
    result$ddmass[i] <- mass$ddvalue[1]
    result$masstauint[i] <- mass$tauint[1]
    result$massdtauint[i] <- mass$dtauint[1]
    result$ddamp[i] <- amp$ddvalue[1]
    result$amptauint[i] <- amp$tauint[1]
    result$ampdtauint[i] <- amp$dtauint[1]
    i=i+1
  }
  attr(result, "class") <- c("massfit", "ampfit", "data.frame")
  if(plot == TRUE) {
    X11()
    plot(result[,1],result[,4], xlab = "t", ylab = "A", main = "One state cosh fit for the amplitude", ylim=c(min(result[,4]-2*result[,5]),max(result[,4]+2*result[,5])))
    arrows(result[,1],result[,4]-result[,5],result[,1],result[,4]+result[,5], length=0.01,angle=90,code=3)
    X11()
    plot(result[,1],result[,2], xlab = "t", ylab = "m", main = "One state cosh fit for the mass", ylim=c(min(result[,2]-2*result[,3]),max(result[,2]+2*result[,3])))
    arrows(result[,1],result[,2]-result[,3],result[,1],result[,2]+result[,3], length=0.01,angle=90,code=3)
  }
  rm(psscar)
  rm(Z)
  return(invisible(result))
}


effmass <- function(data, timeextent, t) {
  mass <- invcosh((data[1]+data[4])/(data[2]+data[3]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

effmass2 <- function(data, timeextent, t) {
  mass <- invcosh((data[1])/(data[2]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

effectivemass <- function(from, to, Time, Z, pl=TRUE, S,...) {
  L <- (to-from+1)
  i <- 1
  result <- data.frame(t = array(0.,dim=c(L)), mass = array(0.,dim=c(L)), dmass = array(0.,dim=c(L)),
                       ddmass = array(0.,dim=c(L)), tauint = array(0.,dim=c(L)), dtauint = array(0.,dim=c(L)))

  for(t in from:to) {
    try(mass <- uwerrderived(effmass2, Z[t:(t+1),], S=S, pl=F, timeextent=Time, t=t, ...))

    result$t[i] <- t-1
    result$mass[i] <- mass$value[1]
    result$dmass[i] <- mass$dvalue[1]
    result$ddmass[i] <- mass$ddvalue[1]
    result$tauint[i] <- mass$tauint[1]
    result$dtauint[i] <- mass$dtauint[1]
    i = i+1
  }

  rm(mass)
  rm(Z)
  attr(result, "class") <- c("massfit", "data.frame")
  if(pl == T) {
    X11()
    plot(result)
  }
  return(invisible(result))
}

ppeffectivemass <- function(filename, from, to, skip=0, S=1.5, plotit=F) {
  if(!missing(filename)) {
    psscar <- read.table(file=filename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename is mandatory!")
  }
  T2 <- (max(psscar$t)-min(psscar$t)+1)
  if(missing(to) || to > (T2/2)) {
    to <- (T2/2)
  }
  if(missing(from) || from < 2) {
    from <- 2
  }
  Z <- array(psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  L <- (to-from+1)

  i <- 1
  result <- data.frame(t = array(0.,dim=c(L)), mass = array(0.,dim=c(L)), dmass = array(0.,dim=c(L)),
                       ddmass = array(0.,dim=c(L)), tauint = array(0.,dim=c(L)), dtauint = array(0.,dim=c(L)))
  cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")
  rm(psscar)
                                        # t counted form 1!!

  for(t in from:to) {
    try(mass <- uwerrderived(effmass, data=rbind(Z[t:(t+1),skip:length(Z[1,])], Z[(T2+2-t-1):(T2+2-t),skip:length(Z[1,])]), S, pl=F, timeextent=T2, t=t))

    result$t[i] <- t-1
    result$mass[i] <- mass$value[1]
    result$dmass[i] <- mass$dvalue[1]
    result$ddmass[i] <- mass$ddvalue[1]
    result$tauint[i] <- mass$tauint[1]
    result$dtauint[i] <- mass$dtauint[1]
    i = i+1
  }

  rm(mass)
  rm(Z)
  attr(result, "class") <- c("massfit", "data.frame")
  if(plotit == T) {
    plot(result)
  }
  return(invisible(result))
}

ppcorr <-  function(filename, skip = 0, S=1.5, plotit=F) {
  if(!missing(filename)) {
    psscar <- read.table(file=filename, col.names=c("t","ps","mps"), header=F)
#    psscar <- read.table(file=filename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename is mandatory!")
  }
  T2 <- (max(psscar$t)-min(psscar$t)+1)
  Z <- array((psscar$ps+psscar$mps)/2, dim=c(T2,length(psscar$ps)/T2))
#  Z <- array(psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  L <- T2

  i <- 1

  result <- data.frame(t = array(0.,dim=c(T2)), corr = array(0.,dim=c(T2)), dcorr = array(0.,dim=c(T2)),
                       ddcorr = array(0.,dim=c(T2)), tauint = array(0.,dim=c(T2)), dtauint = array(0.,dim=c(T2)))

  cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")
  rm(psscar)
  for(t in 1:T2) {
    corr <- uwerrprimary(x=Z[t,skip:length(Z[1,])], S, pl=F)
    
    result$t[i] <- t-1
    result$corr[i] <- corr$value[1]
    result$dcorr[i] <- corr$dvalue[1]
    result$ddcorr[i] <- corr$ddvalue[1]
    result$tauint[i] <- corr$tauint[1]
    result$dtauint[i] <- corr$dtauint[1]
    i = i+1
  }
  rm(corr)
  rm(Z)
  attr(result, "class") <- c("correlator", "data.frame")
  if(plotit == T) {
    plot(result, log="y")
  }
  return(invisible(result))
}

correlator <-  function(filename, skip = 0, S=1.5, plotit=F, fold=F, symmetric=T, scalefactor=1.) {
  if(!missing(filename)) {
    psscar <- read.table(file=filename, col.names=c("t","ps"), header=F)
  }
  else {
    stop("Error! Option filename is mandatory!")
  }
  T2 <- (max(psscar$t)-min(psscar$t)+1)
  Z <- array(scalefactor*psscar$ps, dim=c(T2,length(psscar$ps)/T2))
  L <- T2

  i <- 1
  if(fold) {
    result <- data.frame(t = array(0.,dim=c(T2/2+1)),
                         corr = array(0.,dim=c(T2/2+1)), dcorr = array(0., dim=c(T2/2+1)),
                         ddcorr = array(0.,dim=c(T2/2+1)), tauint = array(0.,dim=c(T2/2+1)),
                         dtauint = array(0.,dim=c(T2/2+1)))
    
    cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " T2=", T2," \n")
    rm(psscar)
    for(t in 1:(T2/2+1)) {
      if((t!=1) && (t!=T2/2+1)) {
        if(symmetric) {
          corr <- uwerrprimary(x=0.5*(Z[t,skip:length(Z[1,])]+Z[T2-t+1,skip:length(Z[1,])]), S, pl=F)
        }
        else {
          corr <- uwerrprimary(x=0.5*(Z[t,skip:length(Z[1,])]-Z[T2-t+1,skip:length(Z[1,])]), S, pl=F)
        }
      }
      else {
        corr <- uwerrprimary(x=Z[t,skip:length(Z[1,])], S, pl=F)
      }
      
      result$t[i] <- t-1
      result$corr[i] <- corr$value[1]
      result$dcorr[i] <- corr$dvalue[1]
      result$ddcorr[i] <- corr$ddvalue[1]
      result$tauint[i] <- corr$tauint[1]
      result$dtauint[i] <- corr$dtauint[1]
      i = i+1
    }
  }
  else {
    result <- data.frame(t = array(0.,dim=c(T2)), corr = array(0.,dim=c(T2)), dcorr = array(0.,dim=c(T2)),
                         ddcorr = array(0.,dim=c(T2)), tauint = array(0.,dim=c(T2)), dtauint = array(0.,dim=c(T2)))
    
    cat("Found", length(psscar$ps)/T2, "measurements, skipping", skip, " \n")
    rm(psscar)
    for(t in 1:T2) {
      corr <- uwerrprimary(x=Z[t,skip:length(Z[1,])], S, pl=F)
      
      result$t[i] <- t-1
      result$corr[i] <- corr$value[1]
      result$dcorr[i] <- corr$dvalue[1]
      result$ddcorr[i] <- corr$ddvalue[1]
      result$tauint[i] <- corr$tauint[1]
      result$dtauint[i] <- corr$dtauint[1]
      i = i+1
    }
  }
  rm(corr)
  rm(Z)
  attr(result, "class") <- c("correlator", "data.frame")
  if(plotit == T) {
    plot(result)
  }
  return(invisible(result))
}


