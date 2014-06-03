# $Id$
#
# this is a more or less line by line translation of
# U. Wolffs UWerr.m for R: The R Project for Statistical Computing
#
# See
# ``Monte Carlo errors with less errors''
#   by Ulli Wolff, hep-lat/0306017
#   for details about the method
#
# R is free software and can be obtained at
# http://www.r-project.org/
#
# a typical call would be like (plaq is a vector with plaquette values)
#  source("UWerr.R")
#  plaq.res <- uwerr(data=plaq)
# or
#  nrep <- c(1234,878)
#  plaq.res <- uwerr(data=plaq, nrep=nrep, S=2.)
# then try
#  summary(plaq.res)
# and
#  source("plotutils.R")
#  plot(plaq.res)
# try help(read.table) for how to read data from a file
#
# 

uwerr <- function(f, data, nrep, S=1.5, pl=FALSE, ...) {
  # f: scalar function handle, needed for derived quantities
  # data: the matrix of data with dim. (Nalpha x N)
  #              N = total number of measurements
  #              Nalpha = number of (primary) observables
  # nrep: the vector (N1, N2, ...) of replica length N1, N2
  # S: initial guess for tau/tauint
  if(missing(f)) {
    if(missing(nrep)) {
      nrep <- c(length(data))
    }
    return(invisible(uwerrprimary(data=data, nrep=nrep, S=S, pl=pl)))
  }
  else {
    if(missing(nrep)) {
      nrep <- c(length(data[1,]))
    }
    return(invisible(uwerrderived(f=f, data=data, nrep=nrep, S=S, pl=pl, ...)))
  }
}

uwerrprimary <- function(data, nrep, S=1.5, pl=FALSE) {
  
  N = length(data)
  if(missing(nrep)) {
    nrep <- c(N)
  }
  if(any(nrep < 1) || (sum(nrep) != N)) {
    stop("Error, inconsistent N and nrep!")
  }
  R <- length(nrep) # number of replica
  mx <- mean(data)
  mxr <- rep(0., times=R) # replica mean
  i0 <- 1
  for( r in 1:R) {
    i1 <- i0-1+nrep[r]
    mxr[r] <- mean(data[i0:i1])
    i0 <- i0+nrep[r]
  }
  Fb=sum(mxr*nrep)/N # weighted mean of replica means
  delpro = data-mx;
  
  if( S == 0 ) {
    Wmax <- 0
    Wopt <- 0
    flag <- 0
  }
  else {
    Wmax <- floor(min(nrep)/2)
    Gint <- 0.
    flag <- 1
  }

  GammaFbb <- rep(0., times=Wmax)
  GammaFbb[1] = mean(delpro^2);
  if(GammaFbb[1] == 0) {
    stop("Error, no fluctuations!")
  }
  # Compute Gamma[t] and find Wopt
  
  W=1
  while(W<=Wmax) {
    GammaFbb[W+1] <- 0.
    i0 <- 1
    for(r in 1:R) {
      i1 <-  i0-1+nrep[r]
      GammaFbb[W+1] <- GammaFbb[W+1] + sum(delpro[i0:(i1-W)]*delpro[(i0+W):i1])
      i0 <- i0+nrep[r]
    }
    GammaFbb[W+1] <- GammaFbb[W+1]/(N-R*W)
    if(flag) {
      Gint <- Gint+GammaFbb[(W+1)]/GammaFbb[1]
      if(Gint<0) {
        tauW <- 5.e-16
      }
      else {
        tauW <- S/(log((Gint+1)/Gint))
      }
      gW <- exp(-W/tauW)-tauW/sqrt(W*N)
      if(gW < 0) {
        Wopt <- W
        Wmax <- min(Wmax,2*Wopt)
        flag <- 0
      }
    }
    W <- W+1
  } # while loop

  if(flag) {
    warning("Windowing condition failed!")
    Wopt <- Wmax
  }

  CFbbopt <- GammaFbb[1] + 2*sum(GammaFbb[2:(Wopt+1)])
  if(CFbbopt <= 0) {
    stop("Gamma pathological: error^2 <0")
  }
  GammaFbb <- GammaFbb+CFbbopt/N            #Correct for bias
  CFbbopt <- GammaFbb[1] + 2*sum(GammaFbb[2:(Wopt+1)]) #refined estimate
  sigmaF <- sqrt(CFbbopt/N) #error of F
  tauintFbb <- cumsum(GammaFbb)/GammaFbb[1]-0.5 # normalised autoccorelation time

  # bias cancellation for the mean value
  if(R>1) {
    bF <-  (Fb-mx)/(R-1)
    mx <- mx-bF
    if(abs(bF) > sigmaF/4) {
      warning("a %.1f sigma bias of the mean has been cancelled", bF/sigmaF)
    }
    mxr <-  mxr - bF*N/nrep
    Fb <- Fb - bF*R
  }

  value <- mx
  dvalue <- sigmaF
  ddvalue <- dvalue*sqrt((Wopt+0.5)/N)
  dtauintofW <- tauintFbb[1:(Wmax+1)]*sqrt(c(0:Wmax)/N)*2 

  tauint  <- tauintFbb[(Wopt+1)]
  dtauint <- tauint*2*sqrt((Wopt-tauint+0.5)/N)
  
  # Q value for replica distribution if R>=2
  if(R>1) {
    chisqr <- sum((mxr-Fb)^2*nrep)/CFbbopt
    Qval <- 1-pgamma(chisqr/2, (R-1)/2)
  }
  else {
    Qval <- NULL
  }
  
  res <- list(value = value, dvalue = dvalue, ddvalue = ddvalue,
              tauint = tauint, dtauint = dtauint,
              Wopt=Wopt, Wmax=Wmax, tauintofW=tauintFbb[1:(Wmax+1)],
              dtauintofW=dtauintofW[1:(Wmax+1)], Qval=Qval, S=S,
              N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb, primary=1)

  attr(res, "class") <- c("uwerr", "list")
  
  if(pl) {
    plot(res)
  }
  
  return(invisible(res))
}



uwerrderived <- function(f, data, nrep, S=1.5, pl=FALSE, ...) {

  Nalpha <- length(data[,1])
  N <- length(data[1,])
  if(missing(nrep)) {
    nrep <- c(N)
  }
  if(any(nrep < 1) || (sum(nrep) != N)) {
    stop("Error, inconsistent N and nrep!")
  }
  R=length(nrep)
  abb <- c(1:Nalpha) #total mean
  abr <- array(0., dim=c(R,Nalpha)) #replicum mean
  for (i in 1:Nalpha) {
    abb[i] <- mean(data[i,])
  }
  i0=1
  for(r in 1:R) {
    i1 <- i0-1+nrep[r]
    for (i in 1:Nalpha) {
      abr[r,i] <- mean(data[i, i0:i1])
    }
    i0 <- i0+nrep[r]
  }

  Fbb <- f(abb, ...)
  Fbr <- rep(0., times=R)
  for(r in 1:R) {
    Fbr[r] <-  f(abr[r,], ...)
  }
  Fb=sum(Fbr*nrep)/N;  # weighted mean of replica means

  fgrad=rep(0., times=Nalpha)
  h <- c(1:Nalpha)
  for (i in 1:Nalpha) {
    h[i] <- stats::sd(data[(i),])/sqrt(N)
  }

  ainc <- abb

  for (alpha in 1:Nalpha) {
    if (h[alpha] == 0) {
      # Data for this observable do not fluctuate
      fgrad[alpha]=0
    }
    else {
      ainc[alpha] <- abb[alpha]+h[alpha]
      fgrad[alpha] <- f(ainc, ...)
      ainc[alpha] <- abb[alpha]-h[alpha]
      brg <- f(ainc, ...)
      fgrad[alpha] <- fgrad[alpha]-brg
      ainc[alpha]  <- abb[alpha];
      fgrad[alpha] <- fgrad[alpha]/(2*h[alpha])
    }
  }

# projected deviations: 
#  delpro <- (t(data) %*% fgrad) - rep(abb %*% fgrad, times=N)
  delpro <- crossprod(data, fgrad) - rep(abb %*% fgrad, times=N)

  GammaFbb<-numeric()
  GammaFbb[1] <- mean(delpro^2)
  if(!any(is.na(GammaFbb))) {
    if(GammaFbb[1] == 0) {
      warning("no fluctuations!")
    }
  }
  else {
    res <- list(value = Fbb, dvalue = NA, ddvalue = NA,
                tauint = NA, dtauint = NA, Wopt=NA, Wmax=NA,
                tauintofW=NA, dtauintofW=NA,
                Qval=NA, S=S,
                N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb, primary=0)
    
    attr(res, "class") <- c("uwerr", "list")
    options(error = NULL)
    return(invisible(res))
    
  }

  if(S==0) {
    Wopt <- 0
    Wmax <- 0
    flag <- 0
  }
  else {
    Wmax <- floor(min(nrep)/2)
    Gint <- 0.
    flag <- 1
  }
  
  W<-1
  while(W <= Wmax) {
    GammaFbb[W+1] <- 0
    i0 <- 1
    for(r in 1:R) {
      i1 <-  i0-1+nrep[r]
      GammaFbb[W+1] <- GammaFbb[W+1] + sum(delpro[i0:(i1-W)]*delpro[(i0+W):i1])
      i0 <- i0+nrep[r]
    }    
    GammaFbb[W+1] = GammaFbb[W+1]/(N-R*W)
    #GammaFbb[(W+1)] <- sum(delpro[1:(N-W)]*delpro[(1+W):N])/(N-W)
    if(flag) {
      Gint <- Gint+GammaFbb[(W+1)]/GammaFbb[1]
      if(Gint<0) {
        tauW <- 5.e-16
      }
      else {
        tauW <- S/(log((Gint+1)/Gint))
      }
      gW <- exp(-W/tauW)-tauW/sqrt(W*N)
      if(gW < 0) {
        Wopt <- W
        Wmax <- min(Wmax,2*W)
        flag=0
      }
    }
    W=W+1
  }

  if(flag) {
    warning("Windowing condition failed!")
    Wopt <- Wmax
  }
  
  CFbbopt <- GammaFbb[1] + 2*sum(GammaFbb[2:(Wopt+1)])
  if(CFbbopt <= 0) {
    #stop("Gamma pathological, estimated error^2 <0\n")
    warning("Gamma pathological, estimated error^2 <0\n")
  }
  GammaFbb <- GammaFbb+CFbbopt/N # bias in Gamma corrected
  CFbbopt <- GammaFbb[1] + 2*sum(GammaFbb[2:(Wopt+1)]) #refined estimate
# Why can CFbbopt get negative?
  sigmaF <- sqrt(abs(CFbbopt)/N) #error of F
  rho <- GammaFbb/GammaFbb[1]
  tauintFbb <- cumsum(rho)-0.5 #normalised autocorrelation

  # bias cancellation for the mean value

  if(R >= 2) {
    bF <- (Fb-Fbb)/(R-1);
    Fbb <- Fbb - bF;
    if(abs(bF) > sigmaF/4) {
      warning("a %.1f sigma bias of the mean has been cancelled",bF/sigmaF)
    }
    
    Fbr = Fbr - bF*N/nrep;
    Fb  = Fb  - bF*R;
  }

  value <- Fbb
  dvalue <- sigmaF
  ddvalue <- sigmaF*sqrt((Wopt+0.5)/N)
  tauint  <- tauintFbb[(Wopt+1)]
  dtauint = tauint*2*sqrt((Wopt-tauint+0.5)/N)
  dtauintofW <- tauintFbb[1:(Wmax+1)]*sqrt(c(0:Wmax)/N)*2 
  
  # Q value for replica distribution if R>=2
  if(R>1) {
    chisqr <- sum((Fbr-Fb)^2*nrep)/CFbbopt
    Qval <- 1-pgamma(chisqr/2, (R-1)/2)
  }
  else {
    Qval <- NULL
  }
  
  res <- list(value = value, dvalue = dvalue, ddvalue = ddvalue,
              tauint = tauint, dtauint = dtauint, Wopt=Wopt, Wmax=Wmax,
              tauintofW=tauintFbb[1:(Wmax+1)], dtauintofW=dtauintofW[1:(Wmax+1)],
              Qval=Qval, S=S,
              N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb, primary=0)

  attr(res, "class") <- c("uwerr", "list")
  
  if(pl) {
    plot(res)
  }

  options(error = NULL)
  return(invisible(res))
  
}

summary.uwerr <- function(uwerr) {
  cat("Error analysis with Gamma method\n")
  cat("based on", uwerr$N, "measurements\n")
  if(uwerr$R>1) {
    cat("split in", uwerr$R, "replica with (", uwerr$nrep, ") measurements, respectively\n")
  }
  cat("The Gamma function was summed up until Wopt=", uwerr$Wopt, "\n\n")
  cat("value    =", uwerr$value, "\n")
  cat("dvalue   =", uwerr$dvalue, "\n")
  cat("ddvalue  =", uwerr$ddvalue, "\n")
  cat("tauint   =", uwerr$tauint, "\n")
  cat("dtauint  =", uwerr$dtauint, "\n")
  if(uwerr$R>1) {
    cat("Qval     =", uwerr$Qval, "\n")
  }
}

plot.uwerr <- function(uwerr, main="x") {
  if(uwerr$primary) {
    X11()
    hist(uwerr$data, main = paste("Histogram of" , main))
  }
  if(!is.null(uwerr$Gamma)) {
    GammaFbb <- uwerr$Gamma/uwerr$Gamma[1]
    Gamma.err <- gammaerror(Gamma=GammaFbb, N=uwerr$N , W=uwerr$Wmax, Lambda=100)
    X11()
    plotwitherror(c(0:uwerr$Wmax),GammaFbb[1:(uwerr$Wmax+1)],
                  Gamma.err[1:(uwerr$Wmax+1)], ylab="Gamma(t)", xlab="t", main=main)
    abline(v=uwerr$Wopt+1)
    abline(h=0)
  }
  X11()
  tauintplot(uwerr$tauintofW, uwerr$dtauintofW, uwerr$Wmax, uwerr$Wopt, main=main)  
  return(invisible(data.frame(t=c(0:uwerr$Wmax),Gamma=GammaFbb[1:(uwerr$Wmax+1)],dGamma=Gamma.err[1:(uwerr$Wmax+1)])))
}

# compute the error of the autocorrelation function using the approximate formula
# given in appendix E of hep-lat/0409106

gammaerror <- function(Gamma, N, W, Lambda) {
  gamma.err <- rep(0., times=W+1)
  Gamma[(W+2):(2*W+W+1)]=0.;
  for(t in 0:W) {
    k <- c((max(1,(t-W))):(t+W))
    gamma.err[t+1] <- sum((Gamma[(k+t+1)]+Gamma[(abs(k-t)+1)]-2*Gamma[t+1]*Gamma[(k+1)])^2);
    gamma.err[t+1] <- sqrt(gamma.err[t+1]/N)
  }
  gamma.err[1] <- 0.001
  return(invisible(gamma.err))
}

tauintplot <- function(ti, dti, Wmax, Wopt, ...) {
  plot(ti[1:Wmax], ylim=c(0.,2*ti[Wopt]), xlab="W", ylab="tauint(W)", ...)
#  plot(ti[1:Wmax], ylim=NULL)
  arrows(c(2:Wmax),ti[2:Wmax]-dti[2:Wmax],c(2:Wmax),ti[2:Wmax]+dti[2:Wmax], length=0.01,angle=90,code=3)
  abline(v=Wopt+1)
  abline(h=ti[Wopt+1],col="red")
}


