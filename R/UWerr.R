#' Time Series Analysis With Gamma Method
#' 
#' Analyse time series data with the so called gamma method
#' 
#' 
#' @aliases uwerr uwerrprimary uwerrderived
#' @param f function computing the derived quantity. If not given it is assumed
#' that a primary quantity is analysed.
#' 
#' f must have the data vector of length Nalpha as the first argument. Further
#' arguments to f can be passed to uwerr via the \code{...} argument.
#' 
#' f may return a vector object of numeric type.
#' @param data array of data to be analysed. It must be of dimension (N x
#' Nalpha) (i.e. N rows and Nalpha columns), where N is the total number of
#' measurements and Nalpha is the number of primary observables
#' @param nrep the vector (N1, N2, ...) of replica length N1, N2
#' @param S initial guess for the ratio tau/tauint, with tau the exponetial
#' autocorrelation length.
#' @param pl logical: if TRUE, the autocorrelation function, the integrated
#' autocorrelation time as function of the integration cut-off and (for primary
#' quantities) the time history of the observable are plotted with plot.uwerr
#' @param ...  arguments passed to function \code{f}.
#' @return In case of a primary observable (\code{uwerrprimary}), an object of
#' class \code{uwerr} with basis class \code{\link{list}} containing the
#' following objects \item{value}{ the expectation value of the obsevable }
#' \item{dvalue}{ the error estimate } \item{ddvalue}{ estimate of the error on
#' the error } \item{tauint}{ estimate of the integrated autocorrelation time
#' for that quantity } \item{dtauint}{ error of tauint } \item{Qval}{ the
#' p-value of the weighted average in case of several replicas } In case of a
#' derived observable (\code{uwerrderived}), i.e. if a function is specified,
#' the above objects are contained in a list called \code{res}.
#' 
#' \code{uwerrprimary} returns in addition \item{data}{ input data } whereas
#' \code{uwerrderived} returns \item{datamean}{ (vector of) mean(s) of the
#' (vector of) data } and in addition \item{fgrad}{ the estimated gradient of
#' \code{f} } and \item{f}{ the input statistics }
#' 
#' In both cases the return object containes \item{Wopt}{ value of optimal
#' cut-off for the Gamma function integration } \item{Wmax}{ maximal value of
#' the cut-off for the Gamma function integration } \item{tauintofW}{
#' integrated autocorrelation time as a function of the cut-off W }
#' \item{dtauintofW}{ error of the integrated autocorrelation time as a
#' function of the cut-off W } \item{S}{ input parameter S } \item{N}{ total
#' number of observations } \item{R}{ number of replicas } \item{nrep}{ vector
#' of observations per replicum } \item{Gamma}{ normalised autocorrelation
#' function } \item{primary}{ set to 1 for \code{uwerrprimary} and 0 for
#' \code{uwerrderived} }
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{plot.uwerr}}
#' @references ``Monte Carlo errors with less errors'', Ulli Wolff,
#'  Comput.Phys.Commun. 156 (2004) 143-153, Comput.Phys.Commun. 176 (2007) 383 (erratum),
#' hep-lat/0306017
#' @keywords optimize ts
#' @examples
#' 
#' data(plaq.sample)
#' plaq.res <- uwerrprimary(plaq.sample)
#' summary(plaq.res)
#' plot(plaq.res)
#' 
#' @export uwerr
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

#' @export
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
    flag <- FALSE
  }
  else {
    Wmax <- floor(min(nrep)/2)
    Gint <- 0.
    flag <- TRUE
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
        flag <- FALSE
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
  dGamma <- gammaerror(Gamma=GammaFbb, N=N , W=Wmax, Lambda=100)
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
              N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb, dGamma = dGamma,
              primary=1)

  attr(res, "class") <- c("uwerr", "list")
  
  if(pl) {
    plot(res)
  }
  
  return(invisible(res))
}

#' @export
uwerrderived <- function(f, data, nrep, S=1.5, pl=FALSE, ...) {
  Nalpha <- dim(data)[2]
  N <- dim(data)[1]
  if(missing(nrep)) {
    nrep <- c(N)
  }
  if(any(nrep < 1) || (sum(nrep) != N)) {
    stop("Error, inconsistent N and nrep!")
  }
  R <- length(nrep)
  abr <- array(0., dim=c(R,Nalpha)) #replicum mean
  ## total mean
  abb <- apply(data, MARGIN=2L, FUN=mean)
  i0 <- 1
  for(r in 1:R) {
    i1 <- i0-1+nrep[r]
    abr[r,] <- apply(data[i0:i1,], MARGIN=2L, FUN=mean)
    i0 <- i0+nrep[r]
  }

  Fbb <- f(abb, ...)
  ## f evaluated on single replicas
  Fbr <- apply(abr, MARGIN=1L, FUN=f, ...)
  if(is.null(dim(Fbr))) {
    dim(Fbr) <- c(1,length(Fbr))
  }
  Fb <- apply(Fbr, MARGIN=1L, FUN=function(x, nrep, N) sum(x*nrep)/N, nrep=nrep, N=N);  # weighted mean of replica means


  h <- apply(data, MARGIN=2L, FUN=sd)/sqrt(N)

  fgrad <- array(0., dim=c(Nalpha, length(Fbb)))
  ainc <- abb
  for (alpha in 1:Nalpha) {
    if (h[alpha] == 0) {
      ## Data for this observable do not fluctuate
      fgrad[alpha,]=0
    }
    else {
      ainc[alpha] <- abb[alpha]+h[alpha]
      fgrad[alpha,] <- f(ainc, ...)
      ainc[alpha] <- abb[alpha]-h[alpha]
      brg <- f(ainc, ...)
      fgrad[alpha,] <- fgrad[alpha,]-brg
      ainc[alpha]  <- abb[alpha];
      fgrad[alpha,] <- fgrad[alpha,]/(2*h[alpha])
    }
  }

  ## projected deviations: 
  delpro <- t(apply(data, MARGIN=1L, FUN=function(x, abb) x-abb, abb=abb)) %*% fgrad

  if(is.null(dim(delpro))) dim(delpro) <- c(1, length(delpro))
  GammaFbb <- as.list(apply(delpro^2, MARGIN=2L, FUN=mean))
  if(!any(is.na(GammaFbb))) {
    if(any(GammaFbb == 0)) {
      warning("no fluctuations!")
    }
  }
  else {
    res <- list(value = Fbb, dvalue = NA, ddvalue = NA,
                tauint = NA, dtauint = NA, Wopt=NA, Wmax=NA,
                tauintofW=NA, dtauintofW=NA,
                Qval=NA, S=S, fgrad=fgrad,
                N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb, primary=0)
    
    attr(res, "class") <- c("uwerr", "list")
    options(error = NULL)
    return(invisible(res))
  }
  Wopt <- as.list(rep(0, times=length(GammaFbb)))
  Wmax <- Wopt
  tauintofW <- Wopt
  dtauintofW <- Wopt
  res <- data.frame(value=rep(NA, times=length(GammaFbb)),
                    dvalue=rep(NA, times=length(GammaFbb)),
                    ddvalue=rep(NA, times=length(GammaFbb)),
                    tauint=rep(NA, times=length(GammaFbb)),
                    dtauint=rep(NA, times=length(GammaFbb)),
                    Qval=rep(NA, times=length(GammaFbb))
                    )

  for(i in c(1:length(GammaFbb))) {
    if(S==0) {
      Wopt[[i]] <- 0
      Wmax[[i]] <- 0
      flag <- FALSE
    }
    else {
      Wmax[[i]] <- floor(min(nrep)/2)
      Gint <- 0.
      flag <- TRUE
    }
    
    W <- 1
    while(W <= Wmax[[i]]) {
      GammaFbb[[i]][W+1] <- 0
      i0 <- 1
      for(r in 1:R) {
        i1 <-  i0-1+nrep[r]
        GammaFbb[[i]][W+1] <- GammaFbb[[i]][W+1] + sum(delpro[i0:(i1-W), i]*delpro[(i0+W):i1, i])
        i0 <- i0+nrep[r]
      }
      GammaFbb[[i]][W+1] <- GammaFbb[[i]][W+1]/(N-R*W)
      ##GammaFbb[(W+1)] <- sum(delpro[1:(N-W)]*delpro[(1+W):N])/(N-W)
      if(flag) {
        Gint <- Gint+GammaFbb[[i]][(W+1)]/GammaFbb[[i]][1]
        if(Gint<0) {
          tauW <- 5.e-16
        }
        else {
          tauW <- S/(log((Gint+1)/Gint))
        }
        gW <- exp(-W/tauW)-tauW/sqrt(W*N)
        if(gW < 0) {
          Wopt[[i]] <- W
          Wmax[[i]] <- min(Wmax[[i]],2*W)
          flag <- FALSE
        }
      }
      W=W+1
    }
    
    if(flag) {
      warning("Windowing condition failed!")
      Wopt[[i]] <- Wmax[[i]]
    }
    
    CFbbopt <- GammaFbb[[i]][1] + 2*sum(GammaFbb[[i]][2:(Wopt[[i]]+1)])
    if(CFbbopt <= 0) {
      ##stop("Gamma pathological, estimated error^2 <0\n")
      warning("Gamma pathological, estimated error^2 <0\n")
    }
    GammaFbb[[i]] <- GammaFbb[[i]] + CFbbopt/N # bias in Gamma corrected
    CFbbopt <- GammaFbb[[i]][1] + 2*sum(GammaFbb[[i]][2:(Wopt[[i]]+1)]) #refined estimate
    sigmaF <- sqrt(abs(CFbbopt)/N) #error of F
    rho <- GammaFbb[[i]]/GammaFbb[[i]][1]
    tauintFbb <- cumsum(rho)-0.5 #normalised autocorrelation
    
    ## bias cancellation for the mean value
    
    if(R >= 2) {
      bF <- (Fb[i]-Fbb[i])/(R-1);
      Fbb <- Fbb[i] - bF;
      if(abs(bF) > sigmaF/4) {
        warning("a %.1f sigma bias of the mean has been cancelled", bF/sigmaF)
      }
      
      Fbr[i] = Fbr[i] - bF*N/nrep;
      Fb[i]  = Fb[i]  - bF*R;
    }
    
    res$value[i] <- Fbb[i]
    res$dvalue[i] <- sigmaF
    res$ddvalue[i] <- sigmaF*sqrt((Wopt[[i]]+0.5)/N)
    res$tauint[i]  <- tauintFbb[(Wopt[[i]]+1)]
    res$dtauint[i] = res$tauint[i]*2*sqrt((Wopt[[i]]-res$tauint[i]+0.5)/N)

    tauintofW[[i]] <- tauintFbb[1:(Wmax[[i]]+1)]
    dtauintofW[[i]] <- tauintFbb[1:(Wmax[[i]]+1)]*sqrt(c(0:Wmax[[i]])/N)*2 

    ## Q value for replica distribution if R>=2
    if(R>1) {
      chisqr <- sum((Fbr[i]-Fb[i])^2*nrep)/CFbbopt
      res$Qval[i] <- 1-pgamma(chisqr/2, (R-1)/2)
    }
  }
  
  res <- list(value=res$value,
              dvalue=res$dvalue,
              ddvalue=res$ddvalue,
              tauint=res$tauint,
              dtauint=res$dtauint,
              Qval=res$Qval,
              Wopt=Wopt, Wmax=Wmax,
              tauintofW=tauintofW, dtauintofW=dtauintofW,
              S=S, fgrad=fgrad, datamean=abb,
              N=N, R=R, nrep=nrep, data=data, Gamma=GammaFbb,
              f=f, primary=0)
  
  attr(res, "class") <- c("uwerr", "list")
  
  if(pl) {
    plot(res)
  }

  options(error = NULL)
  return(invisible(res))
  
}

#' summary.uwerr
#'
#' @param object Object of type \link{uwerr}
#' @param ... Generic parameters to pass on.
#'
#' @return
#' No return value.
#' 
#' @export
summary.uwerr <- function (object, ...) {
  uwerr <- object

  cat("Error analysis with Gamma method\n")
  cat("based on", uwerr$N, "measurements\n")
  if(uwerr$R>1) {
    cat("split in", uwerr$R, "replica with (", uwerr$nrep, ") measurements, respectively\n")
  }
  if(uwerr$primary == 1) {
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
  else {
    cat("The Gamma function was summed up until Wopt=", unlist(uwerr$Wopt), "for the ", length(uwerr$Wopt), "observable(s)\n\n")
    print(uwerr$res) 
  }
}



#' Plot Command For Class UWerr
#' 
#' Plot Command For Class UWerr
#' 
#' 
#' @param x object of class \code{uwerr}
#' @param ...  generic parameters, not used here.
#' @param main main title of the plots.
#' @param plot.hist whether or not to generate a histogram
#' @param index index of the observable to plot.
#' @param Lambda Cutoff to be used in the error computation for the ACF.
#' @return produces various plots, including a histogram, the
#' autocorrelationfunction and the integrated autocorrelation time, all with
#' error bars.
#' @author Carsten Urbach, \email{carsten.urbach@@liverpool.ac.uk}
#' @seealso \code{\link{uwerr}}
#' @keywords methods hplot
#'
#' @return
#' No return value.
#' 
#' @examples
#' 
#' data(plaq.sample)
#' plaq.res <- uwerrprimary(plaq.sample)
#' plot(plaq.res)
#' 
#' @export 
plot.uwerr <- function(x, ..., main="x", plot.hist=TRUE, index=1, Lambda=100) {

  if(x$primary && plot.hist) {
    new_window_if_appropriate()
    hist(x$data, main = paste("Histogram of" , main))
  }
  if(!is.null(x$Gamma)) {
    GammaFbb <- NULL
    Gamma.err <- NULL
    Wopt <- numeric()
    Wmax <- numeric()
    if(x$primary == 1) {
      GammaFbb <- x$Gamma/x$Gamma[1]
      Gamma.err <- gammaerror(Gamma=GammaFbb, N=x$N , W=x$Wmax, Lambda=Lambda)
      Wopt <- x$Wopt
      Wmax <- x$Wmax
    }
    else {
      GammaFbb <- x$Gamma[[index]]/x$Gamma[[index]][1]
      Gamma.err <- gammaerror(Gamma=GammaFbb, N=x$N , W=x$Wmax[[index]], Lambda=Lambda)
      Wopt <- x$Wopt[[index]]
      Wmax <- x$Wmax[[index]]
    }

    new_window_if_appropriate()
    plotwitherror(x=c(0:Wmax),y=GammaFbb[1:(Wmax+1)],
                  dy=Gamma.err[1:(Wmax+1)], ylab="Gamma(t)", xlab="t", main=main)
    abline(v=Wopt+1)
    abline(h=0)
  }

  new_window_if_appropriate()
  if(x$primary == 1) tauintplot(ti=x$tauintofW, dti=x$dtauintofW, Wmax=Wmax, Wopt=Wopt, main=main)  
  else tauintplot(ti=x$tauintofW[[index]], dti=x$dtauintofW[[index]], Wmax=Wmax, Wopt=Wopt, main=main)  
  return(invisible(data.frame(t=c(0:Wmax),Gamma=GammaFbb[1:(Wmax+1)],dGamma=Gamma.err[1:(Wmax+1)])))
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


