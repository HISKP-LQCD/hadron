#' @title Compute the bootstrap error of the mean
#'
#' @description
#' Compute the bootstrap error of the mean
#'
#' @param data Original data to bootstrap
#' @param R Number of bootstrap replicates.
#' @param l Block length.
#'
#' @return
#' Returns a numeric vector with the estimated standard error
#' of the mean.
#' 
#' @export
bootstrap.meanerror <- function(data, R=400, l=20) {
  bootit <- boot::boot(block.ts(data, l=l), meanindexed, R=R)
  return(apply(bootit$t, 2, sd))
}

#' Correlator matrix model.
#' 
#' @param par Numeric vector: Fit parameters of the model. In an 
#'   object of type \code{matrixfit}, this should be located at 
#'   \code{$opt.res$par}.
#' @param t integer vector: Time of interest.
#' @param Time integer: Time extent of the lattice.
#' @param parind See \code{\link{matrixfit}}.
#' @param sign.vec Numeric vector: Relative sign between forward and
#'   backwards propagating part. A plus makes it cosh, a minus makes it sinh.
#' @param ov.sign.vec Numeric vector: Overal sign.
#' @param deltat Numeric: time shift.
#'
#' @return
#' Returns a numeric vector with the same length as the input vector `t`
#' containing the model evaluation for these t-values.
#' 
#' @seealso \code{\link{matrixfit}}
matrixModel <- function(par, t, Time, parind, sign.vec, ov.sign.vec, deltat=0) {
  return(ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*
         (exp(-par[1]*(t-deltat/2)) + sign.vec*exp(-par[1]*(Time-(t-deltat/2))))
         )
}

#' Principal correlator two state model.
#' 
#' @param par Numeric vector: Fit parameters of the model. In an 
#'   object of type \code{matrixfit}, this should be located at 
#'   \code{$opt.res$par}.
#' @param t Numeric vector: Time of interest.
#' @param Time Numeric: Time extent of the lattice.
#' @param reference_time Numeric: GEVP reference time value in physical time convention
#' @param delta1 dummy parameter for compatibility
#' 
#' @return
#' Returns a numeric vector with the same length as the input vector `t`
#' containing the model evaluation for these t-values.
#' 
#' @seealso \code{\link{matrixfit}}
pcModel <- function(par, t, Time, delta1=1, reference_time) {
  return( exp(-abs(par[1])*(t-reference_time))*( par[3] + (1-par[3])*exp(-(abs(par[2]))*(t-reference_time)) ) )
}

matrixChisqr <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(Time-t))))
  return( sum(z %*% M %*% z) )
}

matrixChi <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(Time-t))))
  return( L %*% z )
}

## deltat and reference_time are dummy variable here
dmatrixChi <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  zp <- -ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(-t*exp(-par[1]*t) -(Time-t)*sign.vec*exp(-par[1]*(Time-t)))
  res <- L %*% zp
  for(i in 2:length(par)) {
    zp1 <- rep(0, length(zp))
    j <- which(parind[,1]==i)
    zp1[j] <- -ov.sign.vec*0.5*par[parind[j,2]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(Time-t[j])))
    zp2 <- rep(0, length(zp))
    j <- which(parind[,2]==i)
    zp2[j] <- -ov.sign.vec*0.5*par[parind[j,1]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(Time-t[j])))
    res <- c(res, L %*% zp1 + L %*% zp2)
  }
  return(res)
}

## deltat and reference_time are dummy variable here
dmatrixChisqr <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  res <- rep(0., times=length(par))
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(Time-t))))
  zp <- -ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(-t*exp(-par[1]*t) -(Time-t)*sign.vec*exp(-par[1]*(Time-t)))
  res[1] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  for(i in 2:length(par)) {
    zp <- rep(0, length(z))
    j <- which(parind[,1]==i)
    zp[j] <- -ov.sign.vec*0.5*par[parind[j,2]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(Time-t[j])))
    res[i] <- sum(zp %*% M %*% z + z %*% M %*% zp)
    zp <- rep(0, length(z))
    j <- which(parind[,2]==i)
    zp[j] <- -ov.sign.vec*0.5*par[parind[j,1]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(Time-t[j])))
    res[i] <- res[i] + sum(zp %*% M %*% z + z %*% M %*% zp)
  }
  return(res)
}

## A Chi and Chisqr function for a two-state principal correlator specific model
## for a single correlator only
## The model reads
## A*exp(-E(t-reference_time)) + (1-A)*exp(-E2(t-reference_time))
## = exp(-E(t-reference_time))*(A + (1-A)*exp(-DeltaE(t-reference_time)))
##
## the respective gradient functions follow
##
## deltat is a dummy variable here
pcChi <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time) {
  return( (y - exp(-abs(par[1])*(t-reference_time))*( par[3] + (1-par[3])*exp(-(abs(par[2]))*(t-reference_time)) )) %*% L )
}

## deltat is a dummy variable here
pcChisqr <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  z <- (y - exp(-abs(par[1])*(t-reference_time))*(par[3]+(1-par[3])*exp(-(abs(par[2]))*(t-reference_time))))
  return( sum(z %*% M %*% z) )
}

## deltat is a dummy variable here
dpcChi <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time) {
  zp <- (t-reference_time)*exp(-abs(par[1])*(t-reference_time))*(par[3]+(1-par[3])*exp(-(abs(par[2]))*(t-reference_time)))
  res <- L %*% zp
  zp <- exp(-abs(par[1])*(t-reference_time))*(1-par[3])*(t-reference_time)*exp(-(abs(par[2]))*(t-reference_time))
  res <- c(res, L %*% zp)
  zp <- -exp(-abs(par[1])*(t-reference_time))*(1-exp(-(abs(par[2]))*(t-reference_time)))
  res <- c(res, L %*% zp)
  return(res)
}

## deltat is a dummy variable here
dpcChisqr <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time) {
  res <- rep(0., times=length(par))
  z <- (y - exp(-abs(par[1])*(t-reference_time))*(par[3]+(1-par[3])*exp(-(abs(par[2]))*(t-reference_time))))
  zp <- (t-reference_time)*exp(-abs(par[1])*(t-reference_time))*(par[3]+(1-par[3])*exp(-(abs(par[2]))*(t-reference_time)))
  res[1] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  zp <- exp(-abs(par[1])*(t-reference_time))*(1-par[3])*(t-reference_time)*exp(-(abs(par[2]))*(t-reference_time))
  res[2] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  zp <- -exp(-abs(par[1])*(t-reference_time))*(1-exp(-(abs(par[2]))*(t-reference_time)))
  res[3] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  return(res)
}

## reference_time is a dummy variable here
matrixChisqr.shifted <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(Time-(t-deltat/2)))))
  return( sum(z %*% M %*% z ) )
}

## reference_time is a dummy variable here
dmatrixChisqr.shifted <- function(par, t, y, M, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  res <- rep(0., times=length(par))
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(Time-(t-deltat/2)))))
  zp <- -ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) + (Time-t+deltat/2)*sign.vec*exp(-par[1]*(Time-(t-deltat/2))))
  res[1] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  for(i in 2:length(par)) {
    zp <- rep(0, length(z))
    j <- which(parind[,1]==i)
    zp[j] <- -ov.sign.vec[j]*par[parind[j,2]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(Time-(t[j]-deltat/2))))
    res[i] <- sum(zp %*% M %*% z + z %*% M %*% zp)
    zp <- rep(0, length(z))
    j <- which(parind[,2]==i)
    zp[j] <- -ov.sign.vec[j]*par[parind[j,1]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(Time-(t[j]-deltat/2))))
    res[i] <- res[i] + sum(zp %*% M %*% z + z %*% M %*% zp)
  }
  return(res)
}

## reference_time is a dummy variable here
matrixChi.shifted <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(Time-(t-deltat/2)))))
  return( L %*% z )
}

## reference_time is a dummy variable here
dmatrixChi.shifted <- function(par, t, y, L, Time, parind, sign.vec, ov.sign.vec, deltat=1, reference_time=0) {
  zp <- -ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) +(Time-t+deltat/2)*sign.vec*exp(-par[1]*(Time-(t-deltat/2))))
  res <- L %*% zp
  for(i in 2:length(par)) {
    zp1 <- c(0)
    j <- which(parind[,1]==i)
    zp1[j] <- -ov.sign.vec[j]*par[parind[j,2]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(Time-(t[j]-deltat/2))))
    zp2 <- c(0)
    j <- which(parind[,2]==i)
    zp2[j] <- -ov.sign.vec[j]*par[parind[j,1]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(Time-(t[j]-deltat/2))))
    res <- c(res, L %*% zp1 + L %*% zp2)
  }
  return(res)
}


# the calling code must supply the correct three parameters 
deriv.CExp <- function(par, t, Time, sign) {
  res <- array(0.,dim=c(length(par),length(t)))

  res[1,] <- 0.5*par[2]*par[3]*(-t*exp(-par[1]*t) -(Time-t)*sign*exp(-par[1]*(Time-t)))
  res[2,] <- 0.5*par[3]*(exp(-par[1]*t) + sign*exp(-par[1]*(Time-t)))
  res[3,] <- 0.5*par[2]*(exp(-par[1]*t) + sign*exp(-par[1]*(Time-t)))

  return(res)
}

deriv.CExp.shifted <- function(par, t, Time, sign, deltat=1) {
  res <- array(0.,dim=c(length(par),length(t)))
  
  res[1,] <- par[2]*par[3]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) + (Time-(t-deltat/2))*sign*exp(-par[1]*(Time-(t-deltat/2))))
  res[2,] <- par[3]*(exp(-par[1]*(t-deltat/2)) - sign*exp(-par[1]*(Time-(t-deltat/2))))
  res[3,] <- par[2]*(exp(-par[1]*(t-deltat/2)) - sign*exp(-par[1]*(Time-(t-deltat/2))))

  return(res)
}

deriv.pcModel <- function(par, t, Time, reference_time) {
  res <- array(0.,dim=c(length(par),length(t)))
  res[1, ] <- -(t-reference_time)*exp(-par[1]*(t-reference_time))*(par[3]+(1-par[3])*exp(-par[2]*(t-reference_time)))
  res[2, ] <- -exp(-par[1]*(t-reference_time))*(1-par[3])*(t-reference_time)*exp(-((par[2]))*(t-reference_time))
  res[3, ] <- exp(-par[1]*(t-reference_time))*(1-exp(-((par[2]))*(t-reference_time)))
  return(res)

}



#' Routine For A Factorising Matrix Fit
#' 
#' Performs a factorising fit on a correlation matrix
#' 
#' The routine expects in \code{cf$cf} a set of correlation functions.  The
#' mapping of this linear construct to a matrix or a part of a matrix is
#' achieved via \code{parlist}. The symmetry properties of the individual
#' correlation functions must be encoded in \code{sym.vec}.
#' 
#' \code{matrixfit} will fit to every correlator in \code{cf$cf} a function
#' \eqn{p_i p_j f(t)}. The indices \eqn{i,j} are determined from \code{parlist}
#' and \eqn{f} is either \eqn{cosh}{\cosh} or \eqn{sinh}{\sinh}, depending on
#' \code{sym.vec}.
#' 
#' The inverse covariance matrix is computed using a singular value
#' decomposition. If the sample size N is too small, only sqrt(N) eigenvalues
#' of the matrix are kept exactly, while all others are replaced by the mean of
#' the rest. This helps to reduce instabilities induced by too small
#' eigenvalues of the covariance matrix.
#' 
#' @param cf correlation matrix obtained with a call to \code{extrac.obs}.
#' @param t1 lower bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param t2 upper bound for the fitrange in time (t1,t2). Counting starts with
#' 0.
#' @param parlist a two dimensional array of dimension 2 times number of
#' correlators in cf. Every column assigns a pair of fit parameters to the
#' corresponding correlator in cf. In case this is missing there are defaults
#' provided for certain matrix sizes.
#' @param sym.vec a vector of length number of correlators in cf indicating
#' whether the correlation function is a cosh, a sinh or an exponential.
#' Possible values are \code{"cosh"}, \code{"sinh"} and \code{"exp"}.  In case
#' this is missing there are defaults provided for certain matrix sizes.
#' @param neg.vec a vector of length number of correlators in cf indicating
#' whether the correlation function is to be multiplied globally with a minus
#' sign.  In case this is missing there are defaults provided for certain
#' matrix sizes.
#' @param useCov use correlated or uncorrelated chisquare. Default is
#' \code{useCov=FALSE}.
#' @param boot.fit If set to \code{FALSE}, the fit is not bootstrapped, even if
#' the bootstrapping parameters have been set and the correlation function has
#' been bootstrapped.  This is a useful time-saver if error information is not
#' strictly necessary.  Of course, this affects the return values related to
#' the bootstrap, which are set to \code{NA}.
#' @param fit.method Can be either \code{"optim"} or \code{"lm"}. The latter
#' works only if the library \code{"minpack.lm"} can be loaded. Default and
#' fallback is \code{"optim"}.
#' @param model Sets the fit model to be used in the fit. The default model
#' is\cr \eqn{0.5 p_i p_j (\exp(-Et) \pm c* \exp(-E(Time-t)))}\cr with sign
#' depending on \code{"cosh"} or \code{"sinh"}. c equals one except for the
#' \code{"exp"} functional dependence. When model is set to \code{"shifted"},
#' the fit uses the function\cr \eqn{p_i p_j (\exp(-E(t+1/2)) \mp c*
#' \exp(-E(Time-(t+1/2))))}\cr which is useful when the original correlation
#' function or matrix is shifted, see e.g. \link{bootstrap.gevp}.\cr In case
#' only a single principal correlator from a GEVP is to be fitted the
#' additional model \code{"pc"} is available. It implements\cr
#' \eqn{\exp(-E(t-t_0))(A + (1-A)\exp(-DeltaE(t-t_0))}\cr with \eqn{t_0} the
#' reference timesclice of the GEVP. See \link{bootstrap.gevp} for details.
#' @param autoproceed When the inversion of the variance-covariance matrix
#' fails, the default behaviour is to abort the fit. Setting this to
#' \code{TRUE} means that the fit is instead continued with a diagonal inverse
#' of the variance-covariance matrix.
#' @param every Fit only a part of the data points. Indices that are not
#' multiples of \code{every} are skipped. If no value is provided, all points
#' are taken into account.
#' @return returns an object of class \code{matrixfit} with entries: \item{CF}{
#' object of class cf which contains the mean correlation functions} \item{M}{
#' inverse variance-covariance matrix for weighted Chi squared minimization}
#' \item{L}{squre root of \code{M}.} \item{parind}{indices in the parameter
#' vector used for the different matrix combinations} \item{sign.vec}{vector
#' of signs} \item{ii}{vector of vector indices giving the columns of the
#' correlation function arrays (CF above, say), which are contained in the fit
#' range} \item{opt.res}{return value of the minimization (see ?optim) on the
#' original data.} \item{t0}{Result of the chisqr fit on the original data.
#' \code{t0} is a vector of length npar+1, where \code{npar} the number of fit
#' parameters. The last value is the chisqr value.} \item{t}{Bootstrap
#' samples of the \code{R} Chi squared minimizations of length(par)+1. \code{t}
#' has dimension \eqn{R x (npar+1)}, where \code{R} is the number of bootstrap
#' samples and \code{npar} the number of fit parameters. The last column
#' corresponds to the chisquare values.} \item{se}{Bootstrap estimate of
#' standard error for all parameters. \code{se} is a vector of length
#' \code{npar}, where \code{npar} the number of fit parameters.}
#' \item{useCov}{whether covariances in the data were taken into account}
#' \item{invCovMatrix}{inverse of covariance matrix or inverse variance
#' weighted if useCov=FALSE} \item{Qval}{real number between 0 and 1 giving
#' the "quality" of the fit} \item{chisqr}{total Chi squared of the fit}
#' \item{dof}{fit degrees of freedom} \item{mSize}{integer size of the
#' matrix which was fitted} \item{cf}{object of type cf which contains,
#' amongst other objects, cf$cf which is a concatenated array of raw
#' correlation functions where each row is one of N observations and there are
#' mSize*Time columns (see ?extract.obs)} \item{boot.R}{number of bootstrap
#' samples} \item{boot.l}{block size for blocked bootstrap} \item{t1}{
#' beginning of fit range} \item{t2}{end of fit range} \item{parlist}{array
#' of parameter combinations for the matrix fit} \item{sym.vec}{vector of
#' strings indicating the functional form of correlation functions which were
#' fitted} \item{seed}{RNG seed for bootstrap procedure} \item{model}{see
#' input.} \item{fit.method}{see input.} \item{reference_time}{The GEVP
#' reference time for the principal correlator model}
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{cf}}, \code{\link{bootstrap.cf}}
#' @references C. Michael, `hep-lat/9412087hep-lat/9412087`
#' @keywords optimize ts
#' @examples
#' 
#' data(samplecf)
#' samplecf <- bootstrap.cf(cf=samplecf, boot.R=99, boot.l=2, seed=1442556)
#' fitres <- matrixfit(cf=samplecf, t1=16, t2=24, useCov=FALSE,
#'                     parlist=array(c(1,1), dim=c(2,1)),
#'                     sym.vec=c("cosh"), fit.method="lm")
#' summary(fitres)
#' plot(fitres)
#' 
#' @export matrixfit
matrixfit <- function(cf, t1, t2,
                      parlist,
                      sym.vec,
                      neg.vec,
                      useCov=FALSE,
                      model="single",
                      boot.fit=TRUE,
                      fit.method="optim",
                      autoproceed=FALSE,
                      every) {

  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  pcmodel <- FALSE
  if(model == "pc") pcmodel <- TRUE
  if(pcmodel) {
    stopifnot(inherits(cf, 'cf_principal_correlator'))
  }
  if(cf$symmetrised == FALSE){
    stop('You must symmetrize and bootstrap the function before fitting.')
  }

  t1p1 <- t1 + 1
  t2p1 <- t2 + 1

  reference_time <- numeric()
  if(!inherits(cf, 'cf_principal_correlator')) {
    reference_time <- 1
  } else {
    reference_time <- cf$gevp_reference_time
  }
  
  N <- dim(cf$cf)[1]
  Thalfp1 <- cf$Time/2+1
  t <- c(0:(cf$Time/2))
  deltat <- 1
  if(model == "shifted" && any(names(cf) == "deltat")) {
    deltat <- cf$deltat
  }
  
  ## This is the number of correlators in cf
  if(!is.null(dim(cf$cf)))
    mSize <- dim(cf$cf)[2]/Thalfp1
  else
    mSize <- dim(cf$cf.tsboot$t)[2]/Thalfp1
  
  if(pcmodel && mSize != 1) {
    stop('for model pc only a 1x1 matrix is allowd.')
  }

  if(missing(parlist)) {
    if(mSize == 1) {
      parlist <- array(c(1,1), dim=c(2,1))
      warning("missing parlist, using default for single correlator!\n")
    }
    else if(mSize == 4) {
      parlist <- array(c(1,1,1,2,2,1,2,2), dim=c(2,4))
      warning("missing parlist, using default for four correlators!\n")
    }
    else {
      stop("parlist is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  if(missing(sym.vec)) {
    if(mSize == 1) {
      sym.vec <- c("cosh")
      warning("missing sym.vec, using default for single correlator!\n")
    }
    else if(mSize == 4) {
      sym.vec <- c("cosh","cosh","cosh","cosh")
      warning("missing sym.vec, using default for four correlators!\n")
    }
    else {
      stop("sym.vec is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  if(missing(neg.vec)){
    if(mSize == 1) {
      neg.vec <- c(1)
      warning("missing neg.vec, using default (correlator positive)!\n")
    } 
    else if( mSize == 4 ){
      neg.vec <- c(1,1,1,1)
      warning("missing neg.vec, using default (all correlators positive)!\n")
    }
    else {
      stop("neg.vec is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  
  ## some sanity checks
  if(min(parlist) <= 0) {
    stop("Elements of parlist must be all > 0! Aborting\n")
  }
  for(i in 1:max(parlist)) {
    if(!any(parlist==i)) {
      stop("not all parameters are used in the fit! Aborting\n")
    }
  }
  
  if(dim(parlist)[2] != mSize) {
    warning(mSize, " ", dim(parlist)[2], "\n")
    stop("parlist has not the correct length! Aborting! Use e.g. extractSingleCor.cf or c to bring cf to correct number of observables\n")
  }
  if(length(sym.vec) != mSize) {
    stop("sym.vec does not have the correct length! Aborting\n")
  }
  if(length(neg.vec) != mSize){
    stop("neg.vec does not have the correct length! Aborting\n")
  }

  CF <- data.frame(t=t, Cor=cf$cf0, Err=apply(cf$cf.tsboot$t, 2, cf$error_fn))

  ## index vector for timeslices to be fitted
  ii <- c((t1p1):(t2p1))
  if(mSize > 1) {
    for(j in 2:mSize) {
      ii <- c(ii, (t1p1+(j-1)*Thalfp1):(t2p1+(j-1)*Thalfp1))
    }
  }
  ## for the pc model we have to remove timeslice reference_time, where the error is zero
  if(pcmodel) {
    ii <- ii[which(ii != (reference_time+1))]
  }
  ## use only a part of the time slices for better conditioned cov-matrix
  if(!missing(every)){
    ii <- ii[which(ii%%every == 0)]
  }
  
  ## parind is the index vector for the matrix elements
  ## signvec decides on cosh or sinh
  ## ov.sign.vec indicates the overall sign 
  parind <- array(1, dim=c(length(CF$Cor),2))
  sign.vec <- rep(+1, times=length(CF$Cor))
  ov.sign.vec <- rep(+1, times=length(CF$Cor))
  for(i in 1:mSize) {
    parind[((i-1)*Thalfp1+1):(i*Thalfp1),] <- t(array(parlist[,i]+1, dim=c(2,Thalfp1)))
    if(sym.vec[i] == "sinh") sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
    if(sym.vec[i] == "exp") sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- 0
    if(neg.vec[i] == -1) ov.sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
  }
  
  CovMatrix <- NULL
  # we always use the boostrap samples to estimate the covariance matrix 
  CovMatrix <- cf$cov_fn(cf$cf.tsboot$t[,ii])
  
  ## for uncorrelated chi^2 use diagonal matrix with inverse sd^2
  M <- diag(1/CF$Err[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ## see C. Michael hep-lat/9412087
    M <- try(invertCovMatrix(cf$cf.tsboot$t[,ii], boot.l=cf$boot.l, boot.samples=TRUE, cov_fn=cf$cov_fn), silent=TRUE)
    if(inherits(M, "try-error")) {
      if( autoproceed ){
        M <- diag(1/CF$Err[ii]^2)
        warning("[matrixfit] inversion of variance covariance matrix failed, continuing with uncorrelated chi^2\n")
        useCov <- FALSE
      } else {
        stop("[matrixfit] inversion of variance covariance matrix failed!\n")
      }
    }
  }
  lm.avail <- FALSE
  if(fit.method == "lm") 
    lm.avail <- requireNamespace('minpack.lm')
  LM <- chol(M)
  
  par <- numeric(max(parind))
  ## we get initial guesses for fit parameters from effective masses
  ## first is the mass
  ## (we currently allow for only one)
  if(pcmodel) {
    ## the ground state energy
    par[1] <- log(CF$Cor[reference_time+1]/CF$Cor[reference_time+2])
    ## the deltaE
    par[2] <- log(CF$Cor[reference_time+1]/CF$Cor[reference_time+2]) - par[1]
    par[2] <- 1.
    ## the amplitude
    par[3] <- 1.
  }
  else {
    j <- which(parlist[1,]==1 & parlist[2,]==1)
    par[1] <- invcosh(CF$Cor[t1p1+(j-1)*Thalfp1]/CF$Cor[t1p1+(j-1)*Thalfp1+1], t=t1p1, cf$Time)
    ## catch failure of invcosh
    if(is.na(par[1]) || is.nan(par[1])) par[1] <- 0.2
    ## the amplitudes we estimate from diagonal elements
    for(i in 2:length(par)) {
      j <- which(parlist[1,]==(i-1) & parlist[2,]==(i-1))
      if(length(j) == 0) {
        ##if(full.matrix) warning("one diagonal element does not appear in parlist\n")
        j <- i-1
      }
      par[i] <- sqrt(abs(CF$Cor[t1p1+(j-1)*Thalfp1])/0.5/exp(-par[1]*t1))
    }
  }

  fitfn <- matrixChisqr
  dfitfn <- dmatrixChisqr
  if(model == "shifted") {
    fitfn <- matrixChisqr.shifted
    dfitfn <- dmatrixChisqr.shifted
  }
  if(model == "weighted") {
    if(any(names(cf) == "weighted")) {
      stop('Weighted model is not implemented.')
    }
  }
  if(pcmodel) {
    fitfn <- pcChisqr
    dfitfn <- dpcChisqr
    dfitfn <- NULL
  }
  if(lm.avail) {
    fitfn <- matrixChi
    dfitfn <- dmatrixChi
    if(model == "shifted") {
      fitfn <- matrixChi.shifted
      dfitfn <- dmatrixChi.shifted
    }
    if(pcmodel) {
      fitfn <- pcChi
      dfitfn <- dpcChi
      dfitfn <- NULL
    }
  }
  
  ## check out constrOptim
  ## now perform minimisation
  dof <- (length(CF$t[ii])-length(par))
  opt.res <- NA
  rchisqr <- 0.
  if(lm.avail) {
    opt.res <- minpack.lm::nls.lm(par = par, fn = fitfn, jac=dfitfn, t=CF$t[ii], y=CF$Cor[ii], L=LM, Time=cf$Time, deltat=deltat,
                      parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii], reference_time=reference_time,
                      control = minpack.lm::nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500, maxfev=5000))
    if( !(opt.res$info %in% c(1,2,3) ) ){
      warning(sprintf("Termination reason of nls.lm opt.res$info: %d\n", opt.res$info))
    }
    rchisqr <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fitfn, gr = dfitfn,
                     method="BFGS", control=list(maxit=500, parscale=par, ndeps=rep(1.e-8, times=length(par)), REPORT=50),
                     t=CF$t[ii], y=CF$Cor[ii], M=M, Time=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], reference_time=reference_time,
                     ov.sign.vec=ov.sign.vec[ii], deltat=deltat)
    rchisqr <- opt.res$value
  }
  Qval <- 1-pchisq(rchisqr, dof)
  ## we use absolute values in the fit model
  ## better remove any signs then here!
  if(pcmodel) {
    opt.res$par[1:2] <- abs(opt.res$par[1:2])
  }
  
  opt.tsboot <- NA
  if(boot.fit) {
    opt.tsboot <- apply(X=cf$cf.tsboot$t[,ii], MARGIN=1, FUN=fit.formatrixboot, par=opt.res$par, t=CF$t[ii], deltat=deltat,
                        M=M, Time=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii],
                        L=LM, lm.avail=lm.avail, fitfn=fitfn, dfitfn=dfitfn, reference_time=reference_time)
  }
  N <- length(cf$cf[,1])
  if(is.null(cf$cf)) {
    N <- cf$N
  }
  if(pcmodel) {
    opt.tsboot[c(1:2),] <- abs(opt.tsboot[c(1:2),])
  }
  res <- list(CF=CF, M=M, L=LM, parind=parind, sign.vec=sign.vec, ov.sign.vec=ov.sign.vec, ii=ii, opt.res=opt.res, opt.tsboot=opt.tsboot,
              boot.R=cf$boot.R, boot.l=cf$boot.l, useCov=useCov, CovMatrix=CovMatrix, invCovMatrix=M, seed=cf$seed,
              Qval=Qval, chisqr=rchisqr, dof=dof, mSize=mSize, cf=cf, t1=t1, t2=t2, reference_time=reference_time,
              parlist=parlist, sym.vec=sym.vec, N=N, model=model, fit.method=fit.method, error_fn=cf$error_fn, niter = c(opt.res$niter, opt.tsboot[nrow(opt.tsboot), ]))
  res$t <- t(opt.tsboot)
  res$t0 <- c(opt.res$par, opt.res$value)
  res$se <- apply(opt.tsboot[c(1:(dim(opt.tsboot)[1]-1)),], MARGIN=1, FUN=cf$error_fn)
  attr(res, "class") <- c("matrixfit", "list")
  return(invisible(res))
}

#' Plot a matrixfit
#' 
#' @param x an object of class matrixfit
#' @param plot.errorband Boolean: whether or not to plot an errorband
#' @param ylim Numeric vector: y-limit of the plot
#' @param xlab String: label of x-axis
#' @param ylab String: label of y-axis
#' @param do.qqplot Boolean: whether or not to plot an QQ-plot
#' @param plot.raw Boolean: plot the raw data or multiply out the leading exponetial behaviour
#' @param rep Boolean: whether or not to add to existing plot
#' @param col String vector: vector of colours for the different correlation functions
#' @param every Fit only a part of the data points. Indices that are not multiples of \code{every}
#'    are skipped. If no value is provided, all points are taken into account.
#' @param ... Graphical parameters to be passed on to \link{plot} or \link{plotwitherror}.
#' 
#' @seealso \code{\link{matrixfit}}
#'
#' @return
#' Returns no value, generated only plots.
#' 
#' @export
plot.matrixfit <- function (x, plot.errorband = FALSE, ylim, xlab = "t/a", ylab = "y",
                            do.qqplot = TRUE, plot.raw = TRUE, rep = FALSE, col, every, ...) {
  mfit <- x
  par <- mfit$opt.res$par
  parind <-  mfit$parind
  sign.vec <- mfit$sign.vec
  ov.sign.vec <- mfit$ov.sign.vec
  Time <- mfit$cf$Time
  Thalfp1 <- Time/2+1
  deltat <- 1
  if(mfit$model == "shifted" && any(names(mfit$cf) == "deltat")) {
    deltat <- mfit$cf$deltat
  }
  if(missing(col)){
    col <- c("black",rainbow(n=(mfit$mSize-1)))
  }

  if(missing(ylim)) {
    if(plot.raw) {
      ## prevent stray negative values from ruining the plot
      lbound <- ov.sign.vec*mfit$CF$Cor - 2*mfit$CF$Err
      lbound <- lbound[ lbound > 0 ]
      ylims <- c( min( lbound, na.rm=TRUE ) , max( ov.sign.vec*mfit$CF$Cor + 2*mfit$CF$Err, na.rm=TRUE ) )
    }
    else ylims <- c(0,3)
  }
  else ylims <- ylim
  
  if(!rep) {
    ## generate an empty plot
    if(plot.raw) plot(NA, log="y", ylim=ylims, xlim=range(mfit$CF$t), xlab=xlab, ylab=ylab, ...)
    else plot(NA, ylim=ylims, xlim=range(mfit$CF$t), xlab=xlab, ylab=ylab, ...)
  }
  tx <- seq(mfit$t1, mfit$t2, 0.05)

  for(i in 1:mfit$mSize ) {
    ii <- c(((i-1)*Thalfp1+1):(i*Thalfp1))
    if(!missing(every)){
      ii <- ii[which(ii%%every == 0)]
    }
    tt <- mfit$CF$t[ii]
    
    par.ind <- c(1,parind[(i-1)*Thalfp1+1,1],parind[(i-1)*Thalfp1+1,2])
    if(mfit$model == "pc") par.ind=c(1,2,3)
    pars <- c(par[1],par[par.ind[2]],par[par.ind[3]])
    sgn <- sign.vec[(i-1)*Thalfp1+1]
    
    if(mfit$model == "shifted") y <- pars[2]*pars[3]*( exp(-pars[1]*(tx-deltat/2)) - sgn*exp(-pars[1]*(Time-(tx-deltat/2))))
    else if(mfit$model == "pc") y <- pcModel(par=par[1:3], t=tx, Time=Time, reference_time=mfit$reference_time)
    else y <- 0.5*pars[2]*pars[3]*( exp(-pars[1]*tx) + sgn*exp(-pars[1]*(Time-tx)))
    
    # yp is the physical exponential in case we want to look at the ratio plot
    if(!plot.raw) {
      if(mfit$model == "shifted") {
        yp <- pars[2]*pars[3]*( exp(-pars[1]*(tt-deltat/2)) - sgn*exp(-pars[1]*(Time-(tt-deltat/2))))
      }
      else if(mfit$model == "pc") {
        yp <- exp(-par[1]*(tt - mfit$reference_time))*par[3]
      }
      else {
        yp <- 0.5*pars[2]*pars[3]*( exp(-pars[1]*tt) + sgn*exp(-pars[1]*(Time-tt)))
      }
    }
    else {
      yp <- rep(1, times=length(tt))
    }

    lwd <- c(1.5)
    if(plot.errorband) {

      dummyfn <- function(par, tx, Time, sgn, deltat, reference_time) {
        0.5*pars[2]*pars[3]*( exp(-pars[1]*tx) + sgn*exp(-pars[1]*(Time-tx)))
      }
      if(mfit$model == "shifted") {
        dummyfn <- function(par, tx, Time, sgn, deltat, reference_time) {
          pars[2]*pars[3]*( exp(-pars[1]*(tx - deltat/2)) - sgn*exp(-pars[1]*(Time-(tx-deltat/2))))
        }
      }
      else if(mfit$model == "pc") {
        dummyfn <- function(par, tx, Time, sgn, deltat, reference_time) {
          pcModel(par=par[1:3], t=tx, Time=Time, reference_time=reference_time)
        }        
      }

      se <- apply(X=apply(X=mfit$t[, par.ind], MARGIN=1, FUN=dummyfn, tx=tx, Time=Time, deltat=deltat, sgn=sgn, reference_time=mfit$reference_time), FUN=mfit$cf$error_fn, MARGIN=1)
      
      polyval <- c( (y + se), rev(y - se) )
      ## any of those not on the plot? replace to avoid wrongly drawn band!
      if(any(polyval < ylims[1]) || any(polyval > ylims[2])) {
        polyval[polyval < ylims[1]] <- ylims[1]
        polyval[polyval > ylims[2]] <- ylims[2]
      }
      if(!plot.raw) {
        if(mfit$model == "pc") {
          tmp <- exp(-par[1]*(tx - mfit$reference_time))*par[3]
          polyval <- polyval/c(tmp, rev(tmp))
        }
        else polyval <- polyval/c(y, rev(y))
      }
      polyx <- c(tx,rev(tx))
      polycol <- col2rgb(col[i],alpha=TRUE)/255
      polycol[4] <- 0.65
      ##polygon(x=polyx,y=polyval,col=rgb(red=polycol[1],green=polycol[2],blue=polycol[3],alpha=polycol[4]),border=NA)
      polygon(x=polyx,y=polyval,col="gray",border=NA)
      lwd <- c(1)
    }

    if(plot.raw) lines(tx, y, col=col[i], lwd=lwd)
    else if(mfit$model == "pc") lines(tx, y/(exp(-par[1]*(tx - mfit$reference_time))*par[3]), col=col[i], lwd=lwd)
    else abline(h=1, lwd=lwd, lty=2)

    ## plot data last to have them on top of lines and bands
    plotwitherror(x=mfit$CF$t[ii], y=ov.sign.vec[ii]*mfit$CF$Cor[ii]/yp, dy=mfit$CF$Err[ii]/yp,
                  rep=TRUE, col=col[i])
  }
  if(do.qqplot){
    new_window_if_appropriate()
    s <- seq(0,1,1./nrow(mfit$t))
    x <- qchisq(p=s, df=mfit$dof, ncp=mfit$chisq)
    qqplot(x=x, y=mfit$t[, ncol(mfit$t)-1], xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ-Plot non-central Chi^2 Values")
  }
}

#' summary.matrixfit
#'
#' @param object Object of type \link{matrixfit}
#' @param ... Generic parameters to pass on.
#'
#' @return
#' No return value.
#' 
#' @export
summary.matrixfit <- function (object, ...) {
  mfit <- object
  if(mfit$model == "pc") {
    cat("\n ** Result of two state exponential fit to principal correlator **\n\n")
  }
  else {
    cat("\n ** Result of one state exponential fit **\n\n")
  }
  cat("based on", mfit$N, "measurements\n")
  cat("time range from", mfit$t1, " to ", mfit$t2, "\n")
  if(mfit$model == "pc") cat("GEVP t0\t=\t", mfit$reference_time, "\n")
  cat("\n")
  cat("ground state energy:\n")
  cat("E \t=\t", mfit$t0[1], "\n")
  cat("dE\t=\t", mfit$error_fn(mfit$t[,1]), "\n")
  if(mfit$model == "pc") {
    cat("\nFitted Delta E:\n")
    cat("Delta E \t=\t", mfit$t0[2], "\n")
    cat("dDelta E\t=\t", mfit$error_fn(mfit$t[,2]), "\n")
    cat("\nThis transltes into a second energy level E2 (be careful with interpretation):\n")
    cat("E2 \t=\t", mfit$t0[2]+mfit$t0[1], "\n")
    cat("dE2\t=\t", mfit$error_fn(mfit$t[,2]+mfit$t[,1]), "\n")
    cat("Effective Amplitude:\n")
    cat("A\t=\t", mfit$t0[3], "\n")
    cat("dA\t=\t", mfit$error_fn(mfit$t[,3]), "\n")    
  }
  else {
    cat("\nAmplitudes:\n")
    for(i in 2:length(mfit$opt.res$par)) {
      cat("P",i-1,"\t=\t", mfit$t0[i], "\n")
      cat("dP",i-1,"\t=\t", mfit$error_fn(mfit$t[,i]), "\n")
    }
  }
  cat("\n")
  cat("boot.R\t=\t", mfit$boot.R, " (bootstrap samples)\n")
  cat("boot.l\t=\t", mfit$boot.l, " (block length)\n")
  cat("useCov\t=\t", mfit$useCov, "\n")
  cat("chisqr\t=\t", mfit$chisqr, "\ndof\t=\t", mfit$dof, "\nchisqr/dof=\t",
      mfit$chisqr/mfit$dof, "\n")
  ## probability to find a larger chi^2 value
  ## if the data is generated again with the same statistics
  ## given the model is correct
  cat("Quality of the fit (p-value):", mfit$Qval, "\n")

  if(any(names(mfit) == "fps")) {
    cat("\nDecay Constant (derived quantity):\n")
    cat("mu1 \t=\t", mfit$mu1, "\n")
    cat("mu2 \t=\t", mfit$mu2, "\n")
    if(mfit$normalisation == "cmi") cat("kappa\t=\t", mfit$kappa,"\n")
    cat("fps \t=\t", mfit$fps, "\n")
    cat("dfps\t=\t", mfit$error_fn(mfit$fps.tsboot), "\n")
  }
  if(any(names(mfit) == "fpsOS")) {
    cat("\nOS Decay Constant (derived quantity):\n")
    if(mfit$normalisation == "cmi") cat("kappa\t=\t", mfit$kappa,"\n")
    cat("fps \t=\t", mfit$fpsOS, "\n")
    cat("dfps\t=\t", mfit$error_fn(mfit$fpsOS.tsboot), "\n")
    cat("using\n")
    cat("ZA  \t=\t", mfit$ZA, "\n")
    cat("dZA \t=\t", mfit$error_fn(mfit$ZAboot), "\n")
  }
}

fit.formatrixboot <- function(cf, par, t, M, LM, Time, parind, sign.vec, ov.sign.vec, lm.avail=FALSE, fitfn, dfitfn, deltat=1, reference_time=0) {
  if(lm.avail && !missing(LM)) {
    opt.res <- minpack.lm::nls.lm(par = par, fn = fitfn, jac = dfitfn, t=t, y=cf, L=LM, Time=Time, parind=parind, sign.vec=sign.vec,
                      deltat=deltat, ov.sign.vec=ov.sign.vec, reference_time=reference_time,
                      control = minpack.lm::nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500, maxfev=5000))
    if( !(opt.res$info %in% c(1,2,3) ) ){
      warning(sprintf("Termination reason of nls.lm opt.res$info: %d\n", opt.res$info))
    }
    opt.res$value <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fitfn, gr = dfitfn, reference_time=reference_time,
                     method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                     t=t, y=cf, M=M, Time=Time, parind=parind, sign.vec=sign.vec, deltat=deltat,
                     ov.sign.vec=ov.sign.vec)
  }
  ##opt.res <- optim(opt.res$par, fn = matrixChisqr, gr = dmatrixChisqr,
  ##                 method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
  ##                 t=t, y=apply(cf,2,mean), M=M, Time=Time, parind=parind, sign.vec=sign.vec)
  return(c(opt.res$par, opt.res$value, opt.res$niter))
}


#' Substract excited states.
#'
#' Excited states are subtracted from the given correlation function and
#' matching matrixfit. The fit is usually done on late time slices when the
#' thermal states have decayed so much that they can be neglected. On the early
#' time slices there are contributions which cannot be explained with a single
#' cosh (or sinh) function. These are exactly the contributions that we do not
#' want.
#'
#' The correlation function is altered on the time slices which are earlier than
#' the start of the fit interval. The correlator is replaced by the model
#' function (cosh or sinh or exp) extrapolated until the first time slice. The
#' deviations of the (bootstrap) samples from the mean value are kept.
#'
#' @param cf Correlation function of class `cf`.
#' @param mfit Fit result of class `matrixfit`.
#' @param from.samples Whether to use existing bootstrap samples. If set to
#'   `TRUE`, the same operation will be applied to the bootstrap samples.
#'   Otherwise the result will not contain bootstrap samples, even if the input
#'   correlation function did.
#'
#' @return A correlation function of class `cf` which is computed from the old
#'   correlation function \eqn{C(t)} as \eqn{M(t) + C(t) - \bar{C}(t)}, where
#'   \eqn{M(t)} is the fit model and \eqn{\bar{C}(t)} denotes the average over
#'   the (bootstrap) samples. Only time slices earlier than the fit are altered.
#'
#' @export
subtract.excitedstates <- function(cf, mfit, from.samples=FALSE) {

  if(inherits(cf, "cf") && inherits(mfit, "matrixfit")) {
    ## we only subtract for 0 <= t < t1 (mind the +1 for the index convention)
    t1p1 <- 1
    t2p1 <- mfit$t1
    ii <- c(t1p1:t2p1)
    Thalfp1 <- cf$Time/2+1
    deltat <- 0
    nfac <- 1.
    sign.vec <- mfit$ov.sign.vec
    if("shifted" %in% names(cf)) {
      deltat <- cf$deltat
      nfac <- 2.
      sign.vec <- -sign.vec
    }
    
    if(mfit$mSize > 1) {	
      for(j in 2:mfit$mSize) {
        ii <- c(ii, (t1p1+(j-1)*Thalfp1):(t2p1+(j-1)*Thalfp1))
      }
    }

    tt <- mfit$CF$t[ii]
    ## compute the difference of mean data to model at times smaller than fit range
    dz <- mfit$cf$cf0[ii] - nfac*matrixModel(par=mfit$opt.res$par, t=tt, Time=cf$Time,
                                        parind=mfit$parind[ii,], sign.vec=sign.vec[ii],
                                        ov.sign.vec=mfit$ov.sign.vec[ii], deltat=deltat)
    cf$subtracted.values <- dz
    cf$subtracted.ii <- ii
    for(i in 1:length(cf$cf[,1])) {
      cf$cf[i,ii] <- mfit$cf$cf[i,ii]-dz
    }
    if(from.samples && cf$boot.samples) {
      cf$cf0[ii] <- nfac*matrixModel(par=mfit$opt.res$par, t=tt, Time=cf$Time,
                                parind=mfit$parind[ii,], sign.vec=sign.vec[ii],
                                ov.sign.vec=mfit$ov.sign.vec[ii], deltat=deltat)
      cf$cf.tsboot$t0[ii] <- cf$cf0[ii]
      for(i in 1:cf$boot.R) {
        cf$cf.tsboot$t[i,ii] <- nfac*matrixModel(par=mfit$t[i, c(1:length(mfit$opt.res$par))],
                                            t=tt, Time=cf$Time, parind=mfit$parind[ii,],
                                            sign.vec=mfit$sign.vec[ii],
                                            ov.sign.vec=sign.vec[ii],
                                            deltat=deltat)
      }
    }
    else{
      cf$boot.sample <- FALSE
      cf$boot.R <- NULL
      cf$boot.l <- NULL
    }
    return(cf)
  }
  else {
    stop("subtract.excitedstates: cf must be of class cf and mfit of class matrixfit. Aborting...\n")
  }
}
