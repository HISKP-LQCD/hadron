\name{matrixfit}
\alias{matrixfit}
\title{Routine For A Factorising Matrix Fit}
\description{
  Performs a factorising fit on a correlation matrix
}
\usage{
matrixfit(cf, boot.R=400, boot.l=2, t1=13, t2=22,
          parlist = array(c(1,1,1,2,2,1,2,2), dim=c(2,4)),
          sym.vec=c("cosh","cosh","cosh","cosh"),
          useCov=FALSE, seed=12345,
          boot.fit=TRUE, fit.method="optim")
}
\arguments{
  \item{cf}{
    correlation matrix obtained with a call to \code{extrac.obs}.
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{boot.R}{
    number of bootstrap samples. Default is \code{400}.
  }
  \item{boot.l}{
    block size for autocorrelation analysis. Default is \code{2}.
  }
  \item{parlist}{
    a two dimensional array of dimension 2 times number of correlators
    in cf. Every column assigns a pair of fit parameters to the
    corresponding correlator in cf. Default is \code{parlist =
      array(c(1,1,1,2,2,1,2,2), dim=c(2,4))} for four correlators.
  }
  \item{sym.vec}{
    a vector of length number of correlators in cf indicating whether
    the correlation function is a cosh, a sinh or an
    exponential. Default is 
    \code{sym.vec=c("cosh","cosh","cosh","cosh")}, possible values are
    \code{"cosh"}, \code{"sinh"} and \code{"exp"}.
  }
  \item{useCov}{
    use correlated or uncorrelated chisquare. Default is
    \code{useCov=FALSE}. 
  }
  \item{seed}{
    seed for the random number generation used for boostrapping. Default
    is \code{12345}.
  }
  \item{boot.fit}{
    If set to \code{FALSE}, the fit is not bootstrapped, even if the bootstrapping 
    parameters have been set and the correlation function has been bootstrapped.
    This is a useful time-saver if error information is not strictly necessary. 
    Of course, this affects the return values related to the bootstrap, which are
    set to \code{NA}.
  }
  \item{fit.method}{
    Can be either \code{"optim"} or \code{"lm"}. The latter works only
    if the library \code{"minpack.lm"} can be loaded. Default and
    fallback is \code{"optim"}.
  }
}
\value{
  returns an object of class \code{matrixfit} with entries:
  \item{CF}{
    object of class cf which contains the mean correlation functions
  }
  \item{M}{
    weighting matrix for weighted Chi squared minimization
  }
  \item{parind}{
    indices in the parameter vector used for the different matrix combinations
  }
  \item{sign.vec}{
    vector of signs 
  }
  \item{ii}{
    vector of vector indices giving the columns of the correlation function arrays
    (CF above, say), which are contained in the fit range
  }
  \item{opt.res}{
    return value of the minimization (see ?optim)
  }
  \item{opt.tsboot}{
    bootstrap analysis of Chi squared minimization
  }
  \item{useCov}{
    whether covariances in the data were taken into account
  }
  \item{invCovMatrix}{
    inverse of covariance matrix or inverse variance weighted if useCov=FALSE
  }
  \item{Qval}{
    real number between 0 and 1 giving the "quality" of the fit
  }
  \item{chisqr}{
    total Chi squared of the fit
  }
  \item{dof}{
    fit degrees of freedom
  }
  \item{mSize}{
    integer size of the matrix which was fitted
  }
  \item{cf}{
    object of type cf which contains, amongst other objects, cf$cf which is a
    concatenated array of raw correlation functions where each row is one of N 
    observations and there are mSize*T columns (see ?extract.obs)
  }
  \item{boot.R}{
    number of bootstrap samples
  }
  \item{boot.l}{
    block size for blocked bootstrap
  }
  \item{t1}{
    beginning of fit range
  }
  \item{t2}{
    end of fit range
  }
  \item{parlist}{
    array of parameter combinations for the matrix fit
  }
  \item{sym.vec}{
    vector of strings indicating the functional form of correlation functions which were fitted
  }
  \item{seed}{
    RNG seed for bootstrap procedure
  }
}
\details{
  The routine expects in \code{cf$cf} a set of correlation
  functions.  The
  mapping of this linear construct to a matrix or a part of a matrix is
  achieved via 
  \code{parlist}. The symmetry properties of the individual correlation
  functions must be encoded in \code{sym.vec}.

  \code{matrixfit} will fit to every correlator in \code{cf$cf} a
  function \eqn{p_i p_j f(t)}. The indices \eqn{i,j} are determined from
  \code{parlist} and \eqn{f} is either \eqn{cosh}{\cosh} or
  \eqn{sinh}{\sinh}, depending on \code{sym.vec}.
  
  The inverse covariance matrix is computed using a singular value
  decomposition. If the sample size N is too small, only sqrt(N)
  eigenvalues of the matrix are kept exactly, while all others are
  replaced by the mean of the rest. This helps to reduce instabilities
  induced by too small eigenvalues of the covariance matrix.
}
\references{
  C. Michael, \href{hep-lat/9412087}{hep-lat/9412087}
}
\seealso{
  \code{\link{cf}}, \code{\link{bootstrap.cf}}
}
\examples{
\dontrun{readcmidatafiles(basename="outpr", verbose=TRUE, skip=1,
last.digits=4)}
\dontrun{cf <- extract.obs(cmicor, verbose=T, vec.obs=c(1))}
\dontrun{x <- matrixfit(cf, boot.l=2, t1=13, t2=22, useCov=TRUE)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimize}
\keyword{ts}

