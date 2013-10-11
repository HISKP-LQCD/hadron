\name{matrixfit}
\alias{matrixfit}
\title{Routine For A Factorising Matrix Fit}
\description{
  Performs a factorising fit on a correlation matrix
}
\usage{
matrixfit(cf, boot.l=2, t1=13, t2=22, useCov=TRUE)
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
  \item{symmetrise}{
    Symmetrise the matrix, not implemented yet.
  }
  \item{boot.R}{
    number of bootstrap samples
  }
  \item{boot.l}{
    block size for autocorrelation analysis
  }
  \item{parlist}{
    a two dimensional array of dimension 2 times number of correlators
    in cf. Every column assigns a pair of fit parameters to the
    corresponding correlator in cf. Default is \code{parlist =
      array(c(1,1,1,2,2,1,2,2), dim=c(2,4))} for four correlators.
  }
  \item{sym.vec}{
    a vector of length number of correlators in cf indicating whether
    the correlation function is a cosh or a sinh. Default is
    \code{sym.vec=c("cosh","cosh","cosh","cosh")}. 
  }
  \item{matrix.size}{
    matrix size, currently only 2 implemented
  }
  \item{useCov}{
    use correlated or uncorrelated chisquare. Default is
    \code{useCov=FALSE}. 
  }
  \item{seed}{
    seed for the random number generation used for boostrapping.
  }
}
\value{
  returns an object of class \code{matrixfit}
}
\details{
  The routine expects in \code{cf$cf} a matrix. Every row corresponds to
  one observation of a concatenation of correlation functions. The
  mapping of this linear construct to a matrix is achieved via
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

