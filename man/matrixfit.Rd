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
  \item{noObs}{
    number of gamma matrix combinations in the matrix. Currently only 1
    is implemented
  }
  \item{matrix.size}{
    matrix size, currently only 2 implemented
  }
  \item{useCov}{
    use correlated chisquare
  }
}
\value{
  returns an object of class \code{matrixfit}
}
\details{
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

