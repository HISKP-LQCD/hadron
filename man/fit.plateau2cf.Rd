\name{fit.plateau2cf}
\alias{fit.plateau2cf}
\title{
  fits a plateau to an object of class \code{cf}
}
\description{
  where applicable, a plateau is fitted to the averaged data in
  \code{cf} using a (correlated) chisquare fit.
}
\usage{
fit.plateau2cf(cf, t1, t2, useCov=FALSE)
}
\arguments{
  \item{cf}{
    input object of class \code{cf}
  }
  \item{t1}{
    starting t-value for the fit
  }
  \item{t2}{
    final t-value for the fit.
  }
  \item{boot.samples}{
    use bootstrap samples to estimate the variance-covariance matrix
    instead of the original data
  }
  \item{boot.R}{
    number of bootstrap samples
  }
  \item{boot.l}{
    block length for blocked bootstrap
  }
  \item{useCov}{
    perform a correlated chisquare fit or not.
  }
}
\value{
  Returns a list with elements
  \item{plateau}{
    the fitted plateau value
  }
  \item{dplateau}{
    its error
  }
}
\author{Carsten Urbach \email{curbach@gmx.de}}
\seealso{
  \code{\link{cf}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bootstrap}
\keyword{fit}
