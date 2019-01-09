\name{jackknife.cf}
\alias{jackknife.cf}
\title{jackknife a set of correlation functions}
\description{
  jackknife a set of correlation functions
}
\usage{
jackknife.cf(cf, boot.l = 1)
}
\arguments{
  \item{cf}{
    correlation matrix of class \code{cf} e.g. obtained with a call to \code{extrac.obs}.
  }
  \item{boot.l}{
    block size for autocorrelation analysis
  }
}
\value{
  returns an object of class \code{cf} with blocked jackknife samples
  added for the correlation function called \code{cf.jackknife}.
  Currently, only the moving block jackknife approach is implemented.
  Moreover, the
  original average of \code{cf} is returned as \code{cf0} and the
  bootstrap errors as \code{jackknife.se}. We also copy the input
  parameters over and set \code{jackknife.samples} to \code{TRUE}.
}
\seealso{
  \code{boot::tsboot}, \code{bootstrap.cf}
}
\references{
  H.R. Künsch,
  "The jackknife and the bootstrap for general stationary observations",
  The Annals of Statistics, 1989, Vol. 17, No. 3, 1217-1241

  S.N. Lahiri,
  "On the jackknife-after-bootstrap method for dependent data and its consistency properties",
  Econometric Theory, 2002, Vol. 18, 79-98
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{bootstrap, time series, }
