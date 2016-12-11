\name{jackknife.cf}
\alias{jackknife.cf}
\title{jackknife a set of correlation functions}
\description{
  jackknife a set of correlation functions
}
\usage{
jackknife.cf(cf, boot.l=2)
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
  added for the correlation function called \code{cf.jackknife}. Moreover, the
  original average of \code{cf} is returned as \code{cf0} and the
  bootstrap errors as \code{jackknife.se}. We also copy the input
  parameters over and set \code{jackknife.samples} to \code{TRUE}.
}
\seealso{
  \code{boot::tsboot}, \code{bootstrap.cf}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{bootstrap, time series, }
