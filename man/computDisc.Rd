\name{computeDisc}
\alias{computeDisc}
\title{computes a disconnected correlation function from loops}
\description{
  computes a disconnected correlation function from loops
}
\usage{
cf <- computeDisc(cf, real=TRUE, subtract.vev=TRUE)
}
\arguments{
  \item{cf}{
    loop data as produced by \code{readcmidisc} or \code{readbinarydisc}.
  }
  \item{real}{
    use the real part \code{cf$cf}, if set to \code{TRUE}, otherwise the
    imaginary part \code{cf$icf}. 
  }
  \item{subtract.vev}{
    subtract a vacuum expectation value. It will be estimated as mean
    over all samples, gauges and times available.
  }
}
\value{
  Returns a list with elements \code{cf}, an array of dimension
  \code{dim(N, T)}, where \code{N} is the number of samples and \code{T}
  the time extend, integers \code{Time} for the time extend,
  \code{nrStypes} and \code{nrObs} for the available smearing types and
  operators, and finally \code{nrSamples}, the number of samples used to
  generate \code{cf}. 
}
\details{
  If \code{subtract.vev=TRUE} the vev is estimated as the mean over all
  gauges, samples and times available and subtracted from the original
  loop data.

  The correlation is computed such as to avoid correlation between equal
  samples, unless \code{nrSamples} is equal to 1.
}
\seealso{
  \code{\link{readcmidisc}}, \code{\link{readbinarydisc}}
}
\examples{
\dontrun{Cpi0v4 <- computeDisc(v4data, real=TRUE, subtract.vev=TRUE)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{correlator}
