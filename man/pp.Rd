\name{pp}
\alias{pp}
\alias{coshfit}
\title{Perform A Cosh Fit To 2 Point Correlators}
\description{
  Perform A Cosh Fit To 2 Point Correlators. Mass, amplitude and decay
  constant as well as integrated autocorrelation times will be detemined
  for a moving window in the time extent with fixed upper bound.
}
\usage{
pp(filename, skip=0, from, to, S=1.5, A=0.01, m=0.01, plot=FALSE, debug=FALSE)
}
\arguments{
  \item{filename}{
    filename of the file from which to read the correlator. It is
    supposed to be in GWC code format
  }
  \item{skip}{
    number of measurements to be discarded at the beginning of the time series.
  }
  \item{from}{
    lower bound for the fitrange in time (from,to). Counting starts with
    0. Fits will be done for intervals (from,to), (from+1, to), ...
  }
  \item{to}{
    upper bound for the fitrange in time (from,to). Counting starts with
    0.  Fits will be done for intervals (from,to), (from+1, to), ...
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{A}{
    initial guess for the amplitude
  }
  \item{m}{
    initial guess for the mass
  }
  \item{plot}{
    logical: if set to TRUE the function produces plots
  }
  \item{debug}{
    logical: if set to TRUE some debug output and plots are produced
  }
}
\value{
  returns a \code{\link{data.frame}} containing the results
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{rho}}, \code{\link{pion}}, \code{\link{uwerr}}
}
%\examples{
%}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
