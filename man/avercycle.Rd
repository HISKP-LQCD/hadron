\name{avercycle}
\alias{avercycle}
\title{take the cycle average of correlator data in cmicor format}
\description{
  take the cycle average of correlator data in cmicor format
}
\usage{
avercycle(cmicor, cycle.l)
}
\arguments{
  \item{cmicor}{
    correlator data as obtained with \link{\code{readcmicor}}.
  }
  \item{cycle.l}{
    the cycle length in units of gauge numbering
  }
}
\value{
  returns the averaged correlator data in the same format as provided by
  \link{\code{readcmicor}}.
}
\details{
  the function will determine the average number of measurements per
  cycle, average the data in each cycle and take care of the appropriate
  weight.

  Currently the implementation is very slow. why?

  Averaging the cycles restores translational invariance of the data and
  takes properly care of the correlations in between in within the cycles.
}
\seealso{
  \link{\code{readcmicor}}
}
\examples{
\dontrun{avercycle(cmicor, 10)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{ts}
