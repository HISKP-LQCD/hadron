\name{uwerr}
\alias{uwerrprimary}
\alias{uwerrderived}
\title{Time Series Analysis With Gamma Method}
\description{
  Analyse time series data with the so called gamma method
}
\synopsis{
uwerr(f, data, nrep, S=1.5, pl=FALSE, ...)
}
\usage{
uwerr(f, data, nrep, S=1.5, pl=FALSE, ...)
uwerrprimary(data, nrep, S=1.5, pl=FALSE)
uwerrderived(f, data, nrep, S=1.5, pl=FALSE, ...)
}
\arguments{

  \item{f}{
    function taking the averaged data and returning the derived quantity
    as one numeric object. If not given it is assumed that a primary
    quantity is analysed.
  }
  \item{data}{
    array of data to be analysed. It must be of dimension (Nalpha x N),
    where N is the total number of measurements and Nalpha is the number
    of primary observables
  }
  \item{nrep}{
    the vector (N1, N2, ...) of replica length N1, N2
  }
  \item{S}{
    initial guess for tau/tauint
  }
  \item{pl}{
    boolean indicating whether to produce some plots
  }
}
\value{
  returns a list containing
  \item{value}{
    the expectation value of the obsevable
  }
  \item{dvalue}{
    the error estimate
  }
  \item{ddvalue}{
    estimate of the error on the error
  }
  \item{tauint}{
    estimate of the integrated autocorrelation time for that quantity
  }
  \item{dtauint}{
    error of tauint
  }
  \item{Wopt}{
    value of optimal cut-off for the Gamma function integration
  }
}
%\details{
%}
\references{
  ``Monte Carlo errors with less errors'', Ulli Wolff, hep-lat/0306017
}
%\seealso{
%}
%\examples{
%\dontrun{uwerr(plaqhist)}
%}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
