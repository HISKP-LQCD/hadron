\name{uwerr}
\alias{uwerr}
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
    function computing the derived quantity. If not given it is assumed
    that a primary quantity is analysed.

    f must have the data vector of length Nalpha as the first
    argument. Further arguments to f can be passed to uwerr via the
    \code{...} argument.

    f may return a vector object of numeric type.
  }
  \item{data}{
    array of data to be analysed. It must be of dimension (N x Nalpha)
    (i.e. N rows and Nalpha columns),
    where N is the total number of measurements and Nalpha is the number
    of primary observables
  }
  \item{nrep}{
    the vector (N1, N2, ...) of replica length N1, N2
  }
  \item{S}{
    initial guess for the ratio tau/tauint, with tau the exponetial
    autocorrelation length.
  }
  \item{pl}{
    logical: if TRUE, the autocorrelation function, the integrated
    autocorrelation time as function of the integration cut-off and (for
    primary quantities) the time history of the observable are plotted
    with plot.uwerr
  }
}
\value{
  In case of a primary observable (\code{uwerrprimary}), an object of
  class \code{uwerr} with basis class \code{\link{list}} containing the
  following objects
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

  In case of a derived observable (\code{uwerrderived}), i.e. if a
  function is specified, the above objects are contained in a list
  called \code{res}.

  In both cases the return object containes
  \item{xxx}{
    xxx
  }
}
%\details{
%}
\references{
  ``Monte Carlo errors with less errors'', Ulli Wolff, hep-lat/0306017
}
\seealso{
  \code{\link{plot.uwerr}}
}
\examples{
data(plaq.sample)
plaq.res <- uwerrprimary(plaq.sample)
summary(plaq.res)
plot(plaq.res)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimize}
\keyword{ts}
