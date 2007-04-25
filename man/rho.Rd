\name{rho}
\alias{rho}
\title{Fit Vector Meson}
\description{
  fits one or several cosh and sinh to data for extracting the vector
  meson mass and excited states
}
\usage{
rho(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
    variational=list(ta=4, tb=5, N=6), ind.vec=c(1,3,4,5),
    no.masses=1, matrix.size = 2)
}
\arguments{
  \item{cmicor}{
    data to be fitted to as e.g. the output of \code{\link{readcmicor}}
  }

  \item{mu}{
    the value of the bare quark twisted mass
  }
  \item{kappa}{
    the value of the hopping parameter
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{pl}{
    boolean. if set to TRUE the function produces plots
  }
  \item{skip}{
    number of measurements to be discarded at the beginning of the time series.
  }
  \item{variational}{
    list of parameters used for the variational analysis
  }
  \item{ind.vec}{
    index vector indexing the column numbers in cmicor to be used
  }
  \item{no.masses}{
    number of masses to be extracted. This argument will set the number
    of exponentials to be used in the fit
  }
  \item{matrix.size}{
    matrix size to be used in the fit. Can be currently set to 2,4 and 6.
  }
}
\value{
  returns an object of type cfit with the following items

  \item{fitresult}{
    result from the fit as returned by \code{\link{optim}}
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{N}{
    number of measurements found in the data
  }
  \item{Time}{
    Time extent found in the data
  }
  \item{fitdata}{
    \code{\link{data.frame}} containing the time values used in the fit,
    the averaged correlator and its error and the value of Chi for each
    time value
  }
  \item{uwerrresultmv}{
    the result of the time series analysis for the lowest mass as
    carried out byt \code{\link{uwerr}} 
  }
  \item{uwerrresultmv2}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}} if no.masses larger than 1.
  }
  \item{uwerrresultmv3}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}}, if no.masses larger than 2.
  }
  \item{effmass}{
    ...
  }
  \item{mu}{
    the value of the bare quark twisted mass
  }
  \item{kappa}{
    the value of the hopping parameter
  }
  \item{variational.masses}{
    mass values as determined by the variational analysis
  }
  \item{no.masses}{
    no.masses determined
  }
  \item{matrix.size}{
    size of the data matrix
  }
  
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{pion}}, \code{\link{readcmicor}}, \code{\link{uwerr}}, \code{\link{plot.cfit}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{rhofit <- rho(cmicor, kappa=0.160856, mu=0.004, t1=8, t2=18, no.masses=1, matrix.size=4)}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
