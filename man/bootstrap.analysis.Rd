\name{bootstrap.analysis}
\alias{bootstrap.analysis}}
\title{Performs a Bootstrap and Blocking Analysis of a Timeseries}
\description{
  Performs a Bootstrap and Blocking Analysis of a Timeseries
}
\usage{
bootstrap.analysis <- function(data, skip=0, boot.R=100, tsboot.sim="geom", pl=F)
}
\arguments{

  \item{data}{
    a numerical vector containing the time series
  }
  \item{skip}{
    integer value providing the warm up phase length. 
  }
  \item{boot.R}{
    number of bootstrap samples. See also \code{boot}, and \code{tsboot}
    
  }
  \item{tsboot.sim}{
    The type of simulation required to generate the replicate
    time series. See \code{tsboot} for details.
  }
  \item{pl}{
    logical, indicating whether or not to plot the result.
  }
}
\value{
  returns a data frame containing the mean value, the error
  approximation, the estimate of the error of the error and the value of
  tau int for all block sizes.
}
\details{
  the routine will compute the error, the error of the error and the
  integrated autocorrelation time for different block size using a
  bootstrap analysis. The blocksize is systematically increased starting
  from \code{1} until \code{(length(data)-skip)/blocksize < 20}. 
}

\seealso{
  for an alternative way to analyse such time series see \code{\link{uwerr}}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{ts}
