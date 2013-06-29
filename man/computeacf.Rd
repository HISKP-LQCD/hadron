\name{computeacf}
\alias{}
\title{Computes The ACF and Integrated AC Time}
\description{
  Computes the ACF and integrated autocorrelation time of a time
  series. It also estimates the corresponding standard errors.
}
\usage{
computeacf(tseries, W.max, Lambda=100)
}
\arguments{
  \item{tdata}{
    the time series.
  }
  \item{W.max}{
    maximal time lag to be used.
  }
  \item{Lambda}{
    cut-off needed to estimate the standard error of the ACF.
  }
}
\value{
  It returns a list of class \code{hadronacf} with members
  \item{lags}{
    time lags of the integrated autocorrelation function
  }
  \item{Gamma}{
    normalised autocorrelation function
  }
  \item{dGamma}{
    error of normalised autocorrelation function
  }
  \item{W.max}{
    max time lag used for the call of \code{acf}
  }
  \item{W}{
    the cut-off up to which the ACF is integrated for the integrated
    autocorrelation time
  }
  \item{tdata}{
    the original time series
  }
  \item{tau}{
    the estimated integrated autocorrelation time
  }
  \item{dtau}{
    the estimated error of the integrated autocorrelation time
  }

}
\details{
  The standard error of the ACF is computed using equation (E.11) of
  M. Luescher, hep-lat/0409106. The error of the integrated
  autocorrelation time using the Madras Sokal formula, see also
  hep-lat/0409106. 
}
\references{
  'Monte Carlo errors with less errors', Ulli Wolff,
  \href{http://arxiv.org/abs/hep-lat/0306017}{hep-lat/0306017}
  
  'Schwarz-preconditioned HMC algorithm for two-flavour lattice QCD',
  Martin Luescher, \href{http://arxiv.org/abs/hep-lat/0409106}{hep-lat/0409106}
  
  N. Madras, A. D. Sokal, J. Stat. Phys. 50 (1988) 109
}
\seealso{
  \code{\link{uwerr}}, \code{\link{acf}}, \code{\link{plot.hadronacf}},
  \code{\link{summary.hadronacf}}, \code{\link{bootstrap.analysis}}
}
\examples{
data(plaq.sample)
myacf <- computeacf(tsdata, 300)
plot(myacf)
summary(myacf)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{ts, acf, autocorrelationtime}
