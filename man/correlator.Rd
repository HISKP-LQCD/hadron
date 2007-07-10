\name{correlator}
\alias{correlator}
\title{Computes the correlator estimate plus errors}
\description{
  Computes the correlator estimate plus errors using the Gamma method.
}
\usage{
correlator(filename, psscar, skip = 0, S=1.5, plotit=F, fold=F, symmetric=T, scalefactor=1.)
}
\arguments{
  
  \item{filename}{
    filename of the file containing N measurements of the correlator in
    GWC format, meaning first collumn t, second collumn correlator. If
    it is missing \code{psscar} must contain the data.
  }
  \item{psscar}{
    \code{data.frame} containing the data with first collumn t, second
    collumn correlator. Only used if \code{filename} is missing.
  }
  \item{skip}{
    no of measurements to skip at the beginning of the file
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{plotit}{
    logical. if TRUE a plot i drawn.
  }
  \item{fold}{
    logical. if TRUE positive and negative times are averaged.
  }
  \item{symmetric}{
    logical. if TRUE and \code{fold} is set to TRUE the correlator is
    assumed to be symmetric in time (cosh) and if FALSE anti-symmetric.
  }
  \item{scalefactor}{
    factor to rescale the correlator.
  }
}
\value{
  returns a data.frame with the results. The object is also of type
  \code{correlator} which can be plotted using the generic function
  \code{plot}.
}
%\details{
  
%}
%\references{
%}
\seealso{
  \link{pp}, \link{ppeffectivemass}, \link{ppcorr}
}
\examples{
data(pscor.sample)
ppcorr.res <- correlator(psscar=pscor.sample, from=5, to=21)
plot(ppeff.res)
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{ts}
