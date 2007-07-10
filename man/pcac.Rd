\name{pcac}
\alias{pcac}
\title{Computes the pcac mass}
\description{
  Computes the pcac mass from the PP and the AP (PA) correlators and
  estimates the errors using the gamma method
}
\usage{
pcac(psfilename, apfilename, pafilename, from=3, to=3, fit=F, skip=0, plotit=F, S=1.5)
}
\arguments{
  \item{psfilename}{
    filename of the file from which to read the PP correlator. It is
    supposed to be in GWC code format. mandatory.
  }
  \item{apfilename}{
    filename of the file from which to read the AP correlator. It is
    supposed to be in GWC code format. Either PA or AP correlator (or
    both) must be given. If both are given, both are used by averaging.
  }
  \item{pafilename}{
    filename of the file from which to read the PA correlator. It is
    supposed to be in GWC code format. Either PA or AP correlator (or
    both) must be given. If both are given, both are used by averaging.
  }
  \item{from}{
    the effective mass is computed starting with t=from
  }
  \item{to}{
    the effective mass is computed ending with t=to
  }
  \item{fit}{
    logical. if TRUE a fit is performed to all t-values to determine the
    pcac mass.
  }
  \item{skip}{
    no of measurements to skip at the beginning of the file
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{plotit}{
    logical. if TRUE a plot is drawn.
  }
}
\value{
  returns a data.frame with the results. The object is also of class
  \code{massfit} which can be plotted using the generic function
  \code{plot}.
}
\details{
  the symmetric difference operator is used.
}
%\references{
%}
%\seealso{
%}
%\examples{
%}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}

