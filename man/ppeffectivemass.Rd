\name{ppeffectivemass}
\alias{ppeffectivemass}
\title{Computes the effective masses for a single correlator}
\description{ Computes the effective masses f a single correlator
  which has a cosh form and performs a Gamma method error analysis
}
\usage{
ppeffectivemass(filename, psscar, from, to, skip=0, S=1.5, plotit=F)
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
  \item{from}{
    the effective mass is computed starting with t=from
  }
  \item{to}{
    the effective mass is computed ending with t=to
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
    
}
\value{
  returns a data.frame with the results. The object is also of type
  \code{massfit} which can be plotted using the generic function
  \code{plot}.
}
\details{
  The inverse \code{cosh} is explicitly computed numerically for the
  given ratio at given t. It is averaged over forward and backward in
  time. 
}
\seealso{
  \link{pp}, \link{pion}
}
\examples{
data(pscor.sample)
ppeff.res <- ppeffectivemass(psscar=pscor.sample, from=5, to=21)
plot(ppeff.res)
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{ts}
