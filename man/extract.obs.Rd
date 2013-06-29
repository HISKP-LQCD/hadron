\name{extract.obs}
\title{Extract One or More Gamma Combinations from am CMI Correlator}
\description{
  Extracts one or more gamma matrix combinations (observables) from a
  correlator stored in cmi format
}
\usage{
extract.obs(cmicor, vec.obs=c(1), ind.vec=c(1,2,3,4,5), sym.vec, verbose=FALSE)
}
\arguments{
  \item{vec.obs}{
    vector containing the numbers of observables to be extracted.
  }
  \item{ind.vec}{
    Index vector indexing the column numbers in cmicor to be used. The
    first must be the observable index, the second the smearing type
    index, the third the time, the fourth C(+t) and the fifth C(-t).
  }
  \item{verbose}{
    Increases verbosity of the function.
  }
  \item{sym.vec}{
    a vector of bools of length equal to the number of observables
    indicating whether C(+t) and C(-t) should be added or subtracted. If
    not given C(+t) and C(-t) will be averaged for all observables.
  }
}
\value{
  returns a list containing
  \item{cf}{
    array containing the correlation function with dimension number of
    files times (nrObs*nrStypes*(T/2+1)). C(t) and C(-t) are averaged
    according to \code{sym.vec}.
  }
  \item{Time}{
    The time extend of the correlation functions.
  }
  \item{nrStypes}{
    The number of smearing combinations.
  }
  \item{nrObs}{
    The number of observables.
  }
}
\details{
  C(t) and C(-t) are averaged as indicated by \code{sym.vec}. This
  object can be further processed by e.g. \code{block.ts}.
}
\references{
}
\seealso{
  \code{\link{readcmicor}}, \code{\link{readcmidatafiles}},
  \code{\link{block.ts}}
}
\examples{
\dontrun{cmicor <- readcmidatafiles("outprc", skip=1)}
\dontrun{cf <- extract.obs(cmicor, vec.obs=c(1,3))
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{ts}
