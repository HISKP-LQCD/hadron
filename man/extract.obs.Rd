\name{extract.obs}
\alias{extract.obs}
\alias{extrac.loop}
\title{Extract One or More Gamma Combinations from am CMI Correlator}
\description{
  Extracts one or more gamma matrix combinations (observables) from a
  correlator stored in cmi format
}
\usage{
cmicor <- extract.obs(cmicor, vec.obs=c(1), ind.vec=c(1,2,3,4,5),
                      sym.vec, verbose=FALSE)

cmiloop <- extrac.loop(cmiloop, obs=9, ind.vec=c(2,3,4,5,6,7,8), L)
}
\arguments{
  \item{vec.obs}{
    vector containing the numbers of observables to be extracted.
  }
  \item{obs}{
    number of a single loop to extract.
  }
  \item{ind.vec}{
    Index vector indexing the column numbers in cmicor to be used. The
    first must be the observable index, the second the smearing type
    index, the third the time, the fourth C(+t) and the fifth C(-t).

    Index vector indexing the column numbers in cmiloop to be used. The
    first must be the observable index, the second the smearing type
    index, the third the time, the fourth ReTL, the fifth ImTL, the
    sixth ReTF and the seventh ImTF.
  }
  \item{verbose}{
    Increases verbosity of the function.
  }
  \item{sym.vec}{
    a vector of bools of length equal to the number of observables
    indicating whether C(+t) and C(-t) should be added or subtracted. If
    not given C(+t) and C(-t) will be averaged for all observables.
  }
  \item{sign.vec}{
    a sign vector of length equal to the number of observables
    indicating whether the corresponding correlation function should be 
    multiplied by +-1.
  }
}
\value{
  returns a list containing
  \item{cf}{
    for \code{extract.obs}: array containing the correlation function
    with dimension number of 
    files times (nrObs*nrStypes*(T/2+1)). C(t) and C(-t) are averaged
    according to \code{sym.vec}.

    for \code{extract.loop}: ReTL
  }
  \item{icf}{
    for \code{extract.loop} only: ImTL
  }
  \item{scf}{
    for \code{extract.loop} only: ReTF
  }
  \item{sicf}{
    for \code{extract.loop} only: ImTF
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
  \item{nrSamples}{
    for \code{extrac.loop} only: the number of samples found in the
    files. 
  }
}
\details{
  C(t) and C(-t) are averaged as indicated by \code{sym.vec}.

  The returned object can be further processed by e.g. \code{block.ts}.
}
\references{
}
\seealso{
  \code{\link{readcmicor}}, \code{\link{readcmidatafiles}},
  \code{\link{block.ts}}
}
\examples{
\dontrun{cmicor <- readcmidatafiles("outprc", skip=1)}
\dontrun{cf <- extract.obs(cmicor, vec.obs=c(1,3))}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{ts}
