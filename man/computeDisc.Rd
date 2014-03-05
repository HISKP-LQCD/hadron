\name{computeDisc}
\alias{computeDisc}
\title{computes a disconnected correlation function from loops}
\description{
  computes a disconnected correlation function from loops
}
\usage{
cf <- computeDisc(cf, cf2, real=TRUE, real2=TURE, subtract.vev=TRUE,
subtract.vev2=TRUE, subtract.equal=TRUE, use.samples, use.samples2,
smeared=FALSE, smeared2=FALSE, verbose=FALSE)
}
\arguments{
  \item{cf}{
    loop data as produced by \code{readcmidisc} or \code{readbinarydisc}.
  }
  \item{cf2}{
    second set of loop data as produced by \code{readcmidisc} or
    \code{readbinarydisc}. This is needed for cross-correlators
  }
  \item{real}{
    use the real part \code{cf$cf}, if set to \code{TRUE}, otherwise the
    imaginary part \code{cf$icf}. 
  }
  \item{real2}{
    use the real part \code{cf2$cf}, if set to \code{TRUE}, otherwise the
    imaginary part \code{cf2$icf}. 
  }
  \item{subtract.vev}{
    subtract a vacuum expectation value. It will be estimated as mean
    over all samples, gauges and times available.
  }
  \item{subtract.vev2}{
    subtract a vacuum expectation value for the second set of loops. It will be estimated as mean
    over all samples, gauges and times available.
  }
  \item{subtract.equal}{
    subtract contributions of products computed on identical
    samples. This will introduce a bias, if
    set to FALSE for missing cf2 or if cf and cf2 are computed on the same
    set of random sources.
  }
  \item{use.samples}{
    If set to an integer, only the specified number of samples will be
    used for \code{cf}, instead of all samples.
  }
  \item{use.samples2}{
    Same like \code{use.samples}, but for \code{cf2}.
  }
  \item{smeared}{
    use the loops instead of the local ones for \code{cf}.
  }
  \item{smeared2}{
    use the loops instead of the local ones for \code{cf2}.
  }
  \item{verbose}{
    Print some debug output, like the VEVs of the loops.
  }
}
\value{
  Returns an object of type \code{cf} derived from a \code{list}
  with elements \code{cf}, an array of dimension
  \code{dim(N, T)}, where \code{N} is the number of samples and \code{T}
  the time extend, integers \code{Time} for the time extend,
  \code{nrStypes} and \code{nrObs} for the available smearing types and
  operators, and finally \code{nrSamples}, the number of samples used to
  generate the correlation function \code{cf}. 
}
\details{
  The dimension of \code{cf$cf} and \code{cf$icf} must be \code{dim(T,
    S, N)}, where \code{T} is the time extend, \code{S} is the number of
  samples and \code{N} the number of measurements (gauges). \code{cf2}
  is the same, but needed only for cross-correlators.
  
  If \code{subtract.vev=TRUE} the vev is estimated as the mean over all
  gauges, samples and times available and subtracted from the original
  loop data. (Same for \code{subtrac.vev2}.

  The correlation is computed such as to avoid correlation between equal
  samples, unless \code{nrSamples} is equal to 1.

  \code{cf} and \code{cf2} must agree in \code{T}, number of gauges and
  number of samples. Matching of gauges is assumed. If this is not the
  case results are wrong.
}
\seealso{
  \code{\link{readcmidisc}}, \code{\link{readbinarydisc}},
  \code{\link{bootstrap.cf}}, \code{\link{add.cf}},
  \code{\link{c.cf}}
}
\examples{
\dontrun{files <- getorderedfilelist(basename="disc.0.13872.0.003.k0v.", last.digits=3)}
\dontrun{vdata <- readcmiloopfiles(files)}
\dontrun{loops <- extract.loop(vdata, obs=9)}
\dontrun{Cpi0v4 <- computeDisc(cf=loops, real=TRUE, subtract.vev=TRUE)}
\dontrun{Cpi0v4 <- bootstrap.cf(Cpi0v4, boot.R=1000, boot.l=20)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{correlator}
