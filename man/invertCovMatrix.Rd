\name{invertCovMatrix}
\alias{invertCovMatrix}
\title{Inverts the covariance matrix for noisy data}
\description{
  The covariance matrix of noisy data is inverted. Special care is taken
  in treating spurious small modes of the matrix, which are likely to
  arise from too much noise in the data.
}
\usage{
invertCovMatrix(cf, boot.l=1, boot.samples=FALSE, cov_fn = cov)
}
\arguments{
  \item{cf}{
    The data for which the covariance matrix is to be computed. It is
    expected to be an array or matrix with dimension RxN, where R is the
    number of observations and N the number of observables.

    \code{cf} can be either real data or bootstrap data. In the latter
    case \code{boot.samples=TRUE} must be set for proper normalisation
    of the inverse matrix.
  }
  \item{boot.l}{
    If set to a value larger than 1 the data will be blocked with
    blocklength \code{boot.l} before the
    covariance matrix is computed.
  }
  \item{boot.samples}{
    If set to \code{TRUE} the data is treated a pseudo data from a
    bootstrap procedure.
  }
  \item{cov_fn}{
    Function that computes the covariance matrix from the given samples.
  }
}
\value{
  Returns the inverse covariance matrix as an object of class
  \code{\link{matrix}}. 
}
\details{
  The inverse covariance matrix is estimated. If the number of
  observations is too small the procedure described in the reference is
  used to remove spuriously small eigenvalues of the covariance matrix.

  We always keep the \eqn{\sqrt{R}}{sqrt(R)} largest eigenvalues exactly
  and replace the remaining smallest ones by their mean.
}
\references{
  C.Michael, A.McKerrell,  Phys.Rev. D51 (1995) 3745-3750, hep-lat/9412087
}
\seealso{
  \code{\link{cov}}, \code{\link{matrix}}
}
\examples{
\dontrun{M <- invertCovMatrix(effMass.tsboot, boot.samples=TRUE)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{covariance matrix}
\keyword{correlated chisqr}
