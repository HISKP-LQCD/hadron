\name{bootstrap.gevp}
\alias{bootstrap.gevp}
\title{perform a bootstrap analysis of a GEVP}
\description{
  perform a bootstrap analysis of a GEVP for a real, symmetric correlator matrix
}
\usage{
bootstrap.gevp(cf, t0, boot.R=400, boot.l=2, matrix.size=2, element.order=c(1,2,3,4), seed=1234)
}
\arguments{
  \item{cf}{
    correlation matrix obtained with a call to \code{extrac.obs}.
  }
  \item{t0}{
    initial time value of the GEVP
  }
  \item{matrix.size}{
    the size of the square matrix
  }
  \item{element.order}{
    specifies how to fit the \code{n} linearly ordered single
    correlators into the correlator
    matrix. \code{element.order=c(1,2,3,4)} leads to a matrix
    \code{matrix(cf[element.order], nrow=2)}.
  }
  \item{boot.R}{
    number of bootstrap samples. Used only if cf was not yet
    bootstrapped using \code{bootstrap.cf}.
  }
  \item{boot.l}{
    block size for autocorrelation analysis
  }
  \item{seed}{
    seed for the random number generation used for boostrapping.
  }

}
\value{
  xxx
}
\details{
  xxx
}
%\references{
%}
\seealso{
  \code{gevp}, \code{extract.obs}, \code{bootstrap.cf}
}
%\examples{
%}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}

