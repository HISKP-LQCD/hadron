\name{computefpsOS}
\alias{computefpsOS}
\title{Computes the pseudoscalar decay constant for the Osterwalder Seiler
  case from the pseudoscalar amplitude and mass}
\description{
  From a mass and amplitude determination (using \code{\link{matrixfit}})
  the pseudoscalar decay constant is
  determined for the case of Osterwalder Seiler (OS) fermions
  from the AS and SS amplitude (in the twisted basis), ZA and the OS pion mass.
}
\usage{
computefps(mfit, kappa, normalistaion="cmi", ZA=1, ZAboot, dZA, boot.fit=TRUE)
}
\arguments{
  \item{mfit}{
    An object of type \code{matrixfit} 
    generated with \code{\link{matrixfit}}. The correlation matrix (SS,
    SA, AS, AA) must have been analysed, where the correlators are in
    the twisted basis. 
  }
  \item{ZA}{
    The value of the renormalisation constant \eqn{Z_A}{ZA}.
  }
  \item{dZA}{
    The value of the (normally distributed) error of the renormalisation
    constant \eqn{Z_A}{ZA}.
  }
  \item{ZAboot}{
    Bootstrap samples for \eqn{Z_A}{ZA}. If they are provided, they are
    used for computing fps, if not, bootstrap samples are generated from
    \code{dZA}. If both are missing, the error of \eqn{Z_A}{ZA} is not
    taken into account. 
  }
  \item{Kappa}{
    The \eqn{\kappa}{kappa}-value of the run, needed only if \code{normalisation="cmi"}.
  }
  \item{normalisation}{
    normalisation of the correlators. If set to "cmi" the
    \eqn{\kappa}{kappa} value must be specified.
  }
  \item{boot.fit}{
    If set to \code{FALSE}, the computation is not bootstrapped, even if the 
    \code{matrixfit} or \code{gevp.amplitude} contain bootstrap samples.
    This is a useful time-saver if error information is not strictly necessary.                                                                                                                                 
    Of course, this affects the return values related to the bootstrap, which are
    set to \code{NA}.    
  }
}
\value{
  If \code{mfit} is available, this object will be returned but with
  additional objects added: \code{fpsOS}, \code{fpsOS.tsboot},
  \code{normalistaion}, \code{ZA}, \code{ZAboot} and \code{kappa} if applicable. 
}
\details{
  The pseudoscalar decay constant is computed from\cr
  \deqn{f_\mathrm{PS}^\mathrm{OS} = Z_A \sqrt{2}\kappa\frac{\langle 0| A|\pi\rangle}{m_\mathrm{PS}}}{%
    fpsOS = sqrt(2) kappa ZA <0|A|pi>/mps
  }
  for \code{normalisation="cmi"} or
  \deqn{f_\mathrm{PS}^\mathrm{OS} = Z_A \frac{\langle 0| A|\pi\rangle}{m_\mathrm{PS}}}{%
    fpsOS = ZA <0|A|pi>/mps}
  expecting physical normalisation of the amplitudes.\cr
}
\seealso{
  \code{\link{matrixfit}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimise}
\keyword{ts}
