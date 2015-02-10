\name{computefps}
\alias{computefps}
\title{Computes the pseudoscalar decay constant for the twisted mass
  case from the pseudoscalar amplitude and mass}
\description{
  From a mass and amplitude determination (using \code{\link{matrixfit}}
  or \code{\link{fit.effectivemass}}, \code{\link{bootstrap.gevp}} and
  \code{\link{gevp2amplitude}} the pseudoscalar decay constant is
  determined for the case of Wilson twisted mass fermions from the
  pseudoscalar amplitude and mass
}
\usage{
computefps(mfit, PP, mass, mu1, mu2, Kappa, normalistaion="cmi")
}
\arguments{
  \item{mfit}{
    An object of type \code{matrixfit} or \code{gevp.amplitude}
    generated with \code{\link{matrixfit}} or
    \code{\link{gevp2amplitude}}, respectively.
  }
  \item{PP}{
    If \code{mfit} is missing this must contain the value for the
    pseudoscalar amplitude.
  }
  \item{mass}{
    If \code{mfit} is missing this must contain the value for the
    pseudoscalar mass.
  }
  \item{mu1, mu2}{
    The values for the twisted quark masses involved in the pseudoscalar
    meson. If \code{mu2} is missing it will be assumed to be equal to
    \code{mu1}. 
  }
  \item{Kappa}{
    The \eqn{\kappa}{kappa}-value of the run, needed only if \code{normalisation="cmi"}.
  }
  \item{normalisation}{
    normalisation of the correlators. If set to "cmi" the
    \eqn{\kappa}{kappa} value must be specified.
  }
  \item{disprel}{
    One of "continuum" or "lattice". Indicates whether the formula
    for the decay constant should take into account the lattice
    dispersion relation for the meson. Theoretically this can reduce
    lattice artefacts for heavy mesons.
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
  If \code{mfit} ist missing the value of fps will printed to stdout and
  returned as a simple numerical value.

  If \code{mfit} is available, this object will be returned but with
  additional objects added: \code{fps}, \code{fps.tsboot},
  \code{mu1,mu2}, \code{normalistaion} and \code{Kappa} if applicable. 
}
\details{
  The pseudoscalar decay constant is computed from\cr
  \deqn{f_\mathrm{PS} = 2\kappa(\mu_1+\mu_2)\frac{PP}{\sqrt{2}\sqrt{m_\mathrm{PS}}^3}}{%
    fps = 2 kappa (mu1+mu2) PP/sqrt(2)/sqrt(mps)^3}
  for \code{normalisation="cmi"} or
  \deqn{f_\mathrm{PS} = (\mu_1+\mu_2)\frac{PP}{\sqrt{2}\sqrt{m_\mathrm{PS}}^3}}{%
    fps = (mu1+mu2) PP/sqrt(2)/sqrt(mps)^3}
  expecting physical normalisation of the amplitudes.\cr
  When \code{disprel="lattice"},\cr
  \deqn{\sqrt{m_{\mathrm{PS}}^3}}{%
    sqrt(mps^3)}
  is replaced with
  \deqn{\sqrt{m_{\mathrm{PS}}} \sinh{m_{\mathrm{PS}}}}{%
    sqrt(mps)*sinh(mps))}
  which can reduce lattice artefacts for heavy meson masses.
}
\seealso{
  \code{\link{matrixfit}}, \code{\link{gevp2amplitude}}, \code{\link{pion}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}
\keyword{optimise}
\keyword{ts}
