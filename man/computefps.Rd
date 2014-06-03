\name{computefps}
\alias{computefps}
\title{Computes the pseudoscalar decay constant for the twisted mass
  case from the pseudoscalar amplitude and mass}
\description{
  From a mass and amplitude determination (using \link{\code{matrixfit}}
  or \link{\code{fit.effectivemass}}, \link{\code{bootstrap.gevp}} and
  \link{\code{gevp2amplitude}} the pseudoscalar decay constant is
  determined for the case of Wilson twisted mass fermions from the
  pseudoscalar amplitude and mass
}
\usage{
computefps(mfit, PP, mass, mu1, mu2, Kappa, normalistaion="cmi")
}
\arguments{
  \item{mfit}{
    An object of type \code{matrixfit} or \code{gevp.amplitude}
    generated with \link{\code{matrixfit}} or
    \link{\code{gevp2amplitude}}, respectively.
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
  \item{reduce.latarts}{
    boolean which when set to \code{TRUE} causes the usage of a definition
    of the decay constant involving mps*sinh(mps) rather
    than mps^2 which reduces the effect of lattice artefacts
    for heavy meson masses
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
  When \code{reduce.latarts=TRUE},\cr
  \deqn{\sqrt{m_{\mathrm{PS}}^3}}{%
    sqrt(mps^3)}
  is replaced with
  \deqn{\sqrt{m_{\mathrm{PS}}} \sinh{m_{\mathrm{PS}}}}{%
    sqrt(mps)*sinh(mps))}
  which reduces lattice artefacts for heavy meson masses.
}
\seealso{
  \link{\code{matrixfit}}, \link{\code{gevp2amplitude}}, \link{\code{pion}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}
\keyword{optimise}
\keyword{ts}
