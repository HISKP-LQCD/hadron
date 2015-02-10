\name{gevp2amplitude}
\alias{gevp2amplitude}
\title{Extracts physical amplitudes from a GEVP}
\description{
  Given a GEVP generated with \code{bootstrap.gevp} and masses
  determined from the principle correlator with given \code{id}, the
  physical amplitudes are extracted and bootstraped. The man amplitude
  is determined from a constant fit to the data in the specified time
  range. 
}
\usage{
gevp2amplitude(gevp, mass, id=1, op.id=1, type="cosh", t1, t2, useCov)
}
\arguments{
  \item{gevp}{
    An object of class \code{gevp} as generated with a call to \code{bootstrap.gevp}.
  }
  \item{mass}{
    Optimally, this is an object either of class \code{effectivemassfit}
    generated using \code{\link{fit.effectivemass}} or of class
    \code{matrixfit} generated with \code{\link{matrixfit}} to the
    principal correlator extracted using \code{\link{gevp2cf}} applied
    to \code{gevp} using the same value of \code{id}.

    It can also be given as a numerical vector with the bootstrap
    samples as entries. The mean will then be computed as the bootstrap
    mean over this vector. The number of samples must agree with the
    number of bootstrap samples in \code{gevp}. 
  }
  \item{id}{
    The index of the principal correlator to extract, i.e. the physical
    state to extract.
  }
  \item{op.id}{
    The index of the operator for which to extract the amplitude.
  }
  \item{type}{
    The symmetry of the pricipal correlator in time, can be either
    "cosh" or "sinh".
  }
  \item{t1,t2}{
    The time range in which to fit the amplitude starting with 0. If not
    given it will be tried to infer these from the \code{mass} object.
  }
  \item{useCov}{
    Use the covariance matrix for fitting the constant to the amplitude
    data. 
  }
}
\value{
  
}
\details{
}
\seealso{
  \code{\link{matrixfit}}, \code{\link{fit.effectivemass}},
  \code{\link{gevp}}, \code{\link{gevp2cf}}, \code{\link{computefps}}
}
\examples{
\dontrun{## apply a GEVP analysis}
\dontrun{pion.cor.gevp <- bootstrap.gevp(pion.cor, t0=1)}
\dontrun{## extract the first principal correlator}
\dontrun{pion.pc1 <- gevp2cf(pion.cor.gevp, id=1)}
\dontrun{pion.pc1.effectivemass <- bootstrap.effectivemass(cf=pion.pc1, type="acosh")}
\dontrun{pion.pc1.effectivemass <- fit.effectivemass(pion.pc1.effectivemass, t1=12, t2=23, useCov=TRUE)}
\dontrun{## now determine the amplitude}
\dontrun{pion.pc1.amplitude <- gevp2amplitude(pion.cor.gevp, pion.pc1.effectivemass)}
\dontrun{## compute also the pion decay constant}
\dontrun{pion.pc1.amplitude <- computefps(pion.pc1.amplitude, Kappa=0.125, mu1=0.003)}
\dontrun{summary(pion.pc1.amplitude)}
\dontrun{plot(pion.pc1.amplitude)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{GEVP}
