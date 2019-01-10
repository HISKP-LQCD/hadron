\name{bootstrap.effectivemass}
\alias{bootstrap.effectivemass}
\title{Computes effective masses with bootstrapping errors}
\description{
  Generates bootstrap samples for effective mass values computed from an
  object of class \code{cf} (a correlation function)
}
\usage{
bootstrap.effectivemass(cf, type="solve")
}
\arguments{
  \item{cf}{
    a correlation function as an object of type \code{cf}, preferably
    after a call to \code{\link{bootstrap.cf}}. If the latter has not
    been called yet, it will be called in this function.
  }
  \item{type}{
    The function to be used to compute the effective mass
    values. Possibilities are "acosh", "solve", "log", "temporal",
    "shifted" and "weighted". While the first three assume normal
    cosh behaviour of 
    the correlation function, "temporal" is desigend to remove an
    additional constant stemming from temporal states in two particle
    correlation functions. The same for "shifted" and "weighted", the
    latter for the case of two particle energies with the two particle
    having different energies. In the latter case only the leading
    polution is removed by \code{removeTemporal.cf} and taken into
    account here.
  }
}
\value{
  An object of class \code{effectivemass} is invisibly returned. It has
  objects:
  \code{effMass}:\cr
  The computed effective mass values as a vector of length
  \code{Time/2}. For \code{type="acosh"} also the first value is
  \code{NA}, because this definition requires three time slices.
  
  \code{deffMass}:\cr
  The computed bootstrap errors for the effective masses of the same
  length as \code{effMass}.
  
  \code{effMass.tsboot}:\cr
  The boostrap samples of the effective masses as an array of dimension
  RxN, where \code{R=boot.R} is the number of bootstrap samples and
  \code{N=(Time/2+1)}.
  
  and \code{boot.R},
  \code{boot.l}, \code{Time}
}
\details{
  A number of types is implemented to compute effective mass values from
  the correlation function:

  "solve": the ratio\cr
  \eqn{C(t+1) / C(t) = \cosh(-m*(t+1)) / \cosh(-m*t)}\cr
  is numerically solved for m.

  "acosh": the effective mass is computed from\cr
  \eqn{m=acosh((C(t-1)+C(t+1)) / (2C(t)))}\cr
  Note that this definition is less tolerant against noise.

  "log": the effective mass is defined via\cr
  \eqn{m=\log(C(t) / C(t+1))}\cr
  which has artifacts of the periodicity at large t-values.

  "temporal": the ratio\cr
  \eqn{[C(t)-C(t+1)] / [C(t-1)-C(t)] = [\cosh(-m*(t))-\cosh(-m*(t+1))] / [\cosh(-m*(t-1))-\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)}.

  "shifted": like "temporal", but the differences \eqn{C(t)-C(t+1)}
  are assumed to be taken already at the correlator matrix level using
  \code{removeTemporal.cf} and hence the ratio\cr
  \eqn{[C(t+1)] / [C(t)] = [\cosh(-m*(t))-\cosh(-m*(t+1))] /
    [\cosh(-m*(t-1))-\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)}.

  "weighted": like "shifted", but now there is an additional weight
  factor \eqn{w} from \code{removeTemporal.cf} to be taken into account,
  such that the ratio\cr
  \eqn{[C(t+1)] / [C(t)] = [\cosh(-m*(t))-w*\cosh(-m*(t+1))] /
    [\cosh(-m*(t-1))-w*\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)} with \eqn{w} as input.
  
}
\references{
  arXiv:1203.6041
}
\seealso{
  \code{\link{fit.effectivemass}}, \code{\link{bootstrap.cf}},
  \code{removeTemporal.cf}
}
\examples{
data(samplecf)
samplecf <- bootstrap.cf(cf=samplecf, boot.R=1500, boot.l=2, seed=1442556)
effmass <- bootstrap.effectivemass(cf=samplecf)
summary(effmass)
plot(effmass, ylim=c(0.14,0.15))
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
