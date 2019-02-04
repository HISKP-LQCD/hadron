\name{effectivemass.cf}
\alias{effectivemass.cf}
\title{Computes effective mass values for a correlation function}
\description{
  Computes effective mass values for a correlation function using
  different type of definitions for the effective mass. This function is
  mainly indented for internal usage by \code{\link{bootstrap.effectivemass}}.
}
\usage{
effectivemass.cf(cf, Thalf, type="solve", nrObs=1, replace.inf=TRUE, interval=c(0.000001,2.), weight.factor = NULL, deltat=1, tmax=Thalf)
}
\arguments{
  \item{cf}{
    The correlation function either as a vector of length \code{nrObs*(Thalf+1)}
    or as an array of dimension Nx\code{nrObs*(Thalf+1)}, where N is the number
    of observations. N will be averaged over.
  }
  \item{Thalf}{
    The half time extend of the lattice
  }
  \item{type}{
    The function to be used to compute the effective mass
    values. Possibilities are "acosh", "solve", "log", "temporal",
    "shifted" and "weighted". While the first three assume normal
    cosh behaviour of 
    the correlation function, "temporal" is desigend to remove an
    additional constant stemming from temporal states in two particle
    correlation functions. The same for "subtracted" and "weighted", the
    latter for the case of two particle energies with the two particle
    having different energies. In the latter case only the leading
    polution is removed by \code{removeTemporal.cf} and taken into
    account here.
  }
  \item{nrObs}{
    The number of "observables" included in the correlator
  }
  \item{replace.inf}{
    If set to \code{TRUE}, all \code{\link{Inf}} values will be replaced
    by \code{\link{NA}}. This is needed for instance for
    \code{\link{bootstrap.effectivemass}}. 
  }
  \item{interval}{
    initial interval for the \code{\link{uniroot}} function when
    numerically solving for the effective mass.
  }
  \item{weight.factor}{
    relative weight for type "weighted" only, see details
  }
  \item{deltat}{
    time shift for shifted correlation functions
  }
  \item{tmax}{
    t-value up to which the effectivemass is to be computed
  }
}
\value{
  Returns a vector of length \code{Thalf} with the effective mass values
  for t-values running from 0 to \code{Thalf-1} 
}
\details{
  A number of types is implemented to compute effective mass values from
  the correlation function:

  "solve": the ratio\cr
  \eqn{C(t+1) / C(t) = \cosh(-m*(t+1)) / \cosh(-m*t)}\cr
  is numerically solved for \eqn{m(t)}.

  "acosh": the effective mass is computed from\cr
  \eqn{m(t)=acosh((C(t-1)+C(t+1)) / (2C(t)))}\cr
  Note that this definition is less tolerant against noise.

  "log": the effective mass is defined via\cr
  \eqn{m(t)=\log(C(t) / C(t+1))}\cr
  which has artifacts of the periodicity at large t-values.

  "temporal": the ratio\cr
  \eqn{[C(t)-C(t+1)] / [C(t-1)-C(t)] = [\cosh(-m*(t))-\cosh(-m*(t+1))] / [\cosh(-m*(t-1))-\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)}.

  "subtracted": like "temporal", but the differences \eqn{C(t)-C(t+1)}
  are assumed to be taken already at the correlator matrix level using
  \code{removeTemporal.cf} and hence the ratio\cr
  \eqn{[C(t+1)] / [C(t)] = [\cosh(-m*(t))-\cosh(-m*(t+1))] /
    [\cosh(-m*(t-1))-\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)}.

  "weighted": like "subtracted", but now there is an additional weight
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
  \code{\link{bootstrap.effectivemass}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
