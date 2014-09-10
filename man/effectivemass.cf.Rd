\name{effectivemass.cf}
\alias{effectivemass.cf}
\title{Computes effective mass values for a correlation function}
\description{
  Computes effective mass values for a correlation function using
  different type of definitions for the effective mass. This function is
  mainly indented for internal usage by \link{\code{bootstrap.effectivmass}}.
}
\usage{
effectivemass.cf(cf, Thalf, type="solve", nrObs=1, replace.inf=TRUE)
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
    values. Possibilities are "acosh", "solve", "log" and
    "temporal". While the first three assume normal cosh behaviour of
    the correlation function, "temporal" is desigend to remove an
    additional constant stemming from temporal states in two particle
    correlation functions. 
  }
  \item{nrObs}{
    The number of "observables" included in the correlator
  }
  \item{replace.inf}{
    If set to \code{TRUE}, all \link{\code{Inf}} values will be replaced
    by \link{\code{NA}}. This is needed for instance for
    \link{\code{bootstrap.effectivemass}}. 
  }
}
\value{
  Returns a vector of length \code{Thalf} with the effective mass values
  for t-values running from 0 to \code{Thalf-1} 
}
\details{
  Three functions are implemented to compute effective mass values from
  the correlation function:

  "solve": the ratio\cr
  \eqn{C(t+1) / C(t) = \cosh(-m*(t+1)) / \cosh(-m*t)}\cr
  is numerically solved for \eqn{m(t)}.

  "acosh": the effective mass is computed from\cr
  \eqn{m(t)=\acosh((C(t-1)+C(t+1)) / (2C(t)))}\cr
  Note that this definition is less tolerant against noise.

  "log": the effective mass is defined via\cr
  \eqn{m(t)=\log(C(t) / C(t+1))}\cr
  which has artifacts of the periodicity at large t-values.

  "temporal": the ratio\cr
  \eqn{[C(t)-C(t+1)] / [C(t-1)-C(t)] = [\cosh(-m*(t))-\cosh(-m*(t+1))] / [\cosh(-m*(t-1))-\cosh(-m(t))]}\cr
  is numerically solved for \eqn{m(t)}.
}
\seealso{
  \link{\code{bootstrap.effectivemass}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
