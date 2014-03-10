\name{bootstrap.effectivemass}
\alias{bootstrap.effectivemass}
\title{Computes effective masses with bootstrapping errors}
\description{
  Generates bootstrap samples for effective mass values computed from an
  object of class \code{cf} (a correlation function)
}
\usage{
boostrap.effectivemass(cf, boot.R=400, boot.l=20, seed=123456, type="solve")
}
\arguments{
  \item{cf}{
    a correlation function as an object of type \code{cf}, preferably
    after a call to \link{\code{bootstrap.cf}}. If the latter has not
    been called yet, it will be called in this function.
  }
  \item{boot.R}{
    only used if no preceeding call to \link{\code{bootstrap.cf}} has
    happened and it counts the number of bootstrap samples.
  }
  \item{boot.l}{
    only used if no preceeding call to \link{\code{bootstrap.cf}} has
    happened and it counts the bootstrap block length.
  }
  \item{seed}{
    the seed used in the bootstrap sample generation (for reproducability).
  }
  \item{type}{
    The function to be used to compute the effective mass
    values. Possibilities are "acosh", "solve" and "log".
  }
}
\value{
  An object of class \code{effectivemass} is invisibly returned. It has
  objects:
  \code{effMass}:\cr
  The computed effective mass values as a vector of length
  \code{Time/2+1}. Note that the value at \code{Time/2+1} is always
  \code{NA}, because the effective mass needs at least to time slices to
  be defined. Moreover, for \code{type="acosh"} also the first value is
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
  Three functions are implemented to compute effective mass values from
  the correlation function:

  "solve": the ratio\cr
  \eqn{C(t+1) / C(t) = \cosh(-m*(t+1)) / \cosh(-m*t)}\cr
  is numerically solved for m.

  "acosh": the effective mass is computed from\cr
  \eqn{m=\acosh((C(t-1)+C(t+1)) / (2C(t)))}\cr
  Note that this definition is less tolerant against noise.

  "log": the effective mass is defined via\cr
  \eqn{m=\log(C(t) / C(t+1))}\cr
  which has artifacts of the periodicity at large t-values.
}
\seealso{
  \link{\code{fit.effectivemass}}, \link{\code{bootstrap.cf}}
}
\examples{
\dontrun{pion.pc1.effectivemass <- bootstrap.effectivemass(cf=pion.pc1, type="acosh")}
\dontrun{pion.pc1.effectivemass <- fit.effectivemass(pion.pc1.effectivemass, t1=12, t2=23, useCov=FALSE)}
\dontrun{summary(pion.pc1.effectivemass)}
\dontrun{plot(pion.pc1.effectivemass)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}