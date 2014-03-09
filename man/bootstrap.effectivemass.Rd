\name{bootstrap.effectivemass}
\alias{bootstrap.effectivemass}
\title{Bootstraped effective mass computation}
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
  \code{effMass}, \code{deffMass}, \code{effMass.tsboot}, \code{boot.R},
  \code{boot.l}, \code{Time}
}
\details{
  Three functions are implemented to compute effective mass values from
  the correlation function:

  "solve": the ratio
  \eqn{C(t+1) / C(t) = \cosh(-m*(t+1)) / \cosh(-m*t)}
  is numerically solved for m.

  "acosh": the effective mass is computed from
  \eqn{m=\acosh((C(t-1)+C(t+1)) / (2C(t)))}
  Note that this definition is less tolerant against noise.

  "log": the effective mass is defined via
  \eqn{m=\log(C(t) / C(t+1))}
  which has artifacts of the periodicity at large t-values.
}
\seealso{
  \link{\code{fit.effectivemass}}, \link{\code{bootstrap.cf}}
}
#\examples{
#}
\author{Carsten Urbach, \email{curbach@gmx.de}}