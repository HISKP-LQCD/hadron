\name{cf}
\alias{c.cf, add.cf,"+.cf", "-.cf", "/.cf", as.cf, is.cf, mul.cf}
\title{cf}
\description{
  The function \code{cf()} creates containers for correlation functions
  of class \code{cf}. This class is particularly designed to deal with
  correlation functions emerging in statistical and quantum field theory
  simulations. Arithmetic operations are defined for this class in
  several ways, as well as concatenation and \code{is.cf} and \code{as.cf}.
}
\usage{
correlator <- cf()

correlator <- c(cf1, cf2)

correlator <- add.cf(cf1, cf2, a=1., b=1.)

correlator <- mul.cf(cf1, a=1.)

correlator <- as.cf(x)

iscf <- is.cf(cf1)

cf1 + cf2
cf1 - cf2
cf1 / cf2
}
\arguments{
  \item{cf1, cf2}{
    objects of class \code{cf}
  }
  \item{a, b}{
    real (or complex) valued scalars
  }
}
\value{
  Apart from \code{is.cf} all functions return an object of class
  \code{cf}. \code{cf()} itself does not fill the container with data
  yet. For this purpose one of the reading functions
  \code{\link{readtextcf}} or \code{\link{readbinarycf}} or \code{as.cf}
  can be used. 
}
\details{
}
\seealso{
  \code{\link{readtextcf}}, \code{\link{readbinarycf}}
}
\examples{
}
\author{Carsten Urbach \email{curbach@gmx.de}}
\keyword{correlator, analysis, bootstrap}
