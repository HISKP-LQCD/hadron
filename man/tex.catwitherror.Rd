\name{tex.catwitherror}
\alias{tex.catwitherror}
\title{paste a number with error in tex-ready format}
\description{
  A number with error is converted to a string in tex-ready format like
  xx(yy) thereby automatically determining the digit at which the error
  applies.
}
\usage{
tex.catwitherror(x, dx, digits=1, with.dollar=TRUE, human.readable=TRUE)
}
\arguments{
  \item{x}{
    either a single numeric value, or a numeric vector, where the first
    element is the value and the second is its error
  }
  \item{dx}{
    the error. If supplied, it will be printed as the error and the
    value is the first element of \code{x}. If \code{dx} is missing, the
    second element of \code{x}, if available, is used as the error. If
    \code{dx} is missing and the length of \code{x} is one, only the
    value is converted to a string without error.
  }
  \item{digits}{
    number of error digits
  }
  \item{with.dollar}{
    include the tex dollar in the return string or not
  }
  \item{human.readable}{
    controls whether or not scientific format will be produced in a
    human readable form, or not at all. The latter might be useful for
    output with summary or print.
  }
}
\value{
  writes a string to standard output
}
\details{
  The value of the first element of \code{x} is properly rounded to its
  significant digits determined by the values of \code{dx} or the second
  element of \code{x} (see above) and \code{digits}. Then a tex-ready
  string is returned.
}
%\references{
%}
%\seealso{
%}
\examples{
tex.catwitherror(x=0.375567, dx=0.001)
tex.catwitherror(x=c(0.375567, 0.001))
## it can be used with apply
x = array(c(0.1187, 0.291, 0.388, 0.011, 0.037, 0.021), dim=c(3,2))
apply(x, 1, tex.catwitherror, digits=2)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
