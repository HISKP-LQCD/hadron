\name{tex.catwitherror}
\alias{tex.catwitherror}
\title{paste a number with error in tex-ready format}
\description{
  A number with error is converted to a string in tex-ready format like
  xx(yy) thereby automatically determining the digit at which the error
  applies.
}
\usage{
tex.catwitherror(x, dx, digits=1, with.dollar=TRUE)
}
\arguments{
  \item{x}{
    the values
  }
  \item{dx}{
    the error
  }
  \item{digits}{
    number of error digits
  }
  \item{with.dollar}{
    include the tex dollar in the return string or not
  }
}
\value{
  writes a string to standard output
}
\details{
  The value \code{x} is properly rounded to its significant digits determined
  by the values of \code{dx} and \code{digits} and then a tex-ready
  string is returned.
}
%\references{
%}
%\seealso{
%}
\examples{
tex.catwitherror(x=0.375567, dx=0.001)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
