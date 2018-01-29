\name{tex.catwitherror}
\alias{tex.catwitherror}
\title{cat a number with error in tex-ready format}
\description{
  A number with error is cat to stdoutput in tex-ready format like
  xx(yy) thereby automatically determining the digit at which the error
  applies.
}
\usage{
tex.catwitherror(x, dx, digits=1)
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
}
\value{
  writes a string to standard output
}
\details{
  The value \code{x} is properly rounded to its significant digits determined
  by the values of \code{dx} and \code{digits} and then the appropriate
  output is cat to standard output.
}
%\references{
%}
%\seealso{
%}
\examples{
tex.catwitherror(x=0.375567, dx=0.001, digits=1)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
