\name{addConfIndex2cf}
\alias{addConfIndex2cf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  add a configuration index to an \code{cf} object
}
\description{
  add a configuration number index and adds it to a \code{cf} object.
}
\usage{
addConfIndex2cf(cf, conf.index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cf}{
    and object of class \code{cf}
}
  \item{conf.index}{
    a configuration index of the same length as \code{cf}.
}
}
\value{
  Returns an object of class \code{cf} equal to the input but  with
  element \code{conf.index} added
}
\author{Carsten Urbach, \email{curbach@gmx.de}}

\seealso{
  \link{cf}
}
\keyword{correlator}
\keyword{analysis}
\keyword{bootstrap}
