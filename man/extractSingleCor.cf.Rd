\name{extractSingleCor.cf}
\alias{extractSingleCor.cf}
\title{
  extract one single correlator object as \code{cf} object from a large
  \code{cf} object.
}
\description{
  \code{cf} objects are capable of storing several correlation functions
  in form of a correlation matrix. \code{extractSingleCor.cf} lets one
  extract one of them.
}
\usage{
extractSingleCor.cf(cf, id = c(1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cf}{
    input object of class \code{cf}
}
  \item{id}{
    id of the correlators in \code{cf} to be extracted
}
}
\value{
  A \code{cf} object containing only the single correlator
\seealso{
  \code{\link{cf}}
}
\author{Carsten Urbach \email{curbach@gmx.de}}
\keyword{bootstrap}% use one of  RShowDoc("KEYWORDS")
\keyword{correlator}% __ONLY ONE__ keyword per line
