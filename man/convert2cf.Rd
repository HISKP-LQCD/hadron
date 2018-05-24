\name{convert2cf}
\alias{convert2cf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  convert data into a \code{cf} object
}
\description{
  Lets one convert data of certain type into an object of class \code{cf}
}
\usage{
convert2cf(data, symmetric = TRUE, symmetrise = TRUE)
}
\arguments{
  \item{data}{
    a \code{data.frame}
}
  \item{symmetric}{
    is the data symmetric or anti-symmetric in t
}
  \item{symmetrise}{
    symmetrise the data
}
}
\value{
  An object of class \code{cf}
}
\author{Carsten Urbach \email{curbach@gmx.de}}
\seealso{
  \code{\link{cf}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlator }% use one of  RShowDoc("KEYWORDS")
