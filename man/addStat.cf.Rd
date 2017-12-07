\name{addStat.cf}
\alias{addStat.cf}
\title{Combine statistics of two cf objects}
\description{
  \code{addStat.cf} takes the raw data of two \code{cf} objects and
  combines them into one
}
\usage{
addStat.cf(cf1=cf1, cf2=cf2)
}
\arguments{
  \item{cf1}{
    the first of the two \code{cf} objects to be combined
  }
  \item{cf2}{
    the second of the two \code{cf} objects to be combined
  }
}
\value{
  an object of class \code{cf} with the statistics of the two input
  \code{cf} objects combined
}
\details{
  Note that the two \code{cf} objects to be combined need to be
  compatible. Otherwise, \code{addStat.cf} will abort with an error. 
}
\seealso{
  \code{\link{cf}}
}
\examples{
\dontrun{cf0 <- addStat(cf1=cf1, cf2=cf2)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{correlation function}
