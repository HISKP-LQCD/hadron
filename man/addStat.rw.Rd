\name{addStat.rw}
\alias{addStat.rw}
\title{Combine statistics of two rw objects}
\description{
  \code{addStat.rw} takes the raw data of two \code{rw} objects and
  combines them into one
}
\usage{
addStat.cf(rw1, rw2, reverse1, reverse2)
}
\arguments{
  \item{rw1}{
    the first of the two \code{rw} objects to be combined
  }
  \item{rw2}{
    the second of the two \code{rw} objects to be combined
  }
}
\value{
  an object of class \code{rw} with the statistics of the two input
  \code{rw} objects combined
}
\details{
  Note that the two \code{rw} objects to be combined need to be
  compatible. Otherwise, \code{addStat.rw} will abort with an error. 
}
\seealso{
  \code{\link{rw}}
}
\examples{
\dontrun{rw0 <- addStat(rw1=rw1, rw2=rw2)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{correlation function}
