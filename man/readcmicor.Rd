\name{readcmicor}
\alias{readcmicor}
\title{Read Data In Chris Michael Format}
\description{
  reads data from file in Chris Michael format
}
\usage{
  readcmicor(filename)
}
\arguments{

  \item{filename}{
    filename of the data file
  }
}
\value{
  returns a data frame containing the data
}
%\details{
%}
%\references{
%}
%\seealso{
%}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
}
\author{Carsten Urbach \email{carsten.urbach@liverpool.ac.uk}}
\keyword{file}
