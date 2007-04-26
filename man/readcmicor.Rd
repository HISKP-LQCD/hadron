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
\details{
  This function read data from a file. The data in the file is assumed
  to have 6 collumns: 1) the observable type number; 2) the operator
  type number; 3) the time difference from source going from 0 to \eqn{T/2}
  for each operator type; 4) \eqn{c_1}{c1} correlator value at time value forward in
  time; 5) \eqn{c_2}{c2} correlator value at time value backward in time; 6) number of
  gauge configuration..

  There are scripts shipped with the package converting the output
  written into seperate files for each gauge configuration into the
  expected format. They are called \code{puttogether.sh} and
  \code{puttogether_reverse.sh} which will sort with increasing and
  with decreasing gauge configuration number, respectively.

  Normalisation of c1, c2 needs mult by factor of \eqn{0.5} (and possible
  \eqn{(2*\kappa)^2}{(2*k)^2} 
  and \eqn{L^3} factors dependent on your conventions).
}
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
