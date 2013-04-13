\name{readoutputdata}
\alias{readoutputdata}
\title{Read Data In output.data Format of tmLQCD}
\description{
  reads data from an output.data file written by tmLQCD
}
\usage{
  readouptutdata(filename)
}
\arguments{

  \item{filename}{
    filename of the data file
  }
}
\value{
  returns a data frame of class \dQuote{outputdata} containing the
data.
}
\details{
  The data can be plotted directly using \dQuote{plot}.
}
%\references{
%}
%\seealso{
%}
\examples{
library(hadron)
\dontrun{plaq <- readcmicor("output.data")}
\dontrun{plot(plaq)}
}
\author{Carsten Urbach \email{curbach@gmx.de}}
\keyword{file}
