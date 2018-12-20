\name{plot.outputdata}
\alias{plot.outputdata}
\title{Plot Command For Class Ouputdata}
\description{
  Generic plot routine for class \dQuote{ouputdata}. Currently it plots
  the plaquette history and the history of \eqn{\Delta H}{Delta H}
}
\usage{
\method{plot}{outputdata}(data, skip=0, ...)
}
\arguments{
  \item{data}{
    object of class \dQuote{outputdata} obtained from a read with \code{readoutputdata}
  }
  \item{skip}{
    number of trajectories to be skipped in analysis for plaquette and
    \eqn{\exp(-\Delta H)}{exp(-Delta H)}.
  }
  \item{...}{
    additional arguments passed to the generic plot function.
  }
}
\value{
  list containing the \dQuote{data}, an object of class \dQuote{uwerr}
  called \dQuote{plaq.res} containing the statisical analysis for the
  plaquette and a second object of type \dQuote{uwerr} called
  \dQuote{dH.res} with  the statisical analysis for \eqn{\exp(-\Delta
    H)}{exp(-Delta H)}. 
}
%\details{
%}
%\references{
%}
\seealso{
  \code{\link{readoutputdata}}, \code{\link{uwerr}}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{methods}
\keyword{hplot}
\examples{
\dontrun{plaq <- readoutputdata("output.data")}
\dontrun{plaq.plot <- plot(plaq, skip=100)}
\dontrun{summary(plaq.plot$plaq.res)}
}
