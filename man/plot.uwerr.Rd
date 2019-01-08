\name{plot.uwerr}
\alias{plot.uwerr}
\title{Plot Command For Class UWerr}
\description{
  Plot Command For Class UWerr
}
\usage{
\method{plot}{uwerr}(x, ..., main="x", x11=TRUE, plot.hist=TRUE, index=1, Lambda=100)
}
\arguments{
  \item{x}{
    object of class \code{uwerr}
  }
  \item{...}{
    generic parameters, not used here.
  }
  \item{main}{
    main title of the plot
  }
  \item{x11}{
    whether or not to generate X11 plots
  }
  \item{plot.hist}{
    whether or not to generate a histogram
  }
}
\value{
  produces various plots, including a histogram, the
  autocorrelationfunction and the integrated autocorrelation time, all
  with error bars.
}
\seealso{
  \code{\link{uwerr}}, \code{\link{pion}}
}
\examples{
data(plaq.sample)
plaq.res <- uwerrprimary(plaq.sample)
plot(plaq.res)
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{methods}
\keyword{hplot}

