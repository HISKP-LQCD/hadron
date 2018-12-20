\name{bootstrap.cf}
\alias{bootstrap.cf}
\title{bootstrap a set of correlation functions}
\description{
  bootstrap a set of correlation functions
}
\usage{
bootstrap.cf(cf, boot.R=400, boot.l=2, seed=1234, sim="geom", endcorr=TRUE)
}
\arguments{
  \item{cf}{
    correlation matrix of class \code{cf} e.g. obtained with a call to \code{extrac.obs}.
  }
  \item{boot.R}{
    number of bootstrap samples.
  }
  \item{boot.l}{
    block size for autocorrelation analysis
  }
  \item{seed}{
    seed for the random number generation used for boostrapping.
  }
  \item{sim}{
    The type of simulation required to generate the replicate
    time series.  The possible input values are ‘"fixed"’ (block
    resampling with fixed block lengths of ‘boot.l’) and ‘"geom"’ (block
    resampling with block lengths having a geometric distribution with
    mean ‘boot.l’). Default is ‘"geom"’. See \link{boot::tsboot} for details.
  }
  \item{endcorr}{
    A logical variable indicating whether end corrections are to
    be applied when ‘sim’ is ‘"fixed"’.  When ‘sim’ is
    ‘"geom"’, ‘endcorr’ is automatically set to ‘TRUE’; ‘endcorr’ is
    not used when ‘sim’ is ‘"model"’ or ‘"scramble"’. See
    \link{boot::tsboot} for details.
  }
}
\value{
  returns an object of class \code{cf} with bootstrap samples added for
  th correlation function called \code{cf.tsboot}. Moreover, the
  original average of \code{cf} is returned as \code{cf0} and the
  bootstrap errors as \code{tsboot.se}. We also copy the input
  parameters over and set \code{bootstrap.samples} to \code{TRUE}.
}
\examples{
data(samplecf)
samplecf <- bootstrap.cf(cf=samplecf, boot.R=1500, boot.l=2, seed=1442556)
plot(samplecf, log=c("y"))
}
\seealso{
  \code{\link[boot:tsboot]{tsboot}}, \code{jackknife.cf}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{bootstrap}
