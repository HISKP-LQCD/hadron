\name{cfunction}
\alias{cfunction}
\title{fit one correlation function}
\description{
  fit one correlation function to a cosh or sinh and
  determine the mass
}
\usage{
cfunction(data, t1, t2, S=1.5, pl=FALSE, skip=0, cformat="cmi",
          itype=1, iobs=1, sign=+1., ind.vec=c(1,3,4,5),
          boot.R=99, boot.l=10, tsboot.sim="geom",
          method="uwerr", fit.routine="optim", mass.guess, par.guess, nrep)
}
\arguments{
  \item{data}{
    data to be fitted to as e.g. the output of
    \code{\link{readcmicor}}. Currently only
    \code{\link{cmicor}} format is supported.
  }
  \item{cformat}{
    the format of the correlation function data. Currently only
    \code{\link{cmicor}} format is supported.
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{S}{
    passed to \code{uwerr}, see documentation of \code{\link{uwerr}}.
  }
  \item{pl}{
    logical: if set to TRUE the function produces plots
  }
  \item{skip}{
    number of measurements to be discarded at the beginning of the time
    series. \code{skip} has no effect if two or more replica are used,
    see argument \code{nrep}.
  }
  \item{itype}{
    if there are several gamma matrix combinations in the data, then
    this enumerates these combinations (for cmi format)
  }
  \item{iobs}{
    if there are several operators available (local-local,
    local-smeared, etc.), then this labels these (for cmi format)
  }
  \item{sign}{
    for sign=+1 fit a cosh, for sign=-1. fit a sinh.
  }
  \item{ind.vec}{
    index vector indexing the column numbers in cmicor to be used
  }
  \item{boot.R}{
    number of bootstrap samples for bootstrap analysis
  }
  \item{boot.l}{
    average block size for blocking analysis with tsboot
  }
  \item{tsboot.sim}{
    The type of simulation required to generate the replicate
    time series. See \code{\link{tsboot}} for details.
  }
  \item{method}{
    the type of error analysis to be used. Can be either \dQuote{uwerr},
    \dQuote{boot}, \dQuote{all} or \dQuote{no}. For \dQuote{no} (or any
    other string) no error analysis is performed. This might be helpful
    for a first impression and also to test different initial values for
    the fitting parameters. The latter is in particular needed for more
    than one state in the fit.
  }
  \item{mass.guess}{
    numerical vector of mass-values to be used as initial values in the
    fit. If given, it must have at least \code{no.masses} entries.
  }
  \item{par.guess}{
    numerical vector of couling parameter values to be used as initial
    values in the fit. The order is \eqn{P_L}, \eqn{P_F}, \eqn{A_L},
    \eqn{A_F}, \eqn{4_L}, \eqn{4_F}, for notation see below. If given,
    it must have at least \code{no.masses} times \code{matrix.size}
    entries.
  }
  \item{nrep}{
    vector (N1, N2, ...) of replica length N1, N2. If missing it is
    assumed that there is only one ensemble. If there are two or more
    replica the parameter \code{skip} has no effect.
  }
  \item{fit.routine}{
    The fit routine to be used. Default is \dQuote{gsl}, which uses the
    gnu scientific library \dQuote{gsl_multifit_fdfsolver} solver to
    minimise the chisquare. All other values lead to the usage of R's
    \link{optim} function. The latter choice might be significantly
    slower.
  }  
}
\value{
  returns an object of \code{class} \code{cfit} with the following items

  \item{fitresult}{
    result from the fit as returned by \code{\link{optim}}
  }
  \item{t1}{
    lower bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{t2}{
    upper bound for the fitrange in time (t1,t2). Counting starts with 0.
  }
  \item{N}{
    number of measurements found in the data
  }
  \item{Time}{
    Time extent found in the data
  }
  \item{fitdata}{
    \code{\link{data.frame}} containing the time values used in the fit,
    the averaged correlator and its error and the value of Chi for each
    time value
  }
  \item{uwerrresultmps}{
    the result of the time series analysis for the lowest mass as
    carried out by \code{\link{uwerr}} 
  }
  \item{effmass}{
    effective masses
  }
  \item{matrix.size}{
    size of the data matrix, copied from input
  }
  \item{boot}{
    object returned by the call to \code{\link{boot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{tsboot}{
    object returned by the call to \code{\link{tsboot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{method}{
    error analysis method as copied from input
  }
  \item{fit.routine}{
    \code{fit.routine} as copied from input
  }
  \item{nrep}{
    \code{nrep} as copied from input
  }
}
\details{
  see e.g. \code{\link{pion}} for details.
}
\references{
}
\seealso{
  \code{\link{rho}}, \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{plot.cfit}}, \code{\link{summary.pionfit}},
  \code{\link{pp}}, \code{\link{variational}}, \code{\link{pion}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{fit <- cfunction(cmicor, t1=10, t2=23)}
\dontrun{summary(fit)}
\dontrun{plot(fit)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimize}
\keyword{ts}
