\name{pion}
\alias{pion}
\title{Fit Pseudo Scalar Sector}
\description{
  fits one or several cosh and sinh to data for extracting the
  pseudo-scalar meson mass, excited states in the pseudo-scalar channel
  and determines decays constant and \eqn{m_\mathrm{PCAC}}{m_paca}
}
\usage{
pion <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
                variational=list(ta=3, tb=4, N=6), ind.vec=c(1,3,4,5),
                no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
                method="uwerr")
}
\arguments{
  \item{cmicor}{
    data to be fitted to as e.g. the output of \code{\link{readcmicor}}
  }

  \item{mu}{
    the value of the bare quark twisted mass
  }
  \item{kappa}{
    the value of the hopping parameter
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
    number of measurements to be discarded at the beginning of the time series.
  }
  \item{variational}{
    list of parameters used for the variational analysis
  }
  \item{ind.vec}{
    index vector indexing the column numbers in cmicor to be used
  }
  \item{no.masses}{
    number of masses to be extracted. This argument will set the number
    of exponentials to be used in the fit.

    Note that this is not yet stable for \code{no.masses > 1}.
  }
  \item{matrix.size}{
    matrix size to be used in the fit. Can be currently set to 2,4 and 6.
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
    \dQuote{boot} or \dQuote{all}.
  }
}
\value{
  returns an object of type cfit with the following items

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
  \item{uwerrresultmps2}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}} if no.masses larger than 1.
  }
  \item{uwerrresultmps3}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}}, if no.masses larger than 2.
  }
  \item{uwerrresultfps}{
    the result of the time series analysis for the pseudo-scalar decay
    constant as  carried out by \code{\link{uwerr}}.
  }
  \item{uwerrresultmpcac}{
    the result of the time series analysis for the PCAC mass as
    carried out by \code{\link{uwerr}} if no.masses larger than 1.
  }
  \item{effmass}{
    effective masses
  }
  \item{mu}{
    the value of the bare quark twisted mass
  }
  \item{kappa}{
    the value of the hopping parameter
  }
  \item{variational.masses}{
    mass values as determined by the variational analysis
  }
  \item{no.masses}{
    no.masses determined
  }
  \item{matrix.size}{
    size of the data matrix
  }
  \item{boot}{
    object returned by the call to \code{\link{boot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{tsboot}{
    object returned by the call to \code{\link{tsboot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
}
\details{
  The data must be ordered as in the output of
  \code{\link{readcmicor}}, see \code{help(readcmicor)} for details.
  
  The expected order of gamma matrices and operators (local-local,
  local-fuzzed, fuzzed-local and fuzzed-fuzzed) (fuzzed = non-local)
  is as follows for all charged mesons:
  
  1) the 4 operators for each type must be sorted like local-local,
  local-fuzzed, fuzzed-local and fuzzed-fuzzed. (fuzzed=non-local)
  
  2) The 20 available types must be in the following order:
  order PP PA AP AA 44 P4 4P A4 4A   for pion like  \eqn{P=\gamma_5}{P=g5} \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}
  order 44 VV AA 4V V4 4A A4 VA AV   for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4} \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}
  order BB SS   - total 20    \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5}  \eqn{S=I}
  (cases with space index "i" are summed over i=1,2,3)
  (best choice is weaker coupling at sink - ie PA rather than AP
  order of magnitudes \eqn{P > 4 > A}  (4 mixes A)
  order of magnitudes \eqn{4\sim A > V}{4 ~ A > V}  (A mixes V))
  
  itype=21 is conserved vector current at sink,  \eqn{\gamma_5}{g5} at source 
  ( iobs is LV  1,  FV  5 )

}
%\references{
%}
\seealso{
  \code{\link{rho}}, \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{plot.cfit}}, \code{\link{summary.cfit}}, \code{\link{pp}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{rhofit <- pion(cmicor, kappa=0.160856, mu=0.004, t1=10, t2=23, no.masses=1, matrix.size=4)}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
