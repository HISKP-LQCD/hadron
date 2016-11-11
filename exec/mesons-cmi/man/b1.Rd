\name{b1}
\alias{b1}
\alias{a0}
\title{Fit B1 or A0 Sector}
\description{
  fits b1 or b0 meson masses.

  This is so far only for twisted mass lattice QCD in the twisted basis.
}
\synopsis{
b1(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
   variational=list(ta=3, tb=4, N=2), ind.vec=c(1,3,4,5),
   no.masses=1, boot.R=99, boot.l=10, tsboot.sim="geom",
   method="uwerr", mass.guess, par.guess, nrep)
}
\usage{
b1(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
   variational=list(ta=3, tb=4, N=2), ind.vec=c(1,3,4,5),
   no.masses=1, boot.R=99, boot.l=10, tsboot.sim="geom",
   method="uwerr", mass.guess, par.guess, nrep)
a0(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
   variational=list(ta=3, tb=4, N=2), ind.vec=c(1,3,4,5),
   no.masses=1, boot.R=99, boot.l=10, tsboot.sim="geom",
   method="uwerr", mass.guess, par.guess, nrep)
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
    number of measurements to be discarded at the beginning of the time
    series. \code{skip} has no effect if two or more replica are used,
    see argument \code{nrep}.
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
}
\value{
  returns an object of \code{class} \code{b1fit} (or \code{a0fit}) with
  the following items 

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
  \item{uwerrresultmb1}{
    the result of the time series analysis for the lowest mass as
    carried out by \code{\link{uwerr}} 
  }
  \item{uwerrresultmb12}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}} if no.masses larger than 1.
  }
  \item{uwerrresultmb13}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}}, if no.masses larger than 2.
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
  \item{nrep}{
    vector \code{nrep} as given as input.
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
  
  2) The 20 available types must be in the following order:\cr
  order PP PA AP AA 44 P4 4P A4 4A   for pion like  \eqn{P=\gamma_5}{P=g5} \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}\cr
  order 44 VV AA 4V V4 4A A4 VA AV   for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4} \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}\cr
  order BB SS   - total 20    \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5}  \eqn{S=I}\cr

  In this routine only BB (or SS for a0)  are used. 

  (cases with space index "i" are summed over i=1,2,3)\cr
  (best choice is weaker coupling at sink - ie PA rather than AP\cr
  order of magnitudes \eqn{P > 4 > A}  (4 mixes A)\cr
  order of magnitudes \eqn{4\sim A > V}{4 ~ A > V}  (A mixes V))\cr
  
  itype=21 is conserved vector current at sink,  \eqn{\gamma_5}{g5} at source 
  ( iobs is LV  1,  FV  5 )

  pion will perform a fit of of the following matrix
  \tabular{lcc}{
    \tab BL \tab BF \cr
    BL \tab p1 p1 cosh \tab p1 p2 cosh \cr
    BF \tab p2 p1 cosh \tab p2 p2 cosh \cr
  }
  for coupling parameter \eqn{p_1}{p1} to \eqn{p_2}{p2} and a mass \eqn{m_\mathrm{b1}}{mb1} entering
  \eqn{\cosh}{cosh} as
  \eqn{\cosh(-m_\mathrm{b1}(T/2-t))}{cosh(-mb1(T/2-t))}. The values of \eqn{t} are
  running from \code{t1} to \code{t2} as specified by the user.

  If \code{no.masses} larger
  than one additional masses and coupling parameters are introduced.

  (For \eqn{m_\mathrm{a0}}{ma0} the same is done)
}
%\references{
%}
\seealso{
  \code{\link{rho}}, \code{\link{pion}}, \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{plot.cfit}}, \code{\link{summary.b1fit}}, \code{\link{pp}}
  \code{\link{summary.a0fit}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{b1fit <- b1(cmicor, kappa=0.160856, mu=0.004, t1=10, t2=23, no.masses=1)}
\dontrun{summary(b1fit)}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
