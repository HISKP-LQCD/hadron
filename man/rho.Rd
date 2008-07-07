\name{rho}
\alias{rho}
\title{Fit Vector Channel}
\description{
  fits one or several cosh and sinh to data for extracting the vector
  meson mass and excited states in the vector channel.

  This is so far only for twisted mass lattice QCD in the twisted basis.
}
\usage{
rho(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
    variational=list(ta=4, tb=5, N=6), ind.vec=c(1,3,4,5),
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
    of exponentials to be used in the fit
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
  \item{nrep}{
     vector (N1, N2, ...) of replica length N1, N2. If missing it is
     assumed that there is only one ensemble. If there are two or more
     replica the parameter \code{skip} has no effect.
  }
}
\value{
  returns an object of type cfit \code{rhofit} the following items

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
  \item{uwerrresultmv}{
    the result of the time series analysis for the lowest mass as
    carried out by \code{\link{uwerr}} 
  }
  \item{uwerrresultmv2}{
    the result of the time series analysis for the second lowest mass as
    carried out by \code{\link{uwerr}} if no.masses larger than 1.
  }
  \item{uwerrresultmv3}{
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
  \item{mv.boot}{
    object returned by the call to \code{\link{boot}} if \code{method}
    was set correspodingly. Otherwise \code{NULL}.
  }
  \item{mv.tsboot}{
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
  
  2) The 20 available types must be in the following order:\cr
  order PP PA AP AA 44 P4 4P A4 4A   for pion like  \eqn{P=\gamma_5}{P=g5} \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}\cr
  order 44 VV AA 4V V4 4A A4 VA AV   for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4} \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}\cr
  order BB SS   - total 20    \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5}
  \eqn{S=I}\cr

  In this routine only  44 VV AA 4V V4 4A A4 VA AV are used. See also \link{cmicor}!
  
  (cases with space index "i" are summed over i=1,2,3)\cr
  (best choice is weaker coupling at sink - ie PA rather than AP\cr
  order of magnitudes \eqn{P > 4 > A}  (4 mixes A)\cr
  order of magnitudes \eqn{4\sim A > V}{4 ~ A > V}  (A mixes V))\cr
  
  itype=21 is conserved vector current at sink,  \eqn{\gamma_5}{g5} at source 
  ( iobs is LV  1,  FV  5 )

  rho will perform a fit of of the following matrix
  \tabular{lcccccc}{
    \tab 4L \tab 4F \tab AL \tab AF \tab VL \tab VF \cr
     4L \tab p1 p1 cosh \tab p1 p2
      cosh \tab p1 p3 sinh \tab p1 p4
      sinh \tab p1 p5 sinh \tab p1 p6 sinh \cr
     4F \tab p2 p1 cosh \tab p2 p2 cosh \tab p2
      p3 sinh \tab p2 p4 sinh \tab p2 p5 sinh \tab p2
      p6 sinh \cr 
     AL \tab p3 p1 sinh \tab p3 p2 sinh \tab p3
      p3 cosh \tab p3 p4 cosh \tab p3 p5 cosh \tab p3
      p6 cosh \cr
     AF \tab p4 p1 sinh \tab p4 p2 sinh \tab p4
      p3 cosh \tab p4 p4 cosh \tab p4 p5 cosh \tab p4
      p6 cosh \cr
     VL \tab p5 p1 sinh \tab p5 p2 sinh \tab p5
      p3 cosh \tab p5 p4 cosh \tab p5 p5 cosh \tab p5
      p6 cosh \cr
     VF \tab p6 p1 sinh \tab p6 p2 sinh \tab p6
      p3 cosh \tab p6 p4 cosh \tab p6 p5 cosh \tab p6
      p6 cosh \cr 
  }
  for coupling parameter \eqn{p_1} to \eqn{p_6} and a mass \eqn{m_\mathrm{V}}{mv} entering
  \eqn{\cosh}{cosh} (and \eqn{\sinh}{sinh} in the same way) as
  \eqn{\cosh(-m_\mathrm{V}(T/2-t))}{cosh(-mv(T/2-t))}. The values of \eqn{t} are
  running from \code{t1} to \code{t2} as specified by the user.

  Corresponding to the value of \code{matrix.size} only a submatrix of
  the matrix given above is fitted. Moreover, if \code{no.masses} larger
  than one additional masses and coupling parameters are introduced.
  
}
%\references{
%}
\seealso{
  \code{\link{pion}}, \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{plot.cfit}}, \code{\link{pp}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{rhofit <- rho(cmicor, kappa=0.160856, mu=0.004, t1=8, t2=18, no.masses=1, matrix.size=4)}
}
\author{Carsten Urbach, \email{carsten.urbach@liverpool.ac.uk}}
\keyword{optimize}
\keyword{ts}
