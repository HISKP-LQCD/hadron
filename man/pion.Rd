\name{pion}
\alias{pion}
\title{Fit Pseudo Scalar Sector}
\description{
  fits one or several cosh and sinh to data for extracting the
  pseudo-scalar meson mass, excited states in the pseudo-scalar channel
  and determines decays constant and \eqn{m_\mathrm{PCAC}}{m_paca}

  This is so far only for twisted mass lattice QCD in the twisted basis.
}
\usage{
  pion(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0,
       variational=list(ta=3, tb=4, N=6), ind.vec=c(1,3,4,5),
       no.masses=1, matrix.size=2, boot.R=99, boot.l=10, tsboot.sim="geom",
       method="uwerr", fit.routine="optim", mass.guess, par.guess, nrep)
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
  returns an object of \code{class} \code{pionfit} with the following items

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
  \item{uwerrresultzv}{
    the gamma method analysis for \eqn{Z_V}{Z_V}, if matrix size equals 6
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
    no.masses determined, copied from input
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
  \item{var.res}{
    the full result of the variational analysis, as returned by a call
    to \code{\link{variational}}
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
  A call to \code{pion} returns an object, for which several generic
  functions are available, such as \code{summary} and \code{plot}. The
  function should always be used by also considering the output of
  \code{plot}. In particular, in case of the gamma matrix method it
  should always be checked whether the autocorrelation function could be
  integrated in a sensible way. For this reason the autocorrelation
  function is plotted as well as the integrated autocorrelation time as
  a function of the cut-off value for the integral. Obviously, the
  autocorrelation function should be zero at the cut-off and for larger
  values only fluctuate around zero. In other words, the autocorrelation
  time should reach a plateau. Changing the parameter \code{S} (see
  above) allows to influence the automatic determination of the cut-off
  value (see \code{\link{uwerr}} for details).

  Another useful check is also to devide the full ensemble into two (or more)
  sub-ensembles (using the parameter \code{nrep}). The summary function
  will then also compute a Q-value, which should be larger than
  \eqn{0.1} at least for the error to be trustworthy. Otherwise the
  autocorrlation times are presumably too large to determine the error
  reliably.

  Finally it is also usefull to check the value of \eqn{\chi}{chi} per
  data point. Such a plot is produced when \code{plot} is used. 

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
  order BB SS   - total 20    \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5}\cr
  \eqn{S=I}
  
  In this routine only PP PA AP AA 44 P4 4P A4 4A are used. See also \code{cmicor}!

  (cases with space index "i" are summed over i=1,2,3)
  (best choice is weaker coupling at sink - ie PA rather than AP\cr
  order of magnitudes \eqn{P > 4 > A}  (4 mixes A)\cr
  order of magnitudes \eqn{4\sim A > V}{4 ~ A > V}  (A mixes V))\cr
  
  itype=21 is conserved vector current at sink,  \eqn{\gamma_5}{g5} at source 
  ( iobs is LV  1,  FV  5 )

  pion will perform a fit of of the following matrix
  \tabular{lcccccc}{
    \tab PL \tab PF \tab AL \tab AF \tab 4L \tab 4F \cr
     PL \tab p1 p1 cosh \tab p1 p2
      cosh \tab p1 p3 sinh \tab p1 p4
      sinh \tab p1 p5 sinh \tab p1 p6 sinh \cr
     PF \tab p2 p1 cosh \tab p2 p2 cosh \tab p2
      p3 sinh \tab p2 p4 sinh \tab p2 p5 sinh \tab p2
      p6 sinh \cr 
     AL \tab p3 p1 sinh \tab p3 p2 sinh \tab p3
      p3 cosh \tab p3 p4 cosh \tab p3 p5 cosh \tab p3
      p6 cosh \cr
     AF \tab p4 p1 sinh \tab p4 p2 sinh \tab p4
      p3 cosh \tab p4 p4 cosh \tab p4 p5 cosh \tab p4
      p6 cosh \cr
     4L \tab p5 p1 sinh \tab p5 p2 sinh \tab p5
      p3 cosh \tab p5 p4 cosh \tab p5 p5 cosh \tab p5
      p6 cosh \cr
     4F \tab p6 p1 sinh \tab p6 p2 sinh \tab p6
      p3 cosh \tab p6 p4 cosh \tab p6 p5 cosh \tab p6
      p6 cosh \cr 
  }
  for coupling parameter \eqn{p_1} to \eqn{p_6} and a mass \eqn{m_\mathrm{PS}}{mps} entering
  \eqn{\cosh}{cosh} (and \eqn{\sinh}{sinh} in the same way) as
  \eqn{\cosh(-m_\mathrm{PS}(T/2-t))}{cosh(-mps(T/2-t))}. The values of \eqn{t} are
  running from \code{t1} to \code{t2} as specified by the user.

  Corresponding to the value of \code{matrix.size} only a submatrix of
  the matrix given above is fitted. Moreover, if \code{no.masses} larger
  than one additional masses and coupling parameters are introduced.

  The pion decay constant is computed from the relation
  \deqn{f_\mathrm{PS} = 4\kappa\mu\frac{p_1}{\sqrt{m_\mathrm{PS}}^3}}{%
    fps = 4 kappa mu p_1/sqrt(mps)^3}
  and the PCAC mass (for \code{matrix.size} > 2) from
  \deqn{m_\mathrm{PCAC} = m_\mathrm{PS} \frac{p_3}{2p_1}}{%
    mpcac = mps p_3/(2 p_1)}
  when the \code{\link{summary.pionfit}} function is called.
  For \code{matrix.size} > 4 also \eqn{Z_V} is computed as
  \deqn{Z_\mathrm{V}=\frac{2\mu}{m_\mathrm{PS}}\frac{p_1}{p_5}}{%
    Z_V = 2 mu p_1/(mps p_5)}
}
%\references{
%}
\seealso{
  \code{\link{readcmicor}}, \code{\link{uwerr}},
  \code{\link{variational}}
}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
\dontrun{pionfit <- pion(cmicor, kappa=0.160856, mu=0.004, t1=10, t2=23,}
\dontrun{         no.masses=1, matrix.size=4)}
\dontrun{summary(pionfit)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{optimize}
\keyword{ts}
