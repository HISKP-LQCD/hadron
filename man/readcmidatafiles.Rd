\name{readcmifiles}
\alias{readcmicor}
\alias{readcmifiles}
\alias{readcmidatafiles}
\alias{readcmiloopfiles}
\title{Read Single Data Files in Chris Michael Format}
\description{
  reads data from single files in Chris Michael format
}
\usage{
readcmicor(filename,
           colClasses=c("integer","integer","integer","numeric","numeric","integer"), skip=0)

readcmifiles(files, excludelist=c(""), skip, verbose=FALSE,
             colClasses, obs=NULL, obs.index, avg=1, stride=1)

readcmidatafiles(files, excludelist=c(""), skip=1, verbose=FALSE,
                 colClasses=c("integer", "integer","integer","numeric","numeric"),
                 obs=NULL, obs.index=1, avg=1, stride=1)

readcmiloopfiles(files, excludelist=c(""), skip=0, verbose=FALSE,
                 colClasses=c("integer", "integer","integer","integer",
                              "numeric","numeric","numeric","numeric"),
                 obs=NULL, obs.index=2)
}
\arguments{
  \item{filename}{
    A single file name to read from.
  }
  \item{files}{
    list of filenames to be read. Can be created using
    \code{getorderedfilelist}. 
  }
  \item{skip}{
    Number of lines to be skipped at the beginning of each file
  }
  \item{excludelist}{
    files to exclude from reading.
  }
  \item{verbose}{
    Increases verbosity of the function.
  }
  \item{colClasses}{
    The column data type classes, the \code{read.table}.
  }
  \item{obs}{
    To reduce memory consumption it is possible to extract only one of
    the observales. The column in which to match \code{obs} is to be
    given with \code{obs.index}. This will only be affective if
    \code{obs} is not \code{NULL}.
  }
  \item{obs.index}{
    The column in which to match \code{obs} is to be
    given with \code{obs.index}.
  }
}
\value{
  \code{readcmicor} returns an object of class \code{cmicor}, read from
  a single file.
  
  \code{readcmidatafiles} returns an object of class \code{cmicor},
  which is an \code{rbind} of all \code{data.frame}s read from the
  single files in the filelist.

  \code{readcmiloopfiles} returns an object of class \code{cmiloop},
  which is an \code{rbind} of all \code{data.frame}s read from the
  single files in the filelist.
}
\details{
  These functions reads data from single data files. It is assumed that
  every file has the same number of columns.

  The cmi (Chris Michael) format for connected correlators comprises 6
  colums per file: 1) the observable type number (itype); 2) the operator
  type number (iobs); 3) the time difference from source going from 0 to \eqn{T/2}
  for each operator type; 4) \eqn{c_1}{c1} correlator value at time value forward in
  time; 5) \eqn{c_2}{c2} correlator value at time value backward in time; 6) number of
  gauge configuration.

  There are scripts shipped with the package converting the output
  written into seperate files for each gauge configuration into the
  expected format. They are called \code{puttogether.sh} and
  \code{puttogether_reverse.sh} which will sort with increasing and
  with decreasing gauge configuration number, respectively.

  Note, that the normalisation of correlators needs multiplication by
  factor of \eqn{0.5} (and possible \eqn{(2*\kappa)^2}{(2*k)^2} 
  and \eqn{L^3} factors dependent on your conventions).

  The values of \code{itype} run from \code{1} to the total number of
  gamma matrix combinations available. \code{iobs} equals \code{1} for
  local-local correlators, \code{3} for local-smeared, \code{5} for
  smeared-local and \code{7} for smeared-smeared
  
  For charged mesons the order of gamma-matrix combinations is as
  follows:\cr
  order PP PA AP AA 44 P4 4P A4 4A   for pion like  \eqn{P=\gamma_5}{P=g5} \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{4=\gamma_4}{4=g4}\cr
  order 44 VV AA 4V V4 4A A4 VA AV   for rho-a1 like \eqn{4=\gamma_i\gamma_4}{4=gig4} \eqn{V=\gamma_i}{V=gi} \eqn{A=\gamma_i\gamma_5}{A=gig5}\cr
  order BB SS   - total 20    \eqn{\gamma_i\gamma_4\gamma_5}{B=gig4g5} 
  \eqn{S=I}\cr
  itype=21 is conserved vector current at sink,  \eqn{\gamma_5}{g5} at source

  
  For neutral mesons the order of gamma-matrix combinations is as
  follows:\cr
  order PP PA AP AA II PI IP AI IA   for pion like  \eqn{P=\gamma_5}{P=g5} \eqn{A=\gamma_4\gamma_5}{A=g4g5} \eqn{I=1}{1=1}\cr
  order 44 VV BB 4V V4 4B B4 VB BV   for rho-b1 like \eqn{4=\gamma_i\gamma_4}{4=gig4} \eqn{V=\gamma_i}{V=gi} \eqn{B=\gamma_i\gamma_4\gamma_5}{B=gig4g5}\cr
  order XX AA   - total 20    for a0-X like \eqn{A=\gamma_i\gamma_5}{A=gig5} \eqn{X=\gamma_4}{X=g4}

  For loops (disconnected contributions to neutral mesons) the
  convention is as follows: files are assumed to have eight columns with
  gauge, gamma, t, sample, ReTL, ImTL, ReTF, ImTF, where gamma is 1 to
  16 as list of (hermitian) gamma matrices:  
  order          g_5  g_1    g_2   g_3\cr
  -ig_4*       g_5  g_1    g_2   g_3\cr
  -ig_5*     i*g_5  g_1    g_2   g_3  ie 1,..\cr
  -ig_5g_4* -i*g_5  g_1    g_2   g_3  ie g_4, g_5*row 2\cr
  (so P is 1; A4 is 5;  S is 9;  A_i is 10,11,12    etc)
  
  t is t-value of trace (here spatial momentum is zero)
  sample is sample number 1,...24 (or 96)
  ReTL is real part of trace at time t, with gamma combination given 
  and Local operator  (F is Fuzzed == non-local) operator).
  
  Normalisation is trace M^-1 with M=1+...

  To make a disconnected correlator, one combines these traces for 
  different t (and different sample number) as a product. Note 
  only Re Gamma=1  and Im Gamma=gamma_5  have VEV's, see
  \code{\link{computeDisc}}
}
\seealso{
  \code{\link{getorderedfilelist}},
    \code{\link{extract.obs}}, \code{\link{readcmidisc}}
}
\examples{
library(hadron)
\dontrun{filelist <- getorderedfilelist("ouptrc")}
\dontrun{cmicor <- readcmidatafiles(filelist, skip=1)}
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
\keyword{file}
