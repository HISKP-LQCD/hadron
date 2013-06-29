\name{readcmicor}
\alias{readcmicor}
\alias{cmicor}
\title{Read Data In Chris Michael Format}
\description{
  reads data from file in Chris Michael format
}
\usage{
  readcmicor(filename)
}
\arguments{

  \item{filename}{
    filename of the data file
  }
}
\value{
  returns a data frame of class \dQuote{cmicor} containing the data
}
\details{
  This function read data from a file. The data in the file is assumed
  to have 6 columns: 1) the observable type number (itype); 2) the operator
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

  
}
%\references{
%}
%\seealso{
%}
\examples{
library(hadron)
\dontrun{cmicor <- readcmicor("pion.dat")}
}
\author{Carsten Urbach \email{carsten.urbach@liverpool.ac.uk}}
\keyword{file}
