\name{cdhnew}
\alias{cdhnew}
\title{finite size corrections a la Colangelo, Duerr, Haefeli, but
  re-expanded as series in the quark mass}
\description{
  finite size corrections a la Colangelo, Duerr, Haefeli, but
  re-expanded as series in the quark mass
}
\usage{
cdhnew(parm = rep(0, times=6), rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
       ampiV, afpiV, aF0, a2B0mu, L, printit=FALSE, use.cimpl = TRUE)
}
\arguments{
  \item{rev}{\eqn{rev=-1} corrects from \eqn{L} to \eqn{L=\infty}{L =
      infinity}, \eqn{rev=+1} the other way around}
  \item{aLamb1}{The four low energy
    \eqn{\Lambda_{1-4}}{Lambda1-Lambda4}constants in lattice units.}
  \item{aLamb2}{
    see \code{aLamb1}.
  }
  \item{aLamb3}{
    see \code{aLamb1}.
  }
  \item{aLamb4}{
    see \code{aLamb1}.
  }
  \item{ampiV}{pseudo scalar mass values to be corrected}
  \item{afpiV}{pseudo scalar decay constant values to be corrected}
  \item{aF0}{\eqn{af_0}{af0} in lattice units}
  \item{a2B0mu}{\eqn{2B_0\mu}{2 B0 mu} in lattice units, where
    \eqn{\mu}{mu} is the quark mass and \eqn{B_0}{B0} a low energy
    constant}
  \item{L}{the lattice spatial extent}
  \item{printit}{if set to TRUE the corrections are printed}
  \item{use.cimpl}{use the four times faster direct c Implementation of
    the correction routine}
  \item{parm}{
    m parameters
  }
}
\value{
  a list with the corrected values for mpi and fpi
}
\details{
  see reference for details. We use the simplyfied formulae for the S
  quantities, see eq. (59) in first reference.
}
\references{
  Gilberto Colangelo, Stephan Durr, Christoph Haefeli,
  Nucl.Phys.B721:136-174,2005. hep-lat/0503014

  and

  R. Frezzotti, V. Lubicz, S. Simula,
  arXiv:0812.4042 [hep-lat] 
}
%\seealso{
%}
\examples{
mu <- c(0.004, 0.006, 0.008, 0.010, 0.004)
L <- c(24, 24, 24, 24, 32)
mps <- c(0.14448, 0.17261, 0.19858, 0.22276, 0.14320)
fps <- c(0.06577, 0.07169, 0.07623, 0.07924, 0.06730)
aLamb1 <- 0.05
aLamb2 <- 0.5
aLamb3 <- 0.38
aLamb4 <- 0.66
aF0    <- 0.051
a2B    <- 5.64
cdhres <- cdhnew(rev=+1, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3,
                 aLamb4=aLamb4, ampiV=mps, afpiV=fps, aF0=aF0,
                 a2B0mu=a2B*mu, L=L, printit=TRUE)
cdhres$mpiFV
cdhres$fpiFV
}
\author{Carsten Urbach <curbach@gmx.de>}
%\keyword{}
