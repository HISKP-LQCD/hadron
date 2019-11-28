\name{cdh}
\alias{cdh}
\title{finite size corrections a la Colangelo, Duerr, Haefeli}
\description{
  finite size corrections a la Colangelo, Duerr, Haefeli
}
\usage{
cdh(parm = rep(0, times=6), rev=-1,
               aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
               ampiV, afpiV, aF0, a_fm, L,
               printit=FALSE, incim6 = FALSE,
               rtilde=c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0),
               use.cimpl=TRUE)
}
\arguments{
  \item{parm}{
    parameters
  }
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
  \item{a_fm}{the value of the lattice spacing in fermi}
  \item{L}{the lattice spatial extent}
  \item{printit}{if set to TRUE the corrections are printed}
  \item{incim6}{in- or exclude the NNNLO correction for the mass}
  \item{rtilde}{the low energy constants \eqn{\tilde{r}}{rtilde}, needed
  only if \code{incim6=TRUE}}
  \item{use.cimpl}{use the four times faster direct c Implementation of
    the correction routine}
}
\value{
  a list with the corrected values for mpi and fpi
}
\details{
  see reference for details. We use the simplyfied formulae for the S
  quantities, see eq. (59) in the reference.
}
\references{
  Gilberto Colangelo, Stephan Durr, Christoph Haefeli, Nucl.Phys.B721:136-174,2005. hep-lat/0503014
}
%\seealso{
%}
\examples{
L <- c(24, 24, 24, 24, 32)
mps <- c(0.14448, 0.17261, 0.19858, 0.22276, 0.14320)
fps <- c(0.06577, 0.07169, 0.07623, 0.07924, 0.06730)
aLamb1 <- 0.05
aLamb2 <- 0.5
aLamb3 <- 0.38
aLamb4 <- 0.66
cdhres <- cdh(rev=+1, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3, aLamb4=aLamb4,
              ampiV=mps, afpiV=fps, aF0=fps, a_fm=0.08, L=L, printit=TRUE,
              incim6=FALSE)
cdhres$mpiFV
cdhres$fpiFV
}
\author{Carsten Urbach <curbach@gmx.de>}
