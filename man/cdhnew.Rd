\name{cdhnew}
\alias{cdhnew}
\title{finite size corrections a la Colangelo, Duerr, Haefeli, but as a
  series in the quark mass}
\description{
  finite size corrections a la Colangelo, Duerr, Haefeli, but as a
  series in the quark mass
}
\usage{
cdhnew(rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
       ampiV, afpiV, aF0, a2B0mu, L, printit=FALSE)
}
\arguments{
  \item{rev}{bla}
}
\value{
  a list with the corrected values for mpi and fpi
}
\details{
  bla
}
\references{
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
