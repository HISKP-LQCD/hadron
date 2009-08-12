\name{cdh}
\alias{cdh}
\title{finite size corrections a la Colangelo, Duerr, Haefeli}
\description{
  finite size corrections a la Colangelo, Duerr, Haefeli
}
\usage{
cdh(parm=NULL, rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
                 ampiV, afpiV, aF0, a_fm, L, printit=F, incim6 = TRUE,
                 rtilde=c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0))
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
L <- c(24, 24, 24, 24, 32)
mps <- c(0.14448, 0.17261, 0.19858, 0.22276, 0.14320)
fps <- c(0.06577, 0.07169, 0.07623, 0.07924, 0.06730)
aLamb1 <- 0.05
aLamb2 <- 0.5
aLamb3 <- 0.38
aLamb4 <- 0.66
aF0    <- 0.051
cdhres <- cdh(rev=+1, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3, aLamb4=aLamb4,
              ampiV=mps, afpiV=fps, aF0=aF0, a_fm=0.08, L=L, printit=TRUE,
              incim6=FALSE)
cdhres$mpiFV
cdhres$fpiFV
}
\author{Carsten Urbach <curbach@gmx.de>}
