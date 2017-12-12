\name{fs.qcotdelta}
\alias{fs.mpia0}
\alias{fs.a0}
\title{Finite Size Corrections to \eqn{q\cot\delta}{qcotdelta} for I=2 \eqn{\pi\pi}{pipi}
  near threshold}
\description{
  \code{fs.qcotdelta} computes the finite size corrections to
  \eqn{q\cot\delta}{qcotdelta} while \code{fs.mpia0} computes the
  corresponding finite size corrections to \eqn{M_\pi a_0}{Mpi a0}
  directly using the Gasser Leutwyler result from \eqn{M_\pi}{Mpi}.
}
\usage{
  fs.qcotdelta(L, mps)
  fs.mpia0(fps, mps, L)
  fs.a0(a0, mps, L)
}
\arguments{
  \item{L}{
    spatial lattice extend as a scalar variable (must not be a vector)
  }
  \item{mps}{
    pion mass as a scalar variable (must not be a vector)
  }
  \item{fps}{
    pion decay constant as a scalar variable (must not be a vector)
  }
  \item{a0}{
    scattering length at finite L
  }
}
\value{
  returns a numeric value representing the finite size correction or in
  case of \code{fs.a0} the corrected value for a0.
}
\references{
  For the original formula see Eq. (31) from hep-lat/0601033
}
\examples{
  fs.qcotdelta(mps=0.123, L=24)
}
\author{Carsten Urbach, \email{curbach@gmx.de}}
