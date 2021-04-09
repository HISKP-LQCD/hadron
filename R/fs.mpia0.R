## finite size correction to q cot(delta) in the I=2 pipi scattering
## case from Ref. hep-lat/0601033
##
## this formula is valid near threshold only!
##
## The difference is defined as Delta = FV - Vinfty

## this is the original aymptotic formula Eq (31) from Ref.
## hep-lat/0601033
## which can be applied to a0 by using the effective range expansion
## q*cot(delta) = 1/(-a0)


#' Finite Size Corrections to \eqn{q\cot\delta}{qcotdelta} for I=2
#' \eqn{\pi\pi}{pipi} near threshold
#' 
#' \code{fs.qcotdelta} computes the finite size corrections to
#' \eqn{q\cot\delta}{qcotdelta} while \code{fs.mpia0} computes the
#' corresponding finite size corrections to \eqn{M_\pi a_0}{Mpi a0} directly
#' using the Gasser Leutwyler result from \eqn{M_\pi}{Mpi}.
#' 
#' 
#' @param L spatial lattice extent as a scalar variable (must not be a vector)
#' @param mps pion mass as a scalar variable (must not be a vector)
#' @return returns a numeric value representing the finite size correction or
#' in case of \code{fs.a0} the corrected value for a0.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @references For the original formula see Eq. (31) from hep-lat/0601033
#' @examples
#' 
#'   fs.qcotdelta(mps=0.123, L=24)
#' 
#' @export fs.qcotdelta
fs.qcotdelta <- function(mps, L) {
  ## the 7th is zero, so we skip sqrt(7)
  cn <- c(6, 12, 8, 6, 24, 24, 12)
  n <- c(1, sqrt(2), sqrt(3), 2, sqrt(5), sqrt(6), sqrt(8))
  return(-mps/sqrt(2*pi)*sum( cn*exp(-n*mps*L)/sqrt(n*mps*L)*(1-227/(24*n*mps*L)) ))
}

#' Finite Size Corrections to \eqn{q\cot\delta}{qcotdelta} for I=2
#' \eqn{\pi\pi}{pipi} near threshold
#' 
#' \code{fs.qcotdelta} computes the finite size corrections to
#' \eqn{q\cot\delta}{qcotdelta} while \code{fs.mpia0} computes the
#' corresponding finite size corrections to \eqn{M_\pi a_0}{Mpi a0} directly
#' using the Gasser Leutwyler result from \eqn{M_\pi}{Mpi}.
#' 
#' 
#' @param L spatial lattice extent as a scalar variable (must not be a vector)
#' @param mps pion mass as a scalar variable (must not be a vector)
#' @param a0 scattering length at finite L
#' @return returns a numeric value representing the finite size correction or
#' in case of \code{fs.a0} the corrected value for a0.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @references For the original formula see Eq. (31) from hep-lat/0601033
#' @examples
#' fs.a0(a0=1., mps=0.123, L=24)
#' @export
fs.a0 <- function(a0, mps, L) {
  delta <- fs.qcotdelta(mps, L)
  return(1./(1./a0 + delta))
}

## this is the formula from
## arXiv:0909.3255
## directly for mpi*a0
#' Finite Size Corrections to \eqn{q\cot\delta}{qcotdelta} for I=2
#' \eqn{\pi\pi}{pipi} near threshold
#' 
#' \code{fs.qcotdelta} computes the finite size corrections to
#' \eqn{q\cot\delta}{qcotdelta} while \code{fs.mpia0} computes the
#' corresponding finite size corrections to \eqn{M_\pi a_0}{Mpi a0} directly
#' using the Gasser Leutwyler result from \eqn{M_\pi}{Mpi}.
#' 
#' 
#' @param L spatial lattice extent as a scalar variable (must not be a vector)
#' @param mps pion mass as a scalar variable (must not be a vector)
#' @param fps pion decay constant as a scalar variable (must not be a vector)
#' @return returns a numeric value representing the finite size correction or
#' in case of \code{fs.a0} the corrected value for a0.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @references For the original formula see Eq. (31) from hep-lat/0601033
#' @examples
#' fs.mpia0(mps=0.123, fps=0.2, L=24)
#' @export
fs.mpia0 <- function(mps, fps, L) {

  fn <- function(n, cn, mpsL) {
    cn*exp(-n*mpsL)/sqrt(n*mpsL)*(1-17/(8*n*mpsL))
  }
  cn <- c(6, 12, 8, 6, 24, 24, 0, 12)
  n <- c(1, sqrt(2), sqrt(3), 2, sqrt(5), sqrt(6), sqrt(7), sqrt(8))
  S <- sum(fn(n, cn, mpsL=mps*L))

  return(mps^4/fps^4/2^(13/2)/pi^(5/2)*S)
}

