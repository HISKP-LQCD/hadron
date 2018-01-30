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
fs.qcotdelta <- function(mps, L) {
  ## the 7th is zero, so we skip sqrt(7)
  cn <- c(6, 12, 8, 6, 24, 24, 12)
  n <- c(1, sqrt(2), sqrt(3), 2, sqrt(5), sqrt(6), sqrt(8))
  return(-mps/sqrt(2*pi)*sum( cn*exp(-n*mps*L)/sqrt(n*mps*L)*(1-227/(24*n*mps*L)) ))
}

## use effective range expansion to correct the scattering length a0
## directly
fs.a0 <- function(a0, mps, L) {
  delta <- fs.qcotdelta(mps, L)
  return(1./(1./a0 + delta))
}

## this is the formula from
## arXiv:0909.3255
## directly for mpi*a0
fs.mpia0 <- function(mps, fps, L) {

  fn <- function(n, cn, mpsL) {
    cn*exp(-n*mpsL)/sqrt(n*mpsL)*(1-17/(8*n*mpsL))
  }
  cn <- c(6, 12, 8, 6, 24, 24, 0, 12)
  n <- c(1, sqrt(2), sqrt(3), 2, sqrt(5), sqrt(6), sqrt(7), sqrt(8))
  S <- sum(fn(n, cn, mpsL=mps*L))

  return(mps^4/fps^4/2^(13/2)/pi^(5/2)*S)
}

