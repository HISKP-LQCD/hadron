kappa <- function(m0) {
  return(1./(2*m0 + 8))
}

m0 <- function(kappa) {
  return(1./2/kappa - 4)
}
