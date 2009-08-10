## Nf = 2 - Schema RI'
zetazp<- function(zp0, alpha0, alpha2, nl = 3) {
  
##  alm <- alpha_s(mu = mu2, nl = nl, lam0 = lam0, Nc = Nc, Nf = Nf)
##  al0 <- alpha_s(mu = mu0, nl = nl, lam0 = lam0, Nc = Nc, Nf = Nf)

  if(nl == 3) {
    cmu <- (alpha2)^(-12./29.) * (1. - 8.55727 * alpha2 - 125.423 * alpha2^2 - 3797.71 * alpha2^3) 
    cm0 <- (alpha0)^(-12./29.) * (1. - 8.55727 * alpha0 - 125.423 * alpha0^2 - 3797.71 * alpha0^3) 
  }
  else if(nl == 2) {
    cmu <- (alpha2)^(-12./29.) * (1. - 8.55727 * alpha2 - 125.423 * alpha2^2)
    cm0 <- (alpha0)^(-12./29.) * (1. - 8.55727 * alpha0 - 125.423 * alpha0^2)
  }
  else if(nl == 1) {
    cmu <- (alpha2)^(-12./29.) * (1. - 8.55727 * alpha2)
    cm0 <- (alpha0)^(-12./29.) * (1. - 8.55727 * alpha0)
  }
  else {
    cmu <- (alpha2)^(-12./29.)
    cm0 <- (alpha0)^(-12./29.)
    if(nl > 3 || nl < 0) {
      warning("zeta_zp used with nl > 3 or nl < 0, using nl=0\n")
    }
  } 
  return(cmu/cm0*zp0)
}
