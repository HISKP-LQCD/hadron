## Nf = 2 - Schema RI'


#' Computes the running of Z_P from scale mu0 to scale mu2
#' 
#' Computes the running of the renomalisation constant \eqn{Z_P} from scale
#' \eqn{\mu_0}{mu0} to scale \eqn{\mu_2}{mu2} in the renomalisation schema RI'
#' for \eqn{N_f=2}{Nf=2} only. The running is done using perturbation theory up
#' to \eqn{\alpha_s**3}{alpha_s^3} order. The corresponding values of
#' \eqn{\alpha_s}{alpha_s} at the scales \eqn{\mu_0}{mu0} and \eqn{\mu_2}{mu2}
#' are needed as input, see \code{\link{alphas}}.
#' 
#' 
#' @param zp0 initial value of \eqn{Z_P}
#' @param alpha0 \eqn{\alpha_s}{alpha_s} at initial scale
#' @param alpha2 \eqn{\alpha_s}{alpha_s} at final scale
#' @param nl order in PT, range 0 to 3
#' @return returns the value of Z_P at scale mu2 in the RI' scheme
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{alphas}}
#' @examples
#' 
#' al2 <- alphas(mu = 3.0, nl = 3, lam0 = 0.250, Nc = 3, Nf = 2)
#' al0 <- alphas(mu = 2.0, nl = 3, lam0 = 0.250, Nc = 3, Nf = 2)
#' zetazp(zp0 = 0.6, alpha0 = al0, alpha2 = al2, nl = 3)
#' 
#' @export zetazp
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
