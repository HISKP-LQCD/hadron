#' compute alpha strong at given scale
#' 
#' compute alpha strong (\eqn{\alpha_s}{alpha_s}) at given scale \eqn{\mu}{mu}
#' up to N3LO in PT in the RI' renormalisation scheme.
#' 
#' 
#' @param mu the renormalisation scale \eqn{\mu}{mu} in GeV
#' @param nl order in PT, range 0 to 3
#' @param lam0 \eqn{\Lambda_\mathrm{QCD}}{Lambda_QCD} in GeV
#' @param Nc number of colours \eqn{N_c}{Nc}, defaults to 3
#' @param Nf number of flavours \eqn{N_f}{Nf}, default is 2
#' @param use.cimpl Use the C implementation instead of the R implementation,
#' which might improve speed.
#' @return returns the value of alpha strong \eqn{\alpha_s}{alpha_s} at scale
#' \eqn{\mu}{mu}
#' @author Carsten Urbach, \email{curbach@@gmx.de}, Vittorio Lubicz (of the original Fortran code)
#' @seealso \code{\link{zetazp}}
#' @examples
#' 
#' alphas(mu=2.0, nl=3)
#' 
#' @export alphas
alphas <- function(mu, nl=3, lam0=0.250, Nc=3., Nf=2., use.cimpl=TRUE) {
  if(!use.cimpl) {
    return(alphas.R( mu, nl, lam0, Nc, Nf))
  }
  return(.Call("alphas", mu, nl, lam0, as.numeric(Nc), as.numeric(Nf)))
}


alphas.R <- function(mu, nl, lam0, Nc, Nf) {

  Cf <- (Nc^2-1.)/2./Nc
  Z3 <- 1.20206

  b0 <- 11./3.*Nc - 2./3.*Nf
  b1 <- 34./3.*Nc^2 - 38./3.*Nf
  b2 <- 2857./54.*Nc^3 + Cf^2*Nf - 205./18.*Cf*Nc*Nf -
    1415./54.*Nc^2*Nf + 11./9.*Cf*Nf^2 + 79./54.*Nc*Nf^2
  b3 <- (150653./486. - 44./9.*Z3)*Nc^4 + 
    (-39143./162. + 68./3.*Z3)*Nc^3*Nf +
      (7073./486. - 328./9.*Z3)*Cf*Nc^2*Nf + 
        (-2102./27. + 176./9.*Z3)*Cf^2*Nc*Nf +
          23.*Cf^3*Nf + (3965./162. + 56./9.*Z3)*Nc^2*Nf^2 + 
            (338./27. - 176./9.*Z3)*Cf^2*Nf^2 +
              (4288./243. + 112./9.*Z3)*Cf*Nc*Nf^2 + 53./243.*Nc*Nf^3 +
                154./243.*Cf*Nf^3 + 
                  (-10./27. + 88./9.*Z3)*Nc^2*(Nc^2+36.) +
                    (32./27. - 104./9.*Z3)*Nc*(Nc^2+6)*Nf +
                      (-22./27. + 16./9.*Z3)*(Nc^4 - 6.*Nc^2 + 18.)/Nc^2*Nf^2
  
  b1 <- b1/b0/4./pi
  b2 <- b2/b0/16./pi^2
  b3 <- b3/b0/64./pi^3

  L2 <- log(mu^2/lam0^2)
  LL2 <- log(L2)

  als0 <- 4.*pi/b0/L2
  als1 <- als0 - als0^2*b1*LL2
  als2 <- als1 + als0^3*(b1^2*(LL2^2 - LL2 -1.) + b2)
  als3 <- als2 + als0^4*(b1^3*(-LL2^3+5./2.*LL2^2+2*LL2-1./2.)-
                         3.*b1*b2*LL2 + b3/2.)
  if(nl == 0) return(als0/(4.*pi))
  if(nl == 1) return(als1/(4.*pi))
  if(nl == 2) return(als2/(4.*pi))
  return(als3/(4.*pi))
}
