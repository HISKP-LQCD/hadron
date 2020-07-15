#' finite size corrections a la Colangelo, Duerr, Haefeli
#' 
#' finite size corrections a la Colangelo, Duerr, Haefeli
#' 
#' see reference for details. We use the simplyfied formulae for the S
#' quantities, see eq. (59) in the reference.
#' 
#' @param parm parameters
#' @param rev \eqn{rev=-1} corrects from \eqn{L} to \eqn{L=\infty}{L =
#' infinity}, \eqn{rev=+1} the other way around
#' @param aLamb1 The four low energy
#' \eqn{\Lambda_{1-4}}{Lambda1-Lambda4}constants in lattice units.
#' @param aLamb2 see \code{aLamb1}.
#' @param aLamb3 see \code{aLamb1}.
#' @param aLamb4 see \code{aLamb1}.
#' @param ampiV pseudo scalar mass values to be corrected
#' @param afpiV pseudo scalar decay constant values to be corrected
#' @param aF0 \eqn{af_0}{af0} in lattice units
#' @param a_fm the value of the lattice spacing in fermi
#' @param L the lattice spatial extent
#' @param printit if set to TRUE the corrections are printed
#' @param incim6 in- or exclude the NNNLO correction for the mass
#' @param rtilde the low energy constants \eqn{\tilde{r}}{rtilde}, needed only
#' if \code{incim6=TRUE}
#' @param use.cimpl use the four times faster direct c Implementation of the
#' correction routine
#' @return a list with the corrected values for mpi and fpi
#' @author Carsten Urbach <curbach@@gmx.de>
#' @references Gilberto Colangelo, Stephan Durr, Christoph Haefeli,
#' Nucl.Phys.B721:136-174,2005. hep-lat/0503014
#' @examples
#' 
#' L <- c(24, 24, 24, 24, 32)
#' mps <- c(0.14448, 0.17261, 0.19858, 0.22276, 0.14320)
#' fps <- c(0.06577, 0.07169, 0.07623, 0.07924, 0.06730)
#' aLamb1 <- 0.05
#' aLamb2 <- 0.5
#' aLamb3 <- 0.38
#' aLamb4 <- 0.66
#' cdhres <- cdh(rev=+1, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3, aLamb4=aLamb4,
#'               ampiV=mps, afpiV=fps, aF0=fps, a_fm=0.08, L=L, printit=TRUE,
#'               incim6=FALSE)
#' cdhres$mpiFV
#' cdhres$fpiFV
#' 
#' @export cdh
cdh <- function(parm = rep(0, times=6), rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
                     ampiV, afpiV, aF0, a_fm, L, printit=FALSE, incim6 = FALSE,
                     rtilde=c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0), use.cimpl=TRUE) {

  if(any(ampiV <= 0.)) {
    warning("ampiV must not be negative!\n")
    return(invisible(list(mpiFV=rep(NA, times=length(ampiV)), fpiFV=rep(NA, times=length(ampiV)))))
  }

  if(!use.cimpl) {
    return(invisible(cdh.R(parm=parm, rev=rev, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3, aLamb4=aLamb4,
                           ampiV=ampiV, afpiV=afpiV, aF0=aF0, a_fm=a_fm, L=L, printit=printit, incim6=incim6,
                           rtilde=rtilde)))
  }
  res <- .Call("cdh_c", rev, aLamb1, aLamb2, aLamb3, aLamb4, aF0, a_fm, L, ampiV, afpiV, as.integer(printit), rtilde, as.integer(incim6))
  return(invisible(list(mpiFV=res[1:length(ampiV)], fpiFV=res[(length(ampiV)+1):(2*length(ampiV))])))
}



#' finite size corrections a la Colangelo, Duerr, Haefeli, but re-expanded as
#' series in the quark mass
#' 
#' finite size corrections a la Colangelo, Duerr, Haefeli, but re-expanded as
#' series in the quark mass
#' 
#' see reference for details. We use the simplyfied formulae for the S
#' quantities, see eq. (59) in first reference.
#' 
#' @param rev \eqn{rev=-1} corrects from \eqn{L} to \eqn{L=\infty}{L =
#' infinity}, \eqn{rev=+1} the other way around
#' @param aLamb1 The four low energy
#' \eqn{\Lambda_{1-4}}{Lambda1-Lambda4}constants in lattice units.
#' @param aLamb2 see \code{aLamb1}.
#' @param aLamb3 see \code{aLamb1}.
#' @param aLamb4 see \code{aLamb1}.
#' @param ampiV pseudo scalar mass values to be corrected
#' @param afpiV pseudo scalar decay constant values to be corrected
#' @param aF0 \eqn{af_0}{af0} in lattice units
#' @param a2B0mu \eqn{2B_0\mu}{2 B0 mu} in lattice units, where \eqn{\mu}{mu}
#' is the quark mass and \eqn{B_0}{B0} a low energy constant
#' @param L the lattice spatial extent
#' @param printit if set to TRUE the corrections are printed
#' @param use.cimpl use the four times faster direct c Implementation of the
#' correction routine
#' @param parm m parameters
#' @return a list with the corrected values for mpi and fpi
#' @author Carsten Urbach <curbach@@gmx.de>
#' @references Gilberto Colangelo, Stephan Durr, Christoph Haefeli,
#' Nucl.Phys.B721:136-174,2005. hep-lat/0503014
#' 
#' and
#' 
#' R. Frezzotti, V. Lubicz, S. Simula, arXiv:0812.4042 hep-lat
#' @examples
#' 
#' mu <- c(0.004, 0.006, 0.008, 0.010, 0.004)
#' L <- c(24, 24, 24, 24, 32)
#' mps <- c(0.14448, 0.17261, 0.19858, 0.22276, 0.14320)
#' fps <- c(0.06577, 0.07169, 0.07623, 0.07924, 0.06730)
#' aLamb1 <- 0.05
#' aLamb2 <- 0.5
#' aLamb3 <- 0.38
#' aLamb4 <- 0.66
#' aF0    <- 0.051
#' a2B    <- 5.64
#' cdhres <- cdhnew(rev=+1, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3,
#'                  aLamb4=aLamb4, ampiV=mps, afpiV=fps, aF0=aF0,
#'                  a2B0mu=a2B*mu, L=L, printit=TRUE)
#' cdhres$mpiFV
#' cdhres$fpiFV
#' 
#' @export cdhnew
cdhnew <- function(parm = rep(0, times=6), rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
                     ampiV, afpiV, aF0, a2B0mu, L, printit=FALSE, use.cimpl = TRUE) {

  if(!use.cimpl) {
    return(invisible(cdhnew.R(rev=rev, aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aLamb3, aLamb4=aLamb4,
                              ampiV=ampiV, afpiV=afpiV, aF0=aF0, a2B0mu=a2B0mu, L=L, printit=printit)))
  }
  if(any(ampiV <= 0.)) {
    warning("ampiV must not be negative!\n")
    return(invisible(list(mpiFV=rep(NA, times=length(ampiV)), fpiFV=rep(NA, times=length(ampiV)))))
  }

  res <- .Call("cdhnew_c", rev, aLamb1, aLamb2, aLamb3, aLamb4, aF0, a2B0mu, L, ampiV, afpiV, as.integer(printit))
  return(invisible(list(mpiFV=res[1:length(ampiV)], fpiFV=res[(length(ampiV)+1):(2*length(ampiV))])))
}

cdh.R <-  function(parm, rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
                 ampiV, afpiV, aF0, a_fm, L, printit=FALSE, incim6 = TRUE,
                 rtilde=c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0)) {

  # aF0 at beta=3.9 =0.0534 for correction a la GL 
  # aF0 = af_ps(mu) for higher orders (=F_pi)
  # aLamb1=0.055 and aLamb2=0.58

  if(length(rtilde)!=6) {
    rtilde=c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0)
    warning("rtilde had not enough entries, using defaults instead")
  }
  ## our normalisation
  aF0 <- aF0/sqrt(2)

  ## Eq.(62)
  gg <- c(2-pi/2, pi/4-0.5, 0.5-pi/8, 3*pi/16 - 0.5)
  ## Tab.1
  mm <- c(6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24, 24)
  mml <- length(mm)
  N <- 16*pi^2
  ## physical mass of the rho in lattice units
  amrho_phys <- a_fm*770/197.3
  ## related to those of tab.2a by Eq.(53)
  lb1 <- log((aLamb1/ampiV)^2)
  lb2 <- log((aLamb2/ampiV)^2)
  lb3 <- log((aLamb3/ampiV)^2)
  lb4 <- log((aLamb4/ampiV)^2)
  lpi <- log((ampiV/amrho_phys)^2)
  ## tab 2b -> rtilde
  
  if(missing(parm)) {
    parm <- rep(0, times=6)
  }

  M_P <- ampiV
  F_P <- afpiV
  ## Eq.(10)
  xi_P <- (M_P/(4*pi*aF0))^2
  mmB0 <- rep(0., times=length(ampiV))
  mmB2 <- mmB0
  for(jj in 1:length(ampiV)) {
    ## Eq.(11)
    lambda_pi <-  ampiV[jj]*L[jj]
    ## argument of functions in Eq.(50). sqrt(n) comes from Eq.(26-27)
    z <-  sqrt(c(1:mml))*lambda_pi
    B0 <- 2*besselK(z,1)/z
    B2 <- 2*besselK(z,2)/z^2
    ## remaining factor from Eq.(26-27) and sum, ...
    mmB0[jj] <- sum(mm*B0)
    ## which I can already do since all the dependence on n is here
    mmB2[jj] <- sum(mm*B2)
  }
  ## simplifyed S's: Eq.(59)
  S4mpi <- (13/3)*gg[1] * mmB0 - (1/3)*(40*gg[1] + 32*gg[2] + 26*gg[3])*mmB2
  S6mpi <- 0
  S4fpi <- (1/6)*(8*gg[1] - 13*gg[2])* mmB0 - (1/3)*(40*gg[1] - 12*gg[2] - 8*gg[3] - 13*gg[4])*mmB2

  ## Eq. (49) and (54)
  I2mpi <- -mmB0
  I4mpi <- mmB0*(-55/18 + 4*lb1 + 8/3*lb2 - 5/2*lb3 -2*lb4) +
    mmB2*(112/9 - (8/3)*lb1 - (32/3)*lb2) + S4mpi
  if(incim6) {
    I6mpi <- mmB0*(10049/1296 - 13/72*N + 20/9*lb1 - 40/27*lb2 - 3/4*lb3 - 110/9*lb4
                   - 5/2*lb3^2 - 5*lb4^2 
                   + lb4*(16*lb1 + 32/3*lb2 - 11*lb3)
                   + lpi*(70/9*lpi + 12*lb1 + 32/9*lb2 - lb3 + lb4 + 47/18)
                   + 5*rtilde[1] + 4*rtilde[2] + 8*rtilde[3] + 8*rtilde[4] + 16*rtilde[5] + 16*rtilde[6]) +
                     mmB2*(3476/81 - 77/288*N + 32/9*lb1 + 464/27*lb2 + 448/9*lb4 
                           - 32/3*lb4*(lb1+4*lb2) + lpi*(100/9*lpi + 8/3*lb1 + 176/9*lb2 - 248/9)
                           - 8*rtilde[3] - 56*rtilde[4] - 48*rtilde[5] + 16*rtilde[6])+
                             S6mpi
  }
  else {
    I6mpi <- 0
  }

  I2fpi <- -2*mmB0;
  I4fpi <- mmB0*(-7/9 + 2*lb1 + (4/3)*lb2 - 3*lb4) + 
    mmB2*(112/9 - (8/3)*lb1 -(32/3)*lb2) + 
      S4fpi
  I6fpi <- 0
  ## Eq. (26-27). The sum over n is already done
  Rmpi <- - (xi_P/2) * (I2mpi + xi_P * I4mpi + xi_P^2 * I6mpi)
  Rfpi <- (xi_P)   * (I2fpi + xi_P * I4fpi + xi_P^2 * I6fpi);

  mpiFV <- ampiV * (1 + rev* Rmpi);
  fpiFV <- afpiV * (1 + rev* Rfpi);


  ## print out some further information (i.e. the contribuition of each order)
  if (0) {
    mlo <- - (xi_P/2) * (I2mpi)
    mnlo <- - (xi_P/2) * (xi_P * I4mpi)
    mnnlo <- - (xi_P/2) * (xi_P^2 * I6mpi)
    flo <- (xi_P) * (I2fpi)
    fnlo= (xi_P) * (xi_P * I4fpi)
    message("mpi: ", mlo, " ", mnlo, " ", mnnlo, "\nfpi: ", flo, " ", fnlo, "\n")
  }
  if(printit) {
    message("Rmpi:", Rmpi, "\n")
    message("Rfpi:", Rfpi, "\n")
  }
  
  if(any(is.na(c(mpiFV, fpiFV)))) {
    warning("NaNs produced in: cdh!\n")
    return(invisible(list(mpiFV=ampiV, fpiFV=afpiV)))
  }
  return(invisible(list(mpiFV=mpiFV, fpiFV=fpiFV)))
}

#function [mpiFV fpiFV report] = CDH(ampiV,afpiV,aF0,aLamb1,aLamb2,aLamb3,aLamb4,a_fm,L,print,rev,parm)
#% When rev=1 (-1), compute the finite (infinite) volume am_pi and af_pi from the infinite (finite) volume ones.
#% rev=1 corresponds to the formulae in hep-lat/0503014. Equations and Table numbers refer to that paper.
#% The expansion parameter is 1/aF0, which must be set to 1/af_pi if these formula are to be used beyond LO.
#% It is possible to change the default value of the parameters \tilda{r}_i (i=1,6) with the variable param. 
#% The variables mpiFV fpiFV ampiV,afpiV,aF0, aLamb_i are in lattice units. 
#% "a_fm" is the lattice spacing in fm and it is used only where is necessary (not at LO).
#% L is the number of points in one spatial direction (we assume the spatial volume V=L^3)
#% Eq. (49) and (54)

cdhnew.R <- function(rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4,
                     ampiV, afpiV, aF0, a2B0mu, L, printit=FALSE)
{

  ## our normalisation need this?
  aF0 <- aF0/sqrt(2)

  ## Eq.(62)
  gg <- c(2-pi/2, pi/4-0.5, 0.5-pi/8, 3*pi/16 - 0.5)
  ## Tab.1
  mm <- c(6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24, 24)
  mml <- length(mm)
  N <- 16*pi^2

  ## related to those of tab.2a by Eq.(53)
  lb1 <- log(aLamb1^2/a2B0mu)
  lb2 <- log(aLamb2^2/a2B0mu)
  lb3 <- log(aLamb3^2/a2B0mu)
  lb4 <- log(aLamb4^2/a2B0mu)

  ## Eq.(10)
  xi_P <- a2B0mu/(4*pi*aF0)^2
  mmB0 <- rep(0., times=length(ampiV))
  mmB1 <- mmB0
  mmB2 <- mmB0
  for(jj in 1:length(ampiV)) {
    ## Eq.(11)
    lambda <-  sqrt(c(1:mml))*sqrt(a2B0mu[jj])*L[jj]
    B0 <- 2*besselK(lambda,0)
    ## moved original B0 -> B1
    B1 <- 2*besselK(lambda,1)/lambda
    B2 <- 2*besselK(lambda,2)/lambda^2
    ## remaining factor from Eq.(26-27) and sum, ...
    mmB0[jj] <- sum(mm*B0)
    mmB1[jj] <- sum(mm*B1)
    ## which I can already do since all the dependence on n is here
    mmB2[jj] <- sum(mm*B2)
  }

  ## DeltaM and DeltaF
  DeltaM <- -1/(2*N)*lb3
  DeltaF <- 2/N*lb4
  
  ## simplifyed S's: Eq.(59)
  S4mpi <- (13/3)*gg[1] * mmB1 - (1/3)*(40*gg[1] + 32*gg[2] + 26*gg[3])*mmB2
  S6mpi <- 0
  S4fpi=(1/6)*(8*gg[1] - 13*gg[2])* mmB1 - (1/3)*(40*gg[1] - 12*gg[2] - 8*gg[3] - 13*gg[4])*mmB2

  ## Eq. (49) and (54)
  I2mpi <- -mmB1
  I4mpi <- mmB1*(-55/18 + 4*lb1 + 8/3*lb2 - 5/2*lb3 -2*lb4) +
    mmB2*(112/9 - (8/3)*lb1 - (32/3)*lb2) + S4mpi +
      N/2*DeltaM*mmB0 + N*DeltaF*mmB1

  I2fpi <- -2*mmB1;
  I4fpi <- mmB1*(-7/9 + 2*lb1 + (4/3)*lb2 - 3*lb4) + 
    mmB2*(112/9 - (8/3)*lb1 -(32/3)*lb2) + 
      S4fpi + N*DeltaM*mmB0 + 2*N*DeltaF*mmB1

  ## Eq. (26-27). The sum over n is already done
  ##Rmpi <- - (xi_P/2) * (I2mpi + xi_P * I4mpi)
  ##Rfpi <- (xi_P)   * (I2fpi + xi_P * I4fpi)

  mpiFV <- ampiV * (1 + rev* ( - (xi_P/2) * (I2mpi + xi_P * I4mpi)));
  fpiFV <- afpiV * (1 + rev* ((xi_P)   * (I2fpi + xi_P * I4fpi)));


  ## print out some further information (i.e. the contribuition of each order)
  if (0) {
    mlo <- - (xi_P/2) * (I2mpi)
    mnlo <- - (xi_P/2) * (xi_P * I4mpi)
    flo <- (xi_P) * (I2fpi)
    fnlo= (xi_P) * (xi_P * I4fpi)
    message("mpi: ", mlo, " ", mnlo, "\n", "fpi: ", flo, " ", fnlo, "\n")
  }
  if(printit) {
    message("Rmpi:", ( - (xi_P/2) * (I2mpi + xi_P * I4mpi)), "\n")
    message("Rfpi:", ((xi_P)   * (I2fpi + xi_P * I4fpi)), "\n")
  }
  
  if(any(is.na(c(mpiFV, fpiFV)))) {
    warning("NaNs produced in: cdhnew!\n")
    return(invisible(list(mpiFV=ampiV, fpiFV=afpiV)))
  }
  return(invisible(list(mpiFV=mpiFV, fpiFV=fpiFV)))
}
