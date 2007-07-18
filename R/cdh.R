cdh <-  function(parm, rev=-1, aLamb1=0.055, aLamb2=0.58, aLamb3, aLamb4, ampiV, afpiV, aF0, a_fm, L, printit=F) {

  # aF0 at beta=3.9 =0.0534 for correction a la GL 
  # aF0 = af_ps(mu) for higher orders (=F_pi)
  # aLamb1=0.055 and aLamb2=0.58

  # our normalisation
  aF0 <- aF0/sqrt(2)

  # Eg.(62)
  gg <- c(2-pi/2, pi/4-0.5, 0.5-pi/8, 3*pi/16 - 0.5)
  # Tab.1
  mm <- c(6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24, 24)
  mml <- length(mm)
  N <- 16*pi^2
  # physical mass of the rho in lattice units
  amrho_phys <- a_fm*770/197.3
  # related to those of tab.2a by Eq.(53)
  lb1 <- 2*log(aLamb1/ampiV)
  lb2 <- 2*log(aLamb2/ampiV)
  lb3 <- 2*log(aLamb3/ampiV)
  lb4 <- 2*log(aLamb4/ampiV)
  lpi <- 2*log(ampiV/amrho_phys)
  # tab 2b
  rtilde <- c(-1.5,  3.2, -4.2, -2.5,  3.8, 1.0)
  rtilder <- c(-1.5,  3.2, -4.2, -2.5,  1.0, 0.1)

  if(missing(parm)) {
    parm <- rep(0, times=6)
  }
  rtilde <- rtilde + parm*rtilder
  #transpose rtilde?
  #rtilde <- t(rtilde)

  M_P <- ampiV
  F_P <- afpiV
  # Eq.(10)
  xi_P <- (M_P/(4*pi*aF0))^2
  mmB0 <- rep(0., times=length(ampiV))
  mmB2 <- mmB0
  for(jj in 1:length(ampiV)) {
    # Eq.(11)
    lambda_pi <-  ampiV[jj]*L[jj]
    # argument of unctions in Eq.(50). sqrt(n) comes from Eq.(26-27)
    z <-  sqrt(c(1:mml))*lambda_pi
    B0 <- 2*besselK(z,1)/z
    B2 <- 2*besselK(z,2)/z^2
    # remaining factor from Eq.(26-27) and sum, ...
    mmB0[jj] <- sum(mm*B0)
    # which I can already do since all the dependence on n is here
    mmB2[jj] <- sum(mm*B2)
  }
  # simplifyed S's: Eq.(59)
  S4mpi <- (13/3)*gg[1] * mmB0 - (1/3)*(40*gg[1] + 32*gg[2] + 26*gg[3])*mmB2
  S6mpi <- 0
  S4fpi=(1/6)*(8*gg[1] - 13*gg[2])* mmB0 - (1/3)*(40*gg[1] - 12*gg[2] - 8*gg[3] - 13*gg[4])*mmB2

  # Eq. (49) and (54)
  I2mpi <- -mmB0
  I4mpi <- mmB0*(-55/18 + 4*lb1 + 8/3*lb2 - 5/2*lb3 -2*lb4) +
    mmB2*(112/9 - (8/3)*lb1 - (32/3)*lb2) + S4mpi
  I6mpi <- mmB0*(10049/1296 - 13/72*N + 20/9*lb1 - 40/27*lb2 - 3/4*lb3 - 110/9*lb4
                 - 5/2*lb3^2 - 5*lb4^2 
                 + lb4*(16*lb1 + 32/3*lb2 - 11*lb3)
                 + lpi*(70/9*lpi + 12*lb1 + 32/9*lb2 - lb3 + lb4 + 47/18)
                 + 5*rtilde[1] + 4*rtilde[2] + 8*rtilde[3] + 8*rtilde[4] + 16*rtilde[5] + 16*rtilde[6]) +
                   mmB2*(3476/81 - 77/288*N + 32/9*lb1 + 464/27*lb2 + 448/9*lb4 
                         - 32/3*lb4*(lb1+4*lb2) + lpi*(100/9*lpi + 8/3*lb1 + 176/9*lb2 - 248/9)
                         - 8*rtilde[3] - 56*rtilde[4] - 48*rtilde[5] + 16*rtilde[6])+
                           S6mpi
  I2fpi <- -2*mmB0;
  I4fpi <- mmB0*(-7/9 + 2*lb1 + (4/3)*lb2 - 3*lb4) + 
    mmB2*(112/9 - (8/3)*lb1 -(32/3)*lb2) + 
      S4fpi
  I6fpi <- 0

  # Eq. (26-27). The sum over n is already done
  Rmpi <- - (xi_P/2) * (ampiV/M_P) * (I2mpi + xi_P * I4mpi + xi_P^2 * I6mpi)
  Rfpi <- (xi_P)   * (afpiV/F_P) * (I2fpi + xi_P * I4fpi + xi_P^2 * I6fpi);

  mpiFV <- ampiV * (1 + rev* Rmpi);
  fpiFV <- afpiV * (1 + rev* Rfpi);


  # print out some further information (i.e. the contribuition of each order)
  if (printit) {
    mlo <- - (xi_P/2) * (ampiV/M_P) * (I2mpi)
    mnlo <- - (xi_P/2) * (ampiV/M_P) * (xi_P * I4mpi)
    mnnlo <- - (xi_P/2) * (ampiV/M_P) * (xi_P^2 * I6mpi)
    flo <- (xi_P) * (afpiV/F_P) * (I2fpi)
    fnlo= (xi_P) * (afpiV/F_P) * (xi_P * I4fpi)
    cat(mlo, mnlo, mnnlo, flo, fnlo, "\n")
                                        #    ['percent fin V corr: Mpi_LO   Mpi_NLO   Mpi_NNLO    Fpi_LO   Fpi_NLO']
                                        #[mlo', mnlo', mnnlo', flo', fnlo']
                                        #report=[(ampiV.*L)', ampiV', Rmpi', sqrt(2)*afpiV', Rfpi'];
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
