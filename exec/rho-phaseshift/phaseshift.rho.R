phaseshift.rho <- function(pcfit, L, Mpi, frame="cmf", Mpiboot) {

  if(missing(L)) {
    stop("L must be provided\n")
  }
  if(missing(pcfit)) {
    stop("pcfit must be provided!\n")
  }
  E <- c()
  Eboot  <- c()
  if(inherits(pcfit, "matrixfit")) {
    E <- pcfit$opt.res$par[1]
    Eboot <- pcfit$opt.tsboot[1,]
  }
  else if(inherits(pcfit, "effectivemassfit")) {
    E <- pcfit$opt.res$par[1]
    Eboot <- pcfit$massfit.tsboot[,1]
  }
  else {
    stop("pcfit is not of either type matrixfit or effectivemassfit\n")
  }
  if(missing(Mpiboot)) {
    Mpiboot <- Mpi
  }
  
  if(frame == "cmf") {
    Pcm <- c(0,0,0)
    gamma <- 1.
    Ecm <- E
    Ecmboot <- Eboot
    q <- (E-2.*Mpi)*L/(2.*pi)
    qboot <- (Eboot - 2*Mpiboot)*L/2./pi
    zeta <- Re(LuescherZeta(q^2, gamma = gamma, dvec = Pcm))
    zetaboot <- Re(LuescherZeta(qboot^2, gamma = gamma, dvec = Pcm))
    tandelta <- pi^(3./2.)*q/zeta
    tandeltaboot <- pi^(3./2.)*qboot/zetaboot
    delta <- atan2(pi^(3./2.)*q,zeta)
    deltaboot <- atan2(pi^(3./2.)*qboot,zetaboot)
  }
  else if(frame == "mf1") {
    Pcm <- c(0, 0, 1)

    Ecm <- sqrt(E^2-sum(Pcm^2)*4*pi^2/L^2)
    qsq <- (Ecm^2/4.-Mpi^2)*L^2/4/pi^2
    gamma <- E/Ecm
    Z00 <- Re(LuescherZeta(qsq, gamma=gamma, dvec = Pcm))
    Z20 <- Re(LuescherZeta(qsq, gamma=gamma, dvec = Pcm, l = 2))
    delta <- atan2(gamma*pi**(3./2.) * sqrt(qsq) , (Z00 + (2./(qsq*sqrt(5)))*Z20))
    tandelta <- gamma*pi**(3./2.) * sqrt(qsq) / (Z00 + (2./(qsq*sqrt(5)))*Z20)
    
    Ecmboot <- sqrt(Eboot^2-sum(Pcm^2)*4*pi^2/L^2)
    qsqboot <- (Ecmboot^2/4.-Mpiboot^2)*L^2/4/pi^2
    gammaboot <- Eboot/Ecmboot
    Z00boot <- Re(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm))
    Z20boot <- Re(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm, l = 2))
    deltaboot <- atan2(gammaboot*pi**(3./2.) * sqrt(qsqboot) , (Z00boot + (2./(qsqboot*sqrt(5)))*Z20boot))
    tandeltaboot <- gammaboot*pi**(3./2.) * sqrt(qsqboot) / (Z00boot + (2./(qsqboot*sqrt(5)))*Z20boot)
  }
  else if(frame == "mf2") {
    Ecm <- sqrt(E^2-sum(Pcm^2)*4*pi^2/L^2)
    qsq <- (Ecm^2/4.-Mpi^2)*L^2/4/pi^2
    gamma <- E/Ecm
    Z00 <- Re(LuescherZeta(qsq, gamma=gamma, dvec = Pcm))
    Z20 <- Re(LuescherZeta(qsq, gamma=gamma, dvec = Pcm, l = 2))
    Z22  <- Im(LuescherZeta(qsq, gamma=gamma, dvec = Pcm, l = 2, m = 2))
    Z2.2 <- Im(LuescherZeta(qsq, gamma=gamma, dvec = Pcm, l = 2, m = -2))
    delta <- atan2(gamma*pi^(3./2.) * sqrt(qsq) , (Z00 - (1./(qsq*sqrt(5)))*Z20 + ((sqrt(3./10.)/(qsq))*(Z22-Z2.2))))
    tandelta <- gamma*pi^(3./2.) * sqrt(qsq) / (Z00 - (1./(qsq*sqrt(5)))*Z20 + ((sqrt(3./10.)/(qsq))*(Z22-Z2.2)))
    
    Ecmboot <- sqrt(Eboot^2-sum(Pcm^2)*4*pi^2/L^2)
    qsqboot <- (Ecmboot^2/4.-Mpiboot^2)*L^2/4/pi^2
    gammaboot <- Eboot/Ecmboot
    Z00boot <- Re(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm))
    Z20boot <- Re(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm, l = 2))
    Z22boot  <- Im(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm, l = 2, m = 2))
    Z2.2boot <- Im(LuescherZeta(qsqboot, gamma=gamma, dvec = Pcm, l = 2, m = -2))
    deltaboot <- atan2(gamma*pi^(3./2.) * sqrt(qsqboot) , (Z00boot - (1./(qsqboot*sqrt(5)))*Z20boot + ((sqrt(3./10.)/(qsqboot))*(Z22boot-Z2.2boot))))
    tandelta <- gammaboot*pi^(3./2.) * sqrt(qsqboot) / (Z00boot - (1./(qsqboot*sqrt(5)))*Z20boot + ((sqrt(3./10.)/(qsqboot))*(Z22boot-Z2.2boot)))
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
  return(invisible(list(Ecm=Ecm, Ecmboot=Ecmboot, tandelta=tandelta, tandeltaboot=tandeltaboot, delta=delta, deltaboot=deltaboot)))
}

