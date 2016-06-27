## the ratio defined in Jansen, Renner, Xu
## of shifted 4 and 2pt functions
compRpipi <- function(c4, c2, Thalf) {
  ## we do not go until Thalf+1, because we take the difference and
  ## it would not be defined
  tt <- c(1:Thalf)
  return((c4[tt] - c4[tt+1])/(c2[tt]^2 - c2[tt+1]^2))
}

## in case c4 is already a principal correlator from a weighted and
## shifted matrix
compRpipi2 <- function(c4, c2, Thalf) {
  ## we do not go until Thalf+1, because we take the difference and
  ## it would not be defined
  tt <- c(1:Thalf)
  ## c4 is from a shifted and weighted gevp -> need it at tt+1
  return(c4[tt+1]/(c2[tt]^2 - c2[tt+1]^2))
}

## in case c4 is already a principal correlator from a weighted and
## shifted matrix
compRpipi3 <- function(c4, c21, c22, Thalf, dE) {
  ## we do not go until Thalf+1, because we take the difference and
  ## it would not be defined
  tt <- c(1:Thalf)
  ## c4 is from a shifted and weighted gevp -> need it at tt+1
  return(c4[tt+1]/(c21[tt]*c22[tt]*exp(dE*tt) - c21[tt+1]*c22[tt+1]*exp(dE*(tt+1))))
}


phaseshift.rho <- function(pcfit, L, Mpi, frame="cmf", irrep="A1", Mpiboot, disp="cont", n=1) {

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
  }
  else if(frame == "mf1") {
    Pcm <- n*c(0, 0, 1)
  }
  else if(frame == "mf2") {
    Pcm <- n*c(0,1,1)
  }
  else if(frame == "mf3") {
    Pcm <- n*c(1,1,1)
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
  
  if(disp != "lat") {
    qtilde <- compute.qtildesq.contdisp(E = E, mpi = Mpi, dvec = Pcm, L = L)
    qtildeboot <- compute.qtildesq.contdisp(E = Eboot, mpi = Mpiboot, dvec = Pcm, L = L)
  }
  else {
    qtilde <- compute.qtildesq(E = E, mpi = Mpi, dvec = Pcm, L = L)
    qtildeboot <- compute.qtildesq(E = Eboot, mpi = Mpiboot, dvec = Pcm, L = L)
  }

  Z00 <- Re(LuescherZeta(qtilde$qtsq, gamma = qtilde$gamma, dvec = Pcm))
  Z00boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma = qtilde$gamma, dvec = Pcm))
  y <- qtilde$gamma*pi**(3./2.) * sqrt(qtilde$qtsq)
  yboot <- qtildeboot$gamma*pi**(3./2.) * sqrt(qtildeboot$qtsq)

  x <- numeric()
  xboot <- numeric()
  if(frame == "cmf") {
    x <- Z00
    xboot <- Z00boot
  }
  else if(frame == "mf1") {
    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 0))
    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 0))
    
    if(irrep == "A1") {
      ## arXiv:1212:0830v2: A_1 irrep [00n] page 17
      x <- (Z00 + (2./(qtilde$qtsq*sqrt(5)))*Z20)
      xboot <- (Z00boot + (2./(qtildeboot$qtsq*sqrt(5)))*Z20boot)
    }
    else if(irrep == "E2") {
      ## arXiv:1212:0830v2: E_2 irrep [00n]
      x <- (Z00 - (1./(qtilde$qtsq*sqrt(5)))*Z20)
      xboot <- (Z00boot - (1./(qtildeboot$qtsq*sqrt(5)))*Z20boot)
    }
  }
  else if(frame == "mf2"( {

    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 0))
    Z22  <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))

    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    Z22boot  <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))

    if(irrep == "B2") {
      ## arXiv:1212:0830v2: B_2 irrep [0nn]
      x <- (Z00 - (1./(qtilde$qtsq*sqrt(5)))*Z20 + ((sqrt(6./5.)/(qtilde$qtsq))*Z22))
      xboot <- (Z00boot - (1./(qtildeboot$qtsq*sqrt(5)))*Z20boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
    }
    else {
      Z21  <- -Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 1))
      Z21boot  <- -Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 1))
      if(irrep == "A1") {
        ## arXiv:1212:0830v2: A_1 irrep [0nn]
        x <- (Z00 + (1./(2*qtilde$qtsq*sqrt(5)))*Z20 + ((sqrt(6./5.)/(qtilde$qtsq))*Z21) - ((sqrt(3./10.)/(qtilde$qtsq))*Z22))
        xboot <- (Z00boot + (1./(2*qtildeboot$qtsq*sqrt(5)))*Z20boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z21boot) - ((sqrt(3./10.)/(qtildeboot$qtsq))*Z22boot))
      }
      else if(irrep == "B1") {
        ## arXiv:1212:0830v2: B_1 irrep [0nn]
        x <- (Z00 + (1./(2*qtilde$qtsq*sqrt(5)))*Z20 - ((sqrt(6./5.)/(qtilde$qtsq))*Z21) - ((sqrt(3./10.)/(qtilde$qtsq))*Z22))
        xboot <- (Z00boot + (1./(2*qtildeboot$qtsq*sqrt(5)))*Z20boot - ((sqrt(6./5.)/(qtildeboot$qtsq))*Z21boot) - ((sqrt(3./10.)/(qtildeboot$qtsq))*Z22boot))        
      }
    }
  }
  else if(frame == "mf3") {
    Z22 <- -Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    Z22boot  <- -Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    
    if(irrep == "E2") {
      ## arXiv:1212:0830v2: E_2 irrep [nnn] 
      x <- (Z00 + ((sqrt(6./5.)/(qtilde$qtsq))*Z22))
      xboot <- (Z00boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
    }
    else if(irrep == "A1") {
      ## arXiv:1212:0830v2: A_1 irrep [nnn] 
      Z21 <- LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 1)
      Z21boot  <- LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 1)

      x <- (Z00 - ((sqrt(8./15.)/(qtilde$qtsq))*(Z22 - Re(Z21) -Im(Z21) )))
      xboot <- (Z00boot - ((sqrt(8./15.)/(qtildeboot$qtsq))*(Z22boot - Re(Z21boot) -Im(Z21boot) )))
    }
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
  delta <- atan2(y,x)
  tandelta <- y/x
  shift <- 0.
  if(x < 0 && y >= 0) shift <- pi
  if(x < 0 && y < 0) shift <- -pi
  
  deltaboot <- atan(yboot/xboot) + shift
  tandeltaboot <- yboot/xboot  
  return(invisible(list(Ecm=qtilde$Ecm, Ecmboot=qtildeboot$Ecm, tandelta=tandelta, tandeltaboot=tandeltaboot, delta=delta, deltaboot=deltaboot)))
}

phaseshift.rho.old <- function(pcfit, L, Mpi, frame="cmf", Mpiboot, disp="cont", n=1) {

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
  Pcm <- c(0,0,0)
  if(frame == "cmf") {
    Pcm <- c(0,0,0)
  }
  else if(frame == "mf1") {
    Pcm <- n*c(0, 0, 1)
  }
  else if(frame == "mf2") {
    Pcm <- n*c(1,1,0)
  }
  else if(frame == "mf3") {
    Pcm <- n*c(1,1,1)
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
  
  if(disp != "lat") {
    qtilde <- compute.qtildesq.contdisp(E=E, mpi=Mpi, dvec=Pcm, L=L)
    qtildeboot <- compute.qtildesq.contdisp(E=Eboot, mpi=Mpiboot, dvec=Pcm, L=L)
  }
  else {
    qtilde <- compute.qtildesq(E=E, mpi=Mpi, dvec=Pcm, L=L)
    qtildeboot <- compute.qtildesq(E=Eboot, mpi=Mpiboot, dvec=Pcm, L=L)
  }

  Z00 <- Re(LuescherZeta(qtilde$qtsq, gamma = qtilde$gamma, dvec = Pcm))
  Z00boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma = qtilde$gamma, dvec = Pcm))
  x <- numeric()
  y <- numeric()
  xboot <- numeric()
  yboot <- numeric()
  if(frame == "cmf") {
    y <- pi^(3./2.)*sqrt(qtilde$qtsq)
    x <- Z00
    yboot <- pi^(3./2.)*sqrt(qtildeboot$qtsq)
    xboot <- Z00boot
  }
  else if(frame == "mf1") {
    ## representation A_1 from arXiv:1206:4141 w00 + 2 w20
    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    y <- qtilde$gamma*pi**(3./2.) * sqrt(qtilde$qtsq)
    x <- (Z00 + (2./(qtilde$qtsq*sqrt(5)))*Z20)
    
    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    yboot <- qtildeboot$gamma*pi**(3./2.) * sqrt(qtildeboot$qtsq)
    xboot <- (Z00boot + (2./(qtildeboot$qtsq*sqrt(5)))*Z20boot)
  }
  else if(frame == "mf2") {
    ## representation A_1 from arXiv:1206:4141 w00 - w20 - i \sqrt(6) w22 ??
    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    Z22  <- Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    y <- qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq)
    x <- (Z00 - (1./(qtilde$qtsq*sqrt(5)))*Z20 + ((sqrt(6./5.)/(qtilde$qtsq))*Z22))

    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    Z22boot  <- Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    yboot <- qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq)
    xboot <- (Z00boot - (1./(qtildeboot$qtsq*sqrt(5)))*Z20boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
  }
  else if(frame == "mf3") {
    ## representation A from arXiv:1206:4141 w00 -2i \sqrt(6) w22 ??
    Z22 <- Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    y <- qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq)
    x <- (Z00 +2* ((sqrt(6./5.)/(qtilde$qtsq))*Z22))
    
    Z22boot  <- Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    yboot <- qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq)
    xboot <- (Z00boot +2* ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
  delta <- atan2(y,x)
  tandelta <- y/x
  shift <- 0.
  if(x < 0 && y >= 0) shift <- pi
  if(x < 0 && y < 0) shift <- -pi
  
  deltaboot <- atan(yboot/xboot) + shift
  tandeltaboot <- yboot/xboot  
  return(invisible(list(Ecm=qtilde$Ecm, Ecmboot=qtildeboot$Ecm, tandelta=tandelta, tandeltaboot=tandeltaboot, delta=delta, deltaboot=deltaboot)))
}


tandeltaovEcm <- function(par, Ecm, Mpi) {
  return(par[1]^2/6/pi*sqrt(Ecm^2/4.-Mpi^2)^3/Ecm/(par[2]^2-Ecm^2))
}

ERchi <- function(par, y, dy, x, dx, Mpi) {
  return(c( (y-tandeltaovEcm(par[1:2], Ecm=par[3:length(par)], Mpi))/dy, (x-par[3:length(par)])/dx))
}

fit.rhoresonance <- function() {
  data <- read.table("data.dat")
  ##kk <- c(1:12)
  kk <- c(1:8,11,13,14, 15, 16)
  Mpi <- 0.14142
  MK <- 0.2567
  jj <- which(data$V1[kk] < 2*MK)
  ii <- kk[jj]
  par <- c(6,.4,data$V1[ii])
  lm.avail <- require(minpack.lm)
  require(tikzDevice)
  if(lm.avail) opt.res <- nls.lm(par = par, fn=ERchi, x=data$V1[ii], dx=data$V2[ii], y=data$V7[ii], dy=data$V8[ii], Mpi=Mpi)

  chisq <- opt.res$rsstrace[length(opt.res$rsstrace)]
  cat("chisq = ", chisq, "\n")
  cat("dof = ", length(ii)-2, "\n")
  cat("Mrho = ", opt.res$par[2], "\n")
  cat("g = ", opt.res$par[1], "\n")
  ##plotwitherror(x=data$V1[kk], dx=data$V2[kk], y=data$V7[kk], dy=data$V8[kk], ylim=c(-5,5), xlim=c(2*Mpi, 4*Mpi))
  x <- seq(2*Mpi,2*MK,0.001)
  y <- tandeltaovEcm(opt.res$par[1:2], x, Mpi=Mpi)
  ay <- atan(y)
  ay[which(ay < 0)] <- ay[which(ay < 0)] + pi
  ##lines(x,y)
  tikz(paste("delta1-A40", ".tex", sep=""), standAlone = TRUE, width=6, height=5)
  par(cex=1.3)

  colorlist <- c("blue", "red", "blue", "red", "darkgreen", "darkgreen", "blue", "red", "blue")
  pchlist <- c(15, 15, 16, 16, 15, 16, 17, 17, 18)
  llist <- c("A40.24 CM", "A40.20 CM", "A40.24 MF1", "A40.20 MF1", "A40.32 CM", "A40.32 MF1", "A40.24 MF2", "A40.20 MF2", "A40.24 MF3")
  ilist <- c(1:9)
  
  plot(type="n", x=data$V1[kk], y=data$V3[kk], ylim=c(0,pi), xlim=c(2*Mpi, 2*MK), xlab="$E_{\\mathrm{CM}}$", ylab="$\\delta_1$")
  lines(x,ay, col="red")
  for(i in ilist) {
    kk <- c(2*(i-1)+1,2*i)
    plotwitherror(x=data$V1[kk], dx=data$V2[kk], y=data$V3[kk], dy=data$V4[kk], rep=TRUE, pch=pchlist[i], col=colorlist[i], bg=colorlist[i])
  }
  legend("topleft", legend=llist[ilist], pch=pchlist[ilist], col=colorlist[ilist], pt.bg=colorlist[ilist], bty="n", cex=1.2)
  legend("bottomright", legend=c(paste("$g_{\\pi\\pi\\rho}=", format(opt.res$par[1], digits=3, scientific=FALSE), "$", sep=""), paste("$aM_\\rho=", format(opt.res$par[2], digits=3, scientific=FALSE), "$")), bty="n", cex=1.2)

  dev.off()
  tools::texi2dvi(paste("delta1-A40", ".tex", sep=""), pdf=T)

  return(opt.res)
}
