phaseshift.rho <- function(pcfit, L, Mpi, frame="cmf", Mpiboot, disp="cont", n=1) {

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

  
  if(frame == "cmf") {
    tandelta <- pi^(3./2.)*sqrt(qtilde$qtsq)/Z00
    tandeltaboot <- pi^(3./2.)*sqrt(qtildeboot$qtsq)/Z00boot
    delta <- atan2(pi^(3./2.)*sqrt(qtilde$qtsq),Z00)
    deltaboot <- atan2(pi^(3./2.)*sqrt(qtildeboot$qtsq),Z00boot)
  }
  else if(frame == "mf1") {
    ## representation A_1 from arXiv:1206:4141 w00 + 2 w20
    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    delta <- atan2(qtilde$gamma*pi**(3./2.) * sqrt(qtilde$qtsq) , (Z00 + (2./(qtilde$qtsq*sqrt(5)))*Z20))
    tandelta <- qtilde$gamma*pi**(3./2.) * sqrt(qtilde$qtsq) / (Z00 + (2./(qtilde$qtsq*sqrt(5)))*Z20)
    
    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    deltaboot <- atan2(qtildeboot$gamma*pi**(3./2.) * sqrt(qtildeboot$qtsq) , (Z00boot + (2./(qtildeboot$qtsq*sqrt(5)))*Z20boot))
    tandeltaboot <- qtildeboot$gamma*pi**(3./2.) * sqrt(qtildeboot$qtsq) / (Z00boot + (2./(qtildeboot$qtsq*sqrt(5)))*Z20boot)
  }
  else if(frame == "mf2") {
    ## representation A_1 from arXiv:1206:4141 w00 - w20 - i \sqrt(6) w22 ??
    Z20 <- Re(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    Z22  <- Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    delta <- atan2(qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq) , (Z00 - (1./(qtilde$qtsq*sqrt(5)))*Z20 + ((sqrt(6./5.)/(qtilde$qtsq))*Z22)))
    tandelta <- qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq) / (Z00 - (1./(qtilde$qtsq*sqrt(5)))*Z20 + ((sqrt(6./5.)/(qtilde$qtsq))*Z22))
    
    Z20boot <- Re(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2))
    Z22boot  <- Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    deltaboot <- atan2(qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq) , (Z00boot - (1./(qtildeboot$qtsq*sqrt(5)))*Z20boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot)))
    tandeltaboot <- qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq) / (Z00boot - (1./(qtildeboot$qtsq*sqrt(5)))*Z20boot + ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
  }
  else if(frame == "mf3") {
    ## representation A from arXiv:1206:4141 w00 -2i \sqrt(6) w22 ??
    Z22 <- Im(LuescherZeta(qtilde$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    delta <- atan2(qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq) , (Z00 +2* ((sqrt(6./5.)/(qtilde$qtsq))*Z22)))
    tandelta <- qtilde$gamma*pi^(3./2.) * sqrt(qtilde$qtsq) / (Z00 +2* ((sqrt(6./5.)/(qtilde$qtsq))*Z22))
    
    Z22boot  <- Im(LuescherZeta(qtildeboot$qtsq, gamma=qtilde$gamma, dvec = Pcm, l = 2, m = 2))
    deltaboot <- atan2(qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq) , (Z00boot +2* ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot)))
    tandeltaboot <- qtildeboot$gamma*pi^(3./2.) * sqrt(qtildeboot$qtsq) / (Z00boot +2* ((sqrt(6./5.)/(qtildeboot$qtsq))*Z22boot))
  }
  else {
    stop(paste("value of frame ", frame," not recognised\n", sep=""))
  }
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
