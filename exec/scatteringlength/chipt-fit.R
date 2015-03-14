require(tikzDevice)
lm.avail <- require(minpack.lm)

## x = Mpi^2/fpi^2
chipt.nlo.pipi <- function(x, ell) {
  return( -x/8./pi*(1+x/16./pi^2*(3*log(x) -1 -ell)) )
}

chipt.lo.pipi <- function(x) {
  return(-x/8./pi)
}

chi.chipt.nlo.pipi <- function(par, x, y, dx, dy) {
  ii <- c(2:length(par))
  return(
      c((y-chipt.nlo.pipi(x=par[ii]^2, ell=par[1]))/dy,
        (x-par[ii])/dx)
      )
}

types <- c("-efm", "-ratio", "")
types <- c("-ratio")
## this contains the data for Mpi/fpi not computed in this project
extdata <- read.table("mpi_fpi.dat")

enslist <- c("A30.32", "A40.32", "A40.24", "A40.20", "A60.24", "A80.24", "A100.24", "B55.32", "D45.32", "B35.32", "B85.24", "D15.48", "B35,48")
## indices of all fitted ensembles
## we do not include A40.20 and A40.24 as small volumes
jj <- c(1,2,5:7,8,9,10,11)
## ranges for A,B and D ensembles in the new data
jjA <- c(1:5)
jjB <- c(6,8,9)
jjD <- c(7)
## cut in Mps/fps
xcut <- 3.5
cargs <- commandArgs(trailingOnly = TRUE)
if(length(cargs) > 0) {
  xcut <- as.numeric(cargs[1])
  cat("setting xcut to", xcut, "from command line arguments\n")
}

for(j in c(1:length(types))) {
  ## the data (res) and bootstrap samples (bootres) for mpia0
  load(paste("res", types[j], "-allens.Rdata", sep=""))
  res <- res[jj,]
  a0 <- rep(0, times=length(res[,1]))

  R <- length(bootres[,1,1])
  bootres <- bootres[,jj,]

  ## generate the needed data from the external data file
  ## in correct order
  L <- c()
  x <- c()
  xboot <- array(0, dim=c(R, length(jj)))
  fsM <- c()
  fsMboot <- array(0, dim=c(R, length(jj)))
  dx <- c()
  Mps <- c()
  dfsM <- c()
  Mpsboot <- array(0, dim=c(R, length(jj)))
  for(k in c(1:length(jj))) {
    L[k] <- extdata[extdata$V1 == enslist[jj[k]],12]
    ## Mpi/fpi*Dfs
    x[k] <- extdata[extdata$V1 == enslist[jj[k]],6] / extdata[extdata$V1 == enslist[jj[k]],8] * extdata[extdata$V1 == enslist[jj[k]],10]
    xboot[,k] <- rnorm(R, mean=extdata[extdata$V1 == enslist[jj[k]],6], sd=extdata[extdata$V1 == enslist[jj[k]],7])/rnorm(R, mean=extdata[extdata$V1 == enslist[jj[k]],8], sd=extdata[extdata$V1 == enslist[jj[k]],9])*rnorm(R, mean=extdata[extdata$V1 == enslist[jj[k]],10], sd= extdata[extdata$V1 == enslist[jj[k]],11])
    xboot[1,k] <- x[k]
    fsM[k] <- 1./extdata[extdata$V1 == enslist[jj[k]],8]
    fsMboot[,k] <- 1./rnorm(R, mean=extdata[extdata$V1 == enslist[jj[k]],8], sd=extdata[extdata$V1 == enslist[jj[k]],9])
    fsMboot[1,k] <- fsM[k]
    Mps[k] <- extdata[extdata$V1 == enslist[jj[k]],2]
    dfsM[k] <- extdata[extdata$V1 == enslist[jj[k]],9]
    Mpsboot[,k] <- rnorm(R, mean=extdata[extdata$V1 == enslist[jj[k]],2], sd=extdata[extdata$V1 == enslist[jj[k]],3])
    Mpsboot[1,k] <- Mps[k]
    ## error of Mpi/fpi*Dfs (dominated by the error on Dfs)
    dx[k] <- sd(xboot[,k])
  }

  ## apply analytical FS corrections to a0
  for(i in c(1:length(jj))) {
    a0[i] <- fs.a0(a0=res[i,5], L=L[i], mps=Mps[i])
  }
  ## finite size correct for Mps
  res[,9] <- res[,9]*fsM/res[,5]*a0

  ## error on Mpia0
  ## include systematic error and error on the FS correction for Mps
  ## add in quadrature?
  dmpia0 <- c()
  for(i in c(1:length(jj))) {
    dmpia0[i] <- sqrt(fsM[i]^2*((res[i,10])^2+(max(abs(c(res[i,11],res[i,12]))))^2) + res[i,9]^2*dfsM[i]^2)
  }

  ## indices of points in fitrange
  iifit <- which(x < xcut)
  jpar <- c(1:(length(iifit)+1))
  Npar <- length(jpar)
  opt.bootres <- array(0, dim=c(R, Npar+1))

  cat("mpia0", res[iifit,9], "\n")
  cat("dmpia0", dmpia0[iifit], "\n")
  cat("Mps/fps", x[iifit], "\n")
  cat("dMps/fps", dx[iifit], "\n\n")

  
  ## parameters include one parameter for each Mpi/fpi value
  for(i in c(1:R)) {
    ## finite size corrections
    bootres[i,,3] <- bootres[i,,3]*fsMboot[i,]/res[,5]*a0
    par <- c(1, x[iifit])
    opt.res <- nls.lm(par, fn=chi.chipt.nlo.pipi, x=xboot[i,iifit], y=bootres[i,iifit,3], dx=dx[iifit], dy=dmpia0[iifit])
    chisq <- opt.res$rsstrace[length(opt.res$rsstrace)]
    opt.bootres[i, jpar] <- opt.res$par
    opt.bootres[i, Npar+1] <- chisq
  }

  xs <- seq(0.9,3.0,0.01)
  y <- array(0, dim=c(R, length(xs)))
  for(i in c(1:R)) {
    ## NLO ChiPT curve with fitted parameters
    y[i,] <- chipt.nlo.pipi(x=xs^2, ell=opt.bootres[i,1])
  }
  ## error band
  dy <- apply(y, 2, sd)
  ## LO ChiPT prediction
  ys2 <- chipt.lo.pipi(x=xs^2)

  save(opt.bootres, ensembles, file=paste("mpia0-chiptfit", types[j], ".Rdata", sep=""))
  hist(opt.bootres[,1])
  cat("\n")
  cat("fitted all Mpi/fps <", xcut, "\n")
  cat("l_pipi = ", opt.bootres[1,1], "+-", sd(opt.bootres[,1]), "\n")
  cat("chisq = ", opt.bootres[1,Npar+1], "\n")
  cat("dof = ", length(iifit)-1, "\n")
  cat("Qval = ", 1-pchisq(opt.bootres[1,Npar+1], length(iifit)-1), "\n")
  xphys = 0.135/0.1304
  cat("With neutral pion mass:\n")
  cat("Mpi a0_phys = ", chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1]), "+-", sd(chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[,1])), "\n")
  xphys = 0.13957/0.1304
  cat("With charged pion mass:\n")
  cat("Mpi a0_phys = ", chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1]), "+-", sd(chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[,1])), "\n")

  write.table(rbind(format(c(opt.bootres[1,1], sd(opt.bootres[,1]), chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1]), sd(chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[,1])), opt.bootres[1,Npar+1], length(iifit)-1, 1-pchisq(opt.bootres[1,Npar+1], length(iifit)-1), xcut), digits=3, scientific=FALSE)),
              sep=" & ", quote=FALSE, row.names=FALSE, 
              col.names=c('$M_{\\pi}a_0$', '$dM_{\\pi}a_0$', '$\\ell_{\\pi\\pi}$', '$d\\ell_{\\pi\\pi}$', '$\\chi^2$', '$\\mathrm{dof}$', 'Qval', 'xcut'),
              file=paste("chpt-res-", xcut, ".dat", sep=""), eol=" \\\\ \n")
  
  ## now prepare some plots
  for(p in c(1:3)) {
    tikz(paste("mpia0-fs-cpt", types[j], "-", p, "xcut", xcut, ".tex", sep=""), standAlone = TRUE, width=5, height=5)
    par(cex=1.3)
    plot(x=x, y=res[,9], type="n", xlim=c(0.9,3), ylim=c(-0.35,0), xlab=c("$M_{\\pi}/f_{\\pi}$"), ylab=c("$M_{\\pi}a_0$"))
    
    if(p==3) {
      polygon(x=c(xs, rev(xs)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
      lines(xs,y[1,], col="black")
    }
    lines(xs,ys2, col="black", lty=2)

    ## indicate the fit range
    lines(c(xcut,xcut), c(chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])-0.025, chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])+0.025), col="black" )
    lines(c(xcut,xcut-0.05), c(chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])-0.025, chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])-0.025), col="black")
    lines(c(xcut,xcut-0.05), c(chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])+0.025, chipt.nlo.pipi(x=xcut^2, ell=opt.bootres[i,1])+0.025), col="black")
    xphys = 0.13957/0.1304
    if(p==3) points(x=xphys, y=chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1]), pch=c(24), col="orangered", bg="orangered")
    else points(x=xphys, y=chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1]), pch=c(24), col="navy", bg="navy")
    ##arrows(x0=xphys, x1=xphys, y0=chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1])-sd(chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[,1])), y1=chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[1,1])+sd(chipt.nlo.pipi(x=xphys^2, ell=opt.bootres[,1])), length=0.01,angle=90,code=3, col=c("orange"))
    
    plotwitherror(x=x[jjA], dx=dx[jjA], y=res[jjA,9], dy=dmpia0[jjA], col=c("red"), bg=c("red"), pch=21, rep=TRUE)

    plotwitherror(x=x[jjB], dx=dx[jjB], y=res[jjB,9], dy=dmpia0[jjB], col=c("blue"), bg=c("blue"), pch=22, rep=TRUE)

    plotwitherror(x=x[jjD], dx=dx[jjD], y=res[jjD,9], dy=dmpia0[jjD], col=c("darkgreen"), bg=c("darkgreen"), pch=23, rep=TRUE)
    if(p == 2) {
      nf2data <- read.table("scat_length.dat")
      plotwitherror(x=nf2data$V1[1:6], y=nf2data$V2[1:6], dy=nf2data$V3[1:6], col=c("black"), pch=25, rep=TRUE)
    }
    

    if(p==3) {
      legend("topright", legend=c("LO ChiPT", "NLO ChiPT"), lty=c(2,1), bty="n", cex=1.3)
      legend("bottomleft", legend=c("A ensembles", "B ensembles", "D45", "extrapolated"), pt.bg=c("red", "blue", "darkgreen", "orangered"), col=c("red", "blue", "darkgreen", "orangered"), pch=c(21:24), bty="n", cex=1.3)
    }
    else {
      if(p==1) legend("bottomleft", legend=c("A ensembles", "B ensembles", "D45", "NA48/2 + Roy"), pt.bg=c("red", "blue", "darkgreen", "navy"), col=c("red", "blue", "darkgreen", "navy"), pch=c(21:24), bty="n", cex=1.3)
      else legend("bottomleft", legend=c("A ensembles", "B ensembles", "D45", "Nf=2", "NA48/2 + Roy"), pt.bg=c("red", "blue", "darkgreen", "white", "navy"), col=c("red", "blue", "darkgreen", "black", "navy"), pch=c(21:23,25,24), bty="n", cex=1.3)
      legend("topright", legend=c("LO ChiPT"), lty=c(2), bty="n", cex=1.3)
    }
    
    dev.off()
    tools::texi2dvi(paste("mpia0-fs-cpt", types[j], "-", p, "xcut", xcut, ".tex", sep=""), pdf=T)
  }
}
