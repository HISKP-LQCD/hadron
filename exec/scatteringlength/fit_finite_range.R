require(tikzDevice)

ERchi <- function(par, L, m, x, err) {
  return((-x-4*pi*par[1]/(m*L^3)*(1 - 2.837297*par[1]/L + 6.375183*par[1]^2/L^2 - 8.311951*par[1]^3/L^3) - (8*pi^2*par[1]^3*par[2])/(m*L^6))/err)
}

ERchi.qcotdelta <- function(par, x, y, xerr, yerr) {
  ii <- c((length(par)-length(x)+1):length(par))
  return(c((x-par[ii])/xerr, (y-par[1]-par[2]*par[ii])/yerr))
}

ERchisqr <- function(par, L, m, x, err) {
  return(sum((-x-4*pi*par[1]/(m*L^3)*(1 - 2.837297*par[1]/L + 6.375183*par[1]^2/L^2 - 8.311951*par[1]^3/L^3) - (8*pi^2*par[1]^3*par[2])/(m*L^6))^2/err^2))
}

dERchisqr <- function(par, L, m, x, err) {
  res <- rep(0., times=length(par))
  z <- (-x-4*pi*par[1]/(m*L^3)*(1 - 2.837297*par[1]/L + 6.375183*par[1]^2/L^2 - 8.311951*par[1]^3/L^3) - (8*pi^2*par[1]^3*par[2])/(m*L^6))
  zp1 <- -4*pi/(m*L^3)*(1 - 2.837297*par[1]/L + 6.375183*par[1]^2/L^2 - 8.311951*par[1]^3/L^3) - (8*pi^2*3*par[1]^2*par[2])/(m*L^6) +
    -4*pi*par[1]/(m*L^3)*(1 - 2.837297/L + 6.375183*2*par[1]/L^2 - 8.311951*3*par[1]^2/L^3)
  zp2 <- - (8*pi^2*par[1]^3)/(m*L^6)
  M <- diag(1./err^2)
  res[1] <- sum(zp1 %*% M %*% z + z %*%  M %*% zp1)
  res[2] <- sum(zp2 %*% M %*% z + z %*%  M %*% zp2)
  return(res)
}

fn <- function(par, L, m) {
  return( -4*pi*par[1]/(m*L^3)*(1 - 2.837297*par[1]/L + 6.375183*par[1]^2/L^2 - 8.311951*par[1]^3/L^3) - (8*pi^2*par[1]^3*par[2])/(m*L^6) )
}

compute.weights <- function(err, pvalues) {
  return(pvalues^2 * min(err, na.rm=TRUE)^2/err^2)
}

compute.boots <- function(res, index=1, piononly=FALSE) {
  pvalues <- (1-2*abs(res[1,,5]-0.5)) * (1-2*abs(res[1,,6]-0.5))
  if(piononly) pvalues <- (1-2*abs(res[1,,6]-0.5))
  err <- apply(res[,,index], 2, sd, na.rm=TRUE)
  w <- compute.weights(err, pvalues)
  return(apply(res[,,index], 1, weighted.quantile, prob=c(0.5), w=w, na.rm=TRUE))
}

fit.finite.range.qcotdelta <- function(path="./", type="") {
  load(paste(path, "A40.32/res", type, ".A40.32.Rdata", sep=""))
  resL32 <- res

  load(paste(path, "A40.24/res", type, ".A40.24.Rdata", sep=""))
  resL24 <- res

  load(paste(path, "A40.20/res", type, ".A40.20.Rdata", sep=""))
  resL20 <- res

  R <- length(resL32[,1,1])
  bs <- array(0, dim=c(R, 6))

  jj <- c(1:3)
  Mpi <- compute.boots(resL32, index=4, piononly=TRUE)
  ##cat(fs.qcotdelta(Mpi[1], L=32), fs.qcotdelta(Mpi[1], L=24), fs.qcotdelta(Mpi[1], L=20), "\n")
  ## qcotdelta
  bs[,1] <- compute.boots(resL32, index=7)-fs.qcotdelta(Mpi[1], L=32)
  bs[,2] <- compute.boots(resL24, index=7)-fs.qcotdelta(Mpi[1], L=24)
  bs[,3] <- compute.boots(resL20, index=7)-fs.qcotdelta(Mpi[1], L=20)
  ## q^2
  bs[,4] <- compute.boots(resL32, index=9)
  bs[,5] <- compute.boots(resL24, index=9)
  bs[,6] <- compute.boots(resL20, index=9)
  ##cat(bs[1,1], bs[1,2], bs[1,3], "\n")


  lm.avail <- require(minpack.lm)
  jj1 <- c(1:3)
  jj2 <- jj1+3

  err <- apply(bs, 2, sd)
  ##cat("x= ", bs[1,jj2], "\n")
  ##cat("dx=", err[jj2], "\n")
  ##cat("y= ", bs[1,jj1], "\n")
  ##cat("dy=", err[jj1], "\n")
  
  
  par=c(1,700, bs[1,4], bs[1,5], bs[1,6])
  opt.tsboot <- array(0, dim=c(R, length(par)+1))
  for(i in c(1:R)) {
    if(!lm.avail) {
      ##opt.res <- optim(par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS")#, control=list(ndeps=c(1.e-8, 1.e-8)))
      ##opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-8, times=length(par))))
      ##opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-10, times=length(par))))
    }
    else {
      opt.res <- nls.lm(par, fn=ERchi.qcotdelta, y=bs[i,jj1], yerr=err[jj1], x=bs[i,jj2], xerr=err[jj2])
    }
    opt.tsboot[i,c(1:length(par))] <- opt.res$par
    t <- ERchi.qcotdelta(par=opt.res$par, y=bs[i,jj1], yerr=err[jj1], x=bs[i,jj2], xerr=err[jj2])
    opt.tsboot[i,length(par)+1] <- sum(t^2)
  }

  ## now estimate the systematic uncertainty
  Ns <- 500
  il <- array(0, dim=c(Ns, 3))
  il[,1] <- sample.int(length(resL32[1,,1]), Ns, replace=TRUE)
  il[,2] <- sample.int(length(resL24[1,,1]), Ns, replace=TRUE)
  il[,3] <- sample.int(length(resL20[1,,1]), Ns, replace=TRUE)
  opt.sys <- array(0, dim=c(Ns, length(par)+1))
  for(i in c(1:Ns)) {
    tmp <- c(resL32[1,il[i,1],7], resL24[1,il[i,2],7], resL20[1,il[i,3],7], resL32[1,il[i,1],9], resL24[1,il[i,2],9], resL20[1,il[i,3],9])
    if(!lm.avail) {

    }
    else {
      opt.res <- nls.lm(par, fn=ERchi.qcotdelta, y=tmp[jj1], yerr=err[jj1], x=tmp[jj2], xerr=err[jj2])
    }
    opt.sys[i,c(1:length(par))] <- opt.res$par
    t <- ERchi.qcotdelta(par=opt.res$par, y=tmp[jj1], yerr=err[jj1], x=tmp[jj2], xerr=err[jj2])
    opt.sys[i,length(par)+1] <- sum(t^2)

  }

  cat("a0 =", 1./opt.tsboot[1,1], "+-", sd(1./opt.tsboot[,1]), "sys", abs(1./opt.tsboot[1,1])-abs(quantile(1./opt.sys[,1], probs=c(0.1573, 0.8427))), "\n")
  cat("Mpi a0 =", Mpi[1]/opt.tsboot[1,1], "+-", sd(Mpi/opt.tsboot[,1]), "sys", Mpi[1]*(abs(1./opt.tsboot[1,1])-abs(quantile(1./opt.sys[,1], probs=c(0.1573, 0.8427)))), "\n")
  cat("r =", 2*opt.tsboot[1,2], "+-", 2*sd(opt.tsboot[,2]), "sys", 2*(abs(opt.tsboot[1,2]) - abs(quantile(opt.sys[,2], probs=c(0.1573, 0.8427)))),"\n")
  cat("Mpi r =", 2*Mpi[1]*opt.tsboot[1,2], "+-", 2*sd(Mpi*opt.tsboot[,2]), "sys", 2*Mpi[1]*(abs(opt.tsboot[1,2]) - abs(quantile(opt.sys[,2], probs=c(0.1573, 0.8427)))),"\n")
  cat("Mpi^2 a0 r =", 2*Mpi[1]^2*opt.tsboot[1,2]/opt.tsboot[1,1], "+-", 2*sd(Mpi^2*opt.tsboot[,2]/opt.tsboot[,1]), "sys", 2*Mpi[1]^2*(abs(opt.tsboot[1,2]/opt.tsboot[1,1])-abs(quantile(opt.sys[,2]/opt.sys[,1], probs=c(0.1573, 0.8427)))), "\n")
  cat("chisq =", opt.tsboot[1,length(par)+1], "\n")
  
  tikz(paste("qcotdeltavqsq", type,".tex", sep=""), standAlone = TRUE, width=5, height=5)
  par(cex=1.3)
  x <- seq(0,0.003,0.0001)
  y <- array(0, dim=c(R, length(x)))
  for(i in c(1:R)) {
    y[i,] <- opt.tsboot[i,1] + opt.tsboot[i,2]*x
  }
  dy <- apply(y, 2, sd)
  
  plot(x=bs[1,jj2], y=bs[1,jj1], type="n", xlim=c(0,0.003), ylim=c(-1,-0.6), xlab=c("$a^2q^2$"), ylab=c("$aq\\cot\\delta_0$"))

  polygon(x=c(x, rev(x)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
  lines(x, y[1,], col="red")

  plotwitherror(x=bs[1,jj2], y=bs[1,jj1], dx=err[jj2], dy=err[jj1], rep=TRUE, col=c("blue"), bg=c("blue"), pch=c(21,22,23))

  points(0, opt.tsboot[1,1], pch=c(24), col=c("darkgreen"), bg=c("darkgreen"))

  ##arrows(x0=0, y0=quantile(opt.sys[,1], probs=c(0.1573)), x1=0, y1=quantile(opt.sys[,1], probs=c(0.8427)), length=0.01,angle=90,code=3, col=c("darkgreen"))

  legend("bottomright", pch=c(21,22,23,24), col=c("blue", "blue", "blue", "darkgreen"), bty="n", cex=1.1, pt.bg=c("blue", "blue", "blue", "darkgreen"), legend=c("A40.32", "A40.24", "A40.20", "extrapolated"))
  dev.off()
  tools::texi2dvi(paste("qcotdeltavqsq", type,".tex", sep=""), pdf=T)
  
  return(opt.tsboot)
}

fit.finite.range <- function(path="./", type="") {
  load(paste(path, "A40.32/res", type, ".A40.32.Rdata", sep=""))
  resL32 <- res

  load(paste(path, "A40.24/res", type, ".A40.24.Rdata", sep=""))
  resL24 <- res

  load(paste(path, "A40.20/res", type, ".A40.20.Rdata", sep=""))
  resL20 <- res

  load(paste("res", type, "-allens.Rdata", sep=""))
  x <- rbind(c(res[which(ensembles=="A40.32"),c(1:4)], 32), c(res[which(ensembles=="A40.24"),c(1:4)], 24), c(res[which(ensembles=="A40.20"),c(1:4)], 20))


  R <- length(resL32[,1,1])
  bs <- array(0, dim=c(R, 4))

  jj <- c(1:3)
  bs[,1] <- compute.boots(resL32)
  bs[,2] <- compute.boots(resL24)
  bs[,3] <- compute.boots(resL20)
  bs[,4] <- compute.boots(resL32, index=4, piononly=TRUE)
  ## can we use Levenberg Marquard from minpack?
  lm.avail <- require(minpack.lm)

  par=c(1,700)

  opt.tsboot <- array(0, dim=c(R, length(par)+1))
  for(i in c(1:R)) {
    if(!lm.avail) {
      opt.res <- optim(par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS")#, control=list(ndeps=c(1.e-8, 1.e-8)))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-8, times=length(par))))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-10, times=length(par))))
    }
    else {
      opt.res <- nls.lm(par, fn=ERchi, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj])
    }
    opt.tsboot[i,c(1:2)] <- opt.res$par
    opt.tsboot[i,3] <- ERchisqr(opt.res$par, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj])
  }

  ## now estimate the systematic uncertainty
  Ns <- 500
  il <- array(0, dim=c(Ns, 3))
  il[,1] <- sample.int(length(resL32[1,,1]), Ns, replace=TRUE)
  il[,2] <- sample.int(length(resL24[1,,1]), Ns, replace=TRUE)
  il[,3] <- sample.int(length(resL20[1,,1]), Ns, replace=TRUE)
  opt.sys <- array(0, dim=c(Ns, length(par)+1))
  for(i in c(1:Ns)) {
    tmp <- c(resL32[1,il[i,1],1], resL24[1,il[i,2],1], resL20[1,il[i,3],1])
    if(!lm.avail) {
      opt.res <- optim(par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS")#, control=list(ndeps=c(1.e-8, 1.e-8)))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-8, times=length(par))))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-10, times=length(par))))
    }
    else {
      opt.res <- nls.lm(par, fn=ERchi, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp)
    }
    opt.sys[i,c(1:2)] <- opt.res$par
    opt.sys[i,3] <- ERchisqr(opt.res$par, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp)
  }  
  
  cat("a0 =", opt.tsboot[1,1], "+-", sd(opt.tsboot[,1]), "sys", abs(opt.tsboot[1,1])-abs(quantile(opt.sys[,1], probs=c(0.1573, 0.8427))), "\n")
  cat("mpia0 =", opt.tsboot[1,1]*bs[1,4], "+-", sd(opt.tsboot[,1]*bs[,4]), "sys", bs[1,4]*(abs(opt.tsboot[1,1])-abs(quantile(opt.sys[,1], probs=c(0.1573, 0.8427)))), "\n")
  cat("r =", opt.tsboot[1,2], "+-", sd(opt.tsboot[,2]), "\n")
  cat("mpi r =", opt.tsboot[1,2]*bs[1,4], "+-", sd(opt.tsboot[,2]*bs[,4]), "\n")
  cat("mpi^2 a0 r =", opt.tsboot[1,1]*opt.tsboot[1,2]*bs[1,4], "+-", sd(opt.tsboot[,1]*opt.tsboot[,2]*bs[,4]), "\n")
  cat("chisq =", opt.tsboot[1,3], "\n")
  
  tikz(paste("deltaEovL", type,".tex", sep=""), standAlone = TRUE, width=5, height=5)
  par(cex=1.3)

  Ls <- seq(19,1024,1)
  y <- array(0, dim=c(R, length(Ls)))
  for(i in c(1:R)) {
    y[i,] <- fn(opt.tsboot[i,c(1:2)], L=Ls, m=bs[i,4])
  }
  dy <- apply(y, 2, sd)

  plot(x=1/x[,5], y=x[,1], type="n", xlim=c(0,0.06), ylim=c(0,0.02), xlab=c("$a/L$"), ylab=c("$a\\delta E$"))

  invL <- 1./Ls
  polygon(x=c(invL, rev(invL)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
  
  lines(invL, y[1,], col="red")
  
  plotwitherror(x=1/x[,5], y=x[,1], dy=x[,2], rep=TRUE, col="blue", bg="blue", pch=c(21,22,23))
  
  legend("topleft", pch=c(21,22,23), col=c("blue", "blue", "blue"), bty="n", cex=1.1, pt.bg=c("blue", "blue", "blue"), legend=c("A40.32", "A40.24", "A40.20"))
  
  dev.off()
  tools::texi2dvi(paste("deltaEovL", type,".tex", sep=""), pdf=T)
  

  jj <- c(1:2)
  opt.tsboot2 <- array(0, dim=c(R, length(par)+1))
  for(i in c(1:R)) {
    if(lm.avail) {
      opt.res <- nls.lm(par, fn=ERchi, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj])
    }
    else {
      opt.res <- optim(par, fn=ERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS")
      opt.res <- optim(opt.res$par, fn=ERchisqr, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj], method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-6, times=length(par))))
    }
    opt.tsboot2[i,c(1:2)] <- opt.res$par
    opt.tsboot2[i,3] <- ERchisqr(opt.res$par, L=x[jj,5], m=bs[i,4], err=x[jj,2], x=bs[i,jj])
  }

  opt.sys2 <- array(0, dim=c(Ns, length(par)+1))
  for(i in c(1:Ns)) {
    tmp <- c(resL32[1,il[i,1],1], resL24[1,il[i,2],1])
    if(!lm.avail) {
      opt.res <- optim(par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS")#, control=list(ndeps=c(1.e-8, 1.e-8)))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-8, times=length(par))))
      opt.res <- optim(opt.res$par, fn=ERchisqr, gr=dERchisqr, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp, method="BFGS", control=list(parscale=1./opt.res$par, ndeps=rep(1.e-10, times=length(par))))
    }
    else {
      opt.res <- nls.lm(par, fn=ERchi, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp)
    }
    opt.sys2[i,c(1:2)] <- opt.res$par
    opt.sys2[i,3] <- ERchisqr(opt.res$par, L=x[jj,5], m=bs[1,4], err=x[jj,2], x=tmp)
  }  

  cat("a0 =", opt.tsboot2[1,1], "+-", sd(opt.tsboot2[,1]), "sys", abs(opt.tsboot2[1,1])-abs(quantile(opt.sys2[,1], probs=c(0.1573, 0.8427))), "\n")
  cat("mpia0 =", opt.tsboot2[1,1]*bs[1,4], "+-", sd(opt.tsboot2[,1]*bs[,4]), "sys", bs[1,4]*(abs(opt.tsboot2[1,1])-abs(quantile(opt.sys2[,1], probs=c(0.1573, 0.8427)))), "\n")
  cat("r =", opt.tsboot2[1,2], "+-", sd(opt.tsboot2[,2]), "\n")
  cat("mpi r =", opt.tsboot2[1,2]*bs[1,4], "+-", sd(opt.tsboot2[,2]*bs[,4]), "\n")
  cat("mpi^2 a0 r =", opt.tsboot2[1,1]*opt.tsboot2[1,2]*bs[1,4], "+-", sd(opt.tsboot2[,1]*opt.tsboot2[,2]*bs[,4]), "\n")
  cat("chisq =", opt.tsboot2[1,3], "\n")
  
  tikz(paste("deltaEovL2", type,".tex", sep=""), standAlone = TRUE, width=5, height=5)
  par(cex=1.3)

  Ls <- seq(19,1024,1)
  y <- array(0, dim=c(R, length(Ls)))
  for(i in c(1:R)) {
    y[i,] <- fn(opt.tsboot2[i,c(1:2)], L=Ls, m=bs[i,4])
  }
  dy <- apply(y, 2, sd)

  plot(x=1/x[,5], y=x[,1], type="n", xlim=c(0,0.06), ylim=c(0,0.02), xlab=c("$a/L$"), ylab=c("$a\\delta E$"))

  invL <- 1./Ls
  polygon(x=c(invL, rev(invL)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
  
  lines(invL, y[1,], col="red")
  
  plotwitherror(x=1/x[,5], y=x[,1], dy=x[,2], rep=TRUE, col="blue", bg="blue", pch=c(21,22,23))

  legend("topleft", pch=c(21,22,23), col=c("blue", "blue", "blue"), bty="n", cex=1.1, pt.bg=c("blue", "blue", "blue"), legend=c("A40.32", "A40.24", "A40.20"))
  
  dev.off()
  tools::texi2dvi(paste("deltaEovL2", type,".tex", sep=""), pdf=T)
  
  return(list(opt.tsboot=opt.tsboot, opt.tsboot2=opt.tsboot2))
}
