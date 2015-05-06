source("phaseshift.R")

boot.R <- 1500
rr <- c(2:(boot.R+1))

Lvalues <- c(32, 24, 20)
## too coarse in the first part! q <- sqrt(c(seq(0.000001, .1, 0.00015), seq(0.100001, .15, 0.00025)))
q <- sqrt(c(seq(0.000001, .1, 0.00009), seq(0.100001, .15, 0.00025)))
lm.avail <- require(minpack.lm)
framelist <- c(1:4)
irreplist <- c(1:3)
scanfnlist <- list(c(detEqCMA1scan, detEqMF1A1scan, detEqMF2A1scan, detEqMF3A1scan), c(detEqCMEscan), c(detEqCMT2scan))
fnlist <- list(c(detEqCMA1, detEqMF1A1, detEqMF2A1, detEqMF3A1), c(detEqCME), c(detEqCMT2))

##par <- c(-0.77, 13, -0.011)
par <- c(-0.77, 14, -0.0075)
#par <- c(-0.77, 14, -0.02)

doSamples <- 10
start <- 1

## cut in chisq
cut <- 3.
lat.disp <- TRUE

if(file.exists("parameters.R")) {
  source("parameters.R")
}

load("Energies.Rdata")

qsq <- numeric()
## determine the rough q^2 values for applying a cut
for(i in c(1:length(Epipi$E))) {
  dvec=c(0,0,0)
  if(Epipi$frame[i] == 2) dvec=c(0,0,1)
  if(Epipi$frame[i] == 3) dvec=c(0,1,1)
  if(Epipi$frame[i] == 4) dvec=c(1,1,1)
  if(lat.disp) qtsq <- compute.qtildesq(Epipi$E[i], dvec=dvec, L=Epipi$L[i], mpi=Epi$E[i])
  else qtsq <- compute.qtildesq.contdisp(Epipi$E[i], dvec=dvec, L=Epipi$L[i], mpi=Epi$E[i])
  qsq[i] <- qtsq$q^2
}
Epipi$qsq <- qsq
Epi$qsq <- qsq

cat("applying cut at q^2/Mpi^2 =", cut, "\n")
k <- which(Epi$qsq/Epi$E^2 <= cut)
Epipi <- Epipi[k,]
Epi <- Epi[k,]
kk <- rep(0, times=2*length(k))
l <- c(1:length(k))
kk[2*l-1] <- 2*k-1
kk[2*l] <- 2*k
boot.data <- boot.data[,kk]
print(cbind(Epipi, Epi$E))


chisqr <- function(par, M, energydata, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist) {
  xv <- getEvalues(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist) - energydata
  return(xv %*% M %*% xv)
}

chi <- function(par, LM, energydata, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist) {
  return(LM %*% (getEvalues(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist) - energydata))
}

boot.par <- array(0., dim=c(boot.R+1,length(par)+1))

## we generate index arrays and covariance matrix once at the beginning...
energydata <- c()
iiEdata <- c()
r <- 1
## generate the omega_lm for the moving frames
## we store them in a list
Wlist <- list()
j <- 1
for(L in Lvalues) {
  for(f in framelist) {
    if(any(Epipi[Epipi$L==L,]$frame == f)) {
      jj <- which(Epi$frame==f & Epi$L==L & (Epi$irrep %in% irreplist))
      if(f == 1) Wlist[[j]] <- prepdetEqCMscan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
      if(f == 2) Wlist[[j]] <- prepdetEqMF1scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
      if(f == 3) Wlist[[j]] <- prepdetEqMF2scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
      if(f == 4) Wlist[[j]] <- prepdetEqMF3scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
    }
    else {
      Wlist[[j]] <- NULL
    }
    j <- j+1
  }
}
## save for later
Wlist1 <- Wlist
Mlist <- list()
nolist <- list()
k <- 1
m <- 1
for(L in Lvalues) {
  kk <- c()
  for(f in framelist) {
    for(i in irreplist) {
      if(any(Epipi[Epipi$L == L & Epipi$frame == f,]$irrep == i)) {
        jj <- which(Epi$frame == f & Epi$L == L & Epipi$irrep == i)
        nolist[[k]] <- max(Epipi[Epipi$frame == f & Epipi$L == L & Epipi$irrep == i,]$level)
        iiEdata <- c(iiEdata, jj)
        kk <- c(kk, jj)
      }
      else {
        nolist[[k]] <- 0
      }
      k <- k+1
    }
  }
  Mlist[[m]] <- try(invertCovMatrix(boot.data[,(2*kk-1)], boot.samples=TRUE))
  m <- m+1
}
rm(j,k)
bootevalues <- array(0, dim=c(boot.R+1, length(iiEdata)))
bootqcotdelta <- array(0, dim=dim(bootevalues))
bootq <- array(0, dim=dim(bootevalues))

## we build the covariance matrix from the single ones
M <- matrix(0., nrow=dim(Mlist[[1]])[1], ncol=dim(Mlist[[1]])[2])
if(length(Lvalues) == 2) M <- matrix(0., nrow=dim(Mlist[[1]])[1]+dim(Mlist[[2]])[1], ncol=dim(Mlist[[1]])[2]+dim(Mlist[[2]])[2])
if(length(Lvalues) == 3) M <- matrix(0., nrow=dim(Mlist[[1]])[1]+dim(Mlist[[2]])[1]+dim(Mlist[[3]])[1], ncol=dim(Mlist[[1]])[2]+dim(Mlist[[2]])[2]+dim(Mlist[[3]])[2])
nm1 <- 0
nm2 <- 0
for(i in c(1:length(Lvalues))) {
  nm1 <- nm2 + 1
  nm2 <- nm2 + dim(Mlist[[i]])[1]
  M[c(nm1:nm2), c(nm1:nm2)] <- Mlist[[i]]
}
##M <- try(invertCovMatrix(boot.data[,(2*iiEdata-1)], boot.samples=TRUE))
if(inherits(M, "try-error")) M <- diag(1./Epipi$dE[iiEdata]^2)
LM <- chol(M)

dof <- length(iiEdata) - length(par)

if(start > 1) {
  load(file=paste("fit-detEq-result-cut", cut, ".Rdata", sep=""))
}
## now we bootstrap everything
#for(r in c(1, rr)) {
for(r in c(start: doSamples)) {

  ## generate the omega_lm for the moving frames
  ## we store them in a list
  if(r > 1) {
    Wlist <- list()
    j <- 1
    for(L in Lvalues) {
      for(f in framelist) {
        if(any(Epipi[Epipi$L==L,]$frame == f)) {
          jj <- which(Epi$frame==f & Epi$L==L & (Epi$irrep %in% irreplist))
          if(f == 1) Wlist[[j]] <- prepdetEqCMscan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
          if(f == 2) Wlist[[j]] <- prepdetEqMF1scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
          if(f == 3) Wlist[[j]] <- prepdetEqMF2scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
          if(f == 4) Wlist[[j]] <- prepdetEqMF3scan(q=q, L=L, Mpi=boot.data[r,2*jj[1]])
        }
        else {
          Wlist[[j]] <- NULL
        }
        j <- j+1
      }
    }
  }
  
  ##lm.avail <- FALSE
  if(lm.avail) {
    opt.res <- nls.lm(par, fn=chi, upper = c(0, Inf, 0), lower=c(-Inf, 0, -Inf), control = nls.lm.control(ftol=1.e-10, ptol=1.e-8, nprint=10), LM=LM, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist, irreplist=irreplist)
    boot.par[r, length(par)+1] <- opt.res$rsstrace[length(opt.res$rsstrace)]
    if(opt.res$rsstrace[length(opt.res$rsstrace)] > 600) {
      opt.res <- optim(opt.res$par, chisqr, method = c("BFGS"), control = list(parscale=1./par, trace=10, ndeps=rep(1.e-8, times=length(par))), M=M, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist, irreplist=irreplist)
      boot.par[r, length(par)+1] <- opt.res$value
    }
  }
  if(!lm.avail || any(is.na(opt.res$par))) {
    opt.res <- optim(par, chisqr, method = c("BFGS"), control = list(parscale=1./par, trace=10, ndeps=rep(1.e-8, times=length(par))), M=M, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist, irreplist=irreplist)
##    opt.res <- optim(par, chisqr, control = list(parscale=1./par, trace=1, ndeps=rep(1.e-6, times=length(par))), M=M, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist, irreplist=irreplist)
    ##opt.res <- optim(opt.res$par, chisqr, method = c("BFGS"), control = list(parscale=1./opt.res$par, trace=10, ndeps=rep(1.e-8, times=length(par))), M=M, energydata=boot.data[r,(2*iiEdata-1)], Wlist=Wlist, Lvalues=Lvalues, nolist=nolist, framelist=framelist, scanfnlist=scanfnlist, fnlist=fnlist, irreplist=irreplist)
    boot.par[r, length(par)+1] <- opt.res$value
  }
  boot.par[r, c(1:length(par))] <- opt.res$par
  bootevalues[r,] <- getEvalues(boot.par[r,c(1:length(par))], Wlist=Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist)
  cat("r =", r, "fitted parameters:", opt.res$par, "\n")

  j <- 1
  k <- 1
  p <- 1
  for(L in Lvalues) {
    for(f in framelist) {
      for(i in irreplist) {

        if(nolist[[k]] > 0) {
          for(level in c(1:nolist[[k]])) {
            ##cat("level", level, "irrep", i, "frame", f, "L", L, "\n")
            ##cat(p, k, j, L, f, i, level, Epipi$irrep[iiEdata[p]], Epipi$frame[iiEdata[p]], "\n")
            if(Epipi$irrep[iiEdata[p]] == 1) {
              tmp <- try(solveforCotDelta(par=c(opt.res$par[1], 0, opt.res$par[3]), index=1, level=level,
                                          scanfn=scanfnlist[[i]][[f]], fn=fnlist[[i]][[f]], W=Wlist[[j]], E=boot.data[r,2*iiEdata[p]-1]))
              if(inherits(tmp, "try-error")) {
                bootqcotdelta[r,p] <- NA
                bootq[r,p] <- NA
              }
              else {
                bootqcotdelta[r,p] <- tmp
                bootq[r,p] <- getq(par=c(bootqcotdelta[r,p], 0, opt.res$par[3]), level=level, scanfn=scanfnlist[[i]][[f]], dfn=fnlist[[i]][[f]], W=Wlist[[j]])
              }
              ##cat("par1,2", opt.res$par[1] + opt.res$par[2]*Epipi$q[iiEdata[p]]^2/2.,
              ##    getEDiff(x=bootqcotdelta[r,p], index=1, par=c(opt.res$par[1], 0, opt.res$par[3]), level=level, scanfn=scanfnlist[[i]][[f]], dfn=fnlist[[i]][[f]], W=Wlist[[j]], E=0),
              ##    boot.data[r,2*iiEdata[p]-1], bootq[r,p]^2,"\n")

            }
            else if(Epipi$irrep[iiEdata[p]] == 2) {
              tmp <- try(solveforCotDelta(par=opt.res$par, index=3, level=level, interval=c(opt.res$par[3], -0.005),
                                          scanfn=scanfnlist[[i]][[f]], fn=fnlist[[i]][[f]], W=Wlist[[j]], E=boot.data[r,2*iiEdata[p]-1]))
              if(inherits(tmp, "try-error")) {
                bootqcotdelta[r,p] <- NA
                bootq[r,p] <- NA
              }
              else {
                bootqcotdelta[r,p] <- tmp
                bootq[r,p] <- getq(par=c(opt.res$par[c(1:2)], bootqcotdelta[r,p]), level=level, scanfn=scanfnlist[[i]][[f]], dfn=fnlist[[i]][[f]], W=Wlist[[j]])
              }
              
              ##cat("par3", opt.res$par[3],
              ##    getEDiff(x=bootqcotdelta[r,p], index=3, par=c(opt.res$par[1], opt.res$par[2], opt.res$par[3]), level=level, scanfn=scanfnlist[[i]][[f]], dfn=fnlist[[i]][[f]], W=Wlist[[j]], E=0),
              ##    boot.data[r,2*iiEdata[p]-1], bootq[r,p]^2,"\n")

            }
            p <- p+1
          }
        }
        k <- k+1
      }
      j <- j+1
    }
  }
  save(boot.par, bootevalues, iiEdata, bootqcotdelta, bootq, file=paste("fit-detEq-result-cut", cut, ".Rdata", sep=""))
}

#energydata <- boot.data[1,(2*iiEdata-1)]
#cat(evalues, "\n")
#cat(energydata, "\n")
#cat((evalues-energydata)/Epipi$dE[iiEdata], "\n")

cat("a0 = ", 1/boot.par[1,1], "+-", sd(1/boot.par[c(1:doSamples),1]), "\n")
cat("r0 = ", boot.par[1,2], "+-", sd(boot.par[c(1:doSamples),2]), "\n")
cat("a2 = ", 1/boot.par[1,3], "+-", sd(1/boot.par[c(1:doSamples),3]), "\n")
cat("chi^2 = ", boot.par[1,4], "\n")
cat("dof = ", dof, "\n")

save(boot.par, bootevalues, iiEdata, bootqcotdelta, bootq, file=paste("fit-detEq-result-cut", cut, ".Rdata", sep=""))

x <- seq(0.00001,0.08,0.0001)
y <- atan(sqrt(x)/(boot.par[1,1] + boot.par[1,2]*x/2))
y2 <- atan(x^2*sqrt(x)/(boot.par[1,3]))

y <- array(0, dim=c(doSamples, length(x)))
y2 <- array(0, dim=c(doSamples, length(x)))
for(i in c(1:doSamples)) {
  y[i,] <- atan(sqrt(x)/(boot.par[i,1] + boot.par[i,2]*x/2))
  y2[i,] <- atan(x^2*sqrt(x)/(boot.par[i,3]))
}
## error band
dy <- apply(y, 2, sd)
dy2 <- apply(y2, 2, sd)


tikzfiles <- tikz.init(paste("delta-vs-qsq.cut", cut, sep=""),width=6,height=5)

plot(NA, ylim=c(-0.8, 0.05), xlim=c(0, 0.08), xlab=c("$a^2q^2$"), ylab=c("$\\delta_\\ell$"))
polygon(x=c(x, rev(x)), y=c(y[1,]+dy, rev(y[1,]-dy)), col="gray", lty=0, lwd=0.001, border="gray")
lines(x,y[1,])
polygon(x=c(x, rev(x)), y=c(y2[1,]+dy2, rev(y2[1,]-dy2)), col="gray", lty=0, lwd=0.001, border="gray")
lines(x,y2[1,], lty=2)
legend("bottomleft", legend=c("$\\ell = 0$", "$\\ell = 2$"), bty="n", lty=c(1,2))

#legend(x=0.045, y=-0.45, legend=c("$\\ell = 0$"), bty="n")
#legend(x=0.05, y=0.05, legend=c("$\\ell = 2$"), bty="n")

ind1 <- which(Epipi$irrep[iiEdata] == 1)
indn1 <- which(Epipi$irrep[iiEdata] != 1)

plotwitherror(bootq[1,ind1]^2, atan(bootq[1,ind1]/bootqcotdelta[1,ind1]), dx=apply(bootq[c(1:doSamples),ind1]^2, 2, sd, na.rm=TRUE),
              dy=apply(atan(bootq[c(1:doSamples),ind1]/bootqcotdelta[c(1:doSamples),ind1]), 2, sd, na.rm=TRUE), rep=TRUE, col="red", pch=21)
if(any(is.na(bootq[1,indn1])) || any(is.na(bootqcotdelta[1,indn1]))) {
  plotwitherror(mean(bootq[c(1:doSamples),indn1], na.rm=TRUE)^2, atan(mean(bootq[c(1:doSamples),indn1], na.rm=TRUE)^5/mean(bootqcotdelta[c(1:doSamples),indn1], na.rm=TRUE)), dx=sd(bootq[c(1:doSamples),indn1]^2, na.rm=TRUE),
                dy=sd(atan(bootq[c(1:doSamples),indn1]^5/bootqcotdelta[c(1:doSamples),indn1]), na.rm=TRUE), rep=TRUE, col="blue", pch=22)
}
if(!any(is.na(bootq[1,indn1])) && !any(is.na(bootqcotdelta[1,indn1]))) {
  plotwitherror(bootq[1,indn1]^2, atan(bootq[1,indn1]^5/bootqcotdelta[1,indn1]), dx=sd(bootq[c(1:doSamples),indn1]^2, na.rm=TRUE),
                dy=sd(atan(bootq[c(1:doSamples),indn1]^5/bootqcotdelta[c(1:doSamples),indn1]), na.rm=TRUE), rep=TRUE, col="blue", pch=22)
}

legend("topright", legen=c("$\\delta_0$","$\\delta_2$"), bty="n", pch=c(21,22), col=c("red", "blue"))

tikz.finalize(tikzfiles=tikzfiles,clean=TRUE)

##
