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

## fit formel: R(t+1/2) = A*(cosh(dE*t') +sinh(dE*t')*coth(2*Mpi*t'))
## t' = t+1/2-T/2
## for zero momentum only
Rfn <- function(par, tp, m) {
  return(par[1]*(cosh(par[2]*tp) + sinh(par[2]*tp)/tanh(2*m*tp)))
}

getMatrix.pipi <- function(N=5, tp="TP0", irrep="A1", basename="pipi_pipi_",
                      T, path="./", ens, ind.vector=c(2,3)) {
  
  
  ## read data into Cmatrix
  Cmatrix <- cf()
  
  for(i in c(1:N)) {
    for(j in c(1:N)) {
      filename <- paste(basename, irrep, "_corr_", tp, "_", i-1, j-1, ".dat", sep="")
      tmp <- readtextcf(filename, T=T, check.t=1, path=path, ind.vector=ind.vector)
      Cmatrix <- c(Cmatrix, tmp)
    }
  }
  return(invisible(Cmatrix))
}


## solves the GEVP after temporal states have been removed by shiftig and weighting
## if neccessary.
solveGEVP.pipi <- function(Cmatrix, pion.cor, boot.R, boot.l, seed=123456, p1=c(0,0,0), p2=c(0,0,0),
                           element.order, t0=1, matrix.size=1, L, pionmethod="matrixfit", t1, t2, useCov=TRUE,
                           lat.disp=TRUE) {
  if(missing(element.order)) {
    element.order = c(1:(matrix.size^2))
  }
  ## first determine the pion mass
  if(!pion.cor$boot.samples) {
    pion.cor <- bootstrap.cf(pion.cor, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }
  if(pionmethod == "matrixfit") {
    if(!inherits(pion.cor, "matrixfit")) pion <- matrixfit(pion.cor, t1=t1, t2=t2, useCov=useCov)
    else pion <- pion.cor
  }
  else {
    if(!inherits(pion.cor, "effectivemassfit")) pion <- fit.effectivemass(bootstrap.effectivemass(pion.cor, boot.R=boot.R, boot.l=boot.l, type="solve"), t1=t1, t2=t2, useCov=useCov)
    else pion <- pion.cor
  }

  if(!Cmatrix$boot.samples) Cmatrix <- bootstrap.cf(Cmatrix, boot.R=boot.R, boot.l=boot.l, seed=seed)
  Cmatrix.removed <- removeTemporal.cf(cf=Cmatrix, single.cf1=pion, L=L, p1=p1, p2=p2)
  Cmatrix.bootstrap.gevp <- bootstrap.gevp(Cmatrix.removed, matrix.size=matrix.size, t0=t0, element.order=element.order)
  return(invisible(list(Cmatrix.gevp = Cmatrix.bootstrap.gevp, pion=pion)))
}

## determine deltaE from a fit to the ratio
## three ratios are available
## 1) for zero momentum, see Renner, Xu, Jansen
## 2) for the GEVP case, see compRpipi2 -> pipimetho="gevp", ratiomethod="zeromom"
## 3) for GEVP and non-zero momentum, see compRpipi3 -> pipimetho="gevp", ratiomethod="mom"
deltaEfromRatio <- function(pipi.cor, pion, pion1=NULL, pion2=NULL, boot.R, boot.l, Thalf, L,
                            tr1, tr2, seed=123456, useCov=TRUE, lat.disp=TRUE, p1=c(0,0,0), p2=c(0,0,0),
                            pc.id=1, pipimethod="single", mf=FALSE, ratiomethod="zeromom") {

  if(missing(L)) L <- Thalf
  rr <- c(2:(boot.R+1))
  ## first determine the pion mass
  Mps <- rep(0, times=boot.R+1)
  if(inherits(pion, "matrixfit")) {
    Mps[1] <- pion$opt.res$par[1]
    Mps[rr] <- pion$opt.tsboot[1,]
  }
  else if(inherits(pion, "effectivemassfit")) {
    Mps[1] <- pion$opt.res$par[1]
    Mps[rr] <- pion$massfit.tsboot[,1]
  }
  else stop("pion must be a matrixfit or effectivemassfit object\n")
  
  if(pipimethod != "single") {
    ## principal correlator
    pipi.cor <- gevp2cf(pipi.cor, id=pc.id)
  }
  else if(!pipi.cor$boot.samples) {
    pipi.cor <- bootstrap.cf(pipi.cor, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }

  ## the function to compute the ratio from 4pt and 2pt function
  Rfunction <- compRpipi
  if(pipimethod != "single") {
    if(ratiomethod == "zeromom") {
      Rfunction <- compRpipi2
    }
  }
  
  Rpipi.tsboot <- array(NA, dim=c(boot.R+1, Thalf))
  if(pipimethod == "single" || ratiomethod == "zeromom") {
    ## this is R(t+1/2)
    Rpipi.tsboot[1,] <- Rfunction(c4=pipi.cor$cf0, c2=pion$cf$cf0, Thalf=Thalf)
    
    for(i in c(1:boot.R)) {
      Rpipi.tsboot[i+1,] <- Rfunction(c4=pipi.cor$cf.tsboot$t[i,], c2=pion$cf$cf.tsboot$t[i,], Thalf=Thalf)
    }
  }
  else {
    dE <- rep(0, times=boot.R+1)
    if(lat.disp) {
      p1shift <- 2*sum(sin(pi*p1/L)^2)
      p2shift <- 2*sum(sin(pi*p2/L)^2)
      dE <- abs(acosh( cosh(Mps) + p1shift) - acosh( cosh(Mps) + p2shift))
    }
    else {
      p1shift <- sum((2*pi*p1/L)^2)
      p2shift <- sum((2*pi*p2/L)^2)
      dE <- abs(sqrt( Mps^2 + p1shift ) - sqrt( Mps^2 + p2shift ))
    }

    ## this is R(t)
    Rpipi.tsboot[1,] <- compRpipi3(c4=pipi.cor$cf0, c21=pion1$cf0, c22=pion2$cf0, Thalf=Thalf, dE=dE[1])
    
    for(i in c(1:boot.R)) {
      Rpipi.tsboot[i+1,] <- compRpipi3(c4=pipi.cor$cf.tsboot$t[i,], c21=pion1$cf.tsboot$t[i,], c22=pion2$cf.tsboot$t[i,], Thalf=Thalf, dE=dE[i+1])
    }
  }
  ## error of the ratio
  dRpipi <- apply(Rpipi.tsboot, 2, sd)

  ## chi^2 function includes covariance matrix M
  fitfn <- function(par, y, t, m, Thalf, M) {
    tp <- t - Thalf 
    z <- Rfn(par, tp, m)
    return((z-y) %*% M %*% (z-y))
  }
  ## the same, but in the format for nls.lm
  fitfn.lm <- function(par, y, t, m, Thalf, L) {
    z <- Rfn(par, t-Thalf, m)
    return(L %*% (z-y))
  }
  fitfnmf <- function(par, y, t, M) {
    z <- par[1]*exp(-par[2]*t)
    return((z-y) %*% M %*% (z-y))
  }
  fitfnmf.lm <- function(par, y, t, L) {
    z <- par[1]*exp(-par[2]*t)
    return(L %*% (z-y))
  }

  ## fit range
  tt <- c(tr1:tr2)
  ## ratio is displaced by 1/2
  tphys <- tt-0.5
  ## build the covariance matrix from bootstrap samples, if wanted
  M <- diag(1./dRpipi[tt]^2)
  if(useCov) M <- invertCovMatrix(Rpipi.tsboot[,tt], boot.samples=TRUE)
  
  par <- c(1.6,0.01)
  lm.avail <- require(minpack.lm)
  if(lm.avail) ML <- chol(M)
  opt.tsboot <- array(NA, dim=c(boot.R+1, 5))
  for(i in c(1:(boot.R+1))) {
    par <- c(1.6, 2*Mps[i])
    if(lm.avail) {
      if(!mf) opt <- nls.lm(par, fn=fitfn.lm,  L=ML, Thalf=Thalf, t=tphys, y=Rpipi.tsboot[i,tt], m=Mps[i])
      else opt <- nls.lm(par, fn=fitfnmf.lm,  L=ML, t=tt, y=Rpipi.tsboot[i,tt])
      ## chi^2 value
      opt.tsboot[i, 3] <- opt$rsstrace[length(opt$rsstrace)]
    }
    else {
      if(!mf) opt <- optim(par, fn = fitfn,
                           t=tphys, Thalf=Thalf, control=list(ndeps=rep(1.e-6, times=length(par)), parscale=1/par),
                           method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=Mps[i])
      else opt <- optim(par, fn = fitfnmf,
                        t=tt, control=list(ndeps=rep(1.e-6, times=length(par)), parscale=1/par),
                        method="BFGS", M=M, y = Rpipi.tsboot[i,tt])
      
      ## chi^2 value
      opt.tsboot[i, 3] <- opt$value
    }
    ## fit parameters amplitude and deltaE
    opt.tsboot[i, c(1,2)] <- opt$par
    ## pion mass

  }
  opt.tsboot[, 4] <- Mps
  opt.tsboot[, 5] <- 2*Mps + opt.tsboot[, 2]
  if(!mf) {
    x <- seq(tphys[1], tphys[length(tphys)]+1, 0.01)
    rat <- Rfn(opt.tsboot[1,c(1,2)], tp=x-Thalf, m=opt.tsboot[1,4])
  }
  else {
    x <- seq(tt[1], tt[length(tt)]+1, 0.01)
    rat <- opt.tsboot[1,1]*exp(-opt.tsboot[1,2]*x)
  }
  dof <- length(tt)-2
  Qval <- 1-pchisq(opt.tsboot[1,3], dof)

  return(invisible(list(opt.tsboot=opt.tsboot, Rpipi.tsboot=Rpipi.tsboot, dRpipi = dRpipi, x=x, rat=rat, dof=dof, Qval=Qval, pionQval=pion$Qval,
                        mf=mf, pipimethod=pipimethod, ratiomethod=ratiomethod, pion=pion)))
}

mergeEnergies <- function(irrep, ens, pc.id, tp, path="./", R2=FALSE) {

  filelist <- Sys.glob(paste(path, "dEres.", tp, ".", irrep, ".", ens, ".", pc.id, ".*.Rdata", sep=""))
  cat(ens, irrep, tp, pc.id, "averaging over", length(filelist), "files\n")
  load(filelist[1])
  boot.R <- length(dEres$opt.tsboot[,1])
  res <- array(NA, dim=c(boot.R, length(filelist), 3))
  Qval <- array(NA, dim=c(length(filelist), 2))
  for(i in c(1:length(filelist))) {
    load(filelist[i])
    if(tp != "TP0" && R2 && !is.null(dEresR2) ) dEres <- dEresR2
    if(boot.R != length(dEres$opt.tsboot[,1])) stop(paste("boot.R not consistent for file", filelist[i]))
    res[,i,] <- dEres$opt.tsboot[,c(2,4,5)]
    Qval[i,1] <- dEres$Qval
    Qval[i,2] <- dEres$pionQval
  }
  if(R2) save(Qval, res, file=paste(path, "resR2.", tp, ".", irrep, ".", ens, ".", pc.id, ".Rdata", sep=""))
  else save(Qval, res, file=paste(path, "res.", tp, ".", irrep, ".", ens, ".", pc.id, ".Rdata", sep=""))
  return(invisible(list(res=res, Qval=Qval)))
}

summariseEnergies <- function(irrep, ens, pc.id, tp, prob=c(0.1573, 0.8427), path="./", R2=FALSE) {
  if(R2) load(file=paste(path, "resR2.", tp, ".", irrep, ".", ens, ".", pc.id, ".Rdata", sep=""))
  else load(file=paste(path, "res.", tp, ".", irrep, ".", ens, ".", pc.id, ".Rdata", sep=""))
  
  deltaE <- estimate.error(res, index=1, prob=prob, main="deltaE", Qval=Qval)
  deltaEboots <- compute.boots(res, index=1, Qval=Qval)
  Mpi <- estimate.error(res, index=2, prob=prob, main="Mpi", Qval=Qval)
  Mpiboots <- compute.boots(res, index=2, piononly=TRUE, Qval=Qval)
  E <- estimate.error(res, index=3, prob=prob, main="E", Qval=Qval)
  Eboots <- compute.boots(res, index=3, Qval=Qval)
  reslist <- list(deltaE=deltaE, Mpi=Mpi, E=E, irrep=irrep, ens=ens, pc.id=pc.id, tp=tp, deltaEboots=deltaEboots, Mpiboots=Mpiboots, Eboots=Eboots)
  save(reslist, file=paste(path, "reslist.", tp, ".", irrep, ".", ens, ".", pc.id, ".Rdata", sep=""))
  return(invisible(reslist))
}

