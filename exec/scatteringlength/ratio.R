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

## fit formel: R(t+1/2) = A*(cosh(dE*t') +sinh(dE*t')*coth(2*Mpi*t'))
## t' = t+1/2-T/2
Rfn <- function(par, tp, m) {
  return(par[1]*(cosh(par[2]*tp) + sinh(par[2]*tp)/tanh(2*m*tp)))
}

## 
solveGEVP.pipi <- function(Cmatrix, pion.cor, boot.R, boot.l, seed=123456, p1=c(0,0,0), p2=c(0,0,0),
                           element.order, t0=1, matrix.size=1, L, pionmethod="matrixfit", t1, t2, useCov=TRUE) {
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
deltaEfromRatio <- function(pipi.cor, pion, boot.R, boot.l, Thalf,
                            tr1, tr2, seed=123456, useCov=TRUE,
                            pc.id=1, pipimethod="single") {

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
    Rfunction <- compRpipi2
  }
  
  Rpipi.tsboot <- array(NA, dim=c(boot.R+1, Thalf))
  ## this is R(t+1/2)
  Rpipi.tsboot[1,] <- Rfunction(c4=pipi.cor$cf0, c2=pion$cf$cf0, Thalf=Thalf)

  for(i in c(1:boot.R)) {
    Rpipi.tsboot[i+1,] <- Rfunction(pipi.cor$cf.tsboot$t[i,], pion$cf$cf.tsboot$t[i,], Thalf=Thalf)
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
      opt <- nls.lm(par, fn=fitfn.lm,  L=ML, Thalf=Thalf, t=tphys, y=Rpipi.tsboot[i,tt], m=Mps[i])
      ## chi^2 value
      opt.tsboot[i, 3] <- opt$rsstrace[length(opt$rsstrace)]
    }
    else {
      opt <- optim(par, fn = fitfn,
                   t=tphys, Thalf=Thalf, control=list(ndeps=rep(1.e-6, times=length(par)), parscale=1/par),
                   method="BFGS", M=M, y = Rpipi.tsboot[i,tt], m=Mps[i])
      ## chi^2 value
      opt.tsboot[i, 3] <- opt$value
    }
    ## fit parameters amplitude and deltaE
    opt.tsboot[i, c(1,2)] <- opt$par
    ## pion mass

  }
  opt.tsboot[, 4] <- Mps
  opt.tsboot[, 5] <- 2*Mps + opt.tsboot[, 2]
  x <- seq(tphys[1], tphys[length(tphys)]+1, 0.01)
  rat <- Rfn(opt.tsboot[1,c(1,2)], tp=x-Thalf, m=opt.tsboot[1,4])
  dof <- length(tt)-2
  Qval <- 1-pchisq(opt.tsboot[1,3], dof)

  return(invisible(list(opt.tsboot=opt.tsboot, Rpipi.tsboot=Rpipi.tsboot, dRpipi = dRpipi, x=x, rat=rat, dof=dof, Qval=Qval, pionQval=pion$Qval)))
}
