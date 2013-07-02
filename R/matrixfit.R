bootstrap.meanerror <- function(data, R=400, l=20) {
  bootit <- boot(block.ts(data, l=l), meanindexed, R=R)
  return(sd(bootit$t))
}

matrixChisqr <- function(par, t, y, M, T, parind, sign.vec) {
  z <- (y-0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(T-t)) + sign.vec*exp(-par[1]*t)))
  return( z %*% M %*% z )
}

dmatrixChisqr <- function(par, t, y, M, T, parind, sign.vec) {
  res <- rep(0., times=length(par))
  z <- (y-0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(T-t)) + sign.vec*exp(-par[1]*t)))
  zp <- -0.5*par[parind[,1]]*par[parind[,2]]*(-(T-t)*exp(-par[1]*(T-t)) + (-t)*sign.vec*exp(-par[1]*t))
  res[1] <- zp %*% M %*% z + z %*% M %*% zp
  for(i in 2:length(par)) {
    zp <- rep(0, length(z))
    j <- which(parind[,1]==i)
    zp[j] <- -0.5*par[parind[j,2]]*(exp(-par[1]*(T-t[j])) + sign.vec[j]*exp(-par[1]*t[j]))
    res[i] <- zp %*% M %*% z + z %*% M %*% zp
    zp <- rep(0, length(z))
    j <- which(parind[,2]==i)
    zp[j] <- -0.5*par[parind[j,1]]*(exp(-par[1]*(T-t[j])) + sign.vec[j]*exp(-par[1]*t[j]))
    res[i] <- zp %*% M %*% z + z %*% M %*% zp
  }
  return(res)
}

matrixfit <- function(cf, t1, t2, symmetrise=TRUE, boot.R=400, boot.l=20,
                      parlist = array(c(1,1,1,2,2,1,2,2), dim=c(2,4)),
                      sym.vec=c("cosh","cosh","cosh","cosh"),
                      matrix.size=2, useCov=FALSE, seed=12345) {
  ## some sanity checks
  if(min(parlist) <= 0 || max(parlist) > matrix.size) {
    stop("Elements of parlist must be all > 0 and <= matrix.size\n")
  }
  for(i in 1:max(parlist)) {
    if(!any(parlist==i)) {
      stop("not all parameters are used in the fit\n")
    }
  }

  Cor <- apply(cf$cf, 2, mean)
  Err <- apply(cf$cf, 2, bootstrap.meanerror, R=boot.R, l=boot.l)

  t1p1 <- t1+1
  t2p1 <- t2+1
  N <- length(cf$cf[,1])
  Thalfp1 <- cf$T/2+1
  t <- c(0:(cf$T/2))
  CF=data.frame(t=t, Cor=Cor, Err=Err)
  rm(Cor, Err)

  ## This is the number of correlators in cf
  mSize <- length(CF$Cor)/Thalfp1
  if(length(parlist[1,]) < mSize) {
    stop("parlist has not the correct length\n")
  }
  if(length(sym.vec) < mSize) {
    stop("sym.vec has not the correct length\n")
  }

  ## index vector for timeslices to be fitted
  ii <- c((t1p1):(t2p1))
  for(j in 2:mSize) {
    ii <- c(ii, (t1p1+(j-1)*Thalfp1):(t2p1+(j-1)*Thalfp1))
  }
  
  ## parind is the index vector for the matrix elements
  ## signvec decides on cosh or sinh
  parind <- array(1, dim=c(length(CF$Cor),2))
  sign.vec <- rep(+1, times=length(CF$Cor))
  for(i in 1:mSize) {
    parind[((i-1)*Thalfp1+1):(i*Thalfp1),] <- t(array(parlist[,i]+1, dim=c(2,Thalfp1)))
    if(sym.vec[i] != "cosh") sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
  }
  
  ## for uncorrelated chi^2 use diagonal matrix with inverse sd^2
  M <- diag(1/CF$Err[ii]^2)
  CovMatrix <- numeric()
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ## see C. Michael hep-lat/9412087
    ## block data first
    ncf <- block.ts(cf$cf, l=boot.l)
    ## compute covariance matrix and invert
    CovMatrix <- cov(ncf[,ii])
    cov.svd <- svd(CovMatrix)
    ## replace smallest eigenvalues by their mean, if needed
    if(floor(sqrt(length(cov.svd$d))) < length(cf$cf[,1])) {
      cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)] <-
        mean(cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)])
    }
    D <- diag(1/cov.svd$d)
    M <- floor(N/boot.l)*cov.svd$v %*% D %*% t(cov.svd$u)
  }

  par <- numeric(max(parind))
  ## we get initial guesses for fit parameters from effective masses
  ## first is the mass
  j <- which(parlist[1,]==1 & parlist[2,]==1)
  par[1] <- -log(CF$Cor[t1p1+(j-1)*Thalfp1+1]/CF$Cor[t1p1+(j-1)*Thalfp1])
  ## the amplitudes we estimate from diagonal elements
  for(i in 2:length(par)) {
    j <- which(parlist[1,]==(i-1) & parlist[2,]==(i-1))
    if(length(j) == 0) {
      warning("one diagonal element does not appear in parlist\n")
      j <- 1
    }
    par[i] <- sqrt(CF$Cor[t1p1+(j-1)*Thalfp1]/0.5/exp(-par[1]*t1))
  }

  ## now perform minimisation
  opt.res <- optim(par, fn = matrixChisqr, gr = dmatrixChisqr,
                   method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                   t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$T, parind=parind[ii,], sign.vec=sign.vec[ii])
  ##opt.res <- optim(opt.res$par, matrixChisqr, gr = dmatrixChisqr,
  ##                 method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
  ##                 t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$T, parind=parind[ii,], sign.vec=sign.vec[ii])

  dof <- (length(CF$t[ii])-length(par))
  Qval <- 1-pchisq(opt.res$value, dof)

  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## and bootstrap the fit
  opt.tsboot <- tsboot(tseries=cf$cf[,ii], statistic=fit.formatrixboot, R=boot.R, l=boot.l, sim="geom",
                       par=opt.res$par, t=CF$t[ii], M=M, T=cf$T, parind=parind[ii,], sign.vec=sign.vec[ii])
  
  res <- list(CF=CF, M=M, parind=parind, ii=ii, opt.res=opt.res, opt.tsboot=opt.tsboot,
              boot.R=boot.R, boot.l=boot.l, useCov=useCov, CovMatrix=CovMatrix,
              Qval=Qval, chisqr=opt.res$value, dof=dof, mSize=mSize, cf=cf, t1=t1, t2=t2,
              parlist=parlist, sym.vec=sym.vec, seed=seed)
  attr(res, "class") <- c("matrixfit", "list")
  return(invisible(res))
}

plot.matrixfit <- function(mfit, ...) {
  par <- mfit$opt.res$par
  parind <-  mfit$parind
  T <- mfit$cf$T
  Thalfp1 <- T/2+1
  plotwitherror(mfit$CF$t, mfit$CF$Cor, mfit$CF$Err, log="y", ...)
  tx <- seq(mfit$t1, mfit$t2, 0.005)
  for(i in c(1:mfit$mSize)) {
    y <- 0.5*par[parind[(i-1)*Thalfp1+1,1]]*par[parind[(i-1)*Thalfp1+1,2]]*(exp(-par[1]*(T-tx)) + exp(-par[1]*tx))
    col=c("red", "brown", "green", "blue")
    lines(tx, y, col=col[i], lwd=c(3))
  }
}

summary.matrixfit <- function(mfit) {
  cat("\n ** Result **\n\n")
  cat("based on", length(mfit$cf$cf[,1]), "measurements\n")
  cat("mass:\n")
  cat("m \t=\t", mfit$opt.res$par[1], "\n")
  cat("dm\t=\t", sd(mfit$opt.tsboot$t[,1]), "\n")
  cat("\nAmplitudes:\n")
  for(i in 2:length(mfit$opt.res$par)) {
    cat("P",i-1,"\t=\t", mfit$opt.res$par[i], "\n")
    cat("dP",i-1,"\t=\t", sd(mfit$opt.tsboot$t[,i]), "\n")
  }
  cat("\n")
  cat("boot.R\t=\t", mfit$boot.R, "\n")
  cat("boot.l\t=\t", mfit$boot.l, "\n")
  cat("useCov\t=\t", mfit$useCov, "\n")
  cat("chisqr\t=\t", mfit$opt.res$value, "\ndof\t=\t", mfit$dof, "\nchisqr/dof=\t",
      mfit$opt.res$value/mfit$dof, "\n")
  ## probability to find a larger chi^2 value
  ## if the data is generated again with the same statistics
  ## given the model is correct
  cat("Quality of the fit (p-value):", mfit$Qval, "\n")

}

fit.formatrixboot <- function(cf, par, t, M, T, parind, sign.vec) {
  opt.res <- optim(par, fn = matrixChisqr, gr = dmatrixChisqr,
                   method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                   t=t, y=apply(cf,2,mean), M=M, T=T, parind=parind, sign.vec=sign.vec)
  ##opt.res <- optim(opt.res$par, fn = matrixChisqr, gr = dmatrixChisqr,
  ##                 method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
  ##                 t=t, y=apply(cf,2,mean), M=M, T=T, parind=parind, sign.vec=sign.vec)
  return(c(opt.res$par, opt.res$value))
}
