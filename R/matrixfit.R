bootstrap.meanerror <- function(data, R=400, l=20) {
  bootit <- boot(block.ts(data, l=l), meanindexed, R=R)
  return(sd(bootit$t))
}

matrixChisqr <- function(par, t, y, M, T, parind) {
  z <- (y-0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[3]*(T-t)) + exp(-par[3]*t)))
  return( z %*% M %*% z )
}


matrixfit <- function(cf, t1, t2, symmetrise=TRUE, boot.R=400, boot.l=20, noObs=1,
                      matrix.size=2, useCov=FALSE) {
  ##if(symmetrise) {
  ##  
  ##}
  matrix.size <- 2
  Cor <- apply(cf$cf, 2, mean)
  Err <- apply(cf$cf, 2, bootstrap.meanerror, R=boot.R, l=boot.l)

  t1p1 <- t1+1
  t2p1 <- t2+1
  N <- length(cf$cf[,1])
  Thalfp1 <- cf$T/2+1
  t <- c(0:(cf$T/2))
  CF=data.frame(t=t, Cor=Cor, Err=Err)
  rm(Cor, Err)
  ii <- c((t1p1):(t2p1), (t1p1+Thalfp1):(t2p1+Thalfp1), (t1p1+2*Thalfp1):(t2p1+2*Thalfp1), (t1p1+3*Thalfp1):(t2p1+3*Thalfp1))
  parind <- array(1, dim=c(length(CF$Cor),2))
  parind[(Thalfp1+1):(2*Thalfp1),2] <- 2
  parind[(2*Thalfp1+1):(4*Thalfp1),1] <- 2
  parind[(3*Thalfp1+1):(4*Thalfp1),2] <- 2
  mSize <- length(CF$Cor)/length(t)
  if(mSize != matrix.size^2) {
    stop("matrix is not a square matrix, but we expect a square matrix!")
  }

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
    if(floor(sqrt(length(cov.svd$d))) < length(cf$cf[,1])) {
      cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)] <- mean(cov.svd$d[floor(sqrt(length(cov.svd$d))):length(cov.svd$d)])
    }
    D <- diag(1/cov.svd$d)
    M <- floor(N/boot.l)*cov.svd$v %*% D %*% t(cov.svd$u)
  }

  par <- numeric(3)
  ## we get initial guesses for fit parameters from effective masses
  par[3] <- -log(CF$Cor[t1p1+1]/CF$Cor[t1p1])
  par[1] <- sqrt(CF$Cor[t1p1]/0.5/exp(-par[3]*t1))
  par[2] <- sqrt(CF$Cor[(t1p1+3*Thalfp1)]/0.5/exp(-par[3]*t1))

  opt.res <- optim(par, matrixChisqr, method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                   t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$T, parind=parind[ii,])
  opt.res <- optim(opt.res$par, matrixChisqr, method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
                   t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$T, parind=parind[ii,])

  dof <- (length(CF$t[ii])-length(par))
  Qval <- 1-pchisq(opt.res$value, dof)

  opt.tsboot <- tsboot(tseries=cf$cf[,ii], statistic=fit.formatrixboot, R=boot.R, l=boot.l, sim="geom",
                       par=opt.res$par, t=CF$t[ii], M=M, T=cf$T, parind=parind[ii,])
  
  res <- list(CF=CF, M=M, parind=parind, ii=ii, opt.res=opt.res, opt.tsboot=opt.tsboot,
              boot.R=boot.R, boot.l=boot.l, useCov=useCov, CovMatrix=CovMatrix,
              Qval=Qval, chisqr=opt.res$value, dof=dof, mSize=mSize, cf=cf, t1=t1, t2=t2)
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
    y <- 0.5*par[parind[(i-1)*Thalfp1+1,1]]*par[parind[(i-1)*Thalfp1+1,2]]*(exp(-par[3]*(T-tx)) + exp(-par[3]*tx))
    col=c("red", "brown", "green", "blue")
    lines(tx, y, col=col[i], lwd=c(3))
  }
}

summary.matrixfit <- function(mfit) {
  cat("\n ** Result **\n\n")
  cat("based on", length(mfit$cf$cf[,1]), "measurements\n")
  cat("m \t=\t", mfit$opt.res$par[3], "\n")
  cat("dm\t=\t", sd(mfit$opt.tsboot$t[,3]), "\n")
  cat("P_1 \t=\t", mfit$opt.res$par[1], "\n")
  cat("dP_1\t=\t", sd(mfit$opt.tsboot$t[,1]), "\n\n")
  cat("P_2 \t=\t", mfit$opt.res$par[2], "\n")
  cat("dP_2\t=\t", sd(mfit$opt.tsboot$t[,2]), "\n\n")

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

fit.formatrixboot <- function(cf, par, t, M, T, parind) {
  opt.res <- optim(par, matrixChisqr, method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                   t=t, y=apply(cf,2,mean), M=M, T=T, parind=parind)
  opt.res <- optim(opt.res$par, matrixChisqr, method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
                   t=t, y=apply(cf,2,mean), M=M, T=T, parind=parind)
  return(c(opt.res$par, opt.res$value))
}
