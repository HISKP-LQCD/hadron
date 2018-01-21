bootstrap.meanerror <- function(data, R=400, l=20) {
  bootit <- boot(block.ts(data, l=l), meanindexed, R=R)
  return(apply(bootit$t, 2, sd))
}

matrixModel <- function(par, t, T, parind, sign.vec, ov.sign.vec) {
  return(ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(T-t))))
}

matrixChisqr <- function(par, t, y, M, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(T-t))))
  return( sum(z %*% M %*% z) )
}

matrixChi <- function(par, t, y, L, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(T-t))))
  return( L %*% z )
}

## deltat is dummy variable here
dmatrixChi <- function(par, t, y, L, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  zp <- -ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(-t*exp(-par[1]*t) -(T-t)*sign.vec*exp(-par[1]*(T-t)))
  res <- L %*% zp
  for(i in 2:length(par)) {
    zp1 <- rep(0, length(zp))
    j <- which(parind[,1]==i)
    zp1[j] <- -ov.sign.vec*0.5*par[parind[j,2]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(T-t[j])))
    zp2 <- rep(0, length(zp))
    j <- which(parind[,2]==i)
    zp2[j] <- -ov.sign.vec*0.5*par[parind[j,1]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(T-t[j])))
    res <- c(res, L %*% zp1 + L %*% zp2)
  }
  return(res)
}

## deltat is dummy variable here
dmatrixChisqr <- function(par, t, y, M, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  res <- rep(0., times=length(par))
  z <- (y-ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*t) + sign.vec*exp(-par[1]*(T-t))))
  zp <- -ov.sign.vec*0.5*par[parind[,1]]*par[parind[,2]]*(-t*exp(-par[1]*t) -(T-t)*sign.vec*exp(-par[1]*(T-t)))
  res[1] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  for(i in 2:length(par)) {
    zp <- rep(0, length(z))
    j <- which(parind[,1]==i)
    zp[j] <- -ov.sign.vec*0.5*par[parind[j,2]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(T-t[j])))
    res[i] <- sum(zp %*% M %*% z + z %*% M %*% zp)
    zp <- rep(0, length(z))
    j <- which(parind[,2]==i)
    zp[j] <- -ov.sign.vec*0.5*par[parind[j,1]]*(exp(-par[1]*t[j]) + sign.vec[j]*exp(-par[1]*(T-t[j])))
    res[i] <- res[i] + sum(zp %*% M %*% z + z %*% M %*% zp)
  }
  return(res)
}

matrixChisqr.shifted <- function(par, t, y, M, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(T-(t-deltat/2)))))
  return( sum(z %*% M %*% z ) )
}

dmatrixChisqr.shifted <- function(par, t, y, M, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  res <- rep(0., times=length(par))
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(T-(t-deltat/2)))))
  zp <- -ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) + (T-t+deltat/2)*sign.vec*exp(-par[1]*(T-(t-deltat/2))))
  res[1] <- sum(zp %*% M %*% z + z %*% M %*% zp)
  for(i in 2:length(par)) {
    zp <- rep(0, length(z))
    j <- which(parind[,1]==i)
    zp[j] <- -ov.sign.vec[j]*par[parind[j,2]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(T-(t[j]-deltat/2))))
    res[i] <- sum(zp %*% M %*% z + z %*% M %*% zp)
    zp <- rep(0, length(z))
    j <- which(parind[,2]==i)
    zp[j] <- -ov.sign.vec[j]*par[parind[j,1]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(T-(t[j]-deltat/2))))
    res[i] <- res[i] + sum(zp %*% M %*% z + z %*% M %*% zp)
  }
  return(res)
}

matrixChi.shifted <- function(par, t, y, L, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  z <- (y-ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(exp(-par[1]*(t-deltat/2)) - sign.vec*exp(-par[1]*(T-(t-deltat/2)))))
  return( L %*% z )
}

dmatrixChi.shifted <- function(par, t, y, L, T, parind, sign.vec, ov.sign.vec, deltat=1) {
  zp <- -ov.sign.vec*par[parind[,1]]*par[parind[,2]]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) +(T-t+deltat/2)*sign.vec*exp(-par[1]*(T-(t-deltat/2))))
  res <- L %*% zp
  for(i in 2:length(par)) {
    zp1 <- c(0)
    j <- which(parind[,1]==i)
    zp1[j] <- -ov.sign.vec[j]*par[parind[j,2]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(T-(t[j]-deltat/2))))
    zp2 <- c(0)
    j <- which(parind[,2]==i)
    zp2[j] <- -ov.sign.vec[j]*par[parind[j,1]]*(exp(-par[1]*(t[j]-deltat/2)) - sign.vec[j]*exp(-par[1]*(T-(t[j]-deltat/2))))
    res <- c(res, L %*% zp1 + L %*% zp2)
  }
  return(res)
}


# the calling code must supply the correct three parameters 
deriv.CExp <- function(par, t, T, sign) {
  res <- array(0.,dim=c(length(par),length(t)))

  res[1,] <- 0.5*par[2]*par[3]*(-t*exp(-par[1]*t) -(T-t)*sign*exp(-par[1]*(T-t)))
  res[2,] <- 0.5*par[3]*(exp(-par[1]*t) + sign*exp(-par[1]*(T-t)))
  res[3,] <- 0.5*par[2]*(exp(-par[1]*t) + sign*exp(-par[1]*(T-t)))

  return(res)
}

deriv.CExp.shifted <- function(par, t, T, sign, deltat=1) {
  res <- array(0.,dim=c(length(par),length(t)))
  
  res[1,] <- par[2]*par[3]*(-(t-deltat/2)*exp(-par[1]*(t-deltat/2)) + (T-(t-deltat/2))*sign*exp(-par[1]*(T-(t-deltat/2))))
  res[2,] <- par[3]*(exp(-par[1]*(t-deltat/2)) - sign*exp(-par[1]*(T-(t-deltat/2))))
  res[3,] <- par[2]*(exp(-par[1]*(t-deltat/2)) - sign*exp(-par[1]*(T-(t-deltat/2))))

  return(res)
}

matrixfit <- function(cf, t1, t2, symmetrise=TRUE, boot.R=400, boot.l=20,
                      parlist, sym.vec, neg.vec,
                      useCov=FALSE, seed=12345, model="single",
                      boot.fit=TRUE, fit.method="optim") {
  if(!any(class(cf) == "cf")) {
    stop("matrixfit requires the object to be of class cf! Aborting...!\n")
  }
  t1p1 <- t1+1
  t2p1 <- t2+1
  N <- dim(cf$cf)[1]
  Thalfp1 <- cf$T/2+1
  t <- c(0:(cf$T/2))
  deltat <- 1
  if(model == "shifted" && any(names(cf) == "deltat")) {
    deltat <- cf$deltat
  }
  
  ## This is the number of correlators in cf
  if(!is.null(dim(cf$cf)))
    mSize <- dim(cf$cf)[2]/Thalfp1
  else
    mSize <- dim(cf$cf.tsboot$t)[2]/Thalfp1


  if(missing(parlist)) {
    if(mSize == 1) {
      parlist <- array(c(1,1), dim=c(2,1))
      warning("missing parlist, using default for single correlator!\n")
    }
    else if(mSize == 4) {
      parlist <- array(c(1,1,1,2,2,1,2,2), dim=c(2,4))
      warning("missing parlist, using default for four correlators!\n")
    }
    else {
      stop("parlist is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  if(missing(sym.vec)) {
    if(mSize == 1) {
      sym.vec <- c("cosh")
      warning("missing sym.vec, using default for single correlator!\n")
    }
    else if(mSize == 4) {
      sym.vec <- c("cosh","cosh","cosh","cosh")
      warning("missing sym.vec, using default for four correlators!\n")
    }
    else {
      stop("sym.vec is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  if(missing(neg.vec)){
    if(mSize == 1) {
      neg.vec <- c(1)
      warning("missing neg.vec, using default (correlator positive)!\n")
    } 
    else if( mSize == 4 ){
      neg.vec <- c(1,1,1,1)
      warning("missing neg.vec, using default (all correlators positive)!\n")
    }
    else {
      stop("neg.vec is missing and no default is available for this cf size! Aborting...\n")
    }
  }

  
  ## some sanity checks
  if(min(parlist) <= 0) {
    stop("Elements of parlist must be all > 0! Aborting\n")
  }
  for(i in 1:max(parlist)) {
    if(!any(parlist==i)) {
      stop("not all parameters are used in the fit! Aborting\n")
    }
  }
  
  if(dim(parlist)[2] != mSize) {
    cat(mSize, dim(parlist)[2], "\n")
    stop("parlist has not the correct length! Aborting! Use e.g. extractSingleCor.cf or c to bring cf to correct number of observables\n")
  }
  if(length(sym.vec) != mSize) {
    stop("sym.vec does not have the correct length! Aborting\n")
  }
  if(length(neg.vec) != mSize){
    stop("neg.vec does not have the correct length! Aborting\n")
  }

  ## now we start the real computation
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R, boot.l, seed)
  }
  else {
    boot.R <- cf$boot.R
    boot.l <- cf$boot.l
    seed <- cf$seed
  }

  CF <- data.frame(t=t, Cor=cf$cf0, Err=apply(cf$cf.tsboot$t, 2, sd))

  ## index vector for timeslices to be fitted
  ii <- c((t1p1):(t2p1))
  if(mSize > 1) {
    for(j in 2:mSize) {
      ii <- c(ii, (t1p1+(j-1)*Thalfp1):(t2p1+(j-1)*Thalfp1))
    }
  }

  ## parind is the index vector for the matrix elements
  ## signvec decides on cosh or sinh
  ## ov.sign.vec indicates the overall sign 
  parind <- array(1, dim=c(length(CF$Cor),2))
  sign.vec <- rep(+1, times=length(CF$Cor))
  ov.sign.vec <- rep(+1, times=length(CF$Cor))
  for(i in 1:mSize) {
    parind[((i-1)*Thalfp1+1):(i*Thalfp1),] <- t(array(parlist[,i]+1, dim=c(2,Thalfp1)))
    if(sym.vec[i] == "sinh") sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
    if(sym.vec[i] == "exp") sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- 0
    if(neg.vec[i] == -1) ov.sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
  }
  
  CovMatrix <- NULL
  if(is.null(cf$cf)) {
    CovMatrix <- cov(cf$cf.tsboot$t[,ii])
  }
  else {
    CovMatrix <- cov(cf$cf[,ii])/N
  }
  
  ## for uncorrelated chi^2 use diagonal matrix with inverse sd^2
  M <- diag(1/CF$Err[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ## see C. Michael hep-lat/9412087
    if(is.null(cf$cf)) {
      M <- try(invertCovMatrix(cf$cf.tsboot$t[,ii], boot.l=boot.l, boot.samples=TRUE), silent=TRUE)
    }
    else {
      M <- try(invertCovMatrix(cf$cf[,ii], boot.l=boot.l), silent=TRUE)
    }
    if(inherits(M, "try-error")) {
      M <- diag(1/CF$Err[ii]^2)
      warning("inversion of variance covariance matrix failed in matrixfit, continuing with uncorrelated chi^2\n")
      useCov <- FALSE
    }
  }
  lm.avail=FALSE
  if(fit.method == "lm") 
    lm.avail <- require(minpack.lm)
  LM <- chol(M)
  
  par <- numeric(max(parind))
  ## we get initial guesses for fit parameters from effective masses
  ## first is the mass
  ## (we currently allow for only one)
  j <- which(parlist[1,]==1 & parlist[2,]==1)
  par[1] <- invcosh(CF$Cor[t1p1+(j-1)*Thalfp1]/CF$Cor[t1p1+(j-1)*Thalfp1+1], t=t1p1, cf$T)
  ## catch failure of invcosh
  if(is.na(par[1]) || is.nan(par[1])) par[1] <- 0.2
  ## the amplitudes we estimate from diagonal elements
  for(i in 2:length(par)) {
    j <- which(parlist[1,]==(i-1) & parlist[2,]==(i-1))
    if(length(j) == 0) {
      ##if(full.matrix) warning("one diagonal element does not appear in parlist\n")
      j <- i-1
    }
    par[i] <- sqrt(abs(CF$Cor[t1p1+(j-1)*Thalfp1])/0.5/exp(-par[1]*t1))
  }

  fitfn <- matrixChisqr
  dfitfn <- dmatrixChisqr
  if(model == "shifted") {
    fitfn <- matrixChisqr.shifted
    dfitfn <- dmatrixChisqr.shifted
  }
  if(model == "weighted") {
    if(any(names(cf) == "weighted")) {
      if(cf$weighted) {
        fitfn <- matrixChisqr.weighted
        dfitfn <- dmatrixChisqr.weighted
      }
      else {
        model <- "single"
        warning("model chosen to be weighted, but weighting information not available in cf. Falling back to model=single\n")
      }
    }
  }

  if(lm.avail) {
    fitfn <- matrixChi
    dfitfn <- dmatrixChi
    if(model == "shifted") {
      fitfn <- matrixChi.shifted
      dfitfn <- dmatrixChi.shifted
    }
    if(model == "weighted") {
      fitfn <- matrixChi.weighted
      dfitfn <- dmatrixChi.weighted
    }
  }
  
  ## check out constrOptim
  ## now perform minimisation
  dof <- (length(CF$t[ii])-length(par))
  opt.res <- NA
  rchisqr <- 0.
  if(lm.avail) {
    opt.res <- nls.lm(par = par, fn = fitfn, jac=dfitfn, t=CF$t[ii], y=CF$Cor[ii], L=LM, T=cf$Time, deltat=deltat,
                      parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii],
                      control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500))
    rchisqr <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fitfn, gr = dfitfn,
                     method="BFGS", control=list(maxit=500, parscale=par, ndeps=rep(1.e-8, times=length(par)), REPORT=50),
                     t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], 
                     ov.sign.vec=ov.sign.vec[ii], deltat=deltat)
    rchisqr <- opt.res$value
  }
  Qval <- 1-pchisq(rchisqr, dof)
  
  opt.tsboot <- NA
  if(boot.fit) {
    opt.tsboot <- apply(X=cf$cf.tsboot$t[,ii], MARGIN=1, FUN=fit.formatrixboot, par=opt.res$par, t=CF$t[ii], deltat=deltat,
                        M=M, T=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii],
                        L=LM, lm.avail=lm.avail, fitfn=fitfn, dfitfn=dfitfn)
  }
  N <- length(cf$cf[,1])
  if(is.null(cf$cf)) {
    N <- cf$N
  }

  res <- list(CF=CF, M=M, L=LM, parind=parind, sign.vec=sign.vec, ov.sign.vec=ov.sign.vec, ii=ii, opt.res=opt.res, opt.tsboot=opt.tsboot,
              boot.R=boot.R, boot.l=boot.l, useCov=useCov, CovMatrix=CovMatrix, invCovMatrix=M, seed=seed,
              Qval=Qval, chisqr=rchisqr, dof=dof, mSize=mSize, cf=cf, t1=t1, t2=t2,
              parlist=parlist, sym.vec=sym.vec, seed=seed, N=N, model=model, fit.method=fit.method)
  res$t <- t(opt.tsboot)
  res$t0 <- c(opt.res$par, opt.res$value)
  res$se <- apply(opt.tsboot[c(1:(dim(opt.tsboot)[1]-1)),], MARGIN=1, FUN=sd)
  attr(res, "class") <- c("matrixfit", "list")
  return(invisible(res))
}

plot.matrixfit <- function(mfit, plot.errorband=FALSE, ylim, do.qqplot=TRUE, rep=FALSE, col,...) {
  par <- mfit$opt.res$par
  parind <-  mfit$parind
  sign.vec <- mfit$sign.vec
  ov.sign.vec <- mfit$ov.sign.vec
  T <- mfit$cf$T
  Thalfp1 <- T/2+1
  deltat <- 1
  if(mfit$model == "shifted" && any(names(mfit$cf) == "deltat")) {
    deltat <- mfit$cf$deltat
  }
  # prevent stray negative values from ruining the plot
  lbound <- ov.sign.vec*mfit$CF$Cor - 2*mfit$CF$Err
  lbound <- lbound[ lbound > 0 ]
  if(missing(ylim)) ylims <- c( min( lbound, na.rm=TRUE ) , max( ov.sign.vec*mfit$CF$Cor + 2*mfit$CF$Err, na.rm=TRUE ) )
  else ylims <- ylim

  if(missing(col)){
    col <- c("black",rainbow(n=(mfit$mSize-1)))
  }
  plotwitherror(x=mfit$CF$t, y=ov.sign.vec*mfit$CF$Cor, 
                dy=mfit$CF$Err, log="y", ylim=ylims, rep=rep, col=col,...)
  tx <- seq(mfit$t1, mfit$t2, 0.05)
  for(i in 1:mfit$mSize ) {
    par.ind <- c(1,parind[(i-1)*Thalfp1+1,1],parind[(i-1)*Thalfp1+1,2])
    pars <- c(par[1],par[par.ind[2]],par[par.ind[3]])
    sgn <- sign.vec[(i-1)*Thalfp1+1]

    if(mfit$model == "shifted") y <- pars[2]*pars[3]*( exp(-pars[1]*(tx-deltat/2)) - sgn*exp(-pars[1]*(T-(tx-deltat/2))))
    else y <- 0.5*pars[2]*pars[3]*( exp(-pars[1]*tx) + sgn*exp(-pars[1]*(T-tx)))

    if(plot.errorband) {
      ## in the following the covariance matrix of the paramters and vectors of derivatives
      ## of the correlation function with respect to these parameters will be computed
      ## there is some level of waste because certain combinations of parameters might
      ## occur multiple times, but the overhead is small and this way is the most convenient

      parCov <- cov(mfit$t[,par.ind])

      ## 3 by length(tx) array of derivatives of the model with respect to the three parameters
      ## if any parameters are the same, multiplication with the covariance matrix will give
      ## the correct contributions to the derivative
      if(mfit$model == "shifted") div <- deriv.CExp.shifted(par=pars, t=tx, T=T, sign=sgn, deltat)
      else div <- deriv.CExp(par=pars, t=tx, T=T, sign=sgn)
      yvar <- t(div) %*% parCov %*% div

      ## yvar is a length(tx) by length(tx) matrix of which only the diagonal elements are of
      ## interest here
      ysd <- sqrt(diag(yvar))

      polyval <- c( (y + ysd), rev(y - ysd) )
      polyx <- c(tx,rev(tx))
      polycol <- col2rgb(col[i],alpha=TRUE)/255
      polycol[4] <- 0.65

      polygon(x=polyx,y=polyval,col=rgb(red=polycol[1],green=polycol[2],blue=polycol[3],alpha=polycol[4]),border=NA)
      lines(tx, y, col=col[i], lwd=c(1))
    }
    else {
      lines(tx, y, col=col[i], lwd=c(3))
    }
  }

  if(do.qqplot){
    if(interactive() && (grepl(pattern="X11", x=names(dev.cur()), ignore.case=TRUE) || grepl(pattern="null", x=names(dev.cur()), ignore.case=TRUE))) {
      X11()
    }
    s <- seq(0,1,1./length(mfit$t[,1]))
    x <- qchisq(p=s, df=mfit$dof, ncp=mfit$chisq)
    qqplot(x=x, y=mfit$t[, length(mfit$t[1,])], xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ-Plot non-central Chi^2 Values")
  }
}


summary.matrixfit <- function(mfit) {
  cat("\n ** Result of one state exponential fit **\n\n")
  cat("based on", mfit$N, "measurements\n")
  cat("time range from", mfit$t1, " to ", mfit$t2, "\n")
  cat("mass:\n")
  cat("m \t=\t", mfit$opt.res$par[1], "\n")
  cat("dm\t=\t", sd(mfit$t[,1]), "\n")
  cat("\nAmplitudes:\n")
  for(i in 2:length(mfit$opt.res$par)) {
    cat("P",i-1,"\t=\t", mfit$opt.res$par[i], "\n")
    cat("dP",i-1,"\t=\t", sd(mfit$t[,i]), "\n")
  }
  cat("\n")
  cat("boot.R\t=\t", mfit$boot.R, " (bootstrap samples)\n")
  cat("boot.l\t=\t", mfit$boot.l, " (block length)\n")
  cat("useCov\t=\t", mfit$useCov, "\n")
  cat("chisqr\t=\t", mfit$chisqr, "\ndof\t=\t", mfit$dof, "\nchisqr/dof=\t",
      mfit$chisqr/mfit$dof, "\n")
  ## probability to find a larger chi^2 value
  ## if the data is generated again with the same statistics
  ## given the model is correct
  cat("Quality of the fit (p-value):", mfit$Qval, "\n")

  if(any(names(mfit) == "fps")) {
    cat("\nDecay Constant (derived quantity):\n")
    cat("mu1 \t=\t", mfit$mu1, "\n")
    cat("mu2 \t=\t", mfit$mu2, "\n")
    if(mfit$normalisation == "cmi") cat("kappa\t=\t", mfit$kappa,"\n")
    cat("fps \t=\t", mfit$fps, "\n")
    cat("dfps\t=\t", sd(mfit$fps.tsboot), "\n")
  }
  if(any(names(mfit) == "fpsOS")) {
    cat("\nOS Decay Constant (derived quantity):\n")
    if(mfit$normalisation == "cmi") cat("kappa\t=\t", mfit$kappa,"\n")
    cat("fps \t=\t", mfit$fpsOS, "\n")
    cat("dfps\t=\t", sd(mfit$fpsOS.tsboot), "\n")
    cat("using\n")
    cat("ZA  \t=\t", mfit$ZA, "\n")
    cat("dZA \t=\t", sd(mfit$ZAboot), "\n")
  }
}

fit.formatrixboot <- function(cf, par, t, M, LM, T, parind, sign.vec, ov.sign.vec, lm.avail=FALSE, fitfn, dfitfn, deltat=1) {
  if(lm.avail && !missing(LM)) {
    opt.res <- nls.lm(par = par, fn = fitfn, t=t, y=cf, L=LM, T=T, parind=parind, sign.vec=sign.vec,
                      deltat=deltat, ov.sign.vec=ov.sign.vec,
                      control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500))
    opt.res$value <- opt.res$rsstrac[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fitfn, gr = dfitfn,
                     method="BFGS", control=list(maxit=500, parscale=par, REPORT=50),
                     t=t, y=cf, M=M, T=T, parind=parind, sign.vec=sign.vec, deltat=deltat,
                     ov.sign.vec=ov.sign.vec)
  }
  ##opt.res <- optim(opt.res$par, fn = matrixChisqr, gr = dmatrixChisqr,
  ##                 method="BFGS", control=list(maxit=500, parscale=opt.res$par, REPORT=50),
  ##                 t=t, y=apply(cf,2,mean), M=M, T=T, parind=parind, sign.vec=sign.vec)
  return(c(opt.res$par, opt.res$value))
}


subtract.excitedstates <- function(cf, mfit, from.samples=FALSE) {

  if(inherits(cf, "cf") && inherits(mfit, "matrixfit")) {
    ## we only subtract for 0 <= t < t1 (mind the +1 for the index convention)
    t1p1 <- 1
    t2p1 <- mfit$t1
    ii <- c(t1p1:t2p1)
    Thalfp1 <- cf$Time/2+1
    if(mfit$mSize > 1) {	
      for(j in 2:mfit$mSize) {
        ii <- c(ii, (t1p1+(j-1)*Thalfp1):(t2p1+(j-1)*Thalfp1))
      }
    }

    tt <- mfit$CF$t[ii]
    ## compute the difference of mean data to model at times smaller than fit range
    dz <- mfit$cf$cf0[ii] - matrixModel(mfit$opt.res$par, tt, cf$Time, mfit$parind[ii,], mfit$sign.vec[ii])
    cf$subtracted.values <- dz
    cf$subtracted.ii <- ii
    for(i in 1:length(cf$cf[,1])) {
      cf$cf[i,ii] <- mfit$cf$cf[i,ii]-dz
    }
    if(from.samples && cf$boot.samples) {
      cf$cf0[ii] <- matrixModel(mfit$opt.res$par, tt, cf$Time, mfit$parind[ii,], mfit$sign.vec[ii])
      for(i in 1:cf$boot.R) {
        cf$cf.tsboot$t[i,ii] <- matrixModel(mfit$t[i, c(1:length(mfit$opt.res$par))],
                                            tt, cf$Time, mfit$parind[ii,], mfit$sign.vec[ii])
      }
    }
    else{
      cf$boot.sample <- FALSE
      cf$boot.R <- NULL
      cf$boot.l <- NULL
    }
    return(cf)
  }
  else {
    stop("subtract.excitedstates: cf must be of class cf and mfit of class matrixfit. Aborting...\n")
  }
}
