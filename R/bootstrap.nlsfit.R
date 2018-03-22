bootstrap.nlsfit <- function(fn,
                             par.guess,
                             errormodel="yerrors",
                             sim="parametric",
                             boot.R,
                             y,
                             dy,
                             x,
                             dx=NULL,
                             bsamples,
                             useCov,
                             ...) {

  if(missing(y) || missing(dy) || missing(x) || missing(boot.R) || missing(par.guess) || missing(fn) || missing(useCov)) {
    stop("x, y, dy, par.guess, boot.R, useCov and fn must be provided!")
  }

  ## the chi functions for nls.lm
  fitchi <- function(par, x, y, dy, fitfun) {
    (y - fitfun(par=par, x=x)) %*% dy
  }
  fitchi.xy <- function(par, y, dy, fitfun, nx) {
    ipx <- length(par)-seq(nx-1,0)
    (y - c(fitfun(par=par[-ipx], x=par[ipx]), par[ipx])) %*% dy
  }

  ## the corresponding chisqr functions
  fitchisqr <- function(par, x, y, dy, fitfun) {
    z <- fitchi(par=par, x=x, y=y, dy=dy, fitfun=fitfun)
    return (sum(z * z))
  }
  fitchisqr.xy <- function(par, y, dy, fitfun, nx) {
    z <- fitchi.xy(par=par, y=y, dy=dy, fitfun=fitfun, nx=nx)
    return (sum(z * z))
  }

  ## wrapper functions for apply
  wrapper.lm <- function(y, par, fitfun, dy, x, ...) {
    res <- nls.lm(par=par, fn=fitchi, y=y, fitfun=fitfun, dy=dy, x=x, ...)
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }
  wrapper.lm.xy <- function(y, par, fitfun, dy, nx, ...) {
    res <- nls.lm(par=par, fn=fitchi.xy, y=y, fitfun=fitfun, dy=dy, nx=nx, ...)
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }

  wrapper.optim <- function(y, par, fitfun, dy, x, ...) {
    res <- optim(par=par, fn=fitchisqr, y=y, method=c("BFGS"), fitfun=fitfun, dy=dy, x=x, ...)
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }
  wrapper.optim.xy <- function(y, par, fitfun, dy, nx, ...) {
    res <- optim(par=par, fn=fitchisqr.xy, y=y, method=c("BFGS"), fitfun=fitfun, dy=dy, nx=nx, ...)
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }

  crr <- c(1:(boot.R+1))
  rr <- c(2:(boot.R+1))
  lm.avail <- require(minpack.lm)

  ## cast y and dy to Y and dY, respectively
  Y <- y
  dY <- dy
  if(errormodel == "xyerrors") {
    Y <- c(y, x)
    dY <- c(dy, dx)
  }
  nx <- length(x)
  
  ## the bootstrap samples of the input data
  ## in an array
  if(missing(bsamples)) {
    bsamples <- array(NA, dim=c(boot.R+1, length(Y)))
    bsamples[1,] <- Y
  }
  else {
    ## check consistency
    dbs <- dim(bsamples)
    if(dbs[2] != length(Y)) {
      stop("the provided bootstrap samples do not match the number of data points with errors!")
    }
    if(boot.R != dbs[1]) {
      stop("boot.R inconsistent with dimension one of bsamples!")
    }
    ## add original data as first row
    bsamples <- rbind(Y, bsamples)
  }

  ## generate bootstrap samples if needed
  ## and invert covariance matrix, if applicable
  if(sim == "parametric") {
    for(i in seq_along(Y)) {
      bsamples[rr, i] <- rnorm(n=boot.R, mean = Y[i], sd = dY[i])
    }
    dY <- diag(1./dY)
  }
  else if(useCov) {
    CovMatrix <- cov(bsamples)
    InvCovMatrix <- try(invertCovMatrix(bsamples, boot.l=1, boot.samples=TRUE), silent=TRUE)
    if(inherits(InvCovMatrix, "try-error")) {
      stop("Variance-covariance matrix could not be inverted!")
    }
    dY <- chol(InvCovMatrix)
  }

  if(errormodel == "yerrors") {
    if(lm.avail) {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.lm, x=x, dy=dY, par=par.guess, fitfun=fn)
    }
    else {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.optim, x=x, dy=dY, par=par.guess, fitfun=fn)
    }
  }
  else { ## xyerrors
    if(lm.avail) {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.lm.xy, nx=nx, dy=dY, par=c(par.guess, x), fitfun=fn)
    }
    else {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.optim.xy, nx=nx, dy=dY, par=c(par.guess, x), fitfun=fn)
    }
  }
  res <- list(y=y, dy=dy, x=x, dx=dx, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R, sim=sim,
              bsamples=bsamples,
              errormodel=errormodel,
              t0=boot.res[,1],
              t=t(boot.res),
              useCov=useCov,
              invCovMatrix=dY,
              Qval = 1-pchisq(boot.res[dim(boot.res)[1],1], length(par.guess)),
              chisqr = boot.res[dim(boot.res)[1],1],
              dof = length(y) - length(par.guess))
  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

summary.bootstrapfit <- function(object, digits=2, ...) {

  cat("bootstrap nls fit\n\n")
  cat("model", object$errormodel, "\n")
  errors <- apply(X=object$t[,1:(dim(object$t)[2]-1)], MARGIN=2, FUN=sd)
  values <- object$t[1, 1:(dim(object$t)[2]-1)]
  npar <- length(object$par.guess)
  
  ## parameters with errors as strings
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, with.dollar=FALSE, digits=digits, human.readable=FALSE)
  bias <- object$t0[1:(length(object$t0)-1)]-apply(X=object$t[,1:(dim(object$t)[2]-1)], MARGIN=2, FUN=mean)
  dim(bias) <- c(length(bias), 1)
  bias <- apply(X=bias, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  ci16 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.16), drop=FALSE)
  dim(ci16) <- c(length(ci16), 1)
  ci16 <- apply(X=ci16, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  ci84 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.84), drop=FALSE)
  dim(ci84) <- c(length(ci84), 1)
  ci84 <- apply(X=ci84, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  cat("    best fit parameters with errors, bootstrap bias and 68% confidence interval\n\n")
  print(data.frame(par=tmp[1:npar], bias=bias[1:npar], ci16=ci16[1:npar], ci84=ci84[1:npar]))
  correlation <- cor(object$t[,1:(dim(object$t)[2]-1)], object$t[,1:(dim(object$t)[2]-1)])
  cat("\n   correlation matrix of the fit parameters\n\n")
  print(data.frame(correlation))
  if(object$errormodel != "yerrors") {
    cat("\n estimates for x-values with errors, bootstrap bias and 68% confidence interval\n\n")
    ii <- c((npar+1):length(tmp))
    print(data.frame(x=tmp[ii], bias=bias[ii], ci16=ci16[ii], ci84=ci84[ii]))
  }
  cat("\n   chi^2 and fit quality\n")
  cat("chisqr / dof =", object$chisqr, "/", object$dof, "=", object$chisqr/object$dof, "\n")
  cat("p-value", object$Qval, "\n")
}

print.bootstrapfit <- function(x, digits=2, ...) {
  summary.bootstrapfit(object=x, digits=digits, ...)
}

plot.bootstrapfit <- function(x, ..., xlim, ylim, rep=FALSE, col.line="black", col.band="gray", lty=c(1), lwd=c(1)) {
  rx <- range(x$x)
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/1000)
  npar <- length(x$par.guess)
  Y <- numeric()
  if(!missing(xlim)) {
    my.xlim <- xlim
  }
  if(!missing(ylim)) {
    my.ylim <- ylim
  }
  if(!rep) {
    mylims <- c()
    ## use the xylimits computation of plotwitherror
    if(missing(xlim) || missing(ylim)) {
      if(x$errormodel == "yerrors") {
        mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, ...)
      }
      else {
        mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, ...)
      }
    }
    if(missing(xlim)) {
      my.xlim <- mylims$xlim
    }
    if(missing(ylim)) {
      my.ylim <- mylims$ylim
    }
    
    ## generate empty plot
    
    plot(NA, xlim=my.xlim, ylim=my.ylim, ...)
  }
  
  if(x$errormodel == "yerrors") {
    Y <- x$fn(par=x$t0, x=X)
  }
  else {
    Y <- x$fn(par=x$t0[1:npar], x=X)
  }
  ## error band
  se <- apply(X=apply(X=x$t[, c(1:npar)], MARGIN=1, FUN=x$fn, x=X), MARGIN=1, FUN=sd)
  polygon(x=c(X, rev(X)), y=c(Y+se, rev(Y-se)), col=col.band, lty=0, lwd=0.001, border=col.band)
  ## fitted curve
  lines(x=X, y=Y, col=col.line, lty=lty, lwd=lwd)

  ## plot data with fitted curve on top of everything
  if(x$errormodel == "yerrors") {
    plotwitherror(x=x$x, y=x$y, dy=x$dy, rep=TRUE, ...)
  }
  else {
    plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, rep=TRUE,...)
  }
}

  
