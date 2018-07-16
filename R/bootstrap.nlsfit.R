bootstrap.nlsfit <- function(fn,
                             gr = NULL,
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
  fitchi <- function(par, x, y, dy, fitfun, dfitfun=NULL, ...) {
    (y - fitfun(par=par, x=x, ...)) %*% dy
  }
  ## checked
  dfitchi <- function(par, x, y, dy, dfitfun, fitfun=NULL, ...) {
    -as.vector(t(crossprod(dfitfun(par=par, x=x, ...), dy)))
  }
  fitchi.xy <- function(par, y, dy, fitfun, dfitfun=NULL, nx, ...) {
    ipx <- length(par)-seq(nx-1,0)
    (y - c(fitfun(par=par[-ipx], x=par[ipx], ...), par[ipx])) %*% dy
  }
  ## not checked yet
  dfitchi.xy <- function(par, x, y, dy, dfitfun, fitfun=NULL, nx, ...) {
    ipx <- length(par)-seq(nx-1,0)
    -as.vector(t(crossprod(dfitfun(par=par[-ipx], x=x[ipx], ...), dy)))
    ## pieces missing
  }

  ## the corresponding chisqr functions
  fitchisqr <- function(par, x, y, dy, fitfun, dfitfun, ...) {
    z <- fitchi(par=par, x=x, y=y, dy=dy, fitfun=fitfun, dfitfun=NULL, ...)
    return (sum(z * z))
  }
  dfitchisqr <- function(par, x, y, dy, dfitfun, fitfun, ...) {
    z <- fitchi(par=par, x=x, y=y, dy=dy, fitfun=fitfun, dfitfun=dfitfun, ...)
    dz <- dfitchi(par=par, x=x, y=y, dy=dy, dfitfun=dfitfun, fitfun=fitfun, ...)
    return(2*(z %*% array(dz, dim=c(length(x), length(par)))))
  }
  fitchisqr.xy <- function(par, y, dy, fitfun, nx, ...) {
    z <- fitchi.xy(par=par, y=y, dy=dy, fitfun=fitfun, nx=nx, ...)
    return (sum(z * z))
  }

  ## wrapper functions for apply
  wrapper.lm <- function(y, par, fitfun, dfitfun, dy, x, ...) {
    res <- list()
    if(is.null(dfitfun))  res <- minpack.lm::nls.lm(par=par, fn=fitchi, y=y, fitfun=fitfun, dy=dy, x=x,
                                                    control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500), ...)
    else res <- minpack.lm::nls.lm(par=par, fn=fitchi, y=y, fitfun=fitfun, jac=dfitchi, dy=dy, x=x, dfitfun=dfitfun,
                                   control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500), ...)
    if( !(res$info %in% c(1,2,3) ) ){
      cat(sprintf("Termination wrapper.lm reason of nls.lm res$info: %d\n", res$info))
    }
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }
  wrapper.lm.xy <- function(y, par, fitfun, dy, nx, ...) {
    res <- list()
    res <- minpack.lm::nls.lm(par=par, fn=fitchi.xy, y=y, fitfun=fitfun, dy=dy, nx=nx,
                              control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500),
                              ...)
    if( !(res$info %in% c(1,2,3) ) ){
      cat(sprintf("Termination wrapper.lm.xy reason of nls.lm res$info: %d\n", res$info))
    }
    return (c(res$par, res$rsstrace[length(res$rsstrace)]))
  }

  wrapper.optim <- function(y, par, fitfun, dy, x, dfitfun, ...) {
    res <- NA
    if(is.null(dfitfun)) res <- optim(par=par, fn=fitchisqr, y=y, method=c("BFGS"), fitfun=fitfun, dy=dy, x=x, ...)
    else res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("BFGS"), fitfun=fitfun, dfitfun=dfitfun, dy=dy, x=x, ...)

    return (c(res$par, res$value))
  }
  wrapper.optim.xy <- function(y, par, fitfun, dy, nx, ...) {
    res <- optim(par=par, fn=fitchisqr.xy, y=y, method=c("BFGS"), fitfun=fitfun, dy=dy, nx=nx, ...)
    return (c(res$par, res$value))
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
  else {
    dY <- diag(1./dY)
  }

  if(errormodel == "yerrors") {
    if(lm.avail) {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.lm, x=x, dy=dY, par=par.guess, fitfun=fn, dfitfun=gr, ...)
    }
    else {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.optim, x=x, dy=dY, par=par.guess, fitfun=fn, dfitfun=gr, ...)
    }
  }
  else { ## xyerrors
    if(lm.avail) {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.lm.xy, nx=nx, dy=dY, par=c(par.guess, x), fitfun=fn, ...)
    }
    else {
      boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper.optim.xy, nx=nx, dy=dY, par=c(par.guess, x), fitfun=fn, ...)
    }
  }

  errors <- apply(X=boot.res[1:(dim(boot.res)[1]-1),rr], MARGIN=1, FUN=sd)

  res <- list(y=y, dy=dy, x=x, dx=dx, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R, sim=sim,
              bsamples=bsamples,
              errormodel=errormodel,
              t0=boot.res[,1],
              t=t(boot.res),
              se=errors,
              useCov=useCov,
              invCovMatrix=dY,
              Qval = 1-pchisq(boot.res[dim(boot.res)[1],1], length(y) - length(par.guess)),
              chisqr = boot.res[dim(boot.res)[1],1],
              dof = length(y) - length(par.guess),
              tofn=list(...))
  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

summary.bootstrapfit <- function(object, digits=2, ...) {

  cat("bootstrap nls fit\n\n")
  cat("model", object$errormodel, "\n")
  errors <- object$se
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

plot.bootstrapfit <- function(x, ..., xlim, ylim, rep=FALSE, col.line="black", col.band="gray", lty=c(1), lwd=c(1), xlab="x", ylab="y") {
  rx <- range(x$x)
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/1000)
  npar <- length(x$par.guess)
  Y <- numeric()
  my.xlim <- numeric()
  my.ylim <- numeric()
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
    
    plot(NA, xlim=my.xlim, ylim=my.ylim, xlab=xlab, ylab=ylab, ...)
  }

  ## to include additional parameter to x$fn originally given as ... to bootstrap.nlsfit
  ## requires some pull-ups
  Y <- do.call(what=x$fn, args=c(list(par=x$t0[1:npar], x=X), x$tofn))

  ## error band
  ## define a dummy function to be used in apply
  dummyfn <- function(par, x, object) {
    return(do.call(what=object$fn, args=c(list(par=par, x=x), object$tofn)))
  }
  se <- apply(X=apply(X=x$t[, c(1:npar)], MARGIN=1, FUN=dummyfn, x=X, object=x), MARGIN=1, FUN=sd)

  ## plot it
  polyval <- c(Y+se, rev(Y-se))
  if(any(polyval < my.ylim[1]) || any(polyval > my.ylim[2])) {
    polyval[polyval < my.ylim[1]] <- my.ylim[1]
    polyval[polyval > my.ylim[2]] <- my.ylim[2]
  }

  polygon(x=c(X, rev(X)), y=polyval, col=col.band, lty=0, lwd=0.001, border=col.band)
  ## plot the fitted curve on top
  lines(x=X, y=Y, col=col.line, lty=lty, lwd=lwd)

  ## plot data on top of everything
  if(x$errormodel == "yerrors") {
    plotwitherror(x=x$x, y=x$y, dy=x$dy, rep=TRUE)
  }
  else {
    plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, rep=TRUE)
  }
}

  
