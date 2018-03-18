bootstrap.nlsfit <- function(fn,
                             par.guess,
                             errormodel="yerrors",
                             sim="parametric",
                             boot.R,
                             y,
                             dy,
                             x,
                             dx=NULL,
                             ...) {

  ## the chi functions for nls.lm
  fitchi <- function(par, x, y, dy, fitfun) {
    (y - fitfun(par=par, x=x))/dy
  }
  fitchi.xy <- function(par, y, dy, fitfun, nx) {
    ipx <- length(par)-seq(nx-1,0)
    (y - c(fitfun(par=par[-ipx], x=par[ipx]), par[ipx]))/dy
  }

  ## the corresponding chisqr functions
  fitchisqr <- function(par, x, y, dy, fitfun) {
    z <- fitchi(par=par, x=x, y=y, dy=dy, fitfun=fitfun)
    return (z %*% z)
  }
  fitchisqr.xy <- function(par, y, dy, fitfun, nx) {
    z <- fitchixy(par=par, y=y, dy=dy, fitfun=fitfun, nx=nx)
    return (z %*% z)
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
  bsamples <- array(NA, dim=c(boot.R+1, length(Y)))
  bsamples[1,] <- Y
  
  if(sim == "parametric") {
    for(i in seq_along(Y)) {
      bsamples[rr, i] <- rnorm(n=boot.R, mean = Y[i], sd = dY[i])
    }
  }
  ## else missing currenlty!
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
              ysamples=bsamples,
              errormodel=errormodel,
              t0=boot.res[,1],
              t=boot.res,
              Qval = 1-pchisq(boot.res[dim(boot.res)[1],1], length(par.guess)),
              chisqr = boot.res[dim(boot.res)[1],1],
              dof = length(y) - length(par.guess))
  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

summary.bootstrapfit <- function(x, digits=2) {

  cat("bootstrap nls fit\n\n")
  cat("model", x$errormodel, "\n")
  errors <- apply(X=x$t[1:(dim(x$t)[1]-1), ], MARGIN=1, FUN=sd)
  values <- x$t[1:(dim(x$t)[1]-1), 1]
  npar <- length(x$par.guess)
  
  ## parameters with errors as strings
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, with.dollar=FALSE, digits=2)
  cat("    best fit parameters with errors\n")
  print(data.frame(par=tmp[1:npar]))
  if(x$errormodel != "yerrors") {
    cat("\n estimates for x-values with errors\n")
    print(data.frame(x=tmp[(npar+1):length(tmp)]))
  }
  cat("\n   chi^2 and fit quality\n")
  cat("chisqr / dof =", x$chisqr, "/", x$dof, "=", x$chisqr/x$dof, "\n")
  cat("p-value", x$Qval, "\n")
}
