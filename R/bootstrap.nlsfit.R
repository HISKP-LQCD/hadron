bootstrap.nlsfit <- function(fn,
                             gr = NULL,
                             dfn = NULL,
                             par.guess,
                             errormodel="yerrors",
                             sim="parametric",
                             boot.R,
                             y,
                             dy,
                             x,
                             dx = NULL,
                             bsamples,
                             useCov = FALSE,
                             CovMatrix,
                             use.minpack.lm = TRUE,
                             parallel = FALSE,
                             ...) {

  if(missing(y) || missing(x) || missing(boot.R) || missing(par.guess) || missing(fn)) {
    stop("x, y, par.guess, boot.R and fn must be provided!")
  }
  if(missing(dy)){
    if(!useCov && sim == "parametric"){
      stop("dy has to be provided!")
    }
    if(missing(bsamples) && missing(CovMatrix)){
      stop("dy, bsamples or CovMatrix has to be provided!")
    }
  }

  if(use.minpack.lm){
    lm.avail <- require(minpack.lm)
  }else{
    lm.avail <- FALSE
  }
  if(parallel){
    parallel <- require(parallel)
  }

  crr <- c(1:(boot.R+1))
  rr <- c(2:(boot.R+1))

  ## cast y and dy to Y and dY, respectively
  if(errormodel == "xyerrors") {
    Y <- c(y, x)
    par.Guess <- c(par.guess, x)
  }else{
    Y <- y
    par.Guess <- par.guess
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
  if(!useCov){
    if(missing(dy)){
      dY <- apply(X=bsamples, MARGIN=2, FUN=sd)
      dy <- dY[1:(length(y))]
      if(errormodel == "xyerrors"){
        dx <- dY[(length(y)+1):length(Y)]
      }
    }else{
      ## cast dy to dY
      if(errormodel == "xyerrors") {
        dY <- c(dy, dx)
      }else{
        dY <- dy
      }
    }
  }

  ## generate bootstrap samples if needed
  ## and invert covariance matrix, if applicable
  if(useCov) {
    inversion.worked <- function(InvCovMatrix) {
      if(inherits(InvCovMatrix, "try-error")) {
        stop("Variance-covariance matrix could not be inverted!")
      }
    }
    if(missing(CovMatrix)){
      InvCovMatrix <- try(invertCovMatrix(bsamples, boot.l=1, boot.samples=TRUE), silent=TRUE)
      inversion.worked(InvCovMatrix)
      dY <- chol(InvCovMatrix)
    }else{
      CholCovMatrix <- chol(CovMatrix)
      InvCovMatrix <- try(solve(CholCovMatrix), silent=TRUE)
      inversion.worked(InvCovMatrix)
      dY <- t(InvCovMatrix)
      if(sim == "parametric"){
        std.norm.dist <- matrix(rnorm(boot.R*length(Y)), nrow=boot.R)
        bsamples[rr,] <- matrix(rep(Y,boot.R), nrow=boot.R, byrow=TRUE)
        bsamples[rr,] <- bsamples[rr,] + std.norm.dist %*% CholCovMatrix
      }
    }
    if(missing(dy)){
      dy <- 1./diag(dY)[1:(length(y))]
    }
  }
  else{
    if(sim == "parametric") {
      if(parallel){
        bsamples[rr,] <- mcmapply(rnorm, n=boot.R, mean = Y, sd = dY)
      }else{
        bsamples[rr,] <- mapply(rnorm, n=boot.R, mean = Y, sd = dY)
      }
    }
    dY <- 1./dY
  }

  ## define the chi-vector, the sum of squares of which has to be minimized
  ## the definitions depend on the errormodel and the use of covariance
  ## BUT it always has the same name
  if(errormodel == "yerrors"){
    if(useCov){
      fitchi <- function(y, par) { dY %*% (y - fn(par=par, x=x, ...)) }
    }else{
      fitchi <- function(y, par) { dY * (y - fn(par=par, x=x, ...)) }
    }
  }else{
    ipx <- length(par.Guess)-seq(nx-1,0)
    if(useCov){
      fitchi <- function(y, par) { dY %*% (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx])) }
    }else{
      fitchi <- function(y, par) { dY * (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx])) }
    }
  }

  ## define the derivatives of chi and chi^2
  if(is.null(gr) || (errormodel == "xyerrors" && is.null(dfn))){
    ## in case no derivative is known, the functions are set to NULL
    ## this is the default in the optimization functions anyway
    dfitchi <- NULL
    dfitchisqr <- NULL
  }else{
    ## the format of gr has to be nrows=length(par), ncols=length(Y)
    if(errormodel == "yerrors"){
      if(useCov){
        dfitchi <- function(par, ...) { -dY %*% gr(par=par, x=x, ...) }
      }else{
        dfitchi <- function(par, ...) { -dY * gr(par=par, x=x, ...) }
      }
    }else{
      jacobian <- function(par) {
        df.dpar <- rbind(gr(par=par[-ipx], x=par[ipx], ...), array(0,dim=c(nx,length(par.guess))))
        df.dx <- rbind(diag(dfn(par=par[-ipx], x=par[ipx], ...)), diag(1,nx))
        return(cbind(df.dpar, df.dx))
      }
      if(useCov){
        dfitchi <- function(par, ...) { -dY %*% jacobian(par) }
      }else{
        dfitchi <- function(par, ...) { -dY * jacobian(par) }
      }
    }
    dfitchisqr <- function(y, par) { 2 * crossprod(fitchi(y, par), dfitchi(par)) }
  }

  ## define the wrapper-functions for optimization
  if(lm.avail){
    wrapper <- function(y, par) {
      res <- nls.lm(par=par, fn=fitchi, y=y, jac=dfitchi,
                    control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxfev=10000, maxiter=500))
      if( !(res$info %in% c(1,2,3) ) ){
        cat(sprintf("Termination wrapper reason of nls.lm res$info: %d\n", res$info))
      }
      return (c(res$par, res$rsstrace[length(res$rsstrace)]))
    }
  }else{
    fitchisqr <- function(y, par) { sum(fitchi(y, par)^2) }
    wrapper <- function(y, par) {
      res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("BFGS"))
      return (c(res$par, res$value))
    }
  }

  ## now the actual fit is performed
  first.res <- wrapper(Y, par.Guess)[1:(length(par.Guess))]
  if(parallel){
    boot.list <- mclapply(crr, function(sample) { wrapper(y=bsamples[sample,], par=first.res) })
    boot.res <- do.call(cbind, boot.list)
  }else{
    boot.res <- apply(X=bsamples, MARGIN=1, FUN=wrapper, par=first.res)
  }
  ## the error calculation (this could be parallelized as well, but I don't think this would help)
  errors <- apply(X=boot.res[1:(length(par.Guess)),rr, drop=FALSE], MARGIN=1, FUN=sd)

  res <- list(y=y, dy=dy, x=x, dx=dx, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R, sim=sim,
              bsamples=bsamples[rr,],
              errormodel=errormodel,
              t0=boot.res[,1],
              t=t(boot.res[,rr]),
              se=errors,
              useCov=useCov,
              invCovMatrix=dY,
              Qval = 1-pchisq(boot.res[length(par.Guess)+1,1], length(y) - length(par.guess)),
              chisqr = boot.res[length(par.Guess)+1,1],
              dof = length(y) - length(par.guess),
              tofn=list(...))
  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

summary.bootstrapfit <- function(object, digits=2, print.correlation=TRUE, ...) {

  cat("bootstrap nls fit\n\n")
  cat("model", object$errormodel, "\n")
  errors <- object$se
  values <- object$t0[1:(length(object$t0)-1)]
  npar <- length(object$par.guess)
  
  ## parameters with errors as strings
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, with.dollar=FALSE, digits=digits, human.readable=FALSE)
  bias <- object$t0[1:(length(object$t0)-1)]-apply(X=object$t[,1:(length(object$t0)-1), drop=FALSE], MARGIN=2, FUN=mean)
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
  if(print.correlation){
    correlation <- cor(object$t[,1:(length(object$t0)-1)], object$t[,1:(length(object$t0)-1)])
    cat("\n   correlation matrix of the fit parameters\n\n")
    print(data.frame(correlation))
  }
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

plot.bootstrapfit <- function(x, ..., xlim, ylim, rep=FALSE, col.line="black", col.band="gray", opacity.band=0.65, lty=c(1), lwd=c(1), xlab="x", ylab="y", supports=1000, plot.range) {
  if(missing(plot.range)){
    rx <- range(x$x)
  }else{
    rx <- plot.range
  }
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/supports)
  npar <- length(x$par.guess)

  ## use the xylimits computation of plotwitherror
  if(x$errormodel == "yerrors") {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, rep=rep, xlim=xlim, ylim=ylim, ...)
  }
  else {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, rep=rep, xlim=xlim, ylim=ylim, ...)
  }
  my.xlim <- mylims$xlim
  my.ylim <- mylims$ylim

  ## to include additional parameter to x$fn originally given as ... to bootstrap.nlsfit
  ## requires some pull-ups
  Y <- numeric()
  Y <- do.call(what=x$fn, args=c(list(par=x$t0[1:npar], x=X), x$tofn))

  ## error band
  ## define a dummy function to be used in apply
  dummyfn <- function(par, x, object) {
    return(do.call(what=object$fn, args=c(list(par=par, x=x), object$tofn)))
  }
  se <- apply(X=rbind(apply(X=x$t[, c(1:npar), drop=FALSE], MARGIN=1, FUN=dummyfn, x=X, object=x)), MARGIN=1, FUN=sd)

  ## plot it
  polyval <- c(Y+se, rev(Y-se))
  if(any(polyval < my.ylim[1]) || any(polyval > my.ylim[2])) {
    polyval[polyval < my.ylim[1]] <- my.ylim[1]
    polyval[polyval > my.ylim[2]] <- my.ylim[2]
  }
  pcol <- col2rgb(col.band, alpha=TRUE)/255 
  pcol[4] <- opacity.band
  pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
  polygon(x=c(X, rev(X)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)

  ## plot the fitted curve on top
  lines(x=X, y=Y, col=col.line, lty=lty, lwd=lwd)
}
