#' Parametric bootstrap
#'
#' @param boot.R numeric. Number of bootstrap samples to generate.
#' @param x numeric vector. Actual values for the data.
#' @param dx numeric vector of the same length as `x` or missing. Errors of the
#'   values.
#'
#' @return
#' A matrix with as many columns as there are variables in `x` and as many rows
#' as `boot.R`.
#'
#' @export
#' @family NLS fit functions
#'
#' @examples
#' x <- 1:3
#' dx <- 1:3 * 0.1
#' parametric.bootstrap(5, x, dx)
parametric.bootstrap <- function (boot.R, x, dx) {
  stopifnot(length(x) == length(dx))

  samples.list <- lapply(dx, function (error) rnorm(boot.R, 0, error))
  samples <- do.call(cbind, samples.list)

  samples <- t(t(samples) + x)

  return (samples)
}

#' Parametric bootstrap with covariance
#'
#' @param boot.R numeric. Number of bootstrap samples to generate.
#' @param x numeric vector. Actual values for the data.
#' @param cov numeric matrix, square, length of `x` or missing. Covariance
#'   between the various variables in the vector `x`.
#'
#' @return
#' A matrix with as many columns as there are variables in `x` and as many rows
#' as `boot.R`.
#'
#' @export
#' @family NLS fit functions
#'
#' @examples
#' x <- 1:3
#' cov <- matrix(c(0.1, 0, 0.01,
#'                 0, 0.15, 0.02,
#'                 0.01, 0.02, 0.2), nrow = 3)
#' parametric.bootstrap.cov(5, x, cov)
parametric.bootstrap.cov <- function (boot.R, x, cov) {
  stopifnot(nrow(cov) == length(x))
  stopifnot(ncol(cov) == length(x))
  stopifnot(isSymmetric(cov))

  evv <- eigen(cov)
  evalues <- evv$values
  evectors <- evv$vectors

  ## This is the original matrix now:
  ## evectors %*% diag(evalues) %*% t(evectors) == cov

  errors <- sqrt(evalues)
  stopifnot(!any(is.complex(errors)))
  samples.list <- lapply(errors, function (error) rnorm(boot.R, 0, error))
  samples <- do.call(cbind, samples.list)

  samples <- samples %*% t(evectors)
  samples <- t(t(samples) + x)

  return (samples)
}

#' NLS fit with parametric bootstrap
#'
#' @export
#' @family NLS fit functions
#' 
#' @examples
#'
#' ## Declare some data.
#' value <- c(0.1, 0.2, 0.3)
#' dvalue <- c(0.01, 0.01, 0.015)
#' x <- c(1, 2, 3)
#' dx <- c(0.1, 0.1, 0.1)
#' boot.R <- 1500
#'
#' fn <- function (par, x) par[1] + par[2] * x
#'
#' fit.result <- parametric.nlsfit(fn, c(1, 1), boot.R, value, dvalue, x, dx)
#' summary(fit.result)
parametric.nlsfit <- function (fn, par.guess, boot.R, y, dy, x, dx, ...) {
  stopifnot(length(x) == length(y))
  stopifnot(missing(dx) || length(dx) == length(x))
  stopifnot(missing(dy) || length(dy) == length(y))

  if (missing(dx)) {
    values <- y
    errors <- dy
  } else {
    values <- c(y, x)
    errors <- c(dy, dx)
  }

  bsamples <- parametric.bootstrap(boot.R, values, errors)

  bootstrap.nlsfit(fn, par.guess, y, x, bsamples, ..., dx = dx, dy = dy)
}

#' NLS fit with parametric bootstrap and covariance
#'
#' @export
#' @family NLS fit functions
parametric.nlsfit.cov <- function (fn, par.guess, boot.R, y, x, cov, ...) {
  stopifnot(length(x) == length(y))

  if (ncol(cov) == length(y)) {
    values <- y
  } else if (ncol(cov) == length(y) + length(x)) {
    values <- c(y, x)
  } else {
    stop('The covariance matrix must either be as large as `y` or as `y` and `x` together.')
  }

  bsamples <- parametric.bootstrap.cov(boot.R, values, cov)

  bootstrap.nlsfit(fn, par.guess, y, x, bsamples, cov, ...)
}

#' Bootstrap a non-linear least-squares fit
#'
#' Performs and bootstraps a non-linear least-squares fit to data with y and x
#' errors.
#'
#' @param fn the (non-linear) function to be fitted to the data. Its first
#' argument must be the fit parameters named \code{par}. The second must be
#' \code{x}, the explaining variable.
#' @param gr \code{gr=d(fn) / d(par)} is a function to return the gradient of
#' \code{fn}. It must return an array with \code{length(x)} rows and
#' \code{length(par)} columns.
#' @param dfn \code{dfn=d(fn) / dx} is the canonical derivative of \code{fn} by
#' \code{x} and only relevant if x-errors are provided.
#' @param par.guess initial guess values for the fit parameters.
#' @param y the data as a one-dimensional numerical vector to be described by
#' the fit function. 
#' @param x values of the explaining variable in form of a one-dimensional
#' numerical vector.
#' @param bsamples bootstrap samples of \code{y} (and \code{x}, if applicable).
#' Must be provided as array of dimensions \code{c(boot.R, n)} with \code{n}
#' equals to \code{length(y)} in case of 'yerrors' and For 'xyerrors' to
#' \code{length(y) + length(x)}.
#' @param CovMatrix complete variance-covariance matrix of dimensions
#' \code{c(length(y), length(y))} or \code{c(length(y)+length(x),
#' length(y)+length(x))} depending on the errormodel.
#' @param use.minpack.lm use the \code{minpack.lm} library if available. This
#' is usually faster than the default \code{optim} but somtimes also less
#' stable.
#' @param parallel parallelise over bootstrap samples. The package
#' \code{parallel} is required.
#'
#' @return
#'  returns a list of class 'bootstrapfit'. It returns all input
#'  parameters and adds in addition the following:
#'  \item{t0}{the one dimensional numerical vector of length
#'    \code{npar+1}. \code{npar} is the number of fit parameters. In case
#'    of 'yerrors' this equals \code{length(par.guess)}. For 'xyerrors'
#'    this equals \code{length(par.guess) + length(x)}. \code{t0} contains
#'    the best fit parameters
#'    obtained on the original data. The last element in \code{t0} is the
#'    chisquare value.}
#'  \item{t}{an array of dimensions \code{(npar+1, boot.R)} with
#'    \code{npar} as in \code{t0}. The rows contain the individual
#'    bootstrap observations.}
#'  \item{bsamples}{the bootstrap samples used as an array of dimensions
#'    \code{(length(y), boot.R)} or \code{(length(y)+length(x), boot.R)}
#'    depending on the error model with \code{npar} as in \code{t0}. }
#'  \item{Qval}{the p-value of the fit on the original data}
#'  \item{chisqr}{the residual chisqr value.}
#'  \item{dof}{the residual degrees of freedom of the fit.}
#'  \item{nx}{the number of x-values.}
#'  \item{tofn}{
#'    the original \code{...} list of parameters to be passed on to the
#'    fit function
#'  }
#'
#' @examples
#' ## Declare some data.
#' value <- c(0.1, 0.2, 0.3)
#' dvalue <- c(0.01, 0.01, 0.015)
#' x <- c(1, 2, 3)
#' dx <- c(0.1, 0.1, 0.1)
#' boot.R <- 1500
#'
#' fn <- function (par, x) par[1] + par[2] * x
#'
#' ## Before we can use the fit with this data, we need to create bootstrap
#' ## samples. We do not want to use the correlation matrix here. Note that you
#' ## can simply use the parametric.nlsfit function as a convenient wrapper of
#' ## the two steps.
#' bsamples <- parametric.bootstrap(boot.R, c(value, x), c(dvalue, dx))
#' head(bsamples)
#'
#' fit.result <- bootstrap.nlsfit(fn, c(1, 1), value, x, bsamples)
#' summary(fit.result)
#'
#' @export
#' @family NLS fit functions
bootstrap.nlsfit <- function(fn,
                             par.guess,
                             y,
                             x,
                             bsamples,
                             ...,
                             dy,
                             dx,
                             CovMatrix,
                             gr,
                             dfn,
                             use.minpack.lm = TRUE,
                             parallel = FALSE,
                             error = sd) {
  stopifnot(!missing(y))
  stopifnot(!missing(x))
  stopifnot(!missing(par.guess))
  stopifnot(!missing(fn))
  stopifnot(!missing(bsamples))

  boot.R <- nrow(bsamples)
  useCov <- !missing(CovMatrix)

  if (use.minpack.lm) {
    lm.avail <- require(minpack.lm)
  } else {
    lm.avail <- FALSE
  }

  if (parallel) {
    parallel <- require(parallel)
  }

  crr <- c(1:(boot.R+1))
  rr <- c(2:(boot.R+1))

  ## cast y and dy to Y and dY, respectively
  if (ncol(bsamples) == length(y)) {
    Y <- y
    par.Guess <- par.guess
    errormodel <- "yerrors"
  } else if (ncol(bsamples) == length(y) + length(x)) {
    Y <- c(y, x)
    par.Guess <- c(par.guess, x)
    errormodel <- "xyerrors"
  } else {
    stop("The provided bootstrap samples do not match the number of data points with errors. Make sure that the number of columns is either the length of `y` alone for just y-errors or the length of `y` and `x` for xy-errors.")
  }

  nx <- length(x)
  
  ## generate bootstrap samples if needed
  ## and invert covariance matrix, if applicable
  if (useCov) {
    if (!missing(dx) || !missing(dy)) {
      stop('Specifying a covariance matrix and `dx` and `dy` does not make sense, use either.')
    }

    inversion.worked <- function(InvCovMatrix) {
      if (inherits(InvCovMatrix, "try-error")) {
        stop("Variance-covariance matrix could not be inverted!")
      }
    }

    if (missing(CovMatrix)) {
      InvCovMatrix <- try(invertCovMatrix(bsamples, boot.l = 1, boot.samples = TRUE), silent = TRUE)
      inversion.worked(InvCovMatrix)
      dY <- chol(InvCovMatrix)
    } else {
      CholCovMatrix <- chol(CovMatrix)
      InvCovMatrix <- try(solve(CholCovMatrix), silent = TRUE)
      inversion.worked(InvCovMatrix)
      dY <- t(InvCovMatrix)
    }

    dydx <- 1.0 / diag(dY)

    if (errormodel == 'yerrors') {
      dy <- dydx
    } else {
      dy <- dydx[1:length(y)]
      dx <- dydx[(length(y)+1):length(dydx)]
    }
  }
  else {
    ## The user did not specify the errors, therefore we simply compute them.
    if (missing(dx) && missing(dy)) {
      dydx <- apply(bsamples, 2, error)
      dY <- 1.0 / dydx

      if (errormodel == 'yerrors') {
        dy <- dydx
      } else {
        dy <- dydx[1:length(y)]
        dx <- dydx[(length(y)+1):length(dydx)]
      }
    }
    ## The user has specified either one, so we need to make sure that it is
    ## consistent.
    else {
      if (errormodel == 'yerrors' && ncol(bsamples) == length(dy)) {
        dY <- 1.0 / dy
      } else if (errormodel == 'xyerrors' && ncol(bsamples) == length(dy) + length(dx)) {
        dY <- 1.0 / c(dy, dx)
      } else {
        stop('You have explicitly passed `dy` and/or `dx`, but their combined length does not match the number of columns of the bootstrap samples.')
      }
    }
  }

  ## add original data as first row
  bsamples <- rbind(Y, bsamples)

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
  if(missing(gr) || (errormodel == "xyerrors" && missing(dfn))){
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
  if (lm.avail) {
    wrapper <- function(y, par) {
      res <- nls.lm(par=par, fn=fitchi, y=y, jac=dfitchi,
                    control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxfev=10000, maxiter=500))

      if( !(res$info %in% c(1,2,3) ) ){
        cat(sprintf("Termination wrapper reason of nls.lm res$info: %d\n", res$info))
      }

      list(converged = res$info %in% 1:3,
           par = res$par,
           chisq = res$rsstrace[length(res$rsstrace)])
    }
  } else {
    fitchisqr <- function(y, par) { sum(fitchi(y, par)^2) }
    wrapper <- function(y, par) {
      res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("BFGS"))

      list(converged = res$convergence == 0,
           par = res$par,
           chisq = res$value)
    }
  }

  ## now the actual fit is performed
  first.res <- wrapper(Y, par.Guess)[1:(length(par.Guess))]

  if (parallel)
    my.lapply <- mclapply
  else {
    my.lapply <- lapply
  }

  boot.list <- my.lapply(crr, function(sample) { wrapper(y=bsamples[sample,], par=first.res$par) })

  par.boot <- do.call(rbind, lapply(boot.list, function (elem) elem$par))

  converged <- sapply(boot.list, function (elem) elem$converged)
  par.boot[!converged, ] <- NA

  chisq <- boot.list[[1]]$chisq
  dof = length(y) - length(par.guess)

  errors <- apply(X=par.boot[rr, 1:(length(par.Guess)), drop=FALSE], MARGIN=2, FUN=error)

  res <- list(y=y, dy=dy, x=x, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R,
              bsamples=bsamples[rr, ],
              errormodel=errormodel,
              converged = converged,
              t0=par.boot[1, ],
              t=par.boot[rr, ],
              se=errors,
              useCov=useCov,
              invCovMatrix=dY,
              Qval = 1 - pchisq(chisq, dof),
              chisqr = chisq,
              dof = dof,
              error.function = error,
              tofn=list(...))

  if (errormodel == 'xyerrors') {
    res$dx <- dx
  }

  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

#' Summarize a bootstrap NLS fit
#'
#' @param object object returned by \code{bootstrap.nlsfit}
#' @param digits number of significant digits to print in summary or print.
#'
#' @export
#' @family NLS fit functions
summary.bootstrapfit <- function(object, digits=2, print.correlation=TRUE) {
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

#' Print a bootstrap NLS fit
#'
#' @param x object returned by \code{bootstrap.nlsfit}
#' @param digits number of significant digits to print in summary or print.
#'
#' @family NLS fit functions
print.bootstrapfit <- function(x, digits=2) {
  summary.bootstrapfit(object=x, digits=digits)
}

#' Plot a bootstrap NLS fit
#'
#' @param x object returned by \code{bootstrap.nlsfit}
#' @param xlim x limits of the plot.
#' @param ylim y limits of the plot.
#' @param rep If set to \code{TRUE}, operate like "replot" in gnuplot. Allows
#' adding points with error bars to the current plot. Switches the underlying
#' plotting routine from \code{plot} to \code{points}.
#' @param col.line line colour.
#' @param col.band error band colour.
#' @param opacity.band error band opacity.
#' @param lwd line width for fitted curve.
#' @param lty line type of fitted curve.
#' @param supports number of supporting points for plotting the function.
#' @param plot.range vector with two elements \code{c(min,max)} defining the
#' range in which fitline and errorband are plotted. Default is the range of
#' the data.
#' @param ... Additional parameters passed to the generic `plot` function.
#'
#' @export
#' @family NLS fit functions
plot.bootstrapfit <- function(x, ..., xlim, ylim, rep=FALSE, col.line="black", col.band="gray", opacity.band=0.65, lty=c(1), lwd=c(1), xlab="x", ylab="y", supports=1000, plot.range, error=sd) {
  if(missing(plot.range)){
    rx <- range(x$x)
  }else{
    rx <- plot.range
  }
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/supports)
  npar <- length(x$par.guess)

  ## use the xylimits computation of plotwitherror
  if(x$errormodel == "yerrors") {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, rep=rep, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  }
  else {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, rep=rep, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  }
  my.xlim <- mylims$xlim
  my.ylim <- mylims$ylim

  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  Y <- numeric()
  Y <- do.call(what=x$fn, args=c(list(par=x$t0[1:npar], x=X), x$tofn))

  ## error band
  ## define a dummy function to be used in apply
  dummyfn <- function(par, x, object) {
    return(do.call(what=object$fn, args=c(list(par=par, x=x), object$tofn)))
  }
  se <- apply(X=rbind(apply(X=x$t[, c(1:npar), drop=FALSE], MARGIN=1, FUN=dummyfn, x=X, object=x)), MARGIN=1, FUN=error)

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
