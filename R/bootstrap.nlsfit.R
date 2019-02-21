#' Parametric bootstrap
#'
#' @param boot.R numeric. Number of bootstrap samples to generate.
#' @param x numeric vector. Actual values for the data.
#' @param dx numeric vector of the same length as `x` or missing. Errors of the
#' values.
#' @param seed integer. Seed to use for the random number generation. If it is
#' missing, the seed will not be set to any particular value. If there was a
#' default value, all results would be exactly correlated. So if you want
#' reproducability by fixing the seeds, make sure you choose different seeds
#' for independent variables.
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
parametric.bootstrap <- function (boot.R, x, dx, seed) {
  stopifnot(length(x) == length(dx))

  old_seed <- swap_seed(seed)
  samples.list <- lapply(dx, function (error) rnorm(boot.R, 0, error))
  restore_seed(old_seed)
  samples <- do.call(cbind, samples.list)

  samples <- t(t(samples) + x)


  return (samples)
}

#' Parametric bootstrap with covariance
#'
#' @param cov numeric matrix, square, length of `x` or missing. Covariance
#'   between the various variables in the vector `x`.
#'
#' @inheritParams parametric.bootstrap
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
parametric.bootstrap.cov <- function (boot.R, x, cov, seed) {
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
  old_seed <- swap_seed(seed)
  samples.list <- lapply(errors, function (error) rnorm(boot.R, 0, error))
  restore_seed(old_seed)
  samples <- do.call(cbind, samples.list)

  samples <- samples %*% t(evectors)
  samples <- t(t(samples) + x)

  return (samples)
}

#' NLS fit with parametric bootstrap
#'
#' @inheritParams bootstrap.nlsfit
#' @inheritParams parametric.bootstrap
#'
#' @export
#' @family NLS fit functions
#' 
#' @examples
#' ## Declare some data.
#' value <- c(0.1, 0.2, 0.3)
#' dvalue <- c(0.01, 0.01, 0.015)
#' x <- c(1, 2, 3)
#' dx <- c(0.1, 0.1, 0.1)
#' boot.R <- 1500
#'
#' fn <- function (par, x, ...) par[1] + par[2] * x
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
#' @inheritParams bootstrap.nlsfit
#' @inheritParams parametric.bootstrap.cov
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

get.errors <- function (useCov, y, dy, dx, errormodel, bsamples, cov_fn, error) {
  ## invert covariance matrix, if applicable
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
      InvCovMatrix <- try(invertCovMatrix(bsamples, boot.l = 1, boot.samples = TRUE, cov_fn = cov_fn), silent = TRUE)
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

  if (errormodel == 'yerrors') {
    dx <- NULL
  }
  return(list(dY=dY, dy=dy, dx=dx))
}

set.fitchi <- function (fn, errormodel, useCov, dY, x, ipx) {
  ## define the chi-vector, the sum of squares of which has to be minimized
  ## the definitions depend on the errormodel and the use of covariance
  ## BUT it always has the same name
  if(errormodel == "yerrors"){
    if(useCov){
      fitchi <- function(y, par, ...) { dY %*% (y - fn(par=par, x=x, ...)) }
    }else{
      fitchi <- function(y, par, ...) { dY * (y - fn(par=par, x=x, ...)) }
    }
  }else{
    if(useCov){
      fitchi <- function(y, par, ...) { dY %*% (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx])) }
    }else{
      fitchi <- function(y, par, ...) { dY * (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx])) }
    }
  }

  return(fitchi)
}

set.dfitchi <- function (gr, dfn, errormodel, useCov, dY, x, ipx) {
  ## define the derivatives of chi and chi^2
  if(missing(gr) || (errormodel == "xyerrors" && missing(dfn))){
    ## in case no derivative is known, the functions are set to NULL
    ## this is the default in the optimization functions anyway
    dfitchi <- NULL
  }else{
    ## the format of gr has to be nrows=length(par), ncols=length(Y)
    if(errormodel == "yerrors"){
      if(useCov){
        dfitchi <- function(par, ...) { -dY %*% gr(par=par, x=x, ...) }
      }else{
        dfitchi <- function(par, ...) { -dY * gr(par=par, x=x, ...) }
      }
    }else{
      jacobian <- function(par, ...) {
        df.dpar <- rbind(gr(par=par[-ipx], x=par[ipx], ...), array(0,dim=c(nx,length(par.guess))))
        df.dx <- rbind(diag(dfn(par=par[-ipx], x=par[ipx], ...)), diag(1,nx))
        return(cbind(df.dpar, df.dx))
      }
      if(useCov){
        dfitchi <- function(par, ...) { -dY %*% jacobian(par, ...) }
      }else{
        dfitchi <- function(par, ...) { -dY * jacobian(par, ...) }
      }
    }
  }

  return(dfitchi)
}

set.dfitchisqr <- function (fitchi, dfitchi) {
  if(is.null(dfitchi)){
    dfitchisqr <- NULL
  }else{
    dfitchisqr <- function(y, par, ...) { 2 * crossprod(fitchi(y, par, ...), dfitchi(par, ...)) }
  }
  return(dfitchisqr)
}

set.wrapper <- function (fn, gr, dfn, errormodel, useCov, dY, x, ipx, lm.avail, maxiter) {
  fitchi <- set.fitchi(fn, errormodel, useCov, dY, x, ipx)
  dfitchi <- set.dfitchi(gr, dfn, errormodel, useCov, dY, x, ipx)
  ## define the wrapper-functions for optimization
  if (lm.avail) {
    control = minpack.lm::nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxfev=maxiter*10, maxiter=maxiter)
    wrapper <- function(y, par, ...) {
      suppressWarnings(
        res <- minpack.lm::nls.lm(
          par=par, fn=fitchi, y=y, jac=dfitchi,
          control = control,
          ...))

      list(converged = res$info %in% 1:3,
           info = res$info,
           par = res$par,
           chisq = res$rsstrace[length(res$rsstrace)])
    }
  } else {
    fitchisqr <- function(y, par) { sum(fitchi(y, par)^2) }
    dfitchisqr <- set.dfitchisqr(fitchi, dfitchi)
    wrapper <- function(y, par, ...) {
      res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("BFGS"), control=list(maxit=maxiter), ...)

      list(converged = res$convergence == 0,
           info = NA,
           par = res$par,
           chisq = res$value)
    }
  }

  return(wrapper)
}

#' Bootstrap a non-linear least-squares fit
#'
#' Performs and bootstraps a non-linear least-squares fit to data with y and x
#' errors.
#'
#' @param fn `fn(par, x, ...)`. The (non-linear) function to be fitted to the
#' data. Its first argument must be the fit parameters named \code{par}. The
#' second must be \code{x}, the explaining variable. Additional parameters
#' might be passed to the function. Currently we pass `boot_r` which is `0`
#' for the original data and the ID (1, ...) of the bootstrap sample otherwise.
#' As more parameters might be added in the future it is recommended that the
#' fit function accepts `...` as the last parameter to be forward compatible.
#' @param gr `gr(par, x, ...)`. \code{gr=d(fn) / d(par)} is a function to
#' return the gradient of \code{fn}. It must return an array with
#' \code{length(x)} rows and \code{length(par)} columns.
#' @param dfn `dfn(par, x, ...)`. \code{dfn=d(fn) / dx} is the canonical
#' derivative of \code{fn} by \code{x} and only relevant if x-errors are
#' provided.
#' @param par.guess initial guess values for the fit parameters.
#' @param y the data as a one-dimensional numerical vector to be described by
#' the fit function. 
#' @param x values of the explaining variable in form of a one-dimensional
#' numerical vector.
#' @param bsamples bootstrap samples of \code{y} (and \code{x}, if applicable).
#' Must be provided as array of dimensions \code{c(boot.R, n)} with \code{n}
#' equals to \code{length(y)} in case of 'yerrors' and For 'xyerrors' to
#' \code{length(y) + length(x)}.
#' @param ... Additional parameters passed to `fn`, `gr` and `dfn`.
#' @param dy,dx Numeric vector. Errors of the dependent and independent
#' variable, respectively. These do not need to be specified as they can be
#' computed from the bootstrap samples. In the case of parametric bootstrap it
#' might would lead to a loss of information if they were computed from the
#' pseudo-bootstrap samples. They must not be specified if a covariance matrix
#' is given.
#' @param CovMatrix complete variance-covariance matrix of dimensions
#' \code{c(length(y), length(y))} or \code{c(length(y)+length(x),
#' length(y)+length(x))} depending on the errormodel.
#' @param use.minpack.lm use the \code{minpack.lm} library if available. This
#' is usually faster than the default \code{optim} but somtimes also less
#' stable.
#' @param parallel parallelise over bootstrap samples. The package
#' \code{parallel} is required.
#' @param error Function that takes a sample vector and returns the error
#' estimate. This is a parameter in order to support different resampling
#' methods like jackknife.
#' @param maxiter integer. Maximum number of iterations that can be used in the
#' optimization process.
#' @param relative.weights are the errors on y (and x) to be interpreted as
#' relative weights instead of absolute ones? If TRUE, the covariance martix
#' of the fit parameter results is multiplied by chi^2/dof. This is the default
#' in many fit programs, e.g. gnuplot.
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
#' value <- c(0.1, 0.2, 0.31)
#' dvalue <- c(0.01, 0.01, 0.015)
#' x <- c(1, 2, 3)
#' dx <- c(0.1, 0.1, 0.1)
#' boot.R <- 1500
#'
#' fn <- function (par, x, boot_r, ...) par[1] + par[2] * x
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
#' plot(fit.result)
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
                             error = sd,
                             cov_fn = cov,
                             maxiter = 500,
                             relative.weights = FALSE) {
  stopifnot(!missing(y))
  stopifnot(!missing(x))
  stopifnot(!missing(par.guess))
  stopifnot(!missing(fn))
  stopifnot(!missing(bsamples))

  boot.R <- nrow(bsamples)
  useCov <- !missing(CovMatrix)

  if (use.minpack.lm) {
    lm.avail <- requireNamespace('minpack.lm')
  } else {
    lm.avail <- FALSE
  }

  if (parallel) {
    parallel <- requireNamespace('parallel')
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
  ipx <- length(par.Guess)-seq(nx-1,0)
  
  all.errors <- get.errors(useCov, y, dy, dx, errormodel, bsamples, cov_fn, error)
  dY <- all.errors$dY
  dy <- all.errors$dy
  dx <- all.errors$dx

  ## add original data as first row
  bsamples <- rbind(Y, bsamples)

  wrapper <- set.wrapper(fn, gr, dfn, errormodel, useCov, dY, x, ipx, lm.avail, maxiter)

  ## now the actual fit is performed
  first.res <- wrapper(Y, par.Guess, boot_r = 0, ...)
  if (!first.res$converged) {
    stop(sprintf('The first fit to the original data has failed. The `info` from the algorithm is `%d`', first.res$info))
  }

  if (parallel)
    my.lapply <- parallel::mclapply
  else {
    my.lapply <- lapply
  }

  boot.list <- my.lapply(crr, function(sample) { wrapper(y=bsamples[sample,], par=first.res$par, boot_r = sample - 1, ...) })

  par.boot <- do.call(rbind, lapply(boot.list, function (elem) elem$par))

  converged <- sapply(boot.list, function (elem) elem$converged)
  info <- sapply(boot.list, function (elem) elem$info)
  
  ## The fit on the original data must have converged, otherwise the results
  ## are worthless. We do not check directly after the first fit as the restart
  ## might have helped it to convergence, and we want to give the original data
  ## this second chance.
  if (!converged[1]) {
    stop(sprintf('The second fit to the original data has failed. The `info` from the algorithm is `%d`', info[1]))
  }
  
  if (any(!converged)) {
      warning('There were fits on the samples that did not converge. Check the `converged.boot` and `info.boot` fields of the return value for more information.')
  }
  
  ## We guarantee that the fit on the original data has converged, therefore we can discard this information.
  converged.boot <- converged[2:(boot.R + 1)]
  info.boot <- info[2:(boot.R + 1)]

  ## If most of the bootstrap samples have failed to converged, something else
  ## is clearly wrong. We take 50% as the cutoff.
  stopifnot(mean(converged.boot) >= 0.5)

  par.boot[!converged, ] <- NA

  chisq <- boot.list[[1]]$chisq
  dof = length(y) - length(par.guess)

  errors <- apply(par.boot[rr, , drop=FALSE], 2, error, na.rm = TRUE)
  if(relative.weights){
      errors <- errors * sqrt(chisq/dof)
  }

  res <- list(y=y, dy=dy, x=x, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R,
              bsamples=bsamples[rr, , drop=FALSE],
              errormodel=errormodel,
              converged.boot = converged.boot,
              t0=par.boot[1, ],
              t=par.boot[rr, , drop=FALSE],
              se=errors,
              useCov=useCov,
              invCovMatrix=dY,
              Qval = 1 - pchisq(chisq, dof),
              chisqr = chisq,
              dof = dof,
              error.function = error,
              info.boot = info.boot,
              relative.weights = relative.weights,
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
#' @param ... ignored
#' @param digits number of significant digits to print in summary or print.
#' @param print.correlation Logical. Whether to show the correlation between of
#' the fit parameters.
#'
#' @export
#' @family NLS fit functions
summary.bootstrapfit <- function(object, ..., digits = 2, print.correlation = TRUE) {
  cat("bootstrap nls fit\n\n")
  cat("model", object$errormodel, "\n")
  errors <- object$se
  values <- object$t0
  npar <- length(object$par.guess)
  
  ## parameters with errors as strings
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, with.dollar=FALSE, digits=digits, human.readable=FALSE)
  bias <- object$t0-apply(X=object$t, MARGIN=2, FUN=mean, na.rm=TRUE)
  dim(bias) <- c(length(bias), 1)
  bias <- apply(X=bias, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  ci16 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.16), drop=FALSE, na.rm=TRUE)
  dim(ci16) <- c(length(ci16), 1)
  ci16 <- apply(X=ci16, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  ci84 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.84), drop=FALSE, na.rm=TRUE)
  dim(ci84) <- c(length(ci84), 1)
  ci84 <- apply(X=ci84, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, human.readable=FALSE)
  cat("    best fit parameters with errors, bootstrap bias and 68% confidence interval\n\n")
  print(data.frame(par=tmp[1:npar], bias=bias[1:npar], ci16=ci16[1:npar], ci84=ci84[1:npar]))
  if(print.correlation){
    correlation <- cor(object$t, object$t, use="na.or.complete")
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
  cat('\nRatio of converged fits on samples:', sum(object$converged.boot), '/', length(object$converged.boot), '=', mean(object$converged.boot), '\n')
  if (!all(is.na(object$info.boot))) {
      cat('Table of nls.lm info values (1, 2, 3 are convergence):\n')
      df <- as.data.frame(table(object$info.boot))
      colnames(df) <- c('value', 'frequency')
      print(df)
  }
}

#' Print a bootstrap NLS fit
#'
#' @param x object returned by \code{bootstrap.nlsfit}
#' @param ... Additional parameters passed to the `summary.bootstrapfit` function.
#' @param digits number of significant digits to print in summary or print.
#'
#' @family NLS fit functions
print.bootstrapfit <- function(x, ..., digits = 2) {
  summary.bootstrapfit(object=x, digits=digits)
}

#' Plot a bootstrap NLS fit
#'
#' @param x object returned by \code{bootstrap.nlsfit}
#' @param col.line line colour.
#' @param col.band error band colour.
#' @param opacity.band error band opacity.
#' @param lwd line width for fitted curve.
#' @param lty line type of fitted curve.
#' @param supports number of supporting points for plotting the function.
#' @param plot.range vector with two elements \code{c(min,max)} defining the
#' range in which fitline and errorband are plotted. Default is the range of
#' the data.
#' @param ... Additional parameters passed to the `plotwitherror` function.
#'
#' @export
#' @family NLS fit functions
plot.bootstrapfit <- function(x, ..., col.line="black", col.band="gray", opacity.band=0.65, lty=c(1), lwd=c(1), supports=1000, plot.range, error=sd) {
  if(missing(plot.range)){
    rx <- range(x$x)
  }else{
    rx <- plot.range
  }
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/supports)
  npar <- length(x$par.guess)

  ## use the xylimits computation of plotwitherror
  if(x$errormodel == "yerrors") {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, ...)
  }
  else {
    mylims <- plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, ...)
  }
  my.xlim <- mylims$xlim
  my.ylim <- mylims$ylim

  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  Y <- do.call(x$fn, c(list(par = x$t0[1:npar], x = X, boot_r = 0), x$tofn))

  ## error band
  ## define a dummy function to be used in apply
  prediction_boot_fn <- function (boot_r) {
    par <- x$t[boot_r, 1:npar, drop = FALSE]
    do.call(x$fn, c(list(par = par, x = X, boot_r = boot_r), x$tofn))
  }
  predictions <- do.call(rbind, lapply(1:nrow(x$t), prediction_boot_fn))
  se <- apply(predictions, 2, error, na.rm = TRUE)
  stopifnot(length(se) == length(X))

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
