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
#' @param bootstrap Shall the error calculation be performed using boostrap?
#' If not, the errors are estimated with help of the jacobian (either provided
#' in \code{gr} or calculated using the \code{numDeriv}-package).
#'
#' @return
#' See \link{simple.nlsfit}.
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
parametric.nlsfit <- function (fn, par.guess, boot.R, y, dy, x, dx,
                               lower = rep(x = -Inf, times = length(par.guess)),
                               upper = rep(x = +Inf, times = length(par.guess)),
                               ..., bootstrap=TRUE) {
  stopifnot(length(x) == length(y))
  stopifnot(missing(dx) || length(dx) == length(x))
  stopifnot(missing(dy) || length(dy) == length(y))
  stopifnot(length(lower) == length(par.guess))
  stopifnot(length(upper) == length(par.guess))

  if (missing(dx)) {
    values <- y
    errors <- dy
    errormodel <- "yerrors"
  } else {
    values <- c(y, x)
    errors <- c(dy, dx)
    errormodel <- "xyerrors"
  }

  if (bootstrap) {
    stopifnot(!missing(boot.R))
    bsamples <- parametric.bootstrap(boot.R, values, errors)
    bootstrap.nlsfit(fn, par.guess, y, x, bsamples, ..., lower = lower, upper = upper, dx = dx, dy = dy)
  }else {
    if(missing(boot.R)) {
      boot.R = 0
    }
    simple.nlsfit(fn, par.guess, y, x, errormodel, ..., lower = lower, upper = upper, dx = dx, dy = dy, boot.R = boot.R)
  }
}

#' parametric.nlsfit.cov
#' 
#' @description
#' NLS fit with parametric bootstrap and covariance
#'
#' @inheritParams bootstrap.nlsfit
#' @inheritParams parametric.bootstrap.cov
#' @param bootstrap boolean. If `TRUE`, bootstrap is used.
#' 
#' @return
#' See \link{simple.nlsfit}.
#' 
#' @export
#' @family NLS fit functions
parametric.nlsfit.cov <- function (fn, par.guess, boot.R, y, x, cov,
                                   lower = rep(x = -Inf, times = length(par.guess)),
                                   upper = rep(x = +Inf, times = length(par.guess)),
                                   ..., bootstrap=TRUE, na.rm = FALSE) {
  stopifnot(length(x) == length(y))
  stopifnot(length(lower) == length(par.guess))
  stopifnot(length(upper) == length(par.guess))

  if (ncol(cov) == length(y)) {
    values <- y
    errormodel <- "yerrors"
  } else if (ncol(cov) == length(y) + length(x)) {
    values <- c(y, x)
    errormodel <- "xyerrors"
  } else {
    stop('The covariance matrix must either be as large as `y` or as `y` and `x` together.')
  }

  if (bootstrap) {
    stopifnot(!missing(boot.R))
    bsamples <- parametric.bootstrap.cov(boot.R, values, cov)
    bootstrap.nlsfit(fn, par.guess, y, x, bsamples, ..., CovMatrix = cov, na.rm = na.rm)
  }else {
    if(missing(boot.R)) {
      boot.R = 0
    }
    simple.nlsfit(fn, par.guess, y, x, errormodel, ..., CovMatrix = cov, boot.R = boot.R, na.rm = na.rm)
  }
}

get.errors <- function (useCov, y, dy, dx, CovMatrix, errormodel, bsamples, cov_fn, error) {
  ## invert covariance matrix, if applicable
  if (useCov) {
    inversion.worked <- function(InvCovMatrix) {
      if (inherits(InvCovMatrix, "try-error")) {
        stop("Variance-covariance matrix could not be inverted!")
      }
    }

    if (is.null(CovMatrix)) {
      # no (custom) covariance matrix was passed
      # we want to (potentially) rely on the SV decomposition with replacement of small
      # eigenvalues implemented in `invertCovMatrix`
      CovMatrix <- cov_fn(bsamples)
      InvCovMatrix <- try(invertCovMatrix(bsamples, boot.l = 1, boot.samples = TRUE, cov_fn = cov_fn), silent = TRUE)
      inversion.worked(InvCovMatrix)
      W <- chol(InvCovMatrix)
    } else {
      # a (potentially hand-crafted) covariance matrix was passed
      # we assume that it's cleanly invertible at this stage and simply use `solve`
      CholCovMatrix <- chol(CovMatrix)
      InvCovMatrix <- try(solve(CholCovMatrix), silent = TRUE)
      inversion.worked(InvCovMatrix)
      W <- t(InvCovMatrix)
    }

    dydx <- sqrt(diag(CovMatrix))

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
      W <- 1.0 / dydx

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
        W <- 1.0 / dy
      } else if (errormodel == 'xyerrors' && ncol(bsamples) == length(dy) + length(dx)) {
        W <- 1.0 / c(dy, dx)
      } else {
        stop('You have explicitly passed `dy` and/or `dx`, but their combined length does not match the number of columns of the bootstrap samples.')
      }
    }
  }

  if (errormodel == 'yerrors') {
    dx <- NULL
  }
  return(list(W=W, dy=dy, dx=dx))
}

get.errors.wo.bootstrap <- function (useCov, y, dy, dx, CovMatrix, errormodel) {
  ## invert covariance matrix, if applicable
  if (useCov) {
    inversion.worked <- function(InvCovMatrix) {
      if (inherits(InvCovMatrix, "try-error")) {
        stop("Variance-covariance matrix could not be inverted!")
      }
    }

    if (missing(CovMatrix)) {
        stop("If you want to use covariance, you have to provide the matrix.")
    } else {
      CholCovMatrix <- chol(CovMatrix)
      InvCovMatrix <- try(solve(CholCovMatrix), silent = TRUE)
      inversion.worked(InvCovMatrix)
      W <- t(InvCovMatrix)
    }

    dydx <- 1.0 / diag(W)

    if (errormodel == 'yerrors') {
      dy <- dydx
    } else {
      dy <- dydx[1:length(y)]
      dx <- dydx[(length(y)+1):length(dydx)]
    }
  }
  else {
    if (errormodel == 'yerrors') {
      W <- 1.0 / dy
    } else {
      W <- 1.0 / c(dy, dx)
    }
  }

  if (errormodel == 'yerrors') {
    dx <- NULL
  }
  return(list(W=W, dy=dy, dx=dx))
}

set.fitchi <- function (fn, errormodel, useCov, W, x, ipx, na.rm, priors) {
  ## define the chi-vector, the sum of squares of which has to be minimized
  ## the definitions depend on the errormodel and the use of covariance
  ## BUT it always has the same name
  if(errormodel == "yerrors"){
    if(useCov){
      fitchi <- function(y, par, ...) { W %*% (y - c(fn(par=par, x=x, ...), par[priors$param])) }
    }else{
      fitchi <- function(y, par, ...) { W * (y - c(fn(par=par, x=x, ...), par[priors$param])) }
    }
  }else{
    if(useCov){
      fitchi <- function(y, par, ...) { W %*% (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx], par[priors$param])) }
    }else{
      fitchi <- function(y, par, ...) { W * (y - c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx], par[priors$param])) }
    }
  }

  if(na.rm){
    if(useCov){
      W.na = apply(is.na(W), 1, any)
      fitchi.wo.na <- function(y, par, ...) { ifelse(is.na(y) | W.na, 0, fitchi(y, par, ...)) }
    }else{
      fitchi.wo.na <- function(y, par, ...) { ifelse(is.na(y) | is.na(W), 0, fitchi(y, par, ...)) }
    }
    return(fitchi.wo.na)
  }else{
    return(fitchi)
  }
}

set.dfitchi <- function (gr, dfn, par.guess, errormodel, useCov, W, x, ipx, na.rm, priors, priors.avail) {
  ## define the derivatives of chi and chi^2
  if(missing(gr) || (errormodel == "xyerrors" && missing(dfn))){
    ## in case no derivative is known, the functions are set to NULL
    ## this is the default in the optimization functions anyway
    return(NULL)
  }else{
    ## the format of gr has to be nrows=length(par), ncols=length(Y)
    if(errormodel == "yerrors"){
      grpriors <- c()
      if(priors.avail) {
        npriors <- length(priors$param)
        npar <- length(par.guess)
        for (i in 1:npriors) {
          aux <- rep(0, npar)
          aux[priors$param[i]] <- 1
          grpriors <- rbind(grpriors, aux)
        }
      }
      if(useCov){
        dfitchi <- function(par, ...) { -W %*% rbind(gr(par=par, x=x, ...), grpriors) }
      }else{
        dfitchi <- function(par, ...) { -W * rbind(gr(par=par, x=x, ...), grpriors) }
      }
    }else{
      nx <- length(x)
      jacobian <- function(par, ...) {
        df.dpar <- rbind(gr(par=par[-ipx], x=par[ipx], ...), array(0,dim=c(nx,length(par))))
        df.dx <- rbind(diag(dfn(par=par[-ipx], x=par[ipx], ...)), diag(1,nx))
        return(cbind(df.dpar, df.dx))
      }
      if(useCov){
        dfitchi <- function(par, ...) { -W %*% jacobian(par, ...) }
      }else{
        dfitchi <- function(par, ...) { -W * jacobian(par, ...) }
      }
    }

    if(na.rm){
      if(useCov){
        W.na = apply(is.na(W), 1, any)
        dfitchi.wo.na <- function(y, par, ...) { ifelse(is.na(y) | W.na, 0, dfitchi(y, par, ...)) }
      }else{
        dfitchi.wo.na <- function(y, par, ...) { ifelse(is.na(W), 0, dfitchi(y, par, ...)) }
      }
      return(dfitchi.wo.na)
    }else{
      return(dfitchi)
    }
  }
}

set.dfitchisqr <- function (fitchi, dfitchi) {
  if(is.null(dfitchi)){
    dfitchisqr <- NULL
  }else{
    dfitchisqr <- function(y, par, ...) { 2 * crossprod(fitchi(y, par, ...), dfitchi(par, ...)) }
  }
  return(dfitchisqr)
}

set.wrapper <- function (fn, gr, dfn, par.guess, errormodel, useCov, W, x, ipx, lm.avail, maxiter, success.infos, na.rm, priors, priors.avail, lower, upper) {
  fitchi <- set.fitchi(fn, errormodel, useCov, W, x, ipx, na.rm, priors)
  dfitchi <- set.dfitchi(gr, dfn, par.guess, errormodel, useCov, W, x, ipx, na.rm, priors, priors.avail)
  ## define the wrapper-functions for optimization
  if (lm.avail) {
    control = minpack.lm::nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxfev=maxiter*10, maxiter=maxiter)
    wrapper <- function(y, par, ...) {
      suppressWarnings(
        res <- minpack.lm::nls.lm(
          par=par, fn=fitchi, y=y, jac=dfitchi,
          control = control,
          lower = lower, upper = upper,
          ...))

      list(converged = res$info %in% success.infos,
           info = res$info,
           par = res$par,
           chisq = res$rsstrace[length(res$rsstrace)],
           niter = res$niter)
    }
  } else {
    fitchisqr <- function(y, par, ...) { sum(fitchi(y, par, ...)^2) }
    dfitchisqr <- set.dfitchisqr(fitchi, dfitchi)
      wrapper <- function(y, par, ...) {
        if( any(upper != +Inf) | any(lower != -Inf) ){
          res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("L-BFGS-B"), 
                       lower = lower, upper = upper,
                       control=list(maxit=maxiter), ...)
        } else {
          res <- optim(par=par, fn=fitchisqr, gr=dfitchisqr, y=y, method=c("BFGS"), 
                       control=list(maxit=maxiter), ...)
        }
      
      list(converged = res$convergence == 0,
           info = NA,
           par = res$par,
           chisq = res$value,
           niter = res$counts[1])
    }
  }

  return(wrapper)
}

#' NLS fit with without bootstrap
#'
#' @inheritParams bootstrap.nlsfit
#' @param errormodel Either "yerror" or "xyerror", depending on the x-values having
#' errors or not.
#' @param boot.R If larger than 0, \code{boot.R} paramtetric bootstrap samples are
#' generated on the fit results after fit and error calculation are finished.
#' The original data is never boostraped in this function.
#' @param priors List possessing the elements `param`, `p` and `psamples`.
#' The vector `param` includes the indices of all fit parameters that are
#' to be constrained and the vector `p` the corresponding paramater values
#' (e.g. known from a previous fit). The list element `psamples` is a matrix of
#' dimensions \code{(boot.R, length(param))} and contains the corresponding
#' bootstrap samples. If this list is not specified priors are omitted
#' within the fit.
#'
#' @return
#' Returns an object of class `bootstrapfit`, see \link{bootstrap.nlsfit}.
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
#'
#' fn <- function (par, x, ...) par[1] + par[2] * x
#'
#' fit.result <- simple.nlsfit(fn, c(1, 1), value, x, "xyerrors", dy=dvalue, dx=dx)
#' summary(fit.result)
simple.nlsfit <- function(fn,
                          par.guess,
                          y,
                          x,
                          errormodel,
                          priors = list(param = c(), p = c(), psamples = c()),
                          ...,
                          lower = rep(x = -Inf, times = length(par.guess)),
                          upper = rep(x = +Inf, times = length(par.guess)),
                          dy,
                          dx,
                          CovMatrix,
                          boot.R = 0,
                          gr,
                          dfn,
                          mask,
                          use.minpack.lm = TRUE,
                          error = sd,
                          maxiter = 500,
                          success.infos = 1:3,
                          relative.weights = FALSE,
                          na.rm = FALSE) {
  if(!is.null(priors$psamples)) {
    ncolps <- ncol(as.matrix(priors$psamples))
  } else {
    ncolps <- length(priors$psamples)
  }
  stopifnot(!missing(y))
  stopifnot(!missing(dy))
  stopifnot(!missing(x))
  stopifnot(!missing(par.guess))
  stopifnot(!missing(fn))
  if(!is.null(priors$param)){
    stopifnot(is.vector(priors$param))
  }
  stopifnot( length(priors$param) == length(priors$p) &&
               length(priors$param) == ncolps &&
               length(priors$p) == ncolps )
  if( !is.null(c(priors$param, priors$p, priors$psamples)) ){
    stop("Priors are not implemented in simple.nlsfit yet.")
  }
  stopifnot(length(lower) == length(par.guess))
  stopifnot(length(upper) == length(par.guess))

  # The user might have specified a mask which is logical. Since we later use
  # it as integer indices, we need to make sure that it is an integer mask.
  if (!missing(mask) && is.logical(mask)) {
    mask <- which(mask)
  }

  useCov <- !missing(CovMatrix)
  if (useCov && !(missing(dx) && missing(dy))) {
    stop('Specifying a covariance matrix and `dx` and `dy` does not make sense, use either.')
  }

  if (use.minpack.lm) {
    lm.avail <- requireNamespace('minpack.lm')
  } else {
    lm.avail <- FALSE
  }

  ## Apply the mask as in bootstrap.nlsfit.
  if (!missing(mask)) {
    if(errormodel == "xyerrors") {
      ## CovMatrix has twice the length for xy-errors, therefore
      ## need to apply mask twice.
      mask.twice <- c(mask, length(y)+mask)
    }else{
      dx <- NULL
    }

    full <- list(x=x, y=y,
                 dx=dx, dy=dy)

    x <- x[mask]
    y <- y[mask]
    if(errormodel == "xyerrors") dx <- dx[mask]
    dy <- dy[mask]

    if (!missing(CovMatrix)) {
      full$CovMatrix <- CovMatrix
      if(errormodel == "yerrors"){
        CovMatrix <- CovMatrix[mask, mask]
      }else{
        CovMatrix <- CovMatrix[mask.twice, mask.twice]
      }
    }
  }

  all.errors <- get.errors.wo.bootstrap(useCov, y, dy, dx, CovMatrix, errormodel)

  ## cast y and dy to Y and dY, respectively
  if (errormodel == 'yerrors') {
    Y <- y
    par.Guess <- par.guess
  } else {
    Y <- c(y, x)
    par.Guess <- c(par.guess, x)
  }

  nx <- length(x)
  ipx <- length(par.Guess)-seq(nx-1,0)
  
  ## If a list 'priors' is specified, modify the parameters y and func
  ## by adding p, param and psamples, respectively.
  priors.avail <- !any(c(is.null(priors$param), is.null(priors$p), is.null(priors$psamples)))
  if(priors.avail) {
    Yp <- c(Y, priors$p)
    yp <- c(y, priors$p)
  }
  
  if(!priors.avail) {
    dy <- all.errors$dy
    dx <- all.errors$dx
  } else {
    dy <- head(all.errors$dy, -length(priors$p))
    dydp <- all.errors$dy
    dx <- all.errors$dx
  }
  W <- all.errors$W
  
  if (length(upper) < length(par.Guess)) {
    upper <- c(upper, rep(+Inf, times = length(par.Guess) - length(upper)))
  }
  if (length(lower) < length(par.Guess)) {
    lower <- c(lower, rep(-Inf, times = length(par.Guess) - length(lower)))
  }

  wrapper <- set.wrapper(fn, gr, dfn, par.guess, errormodel, useCov, W, x, ipx, lm.avail, maxiter, success.infos, na.rm, priors, priors.avail, lower, upper)

  ## now the actual fit is performed
  if(!priors.avail){
    first.res <- wrapper(Y, par.Guess, ...)
  } else{
    first.res <- wrapper(Yp, par.Guess, ...)
  }
  if (!first.res$converged) {
    stop(sprintf('The fit has failed. The `info` from the algorithm is `%d`', first.res$info))
  }

  chisq <- first.res$chisq
  if(!priors.avail) {
    dof = length(y) - length(par.guess)
  } else {
    dof = length(yp) - length(par.guess)
  }

  if (missing(gr) || (errormodel == "xyerrors" && missing(dfn))) {
    if (!requireNamespace("numDeriv")) {
      stop("Errors on the fit results cannot be computed. Either provide the jacobian yourself, or install the package numDeriv.")
    }
    if(errormodel == "yerrors"){
      ## We have to hide the argument x from the function jacobian
      hidden <- function(par, my.x, ...) { fn(par=par, x=my.x, ...) }
      jac <- numDeriv::jacobian(hidden, first.res$par, my.x=x, ...)
    }else{
      composition <- function(par, ...) { c(fn(par=par[-ipx], x=par[ipx], ...), par[ipx]) }
      jac <- numDeriv::jacobian(composition, first.res$par, ...)
    }
  }else {
    if(errormodel == "yerrors"){
      jac <- gr(par=first.res$par, x=x, ...)
    }else{
      jacobian <- function(par, ...) {
        df.dpar <- rbind(gr(par=par[-ipx], x=par[ipx], ...), array(0,dim=c(nx,length(par.guess))))
        df.dx <- rbind(diag(dfn(par=par[-ipx], x=par[ipx], ...)), diag(1,nx))
        return(cbind(df.dpar, df.dx))
      }
      jac <- jacobian(par, ...)
    }
  }

  ## Normalise the errors
  if(useCov) {
    jac <- W %*% jac
  }else{
    jac <- diag(W) %*% jac
  }

  cov <- solve(t(jac) %*% jac)
  cov <- 0.5*(cov + t(cov))
  if(relative.weights){
      cov <- cov * chisq/dof
  }
  errors <- sqrt(diag(cov))
  
  if(!priors.avail) {
    dydp = NA
  } else {
    dydp=dydp
  }

  res <- list(y=y, dy=dy, dydp=dydp, x=x, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R,
              errormodel=errormodel,
              t0=first.res$par,
              se=errors,
              cov=cov,
              useCov=useCov,
              invCovMatrix=W,
              Qval = 1 - pchisq(chisq, dof),
              chisqr = chisq,
              dof = dof,
              error.function = error,
              relative.weights = relative.weights,
              tofn=list(...),
              lower=lower,
              upper=upper)

  if (errormodel == 'xyerrors') {
    res$dx <- dx
  }

  if(boot.R > 0){
    res$t <- parametric.bootstrap.cov(boot.R, first.res$par, cov)
  }

  # The user might have supplied a mask, therefore we need to restore all the
  # information and un-apply the mask.
  if (!missing(mask)) {
    for (name in names(full)) {
      if (name %in% names(res)) {
        res[[name]] <- full[[name]]
      }
    }
    res$mask <- mask
  }

  attr(res, "class") <- c("bootstrapfit", "list")
  return(invisible(res))
}

#' Bootstrap a non-linear least-squares fit
#'
#' Performs and bootstraps a non-linear least-squares fit to data with y and x
#' errors.
#'
#' @param fn `fn(par, x, ...)`. The (non-linear) function to be fitted to the
#' data. Its first argument must be the fit parameters named \code{par}. The
#' second must be \code{x}, the explaining variable. Additional parameters
#' might be passed to the function. Currently we pass `boot.r` which is `0`
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
#' @param priors List possessing the elements `param`, `p` and `psamples`.
#' The vector `param` includes the indices of all fit parameters that are
#' to be constrained and the vector `p` the corresponding paramater values
#' (e.g. known from a previous fit). The list element `psamples` is a matrix of
#' dimensions \code{(boot.R, length(param))} and contains the corresponding
#' bootstrap samples. If this list is not specified priors are omitted
#' within the fit.
#' @param ... Additional parameters passed to `fn`, `gr` and `dfn`.
#' @param lower Numeric vector of length \code{length(par.guess)}
#' of lower bounds on the fit parameters. If missing, \code{-Inf}
#' will be set for all.
#' @param upper Numeric vector of length \code{length(par.guess)}
#' of upper bounds on the fit parameters. If missing, \code{+Inf}
#' will be set for all.
#' @param dy,dx Numeric vector. Errors of the dependent and independent
#' variable, respectively. These do not need to be specified as they can be
#' computed from the bootstrap samples. In the case of parametric bootstrap it
#' might would lead to a loss of information if they were computed from the
#' pseudo-bootstrap samples. They must not be specified if a covariance matrix
#' is given.
#' @param CovMatrix complete variance-covariance matrix of dimensions
#' \code{c(length(y), length(y))} or \code{c(length(y)+length(x),
#' length(y)+length(x))} depending on the errormodel. Pass `NULL` if the matrix
#' has to be calculated from the `bsamples`. In that case, if the number of
#' boostrap samples is small compared to the number of variables, singular value
#' decomposition with small eigenvalue replacement will be used (see \link{invertCovMatrix})
#' to attempt a clean inversion.
#' In case a variance-covariance matrix is passed, the inversion will simply be attempted
#' using \code{solve} on the Cholesky decomposition.
#' Finally, if `CovMatrix` is missing, an uncorrelated fit will be performed.
#' @param mask logical or integer index vector. The mask is applied to select the observations from the data that are to be used in the fit. It is applied to `x`, `y`, `dx`, `dy`, `bsamples` and `CovMatrix` as applicable.
#' @param use.minpack.lm use the \code{minpack.lm} library if available. This
#' is usually faster than the default \code{optim} but somtimes also less
#' stable.
#' @param parallel parallelise over bootstrap samples. The package
#' \code{parallel} is required.
#' @param error Function that takes a sample vector and returns the error
#' estimate. This is a parameter in order to support different resampling
#' methods like jackknife.
#' @param cov_fn function. Function to compute the covariance
#'   (matrix). Default is \link{cov}.  
#' @param maxiter integer. Maximum number of iterations that can be used in the
#' optimization process.
#' @param success.infos integer vector. When using `minpack.lm` there is the
#' `info` in the return value. Values of 1, 2 or 3 are certain success. A value
#' of 4 could either be a success or a saddle point. If you want to interpret
#' this as a success as well just pass `1:4` instead of the default `1:3`.
#' @param relative.weights are the errors on y (and x) to be interpreted as
#' relative weights instead of absolute ones? If TRUE, the covariance martix
#' of the fit parameter results is multiplied by chi^2/dof. This is the default
#' in many fit programs, e.g. gnuplot.
#' @param na.rm logical. If set to `true`, NAs in `y` and `dy` will be ignored.
#' If x-errors are taken into account, NAs in `x` and `dx` will be ignored, too.
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
#'    fit function}
#'  \item{mask}{original `mask` value}
#'
#' @examples
#' ## Declare some data.
#' value <- c(0.1, 0.2, 0.31)
#' dvalue <- c(0.01, 0.01, 0.015)
#' x <- c(1, 2, 3)
#' dx <- c(0.1, 0.1, 0.1)
#' boot.R <- 1500
#'
#' fn <- function (par, x, boot.r, ...) par[1] + par[2] * x
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
#' plot(fit.result, main = 'Ribbon on top')
#' plot(fit.result, ribbon.on.top = FALSE, main = 'Ribbon below')
#' residual_plot(fit.result, main = 'Residual Plot')
#'
#' @export
#' @family NLS fit functions
bootstrap.nlsfit <- function(fn,
                             par.guess,
                             y,
                             x,
                             bsamples,
                             priors = list(param = c(), p = c(), psamples = c()),
                             ...,
                             lower = rep(x = -Inf, times = length(par.guess)),
                             upper = rep(x = +Inf, times = length(par.guess)),
                             dy,
                             dx,
                             CovMatrix,
                             gr,
                             dfn,
                             mask,
                             use.minpack.lm = TRUE,
                             parallel = FALSE,
                             error = sd,
                             cov_fn = cov,
                             maxiter = 500,
                             success.infos = 1:3,
                             relative.weights = FALSE,
                             na.rm = FALSE) {
  if(!is.null(priors$psamples)) {
    ncolps <- ncol(as.matrix(priors$psamples))
  } else {
    ncolps <- length(priors$psamples)
  }
  stopifnot(!missing(y))
  stopifnot(!missing(x))
  stopifnot(!missing(par.guess))
  stopifnot(!missing(fn))
  stopifnot(!missing(bsamples))
  if(!is.null(priors$param)){
    stopifnot(is.vector(priors$param))
  }
  stopifnot(length(priors$param) == length(priors$p))
  stopifnot(length(priors$param) == ncolps)
  stopifnot(length(priors$p) == ncolps)
  stopifnot(length(lower) == length(par.guess))
  stopifnot(length(upper) == length(par.guess))

  # The user might have specified a mask which is logical. Since we later use
  # it as integer indices, we need to make sure that it is an integer mask.
  if (!missing(mask) && is.logical(mask)) {
    mask <- which(mask)
  }

  boot.R <- nrow(bsamples)
  useCov <- !missing(CovMatrix)
  if (useCov && !(missing(dx) && missing(dy))) {
    stop('Specifying a covariance matrix and `dx` and `dy` does not make sense, use either.')
  }
  
  if (use.minpack.lm) {
    lm.avail <- requireNamespace('minpack.lm')
  } else {
    lm.avail <- FALSE
  }

  if (parallel) {
    parallel <- requireNamespace('parallel')
  }

  ## determine the errormodel
  if (ncol(bsamples) == length(y)) {
    errormodel <- "yerrors"
  } else if (ncol(bsamples) == length(y) + length(x)) {
    errormodel <- "xyerrors"
  } else {
    stop("The provided bootstrap samples do not match the number of data points with errors. Make sure that the number of columns is either the length of `y` alone for just y-errors or the length of `y` and `x` for xy-errors.")
  }

  ## Apply the mask. The user might have specified a mask that is used to
  ## restrict the selection of the points that are to be used in the fit. In
  ## order to make this additional feature a minimal change to the following code
  ## we will *change* the input parameters here and store them with new names.
  ## Then at the very end we switch them back.
  if (!missing(mask)) {
    if(errormodel == "xyerrors") {
      if(missing(dx)) {
        dx <- apply(bsamples[, (length(y)+1):ncol(bsamples)], 2, error)
      }
      ## bsamples and CovMatrix have twice the length for xy-errors, therefore
      ## need to apply mask twice.
      mask.twice <- c(mask, length(y)+mask)
    }else{
      dx <- NULL
    }

    if(missing(dy)) {
      dy <- apply(bsamples[, 1:length(y)], 2, error)
    }

    full <- list(x=x, y=y,
                 dx=dx, dy=dy,
                 bsamples=bsamples)

    x <- x[mask]
    y <- y[mask]
    dy <- dy[mask]
    if(errormodel == "yerrors"){
      bsamples <- bsamples[, mask]
    }else{
      bsamples <- bsamples[, mask.twice]
      dx <- dx[mask]
    }
    
    if (!missing(CovMatrix)) {
      full$CovMatrix <- CovMatrix
      if(errormodel == "yerrors"){
        CovMatrix <- CovMatrix[mask, mask]
      }else{
        CovMatrix <- CovMatrix[mask.twice, mask.twice]
      }
    }
  }

  crr <- c(1:(boot.R+1))
  rr <- c(2:(boot.R+1))

  if(errormodel == "yerrors"){
    Y <- y
    par.Guess <- par.guess
  }else{
    Y <- c(y, x)
    par.Guess <- c(par.guess, x)
  }


  nx <- length(x)
  ipx <- length(par.Guess)-seq(nx-1,0)
  
  ## If a list 'priors' is specified, modify the parameters y, func and bsamples
  ## by adding p, param and psamples, respectively.
  priors.avail <- !any(c(is.null(priors$param), is.null(priors$p), is.null(priors$psamples)))
  if(priors.avail) {
    Yp <- c(Y, priors$p)
    yp <- c(y, priors$p)
    if(errormodel == "yerrors"){
      dy <- c(dy, apply(matrix(priors$psamples, boot.R), 2, error))
    }else{
      dx <- c(dx, apply(matrix(priors$psamples, boot.R), 2, error))
    }
    bsamples <- cbind(bsamples, priors$psamples)
    CovMatrix <- cov_fn(bsamples)
  }
  
  all.errors <- get.errors(useCov, y, dy, dx, CovMatrix, errormodel, bsamples, cov_fn, error)

  if(!priors.avail) {
    dy <- all.errors$dy
    dx <- all.errors$dx
  } else {
    if(errormodel == "yerrors"){
      dydp <- all.errors$dy
      dy <- head(dydp, -length(priors$p))
      dx <- all.errors$dx
    }else{
      dy <- all.errors$dy
      dydp <- c(dy, tail(all.errors$dx, length(priors$p)))
      dx <- head(all.errors$dx, -length(priors$p))
    }
  }
  W <- all.errors$W
  
  ## add original data as first row
  if(!priors.avail){
    bsamples <- rbind(Y, bsamples)
  } else {
    bsamples <- rbind(Yp, bsamples)
  }
  
  if (length(upper) < length(par.Guess)) {
    upper <- c(upper, rep(+Inf, times = length(par.Guess) - length(upper)))
  }
  if (length(lower) < length(par.Guess)) {
    lower <- c(lower, rep(-Inf, times = length(par.Guess) - length(lower)))
  }
  
  wrapper <- set.wrapper(fn, gr, dfn, par.guess, errormodel, useCov, W, x, ipx, lm.avail, maxiter, success.infos, na.rm, priors, priors.avail, lower, upper)
  
  ## now the actual fit is performed
  if(!priors.avail){
    first.res <- wrapper(Y, par.Guess, boot.r = 0, ...)
  } else{
    first.res <- wrapper(Yp, par.Guess, boot.r = 0, ...)
  }
  if (!first.res$converged) {
    stop(sprintf('The first fit to the original data has failed. The `info` from the algorithm is `%d`', first.res$info))
  }

  if (parallel)
    my.lapply <- parallel::mclapply
  else {
    my.lapply <- lapply
  }

  boot.list <- my.lapply(crr, function(sample) { wrapper(y=bsamples[sample,], par=first.res$par, boot.r = sample - 1, ...) })

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

  niter_valboot <- sapply(boot.list, function (elem) elem$niter)

  ## If most of the bootstrap samples have failed to converged, something else
  ## is clearly wrong. We take 50% as the cutoff.
  stopifnot(mean(converged.boot) >= 0.5)

  par.boot[!converged, ] <- NA

  chisq <- boot.list[[1]]$chisq
  if(!priors.avail) {
    dof = length(y) - length(par.guess)
  } else {
    dof = length(yp) - length(par.guess)
  }

  errors <- apply(par.boot[rr, , drop=FALSE], 2, error, na.rm = TRUE)
  if(relative.weights){
      errors <- errors * sqrt(chisq/dof)
  }
  
  if(!priors.avail) {
    dydp = NA
  } else {
    dydp=dydp
  }
  
  res <- list(y=y, dy=dy, dydp=dydp, x=x, nx=nx,
              fn=fn, par.guess=par.guess, boot.R=boot.R,
              bsamples=bsamples[rr, , drop=FALSE],
              errormodel=errormodel,
              converged.boot = converged.boot,
              t0=par.boot[1, ],
              t=par.boot[rr, , drop=FALSE],
              se=errors,
              useCov=useCov,
              invCovMatrix=W,
              Qval = 1 - pchisq(chisq, dof),
              chisqr = chisq,
              dof = dof,
              error.function = error,
              info.boot = info.boot,
              relative.weights = relative.weights,
              tofn=list(...),
              niter = niter_valboot,
              lower=lower,
              upper=upper)

  if (errormodel == 'xyerrors') {
    res$dx <- dx
  }
  
  # The user might have supplied a mask, therefore we need to restore all the
  # information and un-apply the mask.
  if (!missing(mask)) {
    for (name in names(full)) {
      if (name %in% names(res)) {
        res[[name]] <- full[[name]]
      }
    }
    res$mask <- mask
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
#' @return
#' No return value.
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
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
  if(!is.null(object$t)) {
    bias <- object$t0-apply(X=object$t, MARGIN=2, FUN=mean, na.rm=TRUE)
    dim(bias) <- c(length(bias), 1)
    bias <- apply(X=bias, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
    ci16 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.16), drop=FALSE, na.rm=TRUE)
    dim(ci16) <- c(length(ci16), 1)
    ci16 <- apply(X=ci16, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
    ci84 <- apply(X=object$t, MARGIN=2, FUN=quantile, probs=c(0.84), drop=FALSE, na.rm=TRUE)
    dim(ci84) <- c(length(ci84), 1)
    ci84 <- apply(X=ci84, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
    cat("    best fit parameters with errors, bootstrap bias and 68% confidence interval\n\n")
    print(data.frame(par=tmp[1:npar], bias=bias[1:npar], ci16=ci16[1:npar], ci84=ci84[1:npar]))
  }else{
    cat("    best fit parameters with errors\n\n")
    print(data.frame(par=tmp[1:npar]))
  }
  if(print.correlation){
    if(!is.null(object$cov)) {
      cov.to.cor <- diag(1 / errors)
      correlation <- cov.to.cor %*% object$cov %*% cov.to.cor
    }else if(!is.null(object$t)) {
      correlation <- cor(object$t, object$t, use="na.or.complete")
    }else{
      stop("Correlation cannot be computed.")
    }
    cat("\n   correlation matrix of the fit parameters\n\n")
    print(data.frame(correlation))
  }
  if( any(object$upper != +Inf) ){
    cat("Upper bounds on parameter values:\n")
    print(object$upper)
  }
  if( any(object$lower != -Inf) ){
    cat("Lower bounds on parameter values:\n")
    print(object$lower)
  }
  if(!is.null(object$t) && object$errormodel != "yerrors") {
    cat("\n estimates for x-values with errors, bootstrap bias and 68% confidence interval\n\n")
    ii <- c((npar+1):length(tmp))
    print(data.frame(x=tmp[ii], bias=bias[ii], ci16=ci16[ii], ci84=ci84[ii]))
  }
  cat("\n   chi^2 and fit quality\n")
  cat("chisqr / dof =", object$chisqr, "/", object$dof, "=", object$chisqr/object$dof, "\n")
  cat("p-value", object$Qval, "\n")
  if(!is.null(object$converged.boot)){
    cat('\nRatio of converged fits on samples:', sum(object$converged.boot), '/', length(object$converged.boot), '=', mean(object$converged.boot), '\n')
  }
  if (!is.null(object$info.boot) && !all(is.na(object$info.boot))) {
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
#' @return
#' No return value.
#' 
#' @family NLS fit functions
#'
#' @export
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
#'   range in which fitline and errorband are plotted. Default is the range of
#'   the data.
#' @param ... Additional parameters passed to the `plotwitherror` function.
#' @param error Function to compute the standard error in resampling schemes.
#'   Default is \link{sd} for bootstrap. For other resampling schemes this might
#'   need to be changed.
#' @param ribbon.on.top Logical, controls whether the ribbon should be in
#'   front of the data points. This is recommended when there are very many data
#'   points and a highly constrained model.
#'
#' @importFrom stats predict
#' 
#' @return No return value.
#'
#' @export
#' @family NLS fit functions
plot.bootstrapfit <- function(x, ..., col.line="black", col.band="gray", opacity.band=0.65, lty=c(1), lwd=c(1), supports=1000, plot.range, error=x$error.function, ribbon.on.top = TRUE) {
  # The plot object might not have a mask, we want to have one in either case.
  if (is.null(x$mask)) {
    x$mask <- rep(TRUE, length(x$x))
  }
  if(missing(plot.range)){
    rx <- range(x$x[x$mask])
  }else{
    rx <- plot.range
  }
  X <- seq(rx[1], rx[2], (rx[2]-rx[1])/supports)
  prediction <- predict(x, X)
  
  plot_data <- function (rep) {
    if(x$errormodel == "yerrors") {
      limits <- plotwitherror(x=x$x, y=x$y, dy=x$dy, rep = rep, ...)
    }
    else {
      limits <- plotwitherror(x=x$x, y=x$y, dy=x$dy, dx=x$dx, rep = rep, ...)
    }
    limits
  }
  
  # We plot all the data first to get the xlim and ylim right.
  limits <- plot_data(rep = FALSE)
  
  # Plot the ribbon.
  if(!is.null(prediction$err)) {
    xlim <- limits$xlim
    ylim <- limits$ylim
    
    polyval <- c(prediction$val + prediction$err, rev(prediction$val - prediction$err))
    if(any(polyval < ylim[1]) || any(polyval > ylim[2])) {
      polyval[polyval < ylim[1]] <- ylim[1]
      polyval[polyval > ylim[2]] <- ylim[2]
    }
    pcol <- col2rgb(col.band, alpha=TRUE)/255 
    pcol[4] <- opacity.band
    pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
    polygon(x=c(X, rev(X)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)
  }
  
  ## plot the fitted curve on top
  lines(x=X, y=prediction$val, col=col.line, lty=lty, lwd=lwd)
  
  # Optionally plot the data on top one more time.
  if (!ribbon.on.top) {
    plot_data(rep = TRUE)
  }
}

#' residual_plot
#'
#' generic residual_plot method
#'
#' @param x the object to plot
#' @param ... additional parameters to be passed on to specialised functions
#' 
#' @return
#' No return value.
#' 
#' @export
residual_plot <- function (x, ...) {
  UseMethod("residual_plot", x)
}

#' @export
residual_plot.bootstrapfit <- function (x, ..., error_fn = x$error.function, operation = `/`) {
  if (is.null(x$mask)) {
    x$mask <- rep(TRUE, length(x$x))
  }
  if (is.logical(x$mask)) {
    x$mask <- which(x$mask)
  }

  prediction <- predict(x, x$x)

  residual_val <- operation(x$y, prediction$val)
  # We want to subtract or divide (depending on the given `operation`) the
  # samples of the data and the central value of the prediction. Due to
  # major-layout one needs to transpose twice to get the operation applied the
  # right way.
  residual_boot <- t(operation(t(x$bsamples[, 1:length(x$y)]), prediction$val))
  residual_err <- apply(residual_boot, 2, error_fn)
  
  band_val <- operation(prediction$val, prediction$val)
  band_boot <- t(operation(t(prediction$boot), prediction$val))
  band_err <- apply(band_boot, 2, error_fn)
  
  # First we plot all the data points to get the xlim and ylim right.
  plotwitherror(x=x$x, y=residual_val, dy=residual_err, ...)
  
  # Error band and central fit line.
  polygon(x = c(x$x, rev(x$x)),
          y = c(band_val - band_err, rev(band_val + band_err)),
          border = NA,
          col = rgb(0, 0, 0, alpha = 0.08))
  lines(x = x$x,
        y = band_val,
        col = 'gray70')
  
  # Plot points which are not used in the fit.
  plot_args <- list(x=x$x[x$mask], y=residual_val[x$mask], dy=residual_err[x$mask], rep = TRUE, ...)
  if(x$errormodel == "xyerrors") {
    plot_args$dx <- x$dx[x$mask]
  }
  do.call(plotwitherror, plot_args)
  
  # Plot points which are used in the fit.
  plot_args <- list(x=x$x[-x$mask], y=residual_val[-x$mask], dy=residual_err[-x$mask], col = 'gray40', rep = TRUE, ...)
  if(x$errormodel == "xyerrors") {
    plot_args$dx <- x$dx[-x$mask]
  }
  if (length(plot_args$x) > 0) {
    do.call(plotwitherror, plot_args)
  }
}

#' Predict values for bootstrapfit
#'
#' @param object Object of type bootstrapfit.
#' @param x Numeric vector with independent variable.
#' @param error Function to compute error from samples.
#' @param ... additional parameters to be passed on to the prediction function.
#'
#' @return
#' List with independent variable `x`, predicted central value `val`, error
#' estimate `err` and sample matrix `boot`.
#'
#' @export
#' @family NLS fit functions
predict.bootstrapfit <- function (object, x, error = object$error.function, ...) {
  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  npar <- length(object$par.guess)
  val <- do.call(object$fn, c(list(par = object$t0[1:npar], x = x, boot.r = 0), object$tofn))

  prediction <- list(x = x, val = val)

  if(!is.null(object$t)) {
    ## error band
    ## define a dummy function to be used in apply
    prediction_boot_fn <- function (boot.r) {
      par <- object$t[boot.r, 1:npar, drop = FALSE]
      do.call(object$fn, c(list(par = par, x = x, boot.r = boot.r), object$tofn))
    }
    prediction_boot <- do.call(rbind, lapply(1:nrow(object$t), prediction_boot_fn))
    prediction$boot <- prediction_boot

    err <- apply(prediction_boot, 2, error, na.rm = TRUE)
    stopifnot(length(err) == length(x))
    prediction$err <- err
  }

  return (prediction)
}
