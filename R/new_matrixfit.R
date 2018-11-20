library(R6)

Model <- R6Class(
    'Model',
    list(
        prediction = function (par, t) {
            stop('This is an abstract function.')
        },
        prediction_jacobian = function (par, t) {
            stop('This is an abstract function.')
        }
    )
)

SingleModel <- R6Class(
    'SingleModel',
    inherit = Model,
    public = list(
        initialize = function (time_extent, parind, sign_vec, ov_sign_vec) {
            self$time_extent <- time_extent
            self$parind <- parind
            self$sign_vec <- sign_vec
            self$ov_sign_vec <- ov_sign_vec
        },
        prediction = function (par, t) {
            self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
                (exp(-par[1] * t) + self$sign_vec * exp(-par[1] * (self$time_extent - t)))
        },
        prediction_gradient = function (par, t) {
            ## Derivative with respect to the mass, `par[1]`.
            zp <- self$ov_sign_vec * 0.5 * par[parind[, 1]] * par[parind[, 2]] *
                (-t * exp(-par[1] * t) -
                     (self$time_extent-t) * self$sign_vec * exp(-par[1] * (self$time_extent-t)))
            res <- zp

            ## Derivatives with respect to the amplitudes.
            for (i in 2:length(par)) {
                zp1 <- rep(0, length(zp))
                j <- which(parind[, 1] == i)
                zp1[j] <- -self$ov_sign_vec * 0.5 * par[parind[j, 2]] *
                    (exp(-par[1] * t[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-t[j])))

                zp2 <- rep(0, length(zp))
                j <- which(parind[, 2] == i)
                zp2[j] <- -self$ov_sign_vec * 0.5 * par[parind[j, 1]] *
                    (exp(-par[1] * t[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-t[j])))

                res <- c(res, zp1 + zp2)
            }

            return (res)
        }
    ),
    private = list(
        time_extent = NA,
        parind = NA,
        sign_vec = NA,
        ov_sign_vec = NA
    )
)

Phi4Model <- R6Class(
    'Phi4Model',
    inherit = Model,
    public = list(
        initialize = function (time_extent, sign_vec, ov_sign_vec, n_particle) {
            self$time_extent <- time_extent
            self$sign_vec <- sign_vec
            self$ov_sign_vec <- ov_sign_vec
            self$n_particle <- n_particle
        },
        prediction = function (par, t) {
            return (10 + self$n_particle)
        }
    ),
    private = list(
        time_extent = NA,
        sign_vec = NA,
        ov_sign_vec = NA,
        n_particle = NA
    )
)

new_matrixfit <- function(cf,
                          t1, t2,
                          parlist,
                          sym.vec,
                          neg.vec,
                          useCov = FALSE,
                          model = "single",
                          boot.fit = TRUE,
                          fit.method = "optim",
                          autoproceed = FALSE,
                          every) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  if(model == 'pc') {
    stopifnot(inherits(cf, 'cf_principal_correlator'))
  }

  stopifnot(cf$symmetrised == TRUE)

  t1p1 <- t1 + 1
  t2p1 <- t2 + 1

  N <- dim(cf$cf)[1]
  Thalfp1 <- cf$Time/2 + 1
  t <- c(0:(cf$Time/2))

  ## This is the number of correlators in cf
  if(!is.null(dim(cf$cf)))
    mSize <- dim(cf$cf)[2] / Thalfp1
  else
    mSize <- dim(cf$cf.tsboot$t)[2] / Thalfp1

  if (model == 'pc' && mSize != 1) {
    stop('For model pc only a 1x1 matrix is allowed.')
  }

  if (missing(parlist)) {
    if (mSize == 1) {
      parlist <- array(c(1, 1), dim = c(2, 1))
      warning("missing parlist, using default for single correlator!")
    }
    else if (mSize == 4) {
      parlist <- array(c(1, 1, 1, 2, 2, 1, 2, 2), dim = c(2, 4))
      warning("missing parlist, using default for four correlators!")
    }
    else {
      stop("parlist is missing and no default is available for this cf size!")
    }
  }

  if (missing(sym.vec)) {
    if (mSize == 1) {
      sym.vec <- c("cosh")
      warning("missing sym.vec, using default for single correlator!")
    }
    else if(mSize == 4) {
      sym.vec <- c("cosh", "cosh", "cosh", "cosh")
      warning("missing sym.vec, using default for four correlators!")
    }
    else {
      stop("sym.vec is missing and no default is available for this cf size!")
    }
  }

  if (missing(neg.vec)) {
    if (mSize == 1) {
      neg.vec <- c(1)
      warning("missing neg.vec, using default (correlator positive)!")
    }
    else if (mSize == 4) {
      neg.vec <- c(1, 1, 1, 1)
      warning("missing neg.vec, using default (all correlators positive)!")
    }
    else {
      stop("neg.vec is missing and no default is available for this cf size!")
    }
  }

  ## some sanity checks
  if (min(parlist) <= 0) {
    stop("Elements of parlist must be all > 0!")
  }
  for (i in 1:max(parlist)) {
    if (!any(parlist == i)) {
      stop("not all parameters are used in the fit!")
    }
  }

  if (dim(parlist)[2] != mSize) {
    cat(mSize, dim(parlist)[2], "\n")
    stop("parlist has not the correct length! Aborting! Use e.g. extractSingleCor.cf or c to bring cf to correct number of observables\n")
  }
  if (length(sym.vec) != mSize) {
    stop("sym.vec does not have the correct length! Aborting\n")
  }
  if (length(neg.vec) != mSize){
    stop("neg.vec does not have the correct length! Aborting\n")
  }

  CF <- data.frame(t=t, Cor=cf$cf0, Err=apply(cf$cf.tsboot$t, 2, cf$error_fn))

  ## index vector for timeslices to be fitted
  ii <- c((t1p1):(t2p1))
  if (mSize > 1) {
    for (j in 2:mSize) {
      ii <- c(ii, (t1p1 + (j-1) * Thalfp1):(t2p1 + (j-1) * Thalfp1))
    }
  }
  ## for the pc model we have to remove timeslice reference_time, where the error is zero
  if (model == 'pc') {
    ii <- ii[ii != (reference_time+1)]
  }
  ## use only a part of the time slices for better conditioned cov-matrix
  if (!missing(every)) {
    ii <- ii[ii %% every == 0]
  }

  ## parind is the index vector for the matrix elements
  ## signvec decides on cosh or sinh
  ## ov.sign.vec indicates the overall sign
  parind <- array(1, dim = c(length(CF$Cor), 2))
  sign.vec <- rep(1, times = length(CF$Cor))
  ov.sign.vec <- rep(1, times = length(CF$Cor))
  for (i in 1:mSize) {
    parind[((i-1)*Thalfp1+1):(i*Thalfp1),] <- t(array(parlist[,i]+1, dim=c(2,Thalfp1)))
    if (sym.vec[i] == "sinh")
        sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
    if (sym.vec[i] == "exp")
        sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- 0
    if (neg.vec[i] == -1)
        ov.sign.vec[((i-1)*Thalfp1+1):(i*Thalfp1)] <- -1
  }

  if (model == 'single') {
      model_object <- SingleModel$new(cf$Time, parind, sign.vec, ov.sign.vec)
  } else if (model == 'shifted') {
      cf$deltat
  } else if (model == 'pc') {
      cf$gevp_reference_time
  } else if (model == 'n_particles') {
      cf$n_particles
  }

  # we always use the boostrap samples to estimate the covariance matrix
  CovMatrix <- cf$cov_fn(cf$cf.tsboot$t[, ii])

  ## for uncorrelated chi^2 use diagonal matrix with inverse sd^2
  M <- diag(1 / CF$Err[ii]^2)
  if(useCov) {
    ## compute correlation matrix and compute the correctly normalised inverse
    ## see C. Michael hep-lat/9412087
    M <- try(invertCovMatrix(cf$cf.tsboot$t[, ii], boot.l = cf$boot.l, boot.samples = TRUE, cov_fn = cf$cov_fn), silent = TRUE)
    if (inherits(M, "try-error")) {
      if (autoproceed) {
        M <- diag(1/CF$Err[ii]^2)
        warning("[matrixfit] inversion of variance covariance matrix failed, continuing with uncorrelated chi^2")
        useCov <- FALSE
      } else {
        stop("[matrixfit] inversion of variance covariance matrix failed!")
      }
    }
  }

  par <- numeric(max(parind))
  ## we get initial guesses for fit parameters from effective masses
  ## first is the mass
  ## (we currently allow for only one)
  if (model == 'pc') {
    ## the ground state energy
    par[1] <- log(CF$Cor[reference_time+1]/CF$Cor[reference_time+2])
    ## the deltaE
    par[2] <- log(CF$Cor[reference_time+1]/CF$Cor[reference_time+2]) - par[1]
    par[2] <- 1.
    ## the amplitude
    par[3] <- 1.
  }
  else {
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
  }

  ## check out constrOptim
  ## now perform minimisation
  dof <- (length(CF$t[ii])-length(par))
  opt.res <- NA
  rchisqr <- 0.
  if(lm.avail) {
    opt.res <- nls.lm(par = par, fn = fitfn, jac=dfitfn, t=CF$t[ii], y=CF$Cor[ii], L=LM, T=cf$Time, deltat=deltat,
                      parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii], reference_time=reference_time,
                      control = nls.lm.control(ftol=1.e-8, ptol=1.e-8, maxiter=500, maxfev=5000))
    if( !(opt.res$info %in% c(1,2,3) ) ){
      cat(sprintf("Termination reason of nls.lm opt.res$info: %d\n", opt.res$info))
    }
    rchisqr <- opt.res$rsstrace[length(opt.res$rsstrace)]
  }
  else {
    opt.res <- optim(par, fn = fitfn, gr = dfitfn,
                     method="BFGS", control=list(maxit=500, parscale=par, ndeps=rep(1.e-8, times=length(par)), REPORT=50),
                     t=CF$t[ii], y=CF$Cor[ii], M=M, T=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], reference_time=reference_time,
                     ov.sign.vec=ov.sign.vec[ii], deltat=deltat)
    rchisqr <- opt.res$value
  }
  Qval <- 1-pchisq(rchisqr, dof)
  ## we use absolute values in the fit model
  ## better remove any signs then here!
  if(pcmodel) {
    opt.res$par[1:2] <- abs(opt.res$par[1:2])
  }

  opt.tsboot <- NA
  if(boot.fit) {
    opt.tsboot <- apply(X=cf$cf.tsboot$t[,ii], MARGIN=1, FUN=fit.formatrixboot, par=opt.res$par, t=CF$t[ii], deltat=deltat,
                        M=M, T=cf$Time, parind=parind[ii,], sign.vec=sign.vec[ii], ov.sign.vec=ov.sign.vec[ii],
                        L=LM, lm.avail=lm.avail, fitfn=fitfn, dfitfn=dfitfn, reference_time=reference_time)
  }
  N <- length(cf$cf[,1])
  if(is.null(cf$cf)) {
    N <- cf$N
  }
  if(pcmodel) {
    opt.tsboot[c(1:2),] <- abs(opt.tsboot[c(1:2),])
  }
  res <- list(CF=CF, M=M, L=LM, parind=parind, sign.vec=sign.vec, ov.sign.vec=ov.sign.vec, ii=ii, opt.res=opt.res, opt.tsboot=opt.tsboot,
              boot.R=cf$boot.R, boot.l=cf$boot.l, useCov=useCov, CovMatrix=CovMatrix, invCovMatrix=M, seed=cf$seed,
              Qval=Qval, chisqr=rchisqr, dof=dof, mSize=mSize, cf=cf, t1=t1, t2=t2, reference_time=reference_time,
              parlist=parlist, sym.vec=sym.vec, N=N, model=model, fit.method=fit.method, error_fn=cf$error_fn)
  res$t <- t(opt.tsboot)
  res$t0 <- c(opt.res$par, opt.res$value)
  res$se <- apply(opt.tsboot[c(1:(dim(opt.tsboot)[1]-1)),], MARGIN=1, FUN=cf$error_fn)
  attr(res, "class") <- c("matrixfit", "list")
  return(invisible(res))
}
