library(R6)

MatrixModel <- R6Class(
  'MatrixModel',
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec) {
      self$time_extent <- time_extent
      self$parind <- parind
      self$sign_vec <- sign_vec
      self$ov_sign_vec <- ov_sign_vec
    },
    prediction = function (par, x) {
      stop('This is an abstract function.')
    },
    prediction_jacobian = function (par, x) {
      stop('This is an abstract function.')
    },
    initial_guess = function (corr, parlist, t1, t2) {
      t1p1 <- t1 + 1
      t2p1 <- t2 + 1
      Thalfp1 <- self$time_extent / 2 + 1
      
      par <- numeric(max(self$parind))
      j <- which(parlist[1, ] == 1 & parlist[2, ] == 1)
      par[1] <- invcosh(corr[t1p1 + (j-1) * Thalfp1] / corr[t1p1 + (j-1) * Thalfp1 + 1], t = t1p1, self$time_extent)
      ## catch failure of invcosh
      if(is.na(par[1]) || is.nan(par[1]))
        par[1] <- 0.2
      ## the amplitudes we estimate from diagonal elements
      for (i in 2:length(par)) {
        j <- which(parlist[1, ] == (i-1) & parlist[2, ] == (i-1))
        if (length(j) == 0) {
          ##if(full.matrix) warning("one diagonal element does not appear in parlist\n")
          j <- i-1
        }
        par[i] <- sqrt(abs(corr[t1p1 + (j-1) * Thalfp1]) / 0.5 / exp(-par[1] * t1))
      }

      return (par)
    },
    time_extent = NA,
    parind = NA,
    sign_vec = NA,
    ov_sign_vec = NA
  )
)

SingleModel <- R6Class(
  'SingleModel',
  inherit = MatrixModel,
  public = list(
    prediction = function (par, x) {
      self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
        (exp(-par[1] * x) + self$sign_vec * exp(-par[1] * (self$time_extent - x)))
    },
    prediction_gradient = function (par, x) {
      ## Derivative with respect to the mass, `par[1]`.
      zp <- self$ov_sign_vec * 0.5 * par[parind[, 1]] * par[parind[, 2]] *
        (-t * exp(-par[1] * x) -
           (self$time_extent-x) * self$sign_vec * exp(-par[1] * (self$time_extent-x)))
      res <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
        zp1 <- rep(0, length(zp))
        j <- which(parind[, 1] == i)
        zp1[j] <- -self$ov_sign_vec * 0.5 * par[parind[j, 2]] *
          (exp(-par[1] * x[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        zp2 <- rep(0, length(zp))
        j <- which(parind[, 2] == i)
        zp2[j] <- -self$ov_sign_vec * 0.5 * par[parind[j, 1]] *
          (exp(-par[1] * x[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        res <- c(res, zp1 + zp2)
      }
      
      return (res)
    }
  )
)

TwoStateModel <- R6Class(
  'TwoStateModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec, reference_time) {
      super$initialize(time_extent, parind, sign_vec, ov_sign_vec)
      self$reference_time <- reference_time
    },
    initial_guess = function (corr, parlist, t1, t2) {
      par = numeric(3)

      ## the ground state energy
      par[1] <- log(corr[self$reference_time+1]/corr[self$reference_time+2])
      ## the deltaE
      par[2] <- log(corr[self$reference_time+1]/corr[self$reference_time+2]) - par[1]
      par[2] <- 1.0
      ## the amplitude
      par[3] <- 1.0

      return (par)
    },
    reference_time = NA
  )
)

Phi4Model <- R6Class(
  'Phi4Model',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec, n_particle) {
      super$initialize(time_extent, parind, sign_vec, ov_sign_vec)
      self$n_particle <- n_particle
    },
    prediction = function (par, x) {
      return (10 + self$n_particle)
    },
    n_particle = NA
  ),
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
  
  CF <- data.frame(t = t, Cor = cf$cf0, Err = apply(cf$cf.tsboot$t, 2, cf$error_fn))
  
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
    else if (sym.vec[i] == "exp")
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

  args <- list(fn = model_object$prediction,
               par.guess = model_object$initial_guess(CF$Cor, parlist, t1, t2),
               y = CF$Cor[ii],
               x = CF$t[ii],
               bsamples = cf$cf.tsboot$t[, ii],
               use.minpack.lm = fit.method == 'lm',
               error = cf$error_fn,
               cov_fn = cf$cov_fn)
  
  if (useCov) {
    args$CovMatrix <- cf$cov_fn(cf$cf.tsboot$t[, ii])
  }
  
  do.call(bootstrap.nlsfit, args)
}
