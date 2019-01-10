library(R6)

MatrixModel <- R6Class(
  'MatrixModel',
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec) {
      self$time_extent <- time_extent
      self$parind <- parind
      self$sign_vec <- sign_vec
      self$ov_sign_vec <- ov_sign_vec

      stopifnot(!is.na(self$sign_vec))
      stopifnot(!is.na(self$ov_sign_vec))
    },
    prediction = function (par, x, ...) {
      stop('MatrixModel$prediction is an abstract function.')
    },
    prediction_jacobian = function (par, x, ...) {
      stop('MatrixModel$prediction_jacobian is an abstract function.')
    },
    initial_guess = function (corr, parlist, t1, t2) {
      t1p1 <- t1 + 1
      t2p1 <- t2 + 1
      Thalfp1 <- (self$time_extent / 2) + 1
      
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
    post_process_par = function (par) {
      par
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
    prediction = function (par, x, ...) {
      self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
        (exp(-par[1] * x) + self$sign_vec * exp(-par[1] * (self$time_extent - x)))
    },
    prediction_jacobian = function (par, x, ...) {
      ## Derivative with respect to the mass, `par[1]`.
      zp <- self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
        (-x * exp(-par[1] * x) -
           (self$time_extent-x) * self$sign_vec * exp(-par[1] * (self$time_extent-x)))
      res <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
        zp1 <- rep(0, length(zp))
        j <- which(self$parind[, 1] == i)
        zp1[j] <- self$ov_sign_vec * 0.5 * par[self$parind[j, 2]] *
          (exp(-par[1] * x[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        zp2 <- rep(0, length(zp))
        j <- which(self$parind[, 2] == i)
        zp2[j] <- self$ov_sign_vec * 0.5 * par[self$parind[j, 1]] *
          (exp(-par[1] * x[j]) + self$sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        res <- cbind(res, zp1 + zp2)
      }

      stopifnot(nrow(res) == length(x))
      stopifnot(ncol(res) == length(par))

      dimnames(res) <- NULL
      
      return (res)
    }
  )
)

ShiftedModel <- R6Class(
  'ShiftedModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec, delta_t) {
      super$initialize(time_extent, parind, sign_vec, ov_sign_vec)
      self$delta_t <- delta_t
    },
    prediction = function (par, x, ...) {
      xx <- x - self$delta_t/2
      self$ov_sign_vec * par[self$parind[, 1]] * par[self$parind[, 2]] *
        (exp(-par[1] * xx) - self$sign_vec * exp(-par[1] * (self$time_extent - xx)))
    },
    prediction_jacobian = function (par, x, ...) {
      xx <- x - self$delta_t/2
      
      res <- matrix(0.0, nrow = length(x), ncol = length(par))
      
      ## Derivative with respect to the mass, `par[1]`.
      zp <- self$ov_sign_vec * par[self$parind[, 1]] * par[self$parind[, 2]] *
        (-xx * exp(-par[1] * xx) + (self$time_extent - xx) * self$sign_vec *
           exp(-par[1] * (self$time_extent - xx)))
      stopifnot(length(zp) == nrow(res))
      
      res[, 1] <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
          zp1 <- rep(0, length(zp))##
          j <- which(self$parind[, 1] == i)
          zp1[j] <- self$ov_sign_vec * par[self$parind[j, 2]] *
            (exp(-par[1] * xx[j]) - self$sign_vec[j] *
            exp(-par[1] * (self$time_extent - xx[j])))
        
          zp2 <- rep(0, length(zp))##
          j <- which(self$parind[, 2] == i)
          zp2[j] <- self$ov_sign_vec * par[self$parind[j, 1]] *
            (exp(-par[1] * xx[j]) - self$sign_vec[j] *
            exp(-par[1] * (self$time_extent - xx[j])))
          
          stopifnot(length(zp1) == nrow(res))
          stopifnot(length(zp2) == nrow(res))
        
          res[, i] <- (zp1 + zp2)
      }
      
      return (res)
    },
    delta_t = NA
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
    prediction = function (par, x, ...) {
      par[1] <- abs(par[1])
      par[2] <- abs(par[2])
      xx <- x - self$reference_time

      exp(-par[1] * xx) * (par[3] + (1 - par[3]) * exp(-(par[2]) * xx))
    },
    prediction_jacobian = function (par, x, ...) {
      par[1] <- abs(par[1])
      par[2] <- abs(par[2])
      xx <- x - self$reference_time
      
      res <- array(0.0, dim = c(length(x), length(par)))
      res[, 1] <- -xx * exp(-par[1] * xx) * (par[3] + (1 - par[3]) * exp(-par[2] * xx))
      res[, 2] <- -exp(-par[1] * xx) * (1 - par[3]) * xx * exp(-par[2] * xx)
      res[, 3] <- exp(-par[1] * xx) * (1 - exp(-par[2] * xx))
      return(res)
    },
    initial_guess = function (corr, parlist, t1, t2) {
      par = numeric(3)

      ## the ground state energy
      par[1] <- log(corr[self$reference_time + 1] / corr[self$reference_time + 2])
      ## the deltaE
      par[2] <- log(corr[self$reference_time + 1] / corr[self$reference_time + 2]) - par[1]
      par[2] <- 1.0
      ## the amplitude
      par[3] <- 1.0

      return (par)
    },
    post_process_par = function (par) {
      ## The energy and energy difference parameter enters the model only in its
      ## absolute value, therefore it can become negative. We need to fix that
      ## here.
      par[1] <- abs(par[1])
      par[2] <- abs(par[2])
      return (par)
    },
    reference_time = NA
  )
)

#' 3 particle correlation function in phi^4 theory
#'
#' @description
#' 3 particle correlator with thermal pollution term:
#' \deqn{C_3 (t) = 0.5*A_1^(3)*A_2^(3)*(exp(-E_3*(T-t))+exp(-E_3*t)) + 0.5*A_3^(3)*A_4^(3)*exp(-(2/3)*E_3*(T/2)) * (exp(-(1/3)*E_3*(T-t)) + exp(-(1/3)*E_3*t))}
Phi4Model <- R6Class(
  'Phi4Model',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parind, sign_vec, ov_sign_vec, delta_t) {
      super$initialize(time_extent, parind, sign_vec, ov_sign_vec)
      self$delta_t <- delta_t
    },
    prediction = function (par, x) {
      self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
      (exp(-par[1] * (x - self$delta_t/2)) + self$sign_vec * exp(-par[1] * (self$time_extent - (x - self$delta_t/2)))) +
      self$ov_sign_vec*par[self$parind[,3]]*par[self$parind[,4]]*exp(-(par[1]*(2/3))*(self$time_extent/2))*
      (exp(-(par[1]*(1/3))*(x-self$delta_t/2)) + self$sign_vec*exp(-(par[1]*(1/3))*(self$time_extent-(x-self$delta_t/2))))
    },
    prediction_jacobian = function (par, x, ...) {
      xx <- x - self$delta_t/2
      
      res <- matrix(0.0, nrow = length(x), ncol = length(par))
      
      ## Derivative with respect to the mass, `par[1]`.
      zp <- self$ov_sign_vec * 0.5 * par[self$parind[, 1]] * par[self$parind[, 2]] *
            ((-xx) * exp(-par[1] * xx) - (self$time_extent - xx) * self$sign_vec *
            exp(-par[1] * (self$time_extent - xx))) + self$ov_sign_vec * 0.5 * par[self$parind[,3]] * par[self$parind[,4]] *
            (-2/3)*(self$time_extent/2)*exp(-(par[1]*(-2/3))*(self$time_extent/2))*(exp(-(par[1]*(1/3))*xx) +      
            self$sign_vec*exp(-(par[1]*(1/3))*(self$time_extent-xx))) +
            self$ov_sign_vec * 0.5 * par[self$parind[,3]] * par[self$parind[,4]]*exp(-(par[1]*(2/3))*
            (self$time_extent/2))*(-(1/3)*xx*exp(-(par[1]*(1/3))*xx) +
            self$sign_vec*(-1/3)*(self$time_extent-xx)*exp(-(par[1]*(1/3))*(self$time_extent-xx)))
      
      stopifnot(length(zp) == nrow(res))
      res[, 1] <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
        
        zp1 <- rep(0, length(zp))
        j <- which(self$parind[,1]==i)
        zp1 <- -self$ov_sign_vec*0.5*par[self$parind[j,2]] *
               (exp(-par[1]*(x[j]-self$delta_t/2)) + self$sign_vec[j] * exp(-par[1]*(self$time_extent-(x[j]-self$delta_t/2))))
        
        zp2 <- rep(0, length(zp))
        j <- which(self$parind[,2]==i)
        zp2 <- -self$ov_sign_vec*0.5*par[self$parind[j,1]] *
               (exp(-par[1]*(x[j]-self$delta_t/2)) + self$sign_vec[j] * exp(-par[1]*(self$time_extent-(x[j]-self$delta_t/2))))
        
        zp3 <- rep(0, length(zp))
        j <- which(self$parind[,3]==i)
        zp3 <- -self$ov_sign_vec[j]*par[self$parind[j,4]] * 
               exp(-(par[1]*(2/3))*(self$time_extent/2)) * 0.5*(exp(-(par[1]*(1/3))*(x[j]-self$delta_t/2)) +
               self$sign_vec[j]*exp(-(par[1]*(1/3)) * (self$time_extent-(x[j]-self$delta_t/2))))
        
        zp4 <- rep(0, length(zp))
        j <- which(self$parind[,4]==i)
             zp4 <- -self$ov_sign_vec[j]*par[self$parind[j,3]]*
             exp(-(par[1]*(2/3))*(self$time_extent/2)) * 0.5*(exp(-(par[1]*(1/3))*(x[j]-self$delta_t/2)) +
             self$sign_vec[j]*exp(-(par[1]*(1/3))*(self$time_extent-(x[j]-self$delta_t/2))))
        
        stopifnot(length(zp1) == nrow(res))
        stopifnot(length(zp2) == nrow(res))
        stopifnot(length(zp3) == nrow(res))
        stopifnot(length(zp4) == nrow(res))
        
        res[, i] <- (zp1 + zp2 + zp3 + zp4)
      }
      return(res)
    },
    # n_particle = NA
    delta_t = NA
  )
)

#' @export
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
                          par.guess,
                          every) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  
  if (model == 'pc') {
    stopifnot(inherits(cf, 'cf_principal_correlator'))
  }
  
  stopifnot(cf$symmetrised == TRUE)
  
  t1p1 <- t1 + 1
  t2p1 <- t2 + 1
  
  N <- dim(cf$cf)[1]
  Thalfp1 <- cf$Time/2 + 1
  t <- c(0:(cf$Time/2))
  
  deltat <- 1
  if(model == "shifted" && any(names(cf) == "deltat")) {
    deltat <- cf$deltat
  }
  
  ## This is the number of correlators in cf
  if (!is.null(dim(cf$cf)))
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
    ii <- ii[ii != (cf$gevp_reference_time + 1)]
  }
  ## use only a part of the time slices for better conditioned cov-matrix
  if (!missing(every)) {
    ii <- ii[ii %% every == 0]
  }
  
  ## parind is the index vector for the matrix elements
  ## signvec decides on cosh or sinh
  ## ov.sign.vec indicates the overall sign
  len_t <- length(t1:t2)
  parind <- array(1, dim = c(len_t, 2))
  sign.vec <- rep(1, times = len_t)
  ov.sign.vec <- rep(1, times = len_t)
  for (i in 1:mSize) {
    parind[((i-1)*len_t+1):(i*len_t), ] <- t(array(parlist[, i] + 1, dim = c(2, len_t)))

    if (sym.vec[i] == "sinh")
      sign.vec[((i-1)*len_t+1):(i*len_t)] <- -1
    else if (sym.vec[i] == "exp")
      sign.vec[((i-1)*len_t+1):(i*len_t)] <- 0

    if (neg.vec[i] == -1)
      ov.sign.vec[((i-1)*len_t+1):(i*len_t)] <- -1
  }

  ## perform the bootstrap non-linear least-squares fit (NLS fit):
  
  if (model == 'single') {
    model_object <- SingleModel$new(cf$Time, parind, sign.vec, ov.sign.vec)
  } else if (model == 'shifted') {
    stopifnot(inherits(cf, 'cf_shifted'))
    model_object <- ShiftedModel$new(cf$Time, parind, sign.vec, ov.sign.vec, cf$deltat)
  } else if (model == 'pc') {
    model_object <- TwoStateModel$new(cf$Time, parind, sign.vec, ov.sign.vec, cf$gevp_reference_time)
  } else if (model == 'n_particles') {
    model_object <- Phi4Model$new(cf$Time, parind, sign.vec, ov.sign.vec)#cf$n_particles)
  }
  
  if (missing(par.guess)) {
    par.guess <- model_object$initial_guess(CF$Cor, parlist, t1, t2)
  }

  args <- list(fn = model_object$prediction,
               gr = model_object$prediction_jacobian,
               par.guess = par.guess,
               y = CF$Cor[ii],
               x = CF$t[ii],
               bsamples = cf$cf.tsboot$t[, ii],
               use.minpack.lm = fit.method == 'lm',
               error = cf$error_fn,
               cov_fn = cf$cov_fn)
  
  if (useCov) {
    args$CovMatrix <- cf$cov_fn(cf$cf.tsboot$t[, ii])
  }
  
  res <- do.call(bootstrap.nlsfit, args)

  ## Some fit models have parameters in the absolute value. This means that in
  ## `res` they can be negative and we need to let the model fix that.
  res$t0 <- model_object$post_process_par(res$t0)
  
  old_dim = dim(res$t)
  res$t <- t(apply(res$t, 1, model_object$post_process_par))
  stopifnot(all(old_dim == dim(res$t)))
  
  return (res)
}

#' Create a parameter list for `matrixfit`
#'
#' @param corr_matrix_size integer. Number of correlators in the matrix. This
#' must be a the square of an integer.
#'
#' @export
make_parlist <- function (corr_matrix_size) {
  ## Placeholder, this is not correct.
  return (array(NA, dim = c(2, corr_matrix_size)))
}

#' Create a parameter index matrix for `matrixfit`
#'
#' @param parlist integer array. Parameter list generated with `make_parlist`.
#' @param length_time integer. Number of time slices per correlator.
#' @param summands integer. Number of summands in the fit model that shall be
#' fitted. The signal counts as one summand, each explicit pollution term with
#' independent amplitudes counts as its own summand.
#'
#' @export
make_parind <- function (parlist, length_time, summands = 1) {
  ## Placeholder, this is not correct.
  return (array(NA, dim = c(ncol(parlist) * length_time, 2 * summands)))
}
