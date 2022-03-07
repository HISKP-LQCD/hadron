# signvec decides on cosh or sinh
make_sign_vec <- function (sym_vec, len_t, m_size = 1) {
  sign_vec <- rep(1, times = len_t)

  for (i in 1:m_size) {
    if (sym_vec[i] == "sinh")
      sign_vec[((i-1)*len_t+1):(i*len_t)] <- -1
    else if (sym_vec[i] == "exp")
      sign_vec[((i-1)*len_t+1):(i*len_t)] <- 0
  }

  return (sign_vec)
}

# ov.sign.vec indicates the overall sign
make_ov_sign_vec <- function (neg_vec, len_t, m_size = 1) {
  ov_sign_vec <- rep(1, times = len_t)
  
  for (i in 1:m_size) {
    if (neg_vec[i] == -1)
      ov_sign_vec[((i-1)*len_t+1):(i*len_t)] <- -1
  }

  return (ov_sign_vec)
}

MatrixModel <- R6::R6Class(
  'MatrixModel',
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size) {
      self$time_extent <- time_extent
      self$parlist <- parlist 
      self$sym_vec <- sym_vec
      self$neg_vec <- neg_vec
      self$m_size <- m_size
    },
    prediction = function (par, x, ...) {
      stop('MatrixModel$prediction is an abstract function.')
    },
    prediction_jacobian = function (par, x, ...) {
      stop('MatrixModel$prediction_jacobian is an abstract function.')
    },
    initial_guess = function (corr, t1, t2) {
      t1p1 <- t1 + 1
      t2p1 <- t2 + 1
      Thalfp1 <- (self$time_extent / 2) + 1
      len_t <- length(t1:t2)
      parind <- make_parind(parlist = self$parlist, length_time = len_t, summands = 1)
      
      par <- numeric(max(parind))
      j <- which(self$parlist[1, ] == 1 & self$parlist[2, ] == 1)
      par[1] <- invcosh(corr[t1p1 + (j-1) * Thalfp1] / corr[t1p1 + (j-1) * Thalfp1 + 1], t = t1p1, self$time_extent)
      ## catch failure of invcosh
      if(is.na(par[1]) || is.nan(par[1]))
        par[1] <- 0.2
      ## the amplitudes we estimate from diagonal elements
      for (i in 2:length(par)) {
        j <- which(self$parlist[1, ] == (i-1) & self$parlist[2, ] == (i-1))
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
    parlist = NA,
    sym_vec = NA,
    neg_vec = NA,
    m_size = NA
  )
)

SingleModel <- R6::R6Class(
  'SingleModel',
  inherit = MatrixModel,
  public = list(
    prediction = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      ov_sign_vec * 0.5 * par[parind[, 1]] * par[parind[, 2]] *
        (exp(-par[1] * x) + sign_vec * exp(-par[1] * (self$time_extent - x)))
    },
    prediction_jacobian = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      ## Derivative with respect to the mass, `par[1]`.
      zp <- ov_sign_vec * 0.5 * par[parind[, 1]] * par[parind[, 2]] *
        (-x * exp(-par[1] * x) -
           (self$time_extent-x) * sign_vec * exp(-par[1] * (self$time_extent-x)))
      res <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
        zp1 <- rep(0, length(zp))
        j <- which(parind[, 1] == i)
        zp1[j] <- ov_sign_vec * 0.5 * par[parind[j, 2]] *
          (exp(-par[1] * x[j]) + sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        zp2 <- rep(0, length(zp))
        j <- which(parind[, 2] == i)
        zp2[j] <- ov_sign_vec * 0.5 * par[parind[j, 1]] *
          (exp(-par[1] * x[j]) + sign_vec[j] * exp(-par[1] * (self$time_extent-x[j])))
        
        res <- cbind(res, zp1 + zp2)
      }

      stopifnot(nrow(res) == length(x))
      stopifnot(ncol(res) == length(par))

      dimnames(res) <- NULL
      
      return (res)
    }
  )
)

TwoAmplitudesModel <- R6::R6Class(
  'TwoAmplitudesModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
    },
    prediction = function (par, x, ...) {
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      ov_sign_vec * 0.5 * (par[2]^2 * exp(-par[1] * x) + sign_vec * par[3]^2 * exp(par[1] * x))
    },
    prediction_jacobian = function (par, x, ...) {
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      res <- matrix(NA, nrow = length(x), ncol = length(par))
      res[, 1] <- ov_sign_vec * 0.5 * (par[2]^2 * exp(-par[1] * x) * (-x) + sign_vec * par[3]^2 * exp(par[1] * x) * x)
      res[, 2] <- ov_sign_vec * par[2] * exp(-par[1] * x)
      res[, 3] <- ov_sign_vec * par[3] * sign_vec * exp(par[1] * x)
      return (res)
    },
    initial_guess = function (corr, t1, t2) {
      # This model has two amplitudes but only one energy. We will start by
      # just replicating the forward amplitude to the backwards part.
      normal <- super$initial_guess(corr, t1, t2)
      c(normal, normal[2])
    }
  )
)

ShiftedModel <- R6::R6Class(
  'ShiftedModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size, delta_t) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
      self$delta_t <- delta_t
    },
    prediction = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      xx <- x - self$delta_t/2
      ov_sign_vec * par[parind[, 1]] * par[parind[, 2]] *
        (exp(-par[1] * xx) - sign_vec * exp(-par[1] * (self$time_extent - xx)))
    },
    prediction_jacobian = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)

      xx <- x - self$delta_t/2
      
      res <- matrix(0.0, nrow = length(x), ncol = length(par))
      
      ## Derivative with respect to the mass, `par[1]`.
      zp <- ov_sign_vec * par[parind[, 1]] * par[parind[, 2]] *
        (-xx * exp(-par[1] * xx) + (self$time_extent - xx) * sign_vec *
           exp(-par[1] * (self$time_extent - xx)))
      stopifnot(length(zp) == nrow(res))
      
      res[, 1] <- zp
      
      ## Derivatives with respect to the amplitudes.
      for (i in 2:length(par)) {
          zp1 <- rep(0, length(zp))##
          j <- which(parind[, 1] == i)
          zp1[j] <- ov_sign_vec * par[parind[j, 2]] *
            (exp(-par[1] * xx[j]) - sign_vec[j] *
            exp(-par[1] * (self$time_extent - xx[j])))
        
          zp2 <- rep(0, length(zp))##
          j <- which(parind[, 2] == i)
          zp2[j] <- ov_sign_vec * par[parind[j, 1]] *
            (exp(-par[1] * xx[j]) - sign_vec[j] *
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

WeightedModel <- R6::R6Class(
  'WeightedModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size, delta_t, weight_factor) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
      self$delta_t <- delta_t
      self$weight_factor <- weight_factor
    },
    prediction = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)
      
      # Just defined the variables such that the changes to the Mathematica
      # generated expression are minimal.
      Power <- `^`
      a0 <- par[2]
      deltat <- self$delta_t
      e0 <- par[1]
      sign <- sign_vec
      t <- x
      time <- self$time_extent
      w <- self$weight_factor
      
      # Following is the Mathematica code which has been edited slightly to be
      # R syntax.
      pred <- ov_sign_vec * (a0*(-((exp(2*e0*time) + exp(e0*(2*t + time))*sign)*Power(w,2*t)) - sign*(exp(2*e0*time) + exp(e0*(2*t + time))*sign)*Power(w,2*deltat + time) + (exp(e0*(deltat + 2*time)) + exp(e0*(-deltat + 2*t + time))*sign)*Power(w,deltat)* (Power(w,2*t) + sign*Power(w,time)))) / (exp(e0*(t + 2*time))*Power(w,deltat)*(Power(w,2*t) + sign*Power(w,time)))

      stopifnot(length(pred) == length(x))
      return (pred)
    },
    prediction_jacobian = function (par, x, ...) {
      parind <- make_parind(self$parlist, length(x), summands = 1)
      sign_vec <- make_sign_vec(self$sym_vec, length(x), self$m_size)
      ov_sign_vec <- make_ov_sign_vec(self$neg_vec, length(x), self$m_size)
      
      Power <- `^`
      a0 <- par[2]
      deltat <- self$delta_t
      e0 <- par[1]
      sign <- sign_vec
      t <- x
      time <- self$time_extent
      w <- self$weight_factor

      d1 <- ov_sign_vec * (a0*(exp(e0*(deltat + 2*time))*(deltat - t)*Power(w,deltat)* (Power(w,2*t) + sign*Power(w,time)) - exp(e0*(-deltat + 2*t + time))*sign*(deltat - t + time)*Power(w,deltat)* (Power(w,2*t) + sign*Power(w,time)) + exp(2*e0*time)*t*(Power(w,2*t) + sign*Power(w,2*deltat + time)) - exp(e0*(2*t + time))*sign*(t - time)*(Power(w,2*t) + sign*Power(w,2*deltat + time))))/ (exp(e0*(t + 2*time))*Power(w,deltat)*(Power(w,2*t) + sign*Power(w,time)))

      d2 <- ov_sign_vec * (-((exp(2*e0*time) + exp(e0*(2*t + time))*sign)*Power(w,2*t)) - sign*(exp(2*e0*time) + exp(e0*(2*t + time))*sign)*Power(w,2*deltat + time) + (exp(e0*(deltat + 2*time)) + exp(e0*(-deltat + 2*t + time))*sign)*Power(w,deltat)* (Power(w,2*t) + sign*Power(w,time)))/ (exp(e0*(t + 2*time))*Power(w,deltat)*(Power(w,2*t) + sign*Power(w,time))) 

      cbind(d1, d2)
    },
    delta_t = NA,
    weight_factor = NA
  )
)

TwoStateModel <- R6::R6Class(
  'TwoStateModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size, reference_time) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
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
    initial_guess = function (corr, t1, t2) {
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

## Model for the n particle correlator fit
##
## n particle correlator with thermal pollution term(s)
NParticleModel <- R6::R6Class(
  'NParticleModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
      self$sm <- SingleModel$new(time_extent, parlist, sym_vec, neg_vec, m_size)
     },
    prediction = function (par, x, ...) {
      n <- length(par) / 2
      y <- 0
      par_aux <- numeric(2)
      for (i in 1:n) {
        par_aux[1] <- par[2*i - 1]
        par_aux[2] <- par[2*i]
        corr_part <- self$sm$prediction(par_aux, x, ...)
        y <- y + corr_part
      }
      return(y)
    },
    prediction_jacobian = function (par, x, ...) {
      n <- length(par) / 2
      par_aux <- numeric(2)
      par_aux[1] <- par[1]
      par_aux[2] <- par[2]
      y <- self$sm$prediction_jacobian(par_aux, x, ...)
      if (n > 1) {
        for (i in 2:n) {
          par_aux[1] <- par[2*i - 1]
          par_aux[2] <- par[2*i]
          jacobian_part_aux <- self$sm$prediction_jacobian(par_aux, x, ...)
          y <- cbind(y, jacobian_part_aux)
        }
      }
      return(y)
    },
    sm = NA
  )
)

## Model for a single correlator and a simple additive constant.
##
## @description
## This model is just the “single” model plus a simple additive
## constant. In some way it is a specialization of the “n particles” model as
## the energy difference is just set to zero.
SingleConstantModel <- R6::R6Class(
  'SingleConstantModel',
  inherit = MatrixModel,
  public = list(
    initialize = function (time_extent, parlist, sym_vec, neg_vec, m_size) {
      super$initialize(time_extent, parlist, sym_vec, neg_vec, m_size)
      self$sm <- SingleModel$new(time_extent, parlist, sym_vec, neg_vec, m_size)
    },
    prediction = function (par, x, ...) {
      self$sm$prediction(par[1:2], x, ...) + par[3]
    },
    prediction_jacobian = function (par, x, ...) {
      y <- self$sm$prediction_jacobian(par[1:2], x, ...)
      cbind(y, rep(1, length(x)))
    },
    initial_guess = function (corr, t1, t2) {
      c(self$sm$initial_guess(corr, t1, t2), 1e-10)
    },
    sm = NA
  )
)

#' perform a factorising fit of a matrix of correlation functions
#'
#' Modernised and extended implementation of \link{matrixfit}
#'
#' @param cf Object of class `cf` with `cf_meta` and `cf_boot`.
#' @param t1 Integer, start time slice of fit range (inclusive).
#' @param t2 Integer, end time slie of fit range (inclusive).
#' @param parlist Numeric vector, list of parameters for the model function.
#' @param sym.vec Integer, numeric or vectors thereof specifying the symmetry
#'   properties of the correlation functions stored in `cf`. See
#'   \link{matrixfit} for details.
#' @param neg.vec Integer or integer vector of global signs, see
#'   \link{matrixfit} for details.
#' @param useCov Boolean, specifies whether a correlated chi^2 fit should be
#'   performed.
#' @param model String, specifies the type of model to be assumed for the
#'   correlator. See below for details.
#' @param boot.fit Boolean, specifies if the fit should be bootstrapped.
#' @param fit.method String, specifies which minimizer should be used. See
#'   \link{matrixfit} for details.
#' @param autoproceed Boolean, if TRUE, specifies that if inversion of the
#'   covariance matrix fails, the function should proceed anyway assuming no
#'   correlation (diagonal covariance matrix).
#' @param par.guess Numeric vector, initial values for the paramters, should be
#'   of the same length as `parlist`.
#' @param every Integer, specifies a stride length by which the fit range should
#'   be sparsened, using just `every`th time slice in the fit.
#' @param higher_states List with elements `val` and `boot`. Only used in the
#'   `n_particles` fit model. The member `val` must have the central energy
#'   values for all the states that are to be fitted. The `boot` member will be
#'   a matrix that has the various states as columns and the corresponding
#'   bootstrap samples as rows. The length of `val` must be the column number of
#'   `boot`. The row number of `boot` must be the number of samples.
#' @param ... Further parameters.
#'
#' @details There are different fit models available. The models generally
#'   depend on one or multiple energies \eqn{E} and amplitudes \eqn{p_i} which
#'   for a general matrix are row- and column-amplitudes. The relative sign
#'   factor \eqn{c \in \{-1, 0, +1\}} depends on the chosen symmetry of the
#'   correlator. It is a plus for a “cosh” symmetry and a minus for a “sinh”
#'   symmetry. If the back propagating part is to be neglected (just “exp”
#'   model), it will be zero.
#'
#'   When the back propagating part is not taken into account, then the
#'   `single`, `shifted` and `weighted` model become the same except for changes
#'   in the amplitude.
#'
#'   - `single`: The default model for a single state correlator is
#'   \deqn{\frac{1}{2} p_i p_j (\exp(-p_1 t) \pm c \exp(-p_2 (T-t))) \,.}
#'
#'   - `shifted`: If the correlator has been shifted (using
#'   \link{takeTimeDiff.cf}, then the following model is applicable: \deqn{p_i
#'   p_j (\exp(-p_1 (t+1/2)) \mp c \exp(-p_1 (T-(t+1/2)))) \,.}
#'
#'   - `weighted`: Works similarly to the `shifted` model but includes the
#'   effect of the weight factor from \link{removeTemporal.cf}.
#'
#'   - `pc`: In case only a single principal correlator from a GEVP is to be
#'   fitted this model can be used. It implements \deqn{\exp(-p_1(t-t_0))(p_2 +
#'   (1-p_2)\exp(-p_3 (t-t_0))} with \eqn{t_0} the reference timesclice of the
#'   GEVP. See \link{bootstrap.gevp} for details.
#'
#'   - `two_amplitudes`: Should there be a single state but different amplitudes
#'   in the forward and backwards part, the following method is applicable.
#'   \deqn{\frac{1}{2} (p_2 \exp(-p_1 t) \pm c p_3 \exp(p_1 t))} This only
#'   works with a single correlator at the moment.
#'
#'   - `single_constant`: Uses the `single` model and simply adds \eqn{+ p_3} to
#'   the model such that a constant offset can be fitted. In total the model is
#'   \deqn{\mathrm{single}(p_1, p_2) + p_3 \,.}
#'
#'   - `n_particles`: A sum of `single` models with independent energies and
#'   amplitudes: \deqn{\sum_{i = 1}^n \mathrm{single}(p_{2n-1}, p_{2n}) \,.} Use
#'   the `higher_states` parameter to restrict the thermal states with priors to
#'   stabilize the fit.
#'
#' @importFrom R6 R6Class
#' 
#' @return
#' See \link{bootstrap.nlsfit}.
#' 
#' @export
new_matrixfit <- function(cf,
                          t1, t2,
                          parlist,
                          sym.vec = rep(1, cf$nrObs),
                          neg.vec = rep('cosh', cf$nrObs),
                          useCov = FALSE,
                          model = "single",
                          boot.fit = TRUE,
                          fit.method = "optim",
                          autoproceed = FALSE,
                          par.guess,
                          every,
                          higher_states = list(val = numeric(0), boot = matrix(nrow = 0, ncol = 0), ampl = numeric(0)),
                          ...) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))
  
  if (model == 'pc') {
    stopifnot(inherits(cf, 'cf_principal_correlator'))
  }
  
  stopifnot(cf$symmetrised == TRUE)
  stopifnot(length(higher_states$val) == ncol(higher_states$boot))
  stopifnot(length(higher_states$val) == length(higher_states$ampl))
  stopifnot(nrow(higher_states$boot) %in% c(0, cf$boot.R))
  
  t1p1 <- t1 + 1
  t2p1 <- t2 + 1
  
  N <- dim(cf$cf)[1]
  Thalfp1 <- cf$Time/2 + 1
  t <- c(0:(cf$Time/2))
  
  ## This is the number of correlators in cf
  if (!is.null(dim(cf$cf)))
    mSize <- dim(cf$cf)[2] / Thalfp1
  else
    mSize <- dim(cf$cf.tsboot$t)[2] / Thalfp1
  
  if (model == 'pc' && mSize != 1) {
    stop('For model pc only a 1x1 matrix is allowed.')
  }
  
  if (missing(parlist)) {
    parlist <- make_parlist(mSize)
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
    warning(mSize, " ", dim(parlist)[2], "\n")
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
  
  if (model == 'single') {
    model_object <- SingleModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize)
  } else if (model == 'two_amplitudes') {
    stopifnot(cf$nrObs == 1)
    model_object <- TwoAmplitudesModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize)
  } else if (model == 'shifted') {
    stopifnot(inherits(cf, 'cf_shifted'))
    model_object <- ShiftedModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize, cf$deltat)
  } else if (model == 'weighted') {
    stopifnot(inherits(cf, 'cf_weighted'))
    model_object <- WeightedModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize, cf$deltat, cf$weight.factor)
  } else if (model == 'pc') {
    stopifnot(cf$nrObs == 1)
    model_object <- TwoStateModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize, cf$gevp_reference_time)
  } else if (model == 'n_particles') {
    stopifnot(cf$nrObs == 1)
    model_object <- NParticleModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize)
  } else if (model == 'single_constant') {
    stopifnot(cf$nrObs == 1)
    model_object <- SingleConstantModel$new(cf$Time, parlist, sym.vec, neg.vec, mSize)
  } else {
    stop(sprintf('Unknown matrixfit model `%s`.', model))
  }
  
  if (missing(par.guess)) {
    par.guess <- model_object$initial_guess(CF$Cor, t1, t2)
  }

  ## perform the bootstrap non-linear least-squares fit (NLS fit):
  args <- list(fn = model_object$prediction,
               gr = model_object$prediction_jacobian,
               par.guess = par.guess,
               y = CF$Cor,
               x = CF$t,
               mask = ii,
               bsamples = cf$cf.tsboot$t,
               use.minpack.lm = fit.method == 'lm',
               error = cf$error_fn,
               cov_fn = cf$cov_fn)
  
  if (useCov) {
    args$CovMatrix <- cf$cov_fn(cf$cf.tsboot$t)
  }
  
  if (length(higher_states$val) > 0) {
    args$priors <- list(
      param = seq(3, by = 2, length.out = length(higher_states$val)),
      p = higher_states$val,
      psamples = higher_states$boot)
    
    for (i in 1:length(higher_states$val)) {
      par.guess <- c(par.guess, higher_states$val[i], higher_states$ampl[i])
    }
    
    args$par.guess <- par.guess
    
    # We know that neither amplitude nor masses can become negative. We
    # therefore enforce these as a lower bound.
    args$lower <- rep(0, times = length(par.guess))
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
#' @return
#' Returns a square, integer-valued matrix.
#' @export
make_parlist <- function (corr_matrix_size) {
  n <- sqrt(corr_matrix_size)
  row1 <- rep(1:n, each = n)
  row2 <- rep(1:n, times = n)
  parlist_matrix <- matrix(cbind(row1, row2), nrow = 2, ncol = n * n, byrow = TRUE)
  
  return(parlist_matrix)
}

#' Create a parameter index matrix for `matrixfit`
#'
#' @param parlist integer array. Parameter list generated with `make_parlist`.
#' @param length_time integer. Number of time slices per correlator.
#' @param summands integer. Number of summands in the fit model that shall be
#' fitted. The signal counts as one summand, each explicit pollution term with
#' independent amplitudes counts as its own summand.
#'
#' @return
#' Returns an array with the parameter indices.
#' 
#' @export
make_parind <- function (parlist, length_time, summands = 1) {
  corr_matrix_size <- ncol(parlist)
  ## parind is the index vector for the matrix elements
  parind <- array(1, dim = c(corr_matrix_size * length_time, 2))
  for (i in 1:(corr_matrix_size)) {
    parind[((i - 1) * length_time + 1):(i * length_time), ] <- t(array(parlist[, i] + 1, dim = c(2, length_time)))
  }
  parind_aux <- parind
  if (summands > 1) {
    for(j in 1:(summands - 1)){
      parind <- cbind(parind, (parind_aux + sqrt(corr_matrix_size) * j))
    }
  }
  return(parind)
}
