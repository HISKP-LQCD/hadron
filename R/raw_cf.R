#' Correlation function container
#'
#' This function `raw_cf()` creates containers for raw correlation functions
#' of class `raw_cf`. This class is particularly designed to deal with
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.raw_cf} and \link{as.raw_cf}.
#'
#' @details
#'
#' And last but not least, these are the fields that are used somewhere in the library but we have not figured out which mixin these should belong to:
#'
#' - `conf.index`: TODO
#' - `N`: Integer, number of measurements.
#' - `blockind`: TODO
#' - `jack.boot.se`: TODO
#'
#' @family raw_cf constructors
#'
#' @export
raw_cf <- function () {
  cf <- list()
  class(cf) <- append(class(cf), 'raw_cf')
  return (cf)
}

#' RAW_CF metadata mixin constructor
#'
#' @param nrObs Integer, number of different observables assembled in the data field of this container. One can use \link{c.raw_cf} to add multiple observables into one container.
#' @param Time Integer, full time extent.
#' @param nrStypes Integer, number of smearing types.
#' @param dim Integer vector of "internal" dimensions for matrix-valued correlation functions.
#'
#' @family raw_cf constructors
#'
#' @export
raw_cf_meta <- function (cf = raw_cf(), nrObs = 1, Time = NA, nrStypes = 1, dim=c(1) ) {
  stopifnot(inherits(cf, 'raw_cf'))
  
  cf$nrObs <- nrObs
  cf$Time <- Time
  cf$nrStypes <- nrStypes
  cf$dim <- dim

  class(cf) <- append(class(cf), 'raw_cf_meta')
  return (cf)
}

#' Original data RAW_CF mixin constructor
#'
#' @param cf `raw_cf` object to extend.
#' @param data Numeric or complex array, original data for all observables and measurements. 
#'             This should have dimensions c(Nmeas,cf$Time*cf$nrObs*cf$nrStypes,cf$dim).
#'             Having the internal dimensions innermost is not as efficient, but it allows
#'             different transformations to be applied to different observables in the same
#'             container more easily.
#'
#' @family raw_cf constructors
#'
#' @export
raw_cf_data <- function (cf, data) {
  stopifnot(inherits(cf, 'raw_cf'))
  stopifnot(inherits(cf, 'raw_cf_meta'))
  dims <- dim(data)
  if( any( cf$dim != dims[3:length(dims)] ) ){
    stop("Dimensions of 'data' and 'raw_cf$dim' container meta-data must match")
  } 

  # raw correlation functions are expected to be complex-valued 
  if( mode(data) != "complex" ){
    cf$data <- as.complex(data)
  } else {
    cf$data <- data
  }

  class(cf) <- append(class(cf), 'raw_cf_data')
  return (cf)
}

# Gamma method analysis on all time-slices in a 'raw_cf' object
uwerr.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_data'))
  
  dims <- dim(cf$data)
  
  # get the index set for the correlator tensor
  idcs <- index_set.raw_cf(cf)
  stepping <- dims[1]
  nsteps <- as.integer( nrow(idcs) / stepping )

  # prepare some output tensors which are reduced over the first
  # dimension (the measurements)
  value <- array(as.complex(NA), dim=dims[2:length(dims)]) 
  dvalue <- array(as.complex(NA), dim=dims[2:length(dims)]) 
  ddvalue <- array(as.complex(NA), dim=dims[2:length(dims)]) 
  tauint_real <- array(as.numeric(NA), dim=dims[2:length(dims)])
  tauint_imag <- array(as.numeric(NA), dim=dims[2:length(dims)])
  dtauint_real <- array(as.numeric(NA), dim=dims[2:length(dims)])
  dtauint_imag <- array(as.numeric(NA), dim=dims[2:length(dims)])

  for( step in 0:(nsteps-1) ){
    istart <- step*stepping + 1
    iend <- istart + stepping - 1

    out_idcs <- idcs[istart, 2:ncol(idcs)]

    ts <- cf$data[ idcs[istart:iend,] ]
    uw_real <- try(uwerrprimary(data = Re(ts)), silent=TRUE)
    uw_imag <- try(uwerrprimary(data = Im(ts)), silent=TRUE)
    if( (!any(class(uw_real) == 'try-error')) &
        (!any(class(uw_imag) == 'try-error')) ){
      value[out_idcs] <- complex(real = uw_real$value,
                                 imaginary = uw_imag$value)
      dvalue[out_idcs] <- complex(real = uw_real$dvalue,
                                  imaginary = uw_imag$dvalue)
      ddvalue[out_idcs] <- complex(real = uw_real$ddvalue,
                                   imaginary = uw_imag$ddvalue)

      tauint_real[out_idcs] <- uw_real$tauint
      dtauint_real[out_idcs] <- uw_real$dtauint

      tauint_imag[out_idcs] <- uw_imag$tauint
      dtauint_imag[out_idcs] <- uw_imag$dtauint
    } 
  }
  return(list(t = array( rep(0:(cf$Time-1), times=nsteps), dim=dims[2:length(dims)] ),
              value = value,
              dvalue = dvalue,
              ddvalue = ddvalue,
              tauint_real = tauint_real,
              tauint_imag = tauint_imag,
              dtauint_real = dtauint_real,
              dtauint_imag = dtauint_imag)) 
}

#addConfIndex2cf <- function(cf, conf.index) {
#  if(is.null(cf$conf.index)) {
#    cf$conf.index <- conf.index
#  }
#  return(cf)
#}

block.raw_cf <- function(cf){
  stopifnot(inherits(cf1, 'raw_cf'))
  stopifnot(inherits(cf2, 'raw_cf'))
}

addStat.raw_cf <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'raw_cf'))
  stopifnot(inherits(cf2, 'raw_cf'))

  if (is_empty.raw_cf(cf1)) {
    return (invisible(cf2))
  }
  if (is_empty.raw_cf(cf2)) {
    return (invisible(cf1))
  }

  stopifnot(inherits(cf1, 'raw_cf_meta'))
  stopifnot(inherits(cf2, 'raw_cf_meta'))

  stopifnot(cf1$Time == cf2$Time)
  stopifnot(cf1$nrObs == cf2$nrObs )
  stopifnot(cf1$nrStypes == cf2$nrStypes)
  stopifnot( all(dim(cf1$data)[2:length(dim(cf2$data))] == dim(cf2$data)[2:length(dim(cf2$data))] ) )

  cf <- cf1

  cf$data <- abind(cf1$data, cf2$data, along=1)

  return (invisible(cf))
}

## this is intended for instance for adding diconnected diagrams to connected ones
#add.cf <- function(cf1, cf2, a=1.0, b=1.0) {
#  stopifnot(inherits(cf1, 'cf'))
#  stopifnot(inherits(cf2, 'cf'))
#  stopifnot(inherits(cf1, 'cf_orig'))
#  stopifnot(inherits(cf2, 'cf_orig'))
#  stopifnot(inherits(cf1, 'cf_meta'))
#  stopifnot(inherits(cf2, 'cf_meta'))
#  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
#  stopifnot(cf1$Time == cf2$Time)
#
#  cf <- cf1
#  cf$cf <- a*cf1$cf + b*cf2$cf
#  cf <- invalidate.samples.cf(cf)
#  return(cf)
#}
#
#'+.cf' <- function (cf1, cf2) {
#  add.cf(cf1, cf2, a = 1.0, b = 1.0)
#}
#
#'-.cf' <- function(cf1, cf2) {
#  add.cf(cf1, cf2, a = 1.0, b = -1.0)
#}
#
#'/.cf' <- function(cf1, cf2) {
#  stopifnot(inherits(cf1, 'cf_meta'))
#  stopifnot(inherits(cf2, 'cf_meta'))
#  stopifnot(inherits(cf1, 'cf_orig'))
#  stopifnot(inherits(cf2, 'cf_orig'))
#  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
#  stopifnot(cf1$Time == cf2$Time)
#
#  cf <- cf1
#  cf$cf <- cf1$cf / cf2$cf
#  cf <- invalidate.samples.cf(cf)
#  return (cf)
#}
#
#mul.cf <- function(cf, a=1.) {
#  stopifnot(inherits(cf, 'cf_orig'))
#  stopifnot(is.numeric(a))
#
#  cf$cf <- a*cf$cf
#  cf <- invalidate.samples.cf(cf)
#  return (cf)
#}
#
#extractSingleCor.cf <- function(cf, id=c(1)) {
#  stopifnot(inherits(cf, 'cf_meta'))
#  stopifnot(inherits(cf, 'cf_orig'))
#
#  ii <- c()
#  for(i in c(1:length(id))) {
#    ii <- c(ii, c(1:(cf$Time/2+1)) + (id[i]-1)*(cf$Time/2+1))
#  }
#
#  # TODO: This should be done using constructors.
#  cf$cf <- cf$cf[,ii]
#
#  if (inherits(cf, 'cf_boot')) {
#    cf$cf0 <- cf$cf0[ii]
#    cf$tsboot.se <- cf$tsboot.se[ii]
#    cf$cf.tsboot$t0 <- cf$cf.tsboot$t0[ii]
#    cf$cf.tsboot$t <- cf$cf.tsboot$t[,ii]
#    cf$cf.tsboot$data <- cf$cf.tsboot$data[,ii]
#  }
#  cf$nrObs <- 1
#  cf$nsStypes <- 1
#  return (cf)
#}
#
#is.raw_cf <- function(x){
#  inherits(x, "raw_cf")
#}
#

is_empty.raw_cf <- function(.raw_cf){
  setequal(class(.raw_cf), class(raw_cf())) & is.null(names(.raw_cf))
}

#
##' Concatenate correlation function objects
##'
##' @param ... Zero or multiple objects of type `raw_cf`.
#c.raw_cf <- function (...) {
#  rval <- Reduce(concat.raw_cf, list(...), raw_cf())
#  return (invisible(rval))
#}
#
##' Concatenate two correlation function objects
#concat.raw_cf <- function (left, right) {
#  stopifnot(inherits(left, 'raw_cf'))
#  stopifnot(inherits(right, 'raw_cf'))
#
#  # In case that one of them does not contain data, the other one is the
#  # result. This satisfies the neutral element axiom of a monoid.
#  if (is_empty.cf(left)) {
#    return (right)
#  }
#  if (is_empty.cf(right)) {
#    return (left)
#  }
#
#  stopifnot(inherits(left, 'cf_meta'))
#  stopifnot(inherits(right, 'cf_meta'))
#
#  # At this point both `cf` objects given here have original data, therefore we
#  # need to concatenate them.
#
#  # A few checks for compatability.
#  stopifnot(left$Time == right$Time)
#  stopifnot(nrow(left$cf) == nrow(right$cf))
#  stopifnot(left$symmetrised == right$symmetrised)
#  stopifnot(left$nrStypes == right$nrStypes)
#
#  rval <- cf_meta(nrObs = left$nrObs + right$nrObs,
#                  Time = left$Time,
#                  nrStypes = left$nrStypes,
#                  symmetrised = left$symmetrised)
#  rval <- cf_orig(.cf = rval,
#                  cf = cbind(left$cf, right$cf),
#                  icf = cbind(left$icf, right$icf))
#  return (invisible(rval))
#}
#
#plot.cf <- function(cf, neg.vec = rep(1, times = length(cf$cf0)), rep = FALSE, ...) {
#  stopifnot(any(inherits(cf, c('cf_orig', 'cf_boot', 'cf_jackknife'))))
#  stopifnot(inherits(cf, 'cf_meta'))
#
#  if (inherits(cf, 'cf_boot')) {
#    val <- cf$cf0
#    err <- cf$tsboot.se
#  } else if (inherits(cf, 'cf_jackknife')) {
#    val <- cf$cf0
#    err <- cf$jackknife.se
#  } else {
#    stop('A correlation function must be bootstrapped before it can be plotted.')
#  }
#
#  if(!cf$symmetrised){
#    tmax <- cf$Time - 1
#  } else {
#    tmax <- cf$Time / 2
#  }
#
#  df <- data.frame(t = rep(c(0:tmax), times = length(val)/(tmax+1)),
#                   CF = val,
#                   Err = err)
#
#  plotwitherror(x = df$t, y = neg.vec * df$CF, dy = df$Err, rep = rep, ...)
#
#  return(invisible(df))
#}
#

# shift a correlation function by 'places' time-slices
#   C'(t) = C(t+places)
# where places can be positive or negative as required
# this will of course mix smearings and observables
# and must be taken into account externally by
# invalidating the affected time-slices
shift.raw_cf <- function(cf, places) {
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  if(places == 0){
    return(invisible(cf))
  }

  dims <- dim(cf$data)

  for( oidx in 0:(cf$nrObs-1) ){
    for( sidx in 0:(cf$nrStypes-1) ){
      # the 'time' indices of the observables and smearing types
      istart <- cf$Time*cf$nrStypes*oidx + cf$Time*sidx + 1
      iend <- istart + cf$Time - 1

      # construct an argument list for do.call below
      args <- list()
      args[[1]] <- 1:dims[1]
      args[[2]] <- istart:iend
      for( d in cf$dim ){
        args[[length(args)+1]] <- 1:d
      }
      # construct the tensor index set for the output
      out_dof <- as.matrix(do.call(expand.grid, args))
      
      if( places < 0 ){
        ishift <- c( (iend - abs(places) + 1):iend,
                     (istart:(iend-abs(places))) )
      } else {
        ishift <- c( (istart+places):iend,
                      istart:(istart+places-1) )
      }
      args[[2]] <- ishift
      # construct the tensor index set for the input
      in_dof <- as.matrix(do.call(expand.grid, args))
      
      # shift the correlator tensor
      cf$data[out_dof] <- cf$data[in_dof]
    }
  }
  return(invisible(cf))
}

#' construct the tensor index set for the entire raw correlator
#'
#' @param cf 'raw_cf' container with data and meta-data
index_set.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))
 
  dims <- dim(cf$data)  
  args <- list()
  for( d in dims ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(expand.grid, args))
}

#' construct tensor index set for the internal degrees of freedom
int_index_set.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))
  args <- list()
  for( d in cf$dim ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(exapnd.grid, args))
}

#
#summary.cf <- function(cf, ...) {
#  stopifnot(inherits(cf, 'cf_meta'))
#
#  cat("T = ", cf$Time, "\n")
#  cat("observations = ", dim(cf$cf)[1], "\n")
#  cat("Nr Stypes = ", cf$nrStypes, "\n")
#  cat("Nr Obs    = ", cf$nrObs, "\n")
#
#  if (inherits(cf, 'cf_boot')) {
#    cat("R = ", cf$boot.R, "\n")
#
#    if(!cf$symmetrised){
#      tmax <- cf$Time-1
#    } else {
#      tmax <- cf$Time/2
#    }
#    cat("l = ", cf$boot.l, "\n")
#    out <- data.frame(t=c(0:tmax), C=cf$cf0)
#    cat("sim = ", cf$sim, "\n")
#
#    out <- cbind(out, tsboot.se=cf$tsboot.se)
#  }
#
#
#  if (inherits(cf, 'cf_jackknife')) {
#    out <- cbind(out, jackknife.se=cf$jackknife.se)
#    out <- cbind(out, jab.se=cf$jack.boot.se)
#  }
#
#  if(exists("out")) {
#    print(out)
#  }
#}
#
#print.cf <- function(cf, ...) {
#  summary(cf, ...)
#}
