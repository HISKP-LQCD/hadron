#' @title Container for raw correlation functions
#'
#' @description
#' This function \code{raw_cf()} creates containers for raw correlation functions
#' of class \code{raw_cf}. This class is particularly designed to deal with
#' complex and matrix-valued correlation functions emerging in statistical
#' mechanics and quantum field theory simulations. 
#' Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.raw_cf} and \link{as.raw_cf}.
#'
#' @family raw_cf constructors
#'
#' @export
raw_cf <- function () {
  cf <- list()
  class(cf) <- append(class(cf), 'raw_cf')
  return (cf)
}

#' @title RAW_CF metadata mixin constructor
#'
#' @param nrObs Integer, number of different observables assembled in the data field of this container. 
#' @param Time Integer, full time extent.
#' @param nts Integer, number of time separations actually stored in the data field.
#' @param nrStypes Integer, number of smearing types.
#' @param dim Integer vector of "internal" dimensions for matrix-valued correlation functions. 
#'            For a scalar correlation, this should be specified as \code{c(1,1)}.
#'
#' @family raw_cf constructors
#'
#' @export
raw_cf_meta <- function (cf = raw_cf(), nrObs = 1, Time = NA, nrStypes = 1, dim=c(1,1), nts = Time ) {
  stopifnot(inherits(cf, 'raw_cf'))
  
  cf$nrObs <- nrObs
  cf$Time <- Time
  cf$nts <- nts
  cf$nrStypes <- nrStypes
  cf$dim <- dim

  class(cf) <- append(class(cf), 'raw_cf_meta')
  return (cf)
}

#' @title Original data RAW_CF mixin constructor
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

#' @title Add Gamma-method error analysis to existing `raw_cf` container
#' 
#' @param cf Correlation function container of class 'raw_cf'
#' @family raw_cf constructors
#' 
#' @export
raw_cf_uwerr <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_data'))
  stopifnot(inherits(cf, 'raw_cf_meta'))
  cf$uwerr <- uwerr.raw_cf(cf)
  class(cf) <- append(class(cf), 'raw_cf_uwerr')
  return(cf)
}

#' @title Gamma method analysis on all time-slices in a 'raw_cf' object
#'
#' @param cf Correlation function container of class 'raw_cf'  
uwerr.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_data'))
  
  dims <- dim(cf$data)
  # prepare some output tensors which are reduced over the first
  # dimension (the measurements)
  na_array <- array(as.numeric(NA), dim=c(1,dims[2:length(dims)]))

  value <- list(real = na_array, imag = na_array)
  dvalue <- list(real = na_array, imag = na_array)
  ddvalue <- list(real = na_array, imag = na_array)
  tauint <- list(real = na_array, imag = na_array)
  dtauint <- list(real = na_array, imag = na_array)

  # the measurements are the fastest running dimension so we can easily
  # perform the Gamma method analysis by stepping through the index matrix

  # get the index matrix for the correlator tensor
  idcs <- idx_matrix.raw_cf(cf)
  print(idcs)
  stepping <- dims[1]
  nsteps <- as.integer( nrow(idcs) / stepping )
  for( step in 0:(nsteps-1) ){
    istart <- step*stepping + 1
    iend <- istart + stepping - 1

    out_idcs <- idcs[istart, 2:ncol(idcs)]

    ts <- cf$data[ idcs[istart:iend,] ]
    uw_real <- try(uwerrprimary(data = Re(ts)), silent=TRUE)
    uw_imag <- try(uwerrprimary(data = Im(ts)), silent=TRUE)

    if(!inherits(uw_real,'try-error')){
        value$real[out_idcs] <- uw_real$value
       dvalue$real[out_idcs] <- uw_real$dvalue
      ddvalue$real[out_idcs] <- uw_real$ddvalue
       tauint$real[out_idcs] <- uw_real$tauint
      dtauint$real[out_idcs] <- uw_real$dtauint
    }

    if(!inherits(uw_imag, 'try-error')){
        value$imag[out_idcs] <- uw_imag$value
       dvalue$imag[out_idcs] <- uw_imag$dvalue
      ddvalue$imag[out_idcs] <- uw_imag$ddvalue
       tauint$imag[out_idcs] <- uw_imag$tauint
      dtauint$imag[out_idcs] <- uw_imag$dtauint
    }
  }
  return(list(value = value,
              dvalue = dvalue,
              ddvalue = ddvalue,
              tauint = tauint,
              dtauint = dtauint)) 
}

#addConfIndex2cf <- function(cf, conf.index) {
#  if(is.null(cf$conf.index)) {
#    cf$conf.index <- conf.index
#  }
#  return(cf)
#}

#' @title Block average correlation function data
#' @description Block `block_length` measurements of the correlation
#'              function together. This occurs, for example, when multiple
#'              stochastic noise vectors are used per measurement or multiple
#'              source locations. Alternatively, it can also be used to
#'              account for auto-correlations in the data.
#'
#' @param cf `raw_cf` object
#' @param block_length Integer, number of successive measurements to average over.
block.raw_cf <- function(cf, block_length){
  stopifnot(inherits(cf, 'raw_cf'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  idcs <- idx_matrix.raw_cf(cf)
  dims <- dim(cf$data)

  nelem <- as.integer( nrow(idcs) / dims[1] )
  nblocks <- as.integer( dims[1] / block_length )

  for( elem in 0:(nelem-1) ){
    for( block in 0:(nblocks-1) ){ 
      istart <- elem*dims[1] + block*block_length + 1
      iend <- istart + block_length - 1
      # place block average into first slot
      cf$data[ idcs[istart,] ] <- mean( cf$data[ idcs[istart:iend,] ] )
    }
  }
  # subset the block averages
  iselect <- seq(from = 1, to = dims[1], by = block_length)
  args <- list()
  args[[1]] <- iselect
  args[[2]] <- 1:(cf$nts*cf$nrObs*cf$nrStypes)
  for( d in cf$dim ){
    args[[length(args)+1]] <- 1:d
  }
  idcs <- as.matrix(do.call(expand.grid, args))
  # replace the data member by the blocked data
  cf$data <- array(cf$data[idcs],
                   dim=c(length(iselect), cf$nts*cf$nrObs*cf$nrStypes, cf$dim))
  return(cf)
}

#' @title Extend statistics of an existing `raw_cf` container
#' 
#' @param cf1 `raw_cf` container with or without 'data' and 'meta' mixins
#' @param cf2 `raw_cf` container with or without 'data' and 'meta' mixins
#'
#' @details 
#' When either of `cf1` or `cf2` does not contain any data,
#' the other object is returned. (allows empty `raw_cf` to be extended).
#' If the dimensions (except for the measurements) of the `data` fields of the
#' two containers match, they are concatenated along the measurement dimension.
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

  stopifnot(inherits(cf1, 'raw_cf_data'))
  stopifnot(inherits(cf2, 'raw_cf_data'))

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

plot.raw_cf <- function(cf, plot_im = FALSE, plot_im_same = FALSE, return_only = FALSE, ...) {
  stopifnot(any(inherits(cf, c('raw_cf_data'))))
  stopifnot(inherits(cf, 'raw_cf_meta'))

  if( !inherits(cf, 'raw_cf_uwerr') ){
    data <- uwerr.raw_cf(cf)
  } else {
    data <- cf$uwerr
  }

  dims <- dim(data$value$real)
  args <- list()
  for( d in dims ){
    args[[length(args)+1]] <- 1:d
  }
  idcs <- as.matrix(do.call(expand.grid, args))

  tmax <- cf$nts-1
 
  plotdf <- list()

  for( i in 0:(prod(cf$dim)-1) ){
    istart <- i*cf$nrStypes*cf$nrObs*cf$nts + 1
    iend <- istart + cf$nrStypes*cf$nrObs*cf$nts - 1
   
    iselect <- idcs[istart:iend,]

    #clr_real <- rep("black", times=cf$nrStypes*cf$nrObs*cf$nts)
    #clr_real[ which( Re(data$value[iselect]) < 0) ] <- "red"

    #clr_imag <- rep("blue", times=cf$nrStypes*cf$nrObs*cf$nts)
    #clr_imag[ which( Im(data$value[iselect]) < 0) ] <- "purple"

    plotdf[[length(plotdf)+1]] <-
      list(t = rep(x=0:tmax, times=cf$nrStypes*cf$nrObs),
           real = data.frame(CF = data$value$real[iselect],
                             Err = data$dvalue$real[iselect] ),
           imag = data.frame(CF = data$value$imag[iselect],
                             Err = data$dvalue$imag[iselect] )
           )
    print(str(plotdf[[length(plotdf)]]))

    if( ! return_only ){
      if( plot_im & plot_im_same ){
        plotwitherror(x = rep(x=plotdf[[length(plotdf)]]$t, times=2), 
                      y = c(plotdf[[length(plotdf)]]$real$CF, plotdf[[length(plotdf)]]$imag$CF), 
                      dy = c(plotdf[[length(plotdf)]]$real$Err, plotdf[[length(plotdf)]]$imag$Err),
                      col = c(rep("black", cf$nrStypes*cf$nrObs*cf$nts),
                              rep("red", cf$nrStypes*cf$nrObs*cf$nts)),
                      ...)
      } else {
        print(length(plotdf[[length(plotdf)]]$t))
        print(length(plotdf[[length(plotdf)]]$real$CF))
        plotwitherror(x = plotdf[[length(plotdf)]]$t, 
                      y = plotdf[[length(plotdf)]]$real$CF, 
                      dy = plotdf[[length(plotdf)]]$real$Err,
                      ...)
        if( plot_im ){
          plotwitherror(x = plotdf[[length(plotdf)]]$t, 
                        y = plotdf[[length(plotdf)]]$real$CF, 
                        dy = plotdf[[length(plotdf)]]$real$Err,
                        ...)
        }
      }
    }
  }
  return(invisible(plotdf))
}


#' @title shift a \code{raw_cf} correlation function by 'places' time-slices
#' @param cf \code{raw_cf} container
#' @param places Integer, number of time slices that the correlation function
#'               should be shifted by. Can be positive or negative.
#' @details
#' The correlation funtion \eqn{C(t)} is shifted in time to produce:
#'   \deqn{C'(t) = C(t+places)}
#' using periodic boundary conditions in time.
#' @export
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

#' @title Construct the tensor index set for the entire raw correlator
#' @param cf 'raw_cf' container with data and meta-data
idx_matrix.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))
 
  dims <- dim(cf$data)  
  args <- list()
  for( d in dims ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(expand.grid, args))
}

#' @title Construct tensor index set for the internal degrees of freedom
#' @param cf `raw_cf` container
int_idx_matrix.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))
  args <- list()
  for( d in cf$dim ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(expand.grid, args))
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
