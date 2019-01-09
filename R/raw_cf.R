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

#' @title \code{raw_cf} metadata mixin constructor
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

#' @title Original data RAW_val mixin constructor
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

#' @title Gamma method analysis on all time-slices in a 'raw_cf' object
#'
#' @param cf Correlation function container of class 'raw_cf'  
uwerr.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_data'))
  
  dims <- dim(cf$data)
  # prepare some output tensors which are reduced over the first
  # dimension (the measurements)
  na_array <- array(as.numeric(NA), dim=dims[2:length(dims)])

  value <- list(real = na_array, imag = na_array)
  dvalue <- list(real = na_array, imag = na_array)
  ddvalue <- list(real = na_array, imag = na_array)
  tauint <- list(real = na_array, imag = na_array)
  dtauint <- list(real = na_array, imag = na_array)

  # generate an index matrix for this reduced tensor
  args <- list()
  for( d in dim(na_array) ){
    args[[length(args)+1]] <- 1:d
  }
  idx_matrix_uwerr <- as.matrix(do.call(expand.grid, args))

  # the measurements are the fastest running dimension so we can easily
  # perform the Gamma method analysis by stepping through the index matrix

  # get the index matrix for the correlator tensor
  idcs <- idx_matrix.raw_cf(cf)
  stepping <- dims[1]
  nsteps <- as.integer( nrow(idcs) / stepping )
  for( step in 0:(nsteps-1) ){
    istart <- step*stepping + 1
    iend <- istart + stepping - 1

    out_idcs <- idx_matrix_uwerr[step+1,]

    tseries <- cf$data[ idcs[istart:iend,] ]
    uw_real <- try(uwerrprimary(data = Re(tseries)), silent=TRUE)
    uw_imag <- try(uwerrprimary(data = Im(tseries)), silent=TRUE)

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
              dtauint = dtauint,
              idx_matrix = idx_matrix_uwerr)) 
}

#' @title Block average correlation function data
#' @description Block `block_length` sequential measurements of the correlation
#'              function together. This occurs, for example, when multiple
#'              stochastic noise vectors are used per measurement or multiple
#'              source locations. Alternatively, it can also be used to
#'              account for auto-correlations in the data. If the total numbers
#'              of measurements is not divisible by `block_length`, the last
#'              measurements are discarded.
#'
#' @param cf `raw_cf` object
#' @param block_length Integer, number of successive measurements to average over.
#' @return cf `raw_cf` object with the data member reduced in its first dimension
#'        by a factor of `block_length` and restricted (at the end) to the number
#'        of measurements divisible by `block_length`.
block.raw_cf <- function(cf, block_length){
  stopifnot(inherits(cf, 'raw_cf'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  if( block_length == 1 ){
    return(cf)
  }

  idcs <- idx_matrix.raw_cf(cf)
  dims <- dim(cf$data)

  nblocks <- as.integer( dims[1] / block_length )
  nelem <- as.integer( nrow(idcs) / dims[1] )

  # array for blocked data as well as an index matrix for this array
  out <- array(as.complex(NA),
               dim = c(nblocks, cf$nts*cf$nrObs*cf$nrStypes, cf$dim))
  args <- list()
  for( d in dim(out) ){
    args[[length(args)+1]] <- 1:d
  }
  idcs_out <- as.matrix(do.call(expand.grid, args))

  iout <- 1
  for( elem in 0:(nelem-1) ){
    for( block in 0:(nblocks-1) ){ 
      istart <- elem*dims[1] + block*block_length + 1
      iend <- istart + block_length - 1

      # select block of measurements
      idcs_in <- idcs[istart:iend,]
      # place block average into element of 'out' array
      # be mindful of automatic dropping of dimensions!!
      idx_out <- idcs_out[iout:iout,,drop=FALSE]
      out[idx_out] <- mean( cf$data[ idcs_in ] )
      iout <- iout+1
    }
  }
  # replace the data member by the blocked data
  cf$data <- out
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

#' @title add two `raw_cf` objects
#' @return \code{a*cf1$data + b*cf2$data}
add.raw_cf <- function(cf1, cf2, a=1.0, b=1.0) {
  stopifnot(inherits(cf1, 'raw_cf'))
  stopifnot(inherits(cf2, 'raw_cf'))
  stopifnot(inherits(cf1, 'raw_cf_data'))
  stopifnot(inherits(cf2, 'raw_cf_data'))
  stopifnot(inherits(cf1, 'raw_cf_meta'))
  stopifnot(inherits(cf2, 'raw_cf_meta'))
  stopifnot(all(dim(cf1$data) == dim(cf2$data)))
  stopifnot(cf1$Time == cf2$Time)
  stopifnot(cf1$nts == cf2$nts)

  cf <- cf1
  cf$data <- a*cf1$data + b*cf2$data
  return(cf)
}


#' @title add two `raw_cf` objects
#' @return `raw_cf` object with \code{cf$data == cf1$data + cf2$data}
'+.raw_cf' <- function (cf1, cf2) {
  add.raw_cf(cf1, cf2, a = 1.0, b = 1.0)
}

#' @title add two `raw_cf` objects
#' @return `raw_cf` object with \code{cf$data == cf1$data - cf2$data}
'-.raw_cf' <- function(cf1, cf2) {
  add.raw_cf(cf1, cf2, a = 1.0, b = -1.0)
}

#' @title divide two `raw_cf` objects
#' @return `raw_cf` object with \code{cf$data == cf1$data / cf2$data}
'/.raw_cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'raw_cf_meta'))
  stopifnot(inherits(cf2, 'raw_cf_meta'))
  stopifnot(inherits(cf1, 'raw_cf_data'))
  stopifnot(inherits(cf2, 'raw_cf_data'))
  stopifnot(all(dim(cf1$data) == dim(cf2$data)))
  stopifnot(cf1$Time == cf2$Time)
  stopifnot(cf1$nts == cf2$nts)

  cf <- cf1
  cf$data <- cf1$data / cf2$data
  return (cf)
}

#' @title multiply two `raw_cf` objects
#' @return `raw_cf` object with \code{cf$data == cf1$data * cf2$data}
'*.raw_cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'raw_cf_meta'))
  stopifnot(inherits(cf2, 'raw_cf_meta'))
  stopifnot(inherits(cf1, 'raw_cf_data'))
  stopifnot(inherits(cf2, 'raw_cf_data'))
  stopifnot(all(dim(cf1$data) == dim(cf2$data)))
  stopifnot(cf1$Time == cf2$Time)
  stopifnot(cf1$nts == cf2$nts)

  cf <- cf1
  cf$data <- cf1$data * cf2$data
  return(cf)
}

#' @title scale `raw_cf` data
#' @param a Numeric or complex scaling factor, although it could also be
#'          an array of dimensions compatible with \code{cf$data}
#' @return `raw_cf` object with \code{res$data == a*cf$data}
mul.raw_cf <- function(cf, a=1.) {
  stopifnot(inherits(cf, 'raw_cf_data'))
  stopifnot(is.numeric(a) | is.complex(a))

  cf$data <- a*cf$data
  return (cf)
}

#' @title check if an object is of class `raw_cf`
#' @param x object to be checked 
is.raw_cf <- function(x){
  inherits(x, "raw_cf")
}

#' @title check if an obect is of class `raw_cf` and empty otherwise
#' @param .raw_cf object to be checked
is_empty.raw_cf <- function(x){
  .raw_cf <- x
  setequal(class(.raw_cf), class(raw_cf())) & is.null(names(.raw_cf))
}


#' @title Concatenate `raw_cf` correlation function objects
#'
#' @param ... Zero or multiple objects of type `raw_cf`.
c.raw_cf <- function (...) {
  rval <- Reduce(concat.raw_cf, list(...), raw_cf())
  return (rval)
}

#' @title Concatenate two `raw_cf` correlation function objects
#' @description The data of the \code{left} and \code{right} objects
#'              is concatenated along the second array dimension
#'              such that the output contains the tensor slices of \code{right}
#'              after the slices of \code{left}
#' @param left `raw_cf` object to be concatenated with \code{right}
#' @param right `raw_cf` object to be concatenated with \code{left}
concat.raw_cf <- function (left, right) {
  stopifnot(inherits(left, 'raw_cf'))
  stopifnot(inherits(right, 'raw_cf'))

  # In case that one of them does not contain data, the other one is the
  # result. This satisfies the neutral element axiom of a monoid.
  if (is_empty.raw_cf(left)) {
    return (right)
  } else if (is_empty.raw_cf(right)) {
    return (left)
  } else if ( is_empty.raw_cf(left) & is_empty.raw_cf(right) ) {
    stop("concat.raw_cf: both `raw_cf` objects empty, cannot concatenate!")
  }

  stopifnot(inherits(left, 'raw_cf_meta'))
  stopifnot(inherits(right, 'raw_cf_meta'))

  # At this point both `raw_cf` objects given here have original data, therefore we
  # need to concatenate them.

  # A few checks for compatability.
  stopifnot(left$Time == right$Time)
  stopifnot(left$nts == right$nts)
  stopifnot(nrow(left$data) == nrow(right$data))
  stopifnot(left$nrStypes == right$nrStypes)
  stopifnot(all(left$dim == right$dim))

  rval <- raw_cf_meta(nrObs = left$nrObs + right$nrObs,
                      Time = left$Time,
                      nts = left$nts,
                      dim = left$im,
                      nrStypes = left$nrStypes)
  rval <- raw_cf_data(cf = rval,
                      data = abind::abind(left$data, right$data, along=2))
  return (rval)
}

#' @title extract data in format convenient to plot
#' @description When dealing with with tensorial `raw_cf` objects
#'              pre-processing and reshaping is always required to
#'              prepare the data for plotting (or similar). This function
#'              conveniently prepares a named list of prepared data.
#'              The list elements are themselves lists which contain
#'              `val` and `dval` members with the central value and error
#'              of the element in question. These are in turn
#'              arrays of dimension \code{ c( cf$nts, cf$dim ) } and thus
#'              lack the first index compared to \code{cf$data}.
#' @param cf `raw_cf` object with meta-data and data.
#' @param reim String, one of 'real', 'imag' or 'both'. Specifies
#'             whether the real and/or imaginary parts should be extracted.
#' @param tauint Boolean, specifies if the tensor of auto-correlation times
#'               and corresponding errors should be extracted. 
#' @param relerr Boolean, specifies if the return value should also include
#'               estimates of the relative error and its error.
#' @return List of up to six named elements (depending on what was passed for
#'         \code{reim}, \code{tauint}, \code{relerr}) containing the central
#'         values and errors of the real and/or imaginary part of `cf$data`
#'         as well as the corresponding arrays of auto-correlation times and
#'         relative errors. The list elements come in the order
#'         \code{real}, \code{imag}, \code{relerr_real}, \code{relerr_imag},
#'         \code{tauint_real}, \code{tauint_imag} if \code{reim} is `both` and
#'         \code{tauint} and \code{relerr} are \code{TRUE}. The \code{val} and
#'         \code{dval} members of these list elements are arrays of dimension
#'         \cdoe{ c( cf$nts, cf$dim ) } and thus lack the first index compared
#'         to \code{cf$data}.
get_plotdata_raw_cf <- function(cf,
                                reim,
                                tauint,
                                relerr)
{
  stopifnot(any(inherits(cf, c('raw_cf_data'))))
  stopifnot(inherits(cf, 'raw_cf_meta'))
  
  data <- uwerr.raw_cf(cf)

  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' must be one of 'real', 'imag' or 'both'")
  }

  data <- uwerr.raw_cf(cf)
  # construct a return value
  plotdata <- list()
  # loop over the internal dimensions
  for( idim in 0:(prod(cf$dim)-1) ){
    istart <- idim*cf$nrStypes*cf$nrObs*cf$nts + 1
    iend <- istart + cf$nrStypes*cf$nrObs*cf$nts - 1
  
    # subset of the index matrix for the internal dimension given by the  
    iselect <- data$idx_matrix[istart:iend,]
    
    lidx <- length(plotdata)+1
    plotdata[[lidx]] <- list()
    if( reim %in% c('real','both') ){
      # ensure that we append to the end in the right order
      plotdata[[lidx]] <- append(plotdata[[lidx]],
                               list(real = 
                                    list(val = data$value$real[iselect],
                                         dval = data$dvalue$real[iselect])
                                    )
                               )
    }
    if( reim %in% c('imag','both') ){
      plotdata[[lidx]] <- append(plotdata[[lidx]],
                               list(imag = 
                                    list(val = data$value$imag[iselect],
                                         dval = data$dvalue$imag[iselect])
                                    )
                               )
    }
    if( relerr ){
      if( reim %in% c('real','both') ){
        plotdata[[lidx]] <- append(plotdata[[lidx]],
                                 list(relerr_real =
                                      list(val = abs(data$dvalue$real[iselect] / data$value$real[iselect]),
                                           dval = data$ddvalue$real[iselect])
                                      )
                                 )
      }
      if( reim %in% c('imag','both') ){
        plotdata[[lidx]] <- append(plotdata[[lidx]],
                                 list(relerr_imag =
                                      list(val = abs(data$dvalue$imag[iselect] / data$value$imag[iselect]),
                                           dval = data$ddvalue$imag[iselect])
                                      )
                                 )
      }
    }
    if( tauint ){
      if( reim %in% c('real','both') ){
        plotdata[[lidx]] <- append(plotdata[[lidx]],
                                 list(tauint_real =
                                      list(val = data$tauint$real[iselect],
                                           dval = data$dtauint$real[iselect])
                                      )
                                 )
      }
      if( reim %in% c('imag', 'both') ){
        plotdata[[lidx]] <- append(plotdata[[lidx]],
                                 list(tauint_imag =
                                      list(val = data$tauint$imag[iselect],
                                           dval = data$dtauint$imag[iselect])
                                      )
                                 )
      }
    }
  }
  return(plotdata)
}

#' @title plot all correlators in `raw_cf` object
#' @param x Object of class `raw_cf` with data and meta-data.
#' @param reim Character vector, may contain 'real', 'imag' or 'both'. Determines
#'             whether the real and/or imaginary parts of the correlation funtions
#'             should be plotted.
#' @param reim_same Boolean, determines whether the real and imaginary parts, if both
#'                  are to be plotted, will be plotted in the same plot.
#' @param ... Further parameters passed to \link{plotwitherror}.
#' @return Invisibly returns the plotdata, see \link{plotdata.raw_cf}.
plot.raw_cf <- function(x,
                        ...,
                        reim = 'real', 
                        reim_same = FALSE)
{
  cf <- x
  stopifnot(any(inherits(cf, c('raw_cf_data'))))
  stopifnot(inherits(cf, 'raw_cf_meta'))
  
  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' must be one of 'real', 'imag' or 'both'")
  }
  if( reim_same & !(reim == 'both') ){
    stop("'reim_same' can only be true if 'reim' is 'both'")
  }

  plotdata <- get_plotdata_raw_cf(cf, reim=reim, tauint=FALSE, relerr=FALSE) 
  for( lidx in 1:length(plotdata) ){
    if( plot_reim %in% c('real', 'both') & plot_reim_same == FALSE){
      plotwitherror(x = 0:(cf$nts-1),
                    y = plotdata[[lidx]]$real$val,
                    dy = plotdata[[lidx]]$real$dval,
                    ...)
    }
    if( plot_reim %in% c('imag', 'both') & plot_reim_same == FALSE){
      plotwitherror(x = 0:(cf$nts-1),
                    y = plotdata[[lidx]]$imag$val,
                    dy = plotdata[[lidx]]$imag$dval,
                    ...)
    }
    if( plot_reim_same ){
      plotwitherror(x = rep(0:(cf$nts-1), times=2),
                    y = c(plotdata[[lidx]]$real$val, plotdata[[lidx]]$imag$val),
                    dy = c(plotdata[[lidx]]$real$dval, plotdata[[lidx]]$imag$dval),
                    ...)
    }
  }
  return(invisible(plotdata))
}

overview_plot_raw_cf <- function(cf, 
                                 reim = 'real', 
                                 reim_same = FALSE,
                                 relerr = FALSE,
                                 tauint = FALSE,
                                 logplot = '',
                                 title = '')
{
  stopifnot(any(inherits(cf, c('raw_cf_data'))))
  stopifnot(inherits(cf, 'raw_cf_meta'))

  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' must be one of 'real', 'imag' or 'both'")
  }
  if( reim_same & !(reim == 'both') ){
    stop("'reim_same' can only be true if 'reim' is 'both'")
  }

  tmax <- cf$nts-1  
  plotdata <- get_plotdata_raw_cf(cf, reim=reim, relerr=relerr, tauint=tauint)
  ts <- rep(0:tmax, times=cf$nrStypes*cf$nrObs)
  if( reim_same ){
    ts <- rep(ts,times=2)
  }
  step <- 1
  if( reim_same ) step <- 2
  for( lidx in 1:length(plotdata) ){
    for( qidx in seq(1,length(plotdata[[lidx]]),step) ){
      args <- list()
      args$x <- ts
      if( any(names(plotdata[[lidx]][[qidx]]) == "ylim") ){
        args$ylim <- plotdata[[lidx]][[qidx]]$ylim
      }
      args$xlab <- "t/a"
      args$log <- plotdata[[lidx]][[qidx]]$logplot
      args$main <- title
      if( reim_same ){
        args$y <- c(plotdata[[lidx]][[qidx]]$val, plotdata[[lidx]][[qidx+1]]$val)
        args$dy <- c(plotdata[[lidx]][[qidx]]$dval, plotdata[[lidx]][[qidx+1]]$dval)
        args$col <- c(rep("black", cf$nrStypes*cf$nrObs*cf$nts),
                      rep("red", cf$nrStypes*cf$nrObs*cf$nts))
        args$ylab <- plotdata[[lidx]][[qidx]]$same_ylab
      } else {
        args$y <- plotdata[[lidx]][[qidx]]$val
        args$dy <- plotdata[[lidx]][[qidx]]$dval
        args$ylab <- plotdata[[lidx]][[qidx]]$ylab
      }
      do.call(plotwitherror, args)
    }
  }
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
  return(cf)
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
  args <- list()
  for( d in cf$dim ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(expand.grid, args))
}


#' @title Print summary of data contained in `raw_cf` container
#' @param cf `raw_cf` container with data and meta-data
#' @param ... ignored
summary.raw_cf <- function(cf, ...) {
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  cat("T = ", cf$Time, "\n")
  cat("observations = ", dim(cf$data)[1], "\n")
  cat("Nr Stypes = ", cf$nrStypes, "\n")
  cat("Nr Obs    = ", cf$nrObs, "\n")
  cat("dim       = ", cf$dim, "\n")

  uw <- uwerr.raw_cf(cf)

  idcs <- uw$idx_matrix
  out <- NULL
  for( i in 1:nrow(idcs) ){
    out <- rbind(out,
                 data.frame(t = idcs[i,1]-1, 
                            value_real = uw$value$real[idcs[i,drop=FALSE]], 
                            value_imag = uw$value$imag[idcs[i,drop=FALSE]],
                            dvalue_real = uw$dvalue$real[idcs[i,drop=FALSE]],
                            dvalue_imag = uw$dvalue$imag[idcs[i,drop=FALSE]])
                 )
  }
  rownames(out) <- NULL
  print(out)
  return(invisible(out))
}

#' @title Print summary of data contained in `raw_cf` container
#' @param cf `raw_cf` container with data and meta-data
#' @param ... ignored
print.raw_cf <- function(cf, ...) {
  summary(cf, ...)
}
