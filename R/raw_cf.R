#' @title Container for raw correlation functions
#'
#' @description
#' This function \code{raw_cf()} creates containers for raw correlation functions
#' of class \code{raw_cf}. This class is particularly designed to deal with
#' complex and matrix-valued correlation functions emerging in statistical
#' mechanics and quantum field theory simulations. 
#' Arithmetic operations are defined for this class and utility functions
#' such as \code{is.raw_cf} and \code{is_empty.raw_cf}.
#'
#' @family raw_cf constructors
#'
#' @return
#' An object of S3 class `raw_cf`.
#' 
#' @export
raw_cf <- function () {
  cf <- list()
  class(cf) <- append(class(cf), 'raw_cf')
  return (cf)
}

#' @title \code{raw_cf} metadata mixin constructor
#'
#' @param cf initial \link{raw_cf} object
#' @param nrObs Integer, number of different observables assembled in the data field of this container. 
#' @param Time Integer, full time extent.
#' @param nts Integer, number of time separations actually stored in the data field.
#' @param nrStypes Integer, number of smearing types.
#' @param dim Integer vector of "internal" dimensions for matrix-valued correlation functions. 
#'            For a scalar correlation, this should be specified as \code{c(1,1)}.
#'
#' @family raw_cf constructors
#'
#' @return
#' An object of S3 class `raw_cf` with metadat mixing added.
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

#' @title Original data mixin constructor for `raw_cf`
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
#' @return
#' An object of S3 class `raw_cf` with original data mixin added.
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

#' @title Extract a particular internal component of a 'raw_cf' into a 'cf'
#' @param x 'raw_cf' container with 'raw_cf_data' and 'raw_cf_meta'
#' @param component Integer vector of the same length as the internal dimension
#'                  of the 'raw_cf' specifying which component should be extracted.
#' @return 'cf' object
#' @export
raw_cf_to_cf <- function(x, component){
  stopifnot(inherits(x, 'raw_cf_meta'))
  stopifnot(inherits(x, 'raw_cf_data'))

  if( length(component) != length(x$dim) ){
    stop("The component specifier must be of the same length as the vector of internal dimensions of the 'raw_cf'!")
  }

  if( x$Time != x$nts ){
    stop("The 'nts' dimension must be the same as the 'Time' dimension to convert to 'cf'")
  }

  idcs <- idx_matrix.raw_cf(x, component)

  dims <- dim(x$data)

  subset <- array(x$data[idcs], dim=c(dims[1],dims[2]))

  cf <- cf_meta(Time = x$Time,
                nrObs = x$nrObs,
                nrStypes = x$nrStypes,
                symmetrised = FALSE)

  cf <- cf_orig(cf,
                cf = Re(subset),
                icf = Im(subset))
  return(cf)
}

#' @title Gamma method analysis on all time-slices in a 'raw_cf' object
#'
#' @param cf Correlation function container of class 'raw_cf'
#' @return The return value is a list with elements
#'         \describe{
#'           \item{\code{value}}{central value}
#'           \item{\code{dvalue}}{statistical error}
#'           \item{\code{ddvalue}}{error of the statistical error}
#'           \item{\code{tauint}}{auto-correlation time estimate}
#'           \item{\code{dtauint}}{error of auto-correlation time estimate}
#'         }
#'         Each of these is in turn an array of dimension \code{ c( cf$nts, cf$dim ) } and
#'         hance lacks the first dimension index compared for \code{cf$data}.
#'
#' @export
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

    # this is really annoying. When we subset a single row of the index matrix
    # it's turned into a column vector
    out_idcs <- t(idx_matrix_uwerr[step+1,])

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
#'              account for auto-correlations in the data. If the total number
#'              of measurements is not divisible by `block_length`, the last
#'              measurements are discarded.
#'
#' @param cf `raw_cf` object
#' @param block_length Integer, number of successive measurements to average over.
#' @return cf `raw_cf` object with the data member reduced in its first dimension
#'        by a factor of `block_length` and restricted (at the end) to the number
#'        of measurements divisible by `block_length`.
#'
#' @export
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
      # be mindful of conversion to column vector
      idx_out <- t(idcs_out[iout,])
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
#'
#' @return
#' An object of S3 class `raw_cf` identical to the input object but with extended statistics.
#' 
#' @export
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

  cf$data <- abind::abind(cf1$data, cf2$data, along=1)

  return (invisible(cf))
}

#' @title add two `raw_cf` objects
#' @param cf1 first 'raw_cf' container with data and meta-data
#' @param cf2 second 'raw_cf' container with data and meta-data
#' @param a Numeric or complex, scaling factor applied to \code{cf1}.
#' @param b Numeric or complex, scaling factor applied to \code{cf2}.
#' @return \code{a*cf1$data + b*cf2$data}
#'
#' @export
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
#' @param cf1 first 'raw_cf' container to be added
#' @param cf2 second 'raw_cf' container to be added
#' @return `raw_cf` object with \code{cf$data == cf1$data + cf2$data}
#'
#' @export
'+.raw_cf' <- function (cf1, cf2) {
  add.raw_cf(cf1, cf2, a = 1.0, b = 1.0)
}

#' @title add two `raw_cf` objects
#' @param cf1 first 'raw_cf' container to be subtracted
#' @param cf2 second 'raw_cf' container to be subtracted
#' @return `raw_cf` object with \code{cf$data == cf1$data - cf2$data}
#'
#' @export
'-.raw_cf' <- function(cf1, cf2) {
  add.raw_cf(cf1, cf2, a = 1.0, b = -1.0)
}

#' @title divide two `raw_cf` objects
#' @param cf1 'raw_cf' container with data and meta-data to be the dividend
#' @param cf2 'raw_cf' container with data and meta-data to be the divisor
#' @return `raw_cf` object with \code{cf$data == cf1$data / cf2$data}
#'
#' @export
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
#' @param cf1 first 'raw_cf' container with data and meta-data to be multiplied
#' @param cf2 second 'raw_cf' container with data and meta-data to be multiplied
#' @return `raw_cf` object with \code{cf$data == cf1$data * cf2$data}
#'
#' @export
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

#' @title Take the complex conjugate of a `raw_cf` object
#' @param cf `raw_cf` cotnainer with data
#' @return `raw_cf`
conj_raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_data'))
  cf$data <- Conj(cf$data)
  return(cf)
}

#' @title scale `raw_cf` data
#' @param cf 'raw_cf' container with data to be scaled by the factor \code{a}
#' @param a Numeric or complex scaling factor, although it could also be
#'          an array of dimensions compatible with \code{cf$data}
#' @return `raw_cf` object with \code{res$data == a*cf$data}
#'
#' @export
mul.raw_cf <- function(cf, a=1.) {
  stopifnot(inherits(cf, 'raw_cf_data'))
  stopifnot(is.numeric(a) | is.complex(a))

  cf$data <- a*cf$data
  return (cf)
}

#' @title check if an object is of class `raw_cf`
#' @param x object to be checked 
#'
#' @return
#' Returns TRUE if `x` is an object of class `raw_cf`, FALSE otherwise.
#' 
#' @export
is.raw_cf <- function(x){
  inherits(x, "raw_cf")
}

#' @title check if an obect is of class `raw_cf` and empty otherwise
#' @param x object to be checked
#'
#' @return
#' Returns TRUE if `x` is an empty object of class `raw_cf`, FALSE otherwise.
#' 
#' @export
is_empty.raw_cf <- function(x){
  .raw_cf <- x
  setequal(class(.raw_cf), class(raw_cf())) & is.null(names(.raw_cf))
}


#' @title Concatenate `raw_cf` correlation function objects
#'
#' @param ... Zero or multiple objects of type `raw_cf`.
#'
#' @return
#' Returns an object of S3 class `raw_cf`, the concatenation of the
#' input objects.
#' 
#' @export
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
#'
#' @return
#' Returns an object of S3 class `raw_cf`, the concatenation of the
#' two input objects.
#'
#' @export
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

#' @title extract data from 'raw_cf' in format convenient to plot
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
#'         \code{ c( cf$nts, cf$dim ) } and thus lack the first index compared
#'         to \code{cf$data}.
#'
#' @export
get_plotdata_raw_cf <- function(cf,
                                reim,
                                tauint,
                                relerr)
{
  stopifnot(inherits(cf, 'raw_cf_data'))
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
                                           dval = abs(data$ddvalue$real[iselect] / data$value$real[iselect]))
                                      )
                                 )
      }
      if( reim %in% c('imag','both') ){
        plotdata[[lidx]] <- append(plotdata[[lidx]],
                                 list(relerr_imag =
                                      list(val = abs(data$dvalue$imag[iselect] / data$value$imag[iselect]),
                                           dval = abs(data$ddvalue$imag[iselect] / data$value$imag[iselect]))
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
  class(plotdata) <- append(class(plotdata), "plotdata_raw_cf")
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
#' @return Invisibly returns the plotdata, see \link{get_plotdata_raw_cf}.
#'
#' @export
plot.raw_cf <- function(x,
                        ...,
                        reim = 'real', 
                        reim_same = FALSE)
{
  cf <- x
  stopifnot(inherits(cf, 'raw_cf_data'))
  stopifnot(inherits(cf, 'raw_cf_meta'))
  
  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' must be one of 'real', 'imag' or 'both'")
  }
  if( reim_same & !(reim == 'both') ){
    stop("'reim_same' can only be true if 'reim' is 'both'")
  }

  plotdata <- get_plotdata_raw_cf(cf, reim = reim, tauint = FALSE, relerr = FALSE) 
  for( lidx in 1:length(plotdata) ){
    if( reim %in% c('real', 'both') & reim_same == FALSE){
      plotwitherror(x = 0:(cf$nts-1),
                    y = plotdata[[lidx]]$real$val,
                    dy = plotdata[[lidx]]$real$dval,
                    ...)
    }
    if( reim %in% c('imag', 'both') & reim_same == FALSE){
      plotwitherror(x = 0:(cf$nts-1),
                    y = plotdata[[lidx]]$imag$val,
                    dy = plotdata[[lidx]]$imag$dval,
                    ...)
    }
    if( reim_same ){
      plotwitherror(x = rep(0:(cf$nts-1), times=2),
                    y = c(plotdata[[lidx]]$real$val, plotdata[[lidx]]$imag$val),
                    dy = c(plotdata[[lidx]]$real$dval, plotdata[[lidx]]$imag$dval),
                    ...)
    }
  }
  return(invisible(plotdata))
}

#' @title create convenient overview plots for a `raw_cf` object
#' @param cf 'raw_cf' container with data and meta-data
#' @param grid Optional, integer vector which satisfies 
#'             \code{prod(grid) == prod(cf$dim)}. This is passed to \code{par} via
#'             \code{par(mfrow=grid)} to produce a grid of plots as defined by the
#'             components of \code{grid}.
#' @param reim Vector of strings, one of 'real', 'imag' or 'both'. Specified whether
#'             the real or imaginary parts (or both) should be plotted.
#' @param reim_same Boolean, whether real and imaginary parts should be plotted
#'                  on the same plot. If \code{TRUE}, then \code{reim} must
#'                  be 'both'. If this is given, the imaginary part as well as its
#'                  relative error and per-time-slice integrated autocorreation times
#                   are plotted in red.
#' @param relerr Boolean, whether a plot of the relative error per time slice
#'               should be added.
#' @param tauint Boolean, whether a plot of the integrated auto-correlation time
#'               on each time slice should be added.
#' @param value_logplot Boolean, whether the plot of the correlator should be
#'                      on a logarithmic vertical axis. (does not affect \code{tauint}
#'                      and \code{relerr}).
#' @param value_factor Numeric, either of length '1' or as long as the number of
#'                     correlation functions in \code{cf}. The data will be scaled
#'                     by this factor before plotting.
#' @param title Character vector, will be passed as the \code{main} argument to
#'              \link{plotwitherror} which in turn passes it to \link{plot}. Can
#'              be either of length '1' or \code{prod(cf$dim)}
#'
#' @return
#' No return value, only plots are generated.
#' 
#' @export
overview_plot_raw_cf <- function(cf,
                                 grid,
                                 reim = 'real', 
                                 reim_same = FALSE,
                                 relerr = FALSE,
                                 tauint = FALSE,
                                 value_logplot = TRUE,
                                 value_factor = c(1),
                                 title = '')
{
  stopifnot(inherits(cf, 'raw_cf_data'))
  stopifnot(inherits(cf, 'raw_cf_meta'))

  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' must be one of 'real', 'imag' or 'both'")
  }
  if( reim_same & !(reim == 'both') ){
    stop("'reim_same' can only be true if 'reim' is 'both'")
  }
  if(!missing(grid)){
    if( prod(grid) != prod(cf$dim) ){
      stop("'prod(grid)' must be equal to 'prod(cf$dim)'")
    }
    par_save <- par()
    on.exit(par(par_save))
    par(mfrow=grid)
  }

  # produce a vector of plot symbols such that we have one
  # symbol per observable, all smearing types will be plotted
  # with the same symbol and we wrap around when we run out of symbols
  pch <- expand.grid(rep(1:cf$nts, times=cf$nrStypes), c(0:6, 15:18), KEEP.OUT.ATTRS=FALSE )
  pch <- rep(pch[,2], length.out=cf$nts*cf$nrObs*cf$nrStypes) 

  ylabs <- list()
  ylabs[["real"]] <- "Re[ C(t) ]"
  ylabs[["imag"]] <- "Im[ C(t) ]"
  ylabs[["real_imag"]] <- "Re:Im[ C(t) ]"
  ylabs[["tauint_real"]] <- "tau_int{ Re[ C(t) ] }"
  ylabs[["tauint_imag"]] <- "tau_int{ Im[ C(t) ] }"
  ylabs[["tauint_real_imag"]] <- "tau_int{ Re:Im[ C(t) ] }"
  ylabs[["relerr_real"]] <- "Re[ dC(t) ] / Re[ C(t) ]"
  ylabs[["relerr_imag"]] <- "Im[ dC(t) ] / Im[ C(t) ]"
  ylabs[["relerr_real_imag"]] <- "Re:Im[ dC(t) ] / Re:Im[ C(t) ]"

  tmax <- cf$nts-1  
  plotdata <- get_plotdata_raw_cf(cf, reim=reim, relerr=relerr, tauint=tauint)

  if( length(value_factor) != 1 && length(value_factor) != prod(cf$dim) ){
    message <- sprintf(paste("If 'value_factor' is of length > 1,",
                             "it must be of the same length at the number of correlators",
                             "to be plotted [prod(cf$dim)] !.",
                             "In this case, length(value_factor)=%d, prod(cf$dim)=%d"),
                       length(value_factor), prod(cf$dim) )
    stop(message)
  }
  # if only a single factor has been provided, replicate it as many times as there
  # are correlators to be ploted
  if( length(value_factor) == 1 ){
    value_factor <- rep(value_factor, times = prod(cf$dim) )
  }

  if( length(title) != 1 && length(title) != prod(cf$dim) ){
    message <- sprintf(paste("If 'title' is of length > 1,",
                             "it must be of the same length at the number of correlators",
                             "to be plotted [prod(cf$dim)] !.",
                             "In this case, length(title)=%d, prod(cf$dim)=%d"),
                       length(title), prod(cf$dim) )
    stop(message)
  }
  if( length(title) == 1 ){
    title <- rep(title, times=prod(cf$dim))
  }

  ts <- rep(0:tmax, times=cf$nrStypes*cf$nrObs)
  if( reim_same ){
    ts <- rep(ts,times=2)
  }
  step <- 1

  if( reim_same ) step <- 2
  for( lidx in 1:length(plotdata) ){
    onames <- names(plotdata[[lidx]])
    for( oidx in seq(1,length(onames),step) ){
      args <- list()
      args$x <- ts
      args$xlab <- "t/a"
      args$ylab <- ylabs[[onames[oidx]]]
      args$main <- title[lidx]
      args$pch <- pch
      if( (onames[oidx] == "real" | onames[oidx] == "imag") & value_logplot ){
        args$log <- 'y'
      }
      if( reim_same ){
        args$y <- c(plotdata[[lidx]][[ onames[oidx] ]]$val, plotdata[[lidx]][[ onames[oidx+1] ]]$val)
        args$dy <- c(plotdata[[lidx]][[ onames[oidx] ]]$dval, plotdata[[lidx]][[ onames[oidx+1] ]]$dval)
        args$col <- c(rep("black", cf$nrStypes*cf$nrObs*cf$nts),
                      rep("red", cf$nrStypes*cf$nrObs*cf$nts))
        args$ylab <- ylabs[[sprintf("%s_imag",onames[oidx])]]
      } else {
        args$y <- plotdata[[lidx]][[ onames[oidx] ]]$val
        args$dy <- plotdata[[lidx]][[ onames[oidx] ]]$dval
      }
      # apply the scaling factor to teh expectation value of the observable
      if( (onames[oidx] == "real" | onames[oidx] == "imag") ){
        args$y <- value_factor[lidx]*args$y
        args$dy <- value_factor[lidx]*args$dy
      }

      # when plotting relative errors, it is not useful to plot errors larger than
      # 100%
      if( grepl("relerr",onames[oidx]) ){
        if( max(args$y, na.rm=TRUE) > 1.0 ){
          args$ylim <- c(0,1)
        }
      }
      do.call(plotwitherror, args)
    }
  }
}

#' @title shift a \code{raw_cf} correlation function by 'places' time-slices
#' @param cf \code{raw_cf} container
#' @param places Integer (possibly a vector), number of time slices that the correlation function
#'               should be shifted by. Can be positive or negative. This can either
#'               be a single value such that a shift by this many time slices will be
#'               applied to every measurement or it can be a vector of values of the
#'               same length as the number of measurements in \code{cf}. In that case,
#'               a different shift will be applied to each measurement. This is useful
#'               if it is important to preserve the absolute time coordinates of a
#'               correlation function until some time-dependent transformations
#'               have been applied.
#' @details
#' The correlation funtion \eqn{C(t)} is shifted in time to produce:
#'   \deqn{C'(t) = C(t+places)}
#' using periodic boundary conditions in time.
#'
#' @return
#' Returns an object of class `raw_cf`, shifted compared to the input object.
#'
#' @export
shift.raw_cf <- function(cf, places) {
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  if( (length(places) != 1) & (length(places) != dim(cf$data)[1]) ){
    stop("'places' should be either of length '1' or of a length equalling the number of measurements in 'cf'")
  }
  if(length(places) == 1){
    if(places == 0){
      return(cf)
    }
  }

  dims <- dim(cf$data)

  if( length(places) == 1 ){
    step <- dims[1]
  } else {
    step <- 1
  }

  # if places is of length '1', we do all measurements in one step (see index matrix construction below)
  # else we do each measurement separately
  for( imeas in seq(from = 1, to = dims[1], by = step) ){
    # make sure that we don't do any work when places[
    if( places[imeas] == 0 ) next
    for( oidx in 0:(cf$nrObs-1) ){
      for( sidx in 0:(cf$nrStypes-1) ){
        # the 'time' indices of the observables and smearing types
        istart <- cf$Time*cf$nrStypes*oidx + cf$Time*sidx + 1
        iend <- istart + cf$Time - 1

        # construct an argument list for do.call below
        args <- list()
        args[[1]] <- imeas:(imeas+step-1)
        args[[2]] <- istart:iend
        for( d in cf$dim ){
          args[[length(args)+1]] <- 1:d
        }
        # construct the tensor index set for the output
        out_dof <- as.matrix(do.call(expand.grid, args))
        
        if( places[imeas] < 0 ){
          ishift <- c( (iend - abs(places[imeas]) + 1):iend,
                       (istart:(iend-abs(places[imeas]))) )
        } else {
          ishift <- c( (istart+places[imeas]):iend,
                        istart:(istart+places[imeas]-1))
        }
        args[[2]] <- ishift
        # construct the tensor index set for the input
        in_dof <- as.matrix(do.call(expand.grid, args))

        # shift the correlator tensor
        cf$data[out_dof] <- cf$data[in_dof]
      }
    }
  }
  return(cf)
}

#' @title Construct the tensor index set for the entire raw correlator
#' @param cf 'raw_cf' container with data and meta-data
#' @param component Integer vector. Optional argument to obtain a subset of the
#'                  index matrix to access a particular element of the interior
#'                  dimensions. Must of the the same length as cf$dim.
#'
#' @return
#' An object of type matrix is returned containing the tensor index set.
#' 
#' @export
idx_matrix.raw_cf <- function(cf, component){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))
 
  dims <- dim(cf$data)  
  args <- list()
  for( d in dims ){
    args[[length(args)+1]] <- 1:d
  }

  if(!missing(component)){
    if( length(component) != length(cf$dim) ){
      stop("'component' has to be of length length(cf$dim)!")
    }
    cidx <- 1
    for( didx in 3:length(dims) ){
      args[[didx]] <- component[cidx]
      cidx <- cidx + 1
    }
  }
  as.matrix(do.call(expand.grid, args))
}

#' @title Construct tensor index set for the internal degrees of freedom
#' @param cf `raw_cf` container
#'
#' @return
#' Returns a matrix containing the above mentioned index set.
#' 
#' @export
int_idx_matrix.raw_cf <- function(cf){
  stopifnot(inherits(cf, 'raw_cf_meta'))
  args <- list()
  for( d in cf$dim ){
    args[[length(args)+1]] <- 1:d
  }
  as.matrix(do.call(expand.grid, args))
}


#' @title Print summary of data contained in `raw_cf` container
#' @param object `raw_cf` container with data and meta-data
#' @param ... ignored
#' @param statistics Boolean, return central value and error
#'                   for all components of the 'raw_cf'. This can
#'                   be slow so the default is \code{FALSE}.
#'
#' @return
#' The summary is returned invisibly in form of a data frame.
#' 
#' @export
summary.raw_cf <- function(object, ..., statistics = FALSE) {
  cf <- object
  stopifnot(inherits(cf, 'raw_cf_meta'))
  stopifnot(inherits(cf, 'raw_cf_data'))

  cat("T = ", cf$Time, "\n")
  cat("observations = ", dim(cf$data)[1], "\n")
  cat("Nr Stypes = ", cf$nrStypes, "\n")
  cat("Nr Obs    = ", cf$nrObs, "\n")
  cat("dim       = ", cf$dim, "\n")

  if( statistics ){
    uw <- uwerr.raw_cf(cf)

    idcs <- uw$idx_matrix
    out <- NULL
    for( i in 1:nrow(idcs) ){
      # annoying: subsetting a single row of the index matrix turns the
      # result into a column vector, need to transpose  
      i_out <- t(idcs[i,])
      out <- rbind(out,
                   data.frame(t = idcs[i,1], 
                              value_real = uw$value$real[ i_out ], 
                              value_imag = uw$value$imag[ i_out ],
                              dvalue_real = uw$dvalue$real[ i_out ],
                              dvalue_imag = uw$dvalue$imag[ i_out ])
                   )
    }
    rownames(out) <- NULL
    return(invisible(out))
  }
}

#' @title Print summary of data contained in `raw_cf` container
#' @param x `raw_cf` container with data and meta-data
#' @param ... ignored
#'
#' @return
#' See \link{summary.raw_cf}.
#' 
#' @export
print.raw_cf <- function(x, ...) {
  cf <- x
  summary(cf, ...)
}

#' @title Store a 'raw_cf' correlator in an associative array together with a description
#' The object \code{cf} will be stored as an element of \code{cmap} under key \code{out_key}
#' in the member \code{obj} of \code{cmap}. The data frame passed via \code{desc} will be
#' appended as a row to \code{cmap[[out_key]]$map}. If \code{out_key} does not exist
#' as a key in \code{cmap}, a new element will be created. If it already exists,
#' \code{addStat.raw_cf} is called to add statistics to the existing \code{raw_cf}. Requires
#' the 'hash' package.
#' @return Since objects of class \code{hash} are passed and modified by reference, there
#'         is no explicit return value. Instead, the passed \code{cmap} is modified.
#' @param cmap Object of class \code{hash} to act as storage for 'raw_cf' correlators.
#' @param cf Object of class \code{raw_cf} to be stored in \code{cmap}.
#' @param out_key String, key associated with \code{cf} object to be stored in \code{cmap}.
#' @param desc Single row data frame containing some descriptive parameters for \code{cf}.
store_correl <- function(cmap, cf, out_key, desc)
{
  hash_avail <- requireNamespace("hash")
  if( !hash_avail ){
    stop("The 'hash' package is required to use this function!\n")
  }
  stopifnot( "hash" %in% class(cmap) )
  stopifnot( "raw_cf" %in% class(cf) )
  stopifnot( "raw_cf_meta" %in% class(cf) )
  stopifnot( "raw_cf_data" %in% class(cf) )

  if( hash::has.key(out_key, cmap) ){
    cmap[[out_key]]$obj <- addStat.raw_cf(cmap[[out_key]]$obj, cf)
    cmap[[out_key]]$map <- rbind(cmap[[out_key]]$map, desc)
  } else {
    cmap[[out_key]] <- list()
    cmap[[out_key]]$obj <- cf
    cmap[[out_key]]$map <- desc
  }
}
