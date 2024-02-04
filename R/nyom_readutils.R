#' @title Convert correlation function read from nyom HDF5 format to 'raw_cf'
#' @description Given a data frame with two columns 'r' and 'i' of real and imaginary components
#'              of correlation function, creates an object of class 'raw_cf' with
#'              a single measurement, inferring \code{Time} from the passed numeric vector
#'              while the shape of the internal dimensions has to be specified explicitly
#'              if larger than `one by one` (\code{c(1,1)}).
#' @param cf_dat data frame with two columns 'r' and 'i' corresponding to the real and imaginary parts
#'               of a single measurement of a correlation function. 
#'               Ordering of the input should internal dimensions, time (left to right, fastest to slowest).
#' @param dims Integer vector with the sizes of the internal dimensions. For example,
#'             \code{c(4,4)} for spin correlators or \code{c(8,8)} for a spin-flavour correlator.
#' @param mult Numeric or complex vector of length 1 or \code{prod(dims)} to apply an element-wise
#'             scaling or sign change to the correlator along the internal dimensions.    
#' @return `raw_cf` object with a \code{data} member which contains the data (as complex numbers)
#'         in the shape \code{c(1,nts,dims)}, where `nts` is the number of time slices
#'         inferred from the length of \code{cfdat} and the product of the internal dimensions \code{dims}.
#'                 
#' @export
nyom_to_raw_cf <- function(cf_dat, dims = c(1,1), mult = 1.0)
{
  number_of_internal_dims <- prod(dims)
  stopifnot( length(mult) == 1 | length(mult) == prod(dims) )

  nts <- nrow(cf_dat)/number_of_internal_dims

  # turn it into a complex 4D array of with ordering of
  # internal dimensions and time in the 'wrong' order
  cf_dat <- array(complex(real = cf_dat$r,
                          imaginary = cf_dat$i)*mult,
                  dim = c(1,dims,nts))
  cf_dims <- dim(cf_dat)
  # reshape the array, ordering 'measurement' (a single one), 'time',
  # internal dimensions
  cf_dat <- aperm(cf_dat, c(1,length(cf_dims),2:(length(cf_dims)-1)))

  cf <- raw_cf_data(cf = raw_cf_meta(Time=nts, dim=dims),
                    data = cf_dat)
  return(cf)
}

