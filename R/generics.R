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

#' generic function to extract a fitted mass
#'
#' @description
#' One of the main analysis tasks in \link{hadron} is the estimation
#'   of energy levels or masses from correlation functions. The
#'   corresponding analysis functions return objects, typically lists,
#'   containing the masses or energy levels. `extract_mass` is a
#'   generic function to extrac such fitted mass values.
#' 
#' @param object Object to extract the mass from.
#'
#' @return Numeric. The mass value.
#'
#' @export
extract_mass <- function (object) {
  UseMethod('extract_mass')
}

#' generic function to check if resampling samples are compatible
#'
#' @description
#' When binary operations are performed on resampled data, it is
#' necessary to check if the samples are compatible, at the very
#' least they must have the same dimensions.
#' 
#' @details
#' Note that since R's double dispatch doesn't really work with S3
#' classes, the class of \code{x} decides which method is called.
#' 
#' @param x lhs object
#' @param y rhs object
#'
#' @return 
#' list of booleans which correspond to the results of checks
#' for equality of different properties of the resampling samples
#'
#' @export
resampling_is_compatible <- function(x,y){
  UseMethod('resampling_is_compatible', x)
}

#' generic function to check if resampling samples are concatenable
#'
#' @description
#' When we want to combine resampled data along an axis orthogonal
#' to the axis of samples (for example, if we want to turn two
#' vectors with 'boot.R' samples into a matrix with 2 columns),
#' then we need to check if the number of samples for both
#' data are the same.
#' 
#' @details
#' Note that since R's double dispatch doesn't really work with S3
#' classes, the class of \code{x} decides which method is called.
#' 
#' @param x lhs object
#' @param y rhs object
#'
#' @return
#' list of booleans which correspond to the results of checks
#' for equality of different properties of the resampling samples
#'
#' @export
resampling_is_concatenable <- function(x,y){
  UseMethod('resampling_is_concatenable', x)
}

#' generic function to multiply two objects with each other
#'
#' @details
#' Note that since R's double dispatch doesn't really work with S3
#' classes, the class of \code{x} decides which method is called.
#'
#' @description
#' 
#' 
#' @export
mul <- function(x, y, ...){
  UseMethod('mul', x)
}

#' generic function to add two objects to each other
#'
#' This function provides 
#'
#' @export
add <- function(x, y, ...){
  UseMethod('add', x)
}

#' @export
subtract <- function(x, y, ...){
  UseMethod('subtract', x)
}

#' generic function to divide objects by each other
#' @export
div <- function(x, y, ...){
  UseMethod('div', x)
}
