#' @title numerically invert the cosh function for the mass
#'
#' @param ratio Numeric. The value of the ratio.
#' @param timeextent Integer. Time extent of the lattice.
#' @param t Integer. The t-value where the ratio was taken.
#' @param eps Numeric. Precision of the numerical solution
#' @param maxiterations Integer. Maximal number of iterations to be
#'   used in the iterative solver.
#'
#' @useDynLib hadron
#' @importFrom Rcpp evalCpp
#'
#' @return
#' A single numeric value is returned corresponding to the mass.
#' @examples
#'
#' invcosh(1.2, timeextent=24, t=12)
#' @export
invcosh <- function(ratio, timeextent, t, eps=1.e-9, maxiterations=1000) {

  if(ratio < 1 || is.na(ratio)) {
    return(NA);
#    stop("Error: ratio is smaller than 1 in invcosh!")
  }

  return(.Call("invcosh", ratio, timeextent, t, eps, maxiterations))
}
