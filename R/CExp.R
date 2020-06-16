#' Cosh Or Sinh Build Out Of Two Exps
#' 
#' Evaluates \deqn{f(x) = \frac{1}{2}(\exp(-m(T-x))\pm\exp(-m x))}{% f(x) = 1/2
#' (exp(-m(T-x)) +/- exp(-m x))} for given mass \eqn{m}, vector \eqn{x} and
#' time extent \eqn{T}. This form is better usable in \eqn{\chi^2}{chi^2}
#' fitting than cosh or sinh.
#' 
#' 
#' @param m mass value
#' @param Time Time extent
#' @param x vector of values on which to evaluate the function
#' @param sign with sign=1 cosh is evaluated, with sign=-1 sinh
#' @return vector \eqn{f(x)}
#' @author Carsten Urbach \email{carsten.urbach@@liverpool.ac.uk}
#' @keywords math
#' @examples
#'
#' m <- 0.1
#' Time <- 48
#' x <- seq(0, 48, 1)
#' CExp(m=m, Time=Time, x=x)
#' @export CExp
CExp <- function(m, Time, x, sign=1.) {
  return(0.5*(exp(-m*(Time-x)) + sign*exp(-m*x)))
}

dCExpdm <- function(m, Time, x, sign=1.) {
  return(0.5*(-(Time-x)*exp(-m*(Time-x)) -x* sign*exp(-m*x)))
}

