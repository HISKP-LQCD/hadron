#' plot.pionfit
#'
#' @description
#' Generic function to plot an object of type `pionfit`
#'
#' @param x Object of type `pionfit`
#' @param ... Generic graphical parameter, ignored.
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.pionfit <- function(x, ...) {
  plot.cfit(x)
}

#' plot.rhofit
#'
#' @description
#' Generic function to plot an object of type `rhofit`
#'
#' @param x Object of type `rhofit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.rhofit <- function(x, ...) {
  plot.cfit(x)
}

#' plot.b1fit
#'
#' @description
#' Generic function to plot an object of type `b1fit`
#'
#' @param x Object of type `b1fit`
#' @param ... Generic graphical parameter to be passed on to \link{plotwitherror}
#'
#' @return
#' See \link{plot.cfit}
#' 
#' @export
plot.b1fit <- function(x, ...) {
  plot.cfit(x)
}

