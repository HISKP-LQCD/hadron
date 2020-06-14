## The functions in this file date back to the very beginnings of hadron
## and are still used in a number of legacy routines, so they are kept
## around


#' effmass
#'
#' computes the effective mass via the inverse cosh
#'
#' @param data numeric vector. data vector of length 4
#' @param timeextent integer. time extent of the lattice
#' @param t integer. physical time at which to evaluate the cosh
#'
#' @return
#' Returns the effective mass as a single numeric value.
#' 
effmass <- function(data, timeextent, t) {
  mass <- invcosh((data[1]+data[4])/(data[2]+data[3]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

#' effmass2
#'
#' computes the effective mass via the inverse cosh
#'
#' @param data numeric vector. data vector of length 4
#' @param timeextent integer. time extent of the lattice
#' @param t integer. physical time at which to evaluate the cosh
#' 
#' @return
#' Returns the effective mass as a single numeric value.
effmass2 <- function(data, timeextent, t) {
  mass <- invcosh(ratio=(data[1])/(data[2]), timeextent=timeextent, t=t)
  return(invisible(mass))
}

#' effectivemass
#'
#' computes the effective mass with error analysis using UWerr
#'
#' @param from integer. Fit in fitrange (from, to)
#' @param to integer. see from.
#' @param Time integer. time extent of the lattice
#' @param Z data
#' @param pl boolean. plot
#' @param S numeric. see \link{uwerr}
#' @param ... additional parameters passed to \link{uwerr}
#'
#' @seealso \link{uwerr}
#' @return
#' Returns a \link{data.frame} with named columns `t`, `mass`, `dmass`,
#' `ddmass`, `tauint` and `dtauint`.
effectivemass <- function(from, to, Time, Z, pl=TRUE, S,...) {
  L <- (to-from+1)
  i <- 1
  result <- data.frame(t = array(0.,dim=c(L)), mass = array(0.,dim=c(L)), dmass = array(0.,dim=c(L)),
                       ddmass = array(0.,dim=c(L)), tauint = array(0.,dim=c(L)), dtauint = array(0.,dim=c(L)))
  for(t in from:to) {
    try(mass <- uwerrderived(f=effmass2, data=t(Z[t:(t+1),]), S=S, pl=F, timeextent=Time, t=t, ...))

    result$t[i] <- t-1
    # these are NA, rather than single element vectors containing NA, so we need
    # to check for their lengths, otherwise the assignments below fail 
    if( length(mass$value) > 0 ){
      result$mass[i] <- mass$value[1]
    }
    if( length(mass$dvalue) > 0 ){
      result$dmass[i] <- mass$dvalue[1]
    }
    if( length(mass$ddvalue) > 0 ){
      result$ddmass[i] <- mass$ddvalue[1]
    }
    if( length(mass$tauint) > 0 ){
      result$tauint[i] <- mass$tauint[1]
    }
    if( length(mass$dtauint) > 0 ){  
      result$dtauint[i] <- mass$dtauint[1]
    }
    i = i+1
  }
  rm(mass)
  rm(Z)
  attr(result, "class") <- c("massfit", "data.frame")
  if(pl == T) {
    new_window_if_appropriate()
    plot(result)
  }
  return(invisible(result))
}

