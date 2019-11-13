#' string2error
#'
#' @description
#' takes a string of the form "x(dx)", where dx are the error digits
#' and returns a numeric vector c(x, y), where y is dx as a proper
#' numeric value.
#'
#' @param x Input character string.
#'
#' @return a numeric vector with the first element the value and the
#' second the error
#'
#' @details
#' can be used in combination with \link{apply}
#' 
#' @examples
#' string2error("0.35667(25)")
#'
#' s <- c("0.35667(25)", "0.667(50)")
#' apply(array(s, dim=c(1, length(s))), 2, string2error)
#' 
string2error <- function(x) {
  if(is.na(x)) {
    return(c(NA, NA))
  }
  stopifnot(is.character(x))

  X <- strsplit(x=x, split="\\(|\\)")[[1]]
  N <- nchar(X)
  err <- paste0(gsub("[0-9]", "0", substr(X[1], 1, N[1]-N[2])), X[2])
  return(as.numeric(c(X[1], err)))
}
