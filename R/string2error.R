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
#' @export
#' @examples
#' string2error("0.35667(25)")
#'
#' s <- c("0.35667(25)", "0.667(50)")
#' apply(array(s, dim=c(1, length(s))), 2, string2error)
#' 
string2error <- function (x) {
  if (is.na(x)) {
    return (c(NA, NA))
  }
  stopifnot(is.character(x))
  
  # We use a regular expression to match the whole string. This makes sure that
  # we reject strings that do not match the format.
  match <- stringr::str_match(x, '^([+-]?[\\d.]+)\\((\\d+)\\)$')
  stopifnot(all(dim(match) == c(1, 3)))
  first <- match[[1, 2]]
  second <- match[[1, 3]]
  
  value <- as.numeric(first)
  
  # We need to differentiate whether the number actually has a decimal point as
  # this changes the interpretation of the error.
  first_parts <- strsplit(first, '.', fixed = TRUE)[[1]]
  if (length(first_parts) == 1) {
    # There is no decimal point. Therefore the error is just to be taken as is.
    error <- as.numeric(second)
  } else {
    # There is a decimal point, therefore the error is to be scaled down by that
    # many digits.
    error <- as.numeric(second) * 10^(- nchar(first_parts[2]))
  }

  return(c(value, abs(error)))
}
