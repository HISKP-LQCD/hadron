#' paste a number with error in tex-ready format
#' 
#' A number with error is converted to a string in tex-ready format like xx(yy)
#' thereby automatically determining the digit at which the error applies.
#'
#' It is strongly recommended to install the \code{errors}-package. Otherwise 
#' the formatting options are significantly reduced.
#' 
#' The value of the first element of \code{x} is properly rounded to its
#' significant digits determined by the values of \code{dx} or the second
#' element of \code{x} (see above) and \code{digits}. Then a tex-ready string
#' is returned.
#' 
#' @param x either a single numeric value, or a numeric vector, where the first
#' element is the value and the second is its error
#' @param dx the error. If supplied, it will be printed as the error and the
#' value is the first element of \code{x}. If \code{dx} is missing, the second
#' element of \code{x}, if available, is used as the error. If \code{dx} is
#' missing and the length of \code{x} is one, only the value is converted to a
#' string without error.
#' @param digits number of error digits
#' @param with.dollar include the tex dollar in the return string or not
#' @param with.cdot replace the "e" in scientific notation by tex-style "cdot" or not
#' @param ... arguments to be passed to \code{formatC} in case of no error, or 
#' to \code{format.errors} otherwise
#' @return writes a string to standard output
#' @author Carsten Urbach, \email{curbach@@gmx.de} and Johann Ostmeyer
#' @examples
#' 
#' tex.catwitherror(x=0.375567, dx=0.001)
#' tex.catwitherror(x=c(0.375567, 0.001))
#' ## it can be used with apply
#' x = array(c(0.1187, 0.291, 0.388, 0.011, 0.037, 0.021), dim=c(3,2))
#' apply(x, 1, tex.catwitherror, digits=2)
#' 
#' @export tex.catwitherror
tex.catwitherror <- function(x, dx, digits=1, with.dollar=TRUE, with.cdot=TRUE, ...) {
  if(missing(x) || length(x) == 0) {
    stop("x must be a numeric vector with length > 0")
  }
  if(length(x) == 2){
    dx <- x[2]
    x <- x[1]
    have.error <- TRUE
  }else if(!missing(dx)){
    have.error <- TRUE
  }else{
    have.error <- FALSE
  }

  if(!have.error){
    ## just a number without error
    tmp <- formatC(x, digits=digits, ...)
  }
  else if(is.numeric(dx) && dx == 0){
    tmp <- formatC(x, digits=digits, ...)
    if(grepl("e", tmp, fixed=TRUE)){
      tmp <- sub("e", "(0)e", tmp)
    }else{
      tmp <- paste0(tmp, "(0)")
    }
  }
  else{
    if(requireNamespace('errors')){
      tmp <- format(errors::set_errors(x, dx), digits=digits, ...)
    }else{
      warning("The `errors`-package is not installed. The output of `tex.catwitherror` might not be as you want it.")
      tmp <- formatC(x, digits=digits, ...)
      tmp <- paste0(tmp, " +- ", formatC(dx, digits=digits, ...))
    }
  } 

  if(with.cdot){
    if(grepl("e", tmp, fixed=TRUE)){
      tmp <- sub("e", "\\\\cdot 10^{", tmp)
      tmp <- paste0(tmp, "}")
    }
    tmp <- sub("+-", "\\\\pm", tmp, fixed=TRUE)
  }
  if (with.dollar) {
    tmp <- paste0("$", tmp, "$")
  }
  return (tmp)
}


#' @title Escape special LaTeX characters for use in LaTeX labels
#'
#' @param x String or vector of strings.
#' @return String or vector of strings with all occurences of "#", "$", "%",
#'        "&", "~", "_", "^", ">", "<" replaced by escaped
#'        counterparts which should render fine when used in a tikz plot, for
#'        example.
#' @export
#'
#' @references 
#' from https://stackoverflow.com/questions/36338629/escaping-special-latex-characters-in-r
escapeLatexSpecials <- function(x) {
  x <- gsub("\\", "$\\backslash$", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed=TRUE)
  x <- gsub("$", "\\$", x, fixed=TRUE)
  x <- gsub("%", "\\%", x, fixed=TRUE)
  x <- gsub("&", "\\&", x, fixed=TRUE)
  x <- gsub("~", "\\~", x, fixed=TRUE)
  x <- gsub("_", "\\_", x, fixed=TRUE)
  x <- gsub("^", "\\^", x, fixed=TRUE)
  x <- gsub(">", "$>$", x, fixed=TRUE)
  x <- gsub("<", "$<$", x, fixed=TRUE)
  return(x)
}

