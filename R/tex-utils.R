## function to convert the scientific eE into proper 10^XX
## for tex output
convert.scientific <- function(str, errstr) {
  ## check whether or not str is in scientific format
  if(grepl("[+-]?(\\d+(\\.\\d+)?|\\.\\d+)([eE][+-]?\\d+)", str)) {
    ## position of e or E in str
    posE <- regexpr("[eE]", str)
    ## replace by \cdot 10^{
    if(missing(errstr)) {
      str <- paste(substr(str, 1, posE[1]-1), "\\cdot 10^{", substr(str, posE[1]+1, nchar(str)), "}", sep="")
    }
    else {
      str <- paste(substr(str, 1, posE[1]-1), "(", errstr, ")",
                   "\\cdot 10^{", substr(str, posE[1]+1, nchar(str)), "}", sep="")
    }
  }
  else if(!missing(errstr)){
    str <- paste(str, "(", errstr, ")", sep="")
  }
  return(str)
}


## Convert the number of significant digits to the absolute number of digits, 
## given a value.
absolute.number.digits <- function(x, digits){
  if(x == 0 || is.na(x)){
    return(digits-1)
  }else{
    return(ceiling(digits-1-log10(abs(x))))
  }
}

#' paste a number with error in tex-ready format
#' 
#' A number with error is converted to a string in tex-ready format like xx(yy)
#' thereby automatically determining the digit at which the error applies.
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
#' @param human.readable controls whether or not scientific format will be
#' produced in a human readable form, or not at all. The latter might be useful
#' for output with summary or print.
#' @return writes a string to standard output
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @examples
#' 
#' tex.catwitherror(x=0.375567, dx=0.001)
#' tex.catwitherror(x=c(0.375567, 0.001))
#' ## it can be used with apply
#' x = array(c(0.1187, 0.291, 0.388, 0.011, 0.037, 0.021), dim=c(3,2))
#' apply(x, 1, tex.catwitherror, digits=2)
#' 
#' @export tex.catwitherror
tex.catwitherror <- function(x, dx, digits=1, with.dollar=TRUE, human.readable=TRUE) {
  if(missing(x) || length(x) == 0) {
    stop("x must be a numeric vector with length > 0")
  }
  lx <- length(x)
  tmp <- ""

  ## In case the number is very small, printing it with fixed point (`%f`) will
  ## not work. For this case we divide out the exponent from both value and
  ## error and attach it at the end.
  if (is.infinite(x[1]) || is.na(x[1]) || x == 0) {
    base <- 0
  } else {
    base <- floor(log10(abs(x[1])))
  }
  split_base <- abs(base) > 19
  if (split_base) {
    x <- x / 10**base
    if (!missing(dx)) {
      dx <- dx / 10**base
    }
  }

  if(missing(dx) && lx < 2) {
    ## just a number without error
    N <- absolute.number.digits(x, digits)
    tmp <- paste(format(round(x, digits=N), nsmall=N), sep="")
    if(human.readable) tmp <- convert.scientific(str=tmp)
    else tmp <- paste(format(round(x, digits=N), nsmall=N, scientific=FALSE), sep="")
  }
  else {
    ## now we need to typeset the error as well
    err <- 0.
    if(missing(dx)){
      err <- x[2]
    } else {
      err <- dx
    }
    if(lx > 1){
      x <- x[1]
    }

    if(!is.na(err) && err > 0) {
      N <- absolute.number.digits(err, digits)
      # if the error is large it may exceed the number of digits that one actually desires
      # also, the error may be larger or similar in size to the value itself
      # in these cases, we display it in the same format as the value, rounded to the
      # desired number of digits
      displayerr <- paste(round(10^N*err))
      if(is.na(x) |
         nchar(displayerr) > digits |
         ( ceiling(log10(abs(err)/abs(x))) >= 0 && ( abs(err) >= 1.0 ) ) |
         ( abs(err) >= abs(10*x) ) ){
        displayerr <- paste(format(round(err, digits=N)))
      }

      if(human.readable){
        tmp <- convert.scientific(str = format(round(x, digits=N), nsmall=N), 
                                  errstr = displayerr)
      } else {
        tmp <- paste(format(round(x, digits=N), nsmall=N, scientific=FALSE), "(", displayerr, ")", sep="")
      }
    }else {
      N <- absolute.number.digits(x, digits)
      displayerr <- paste(format(err))
      tmp <- paste(format(round(x, digits=N), nsmall=N), sep="")
      if(human.readable) tmp <- convert.scientific(str=tmp, errstr = displayerr)
      else tmp <- paste(format(round(x, digits=N), nsmall=N, scientific=FALSE), "(", displayerr, ")", sep="")
    }
  }
  if (with.dollar) {
    if (split_base) {
      tmp <- paste0(tmp, '\\cdot 10^{', base, '}')
    }
    tmp <- paste0("$", tmp, "$")
  } else {
    if (split_base) {
      tmp <- paste0(tmp, 'e', base)
    }
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

