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

tex.catwitherror <- function(x, dx, digits=1, with.dollar=TRUE, human.readable=TRUE) {
  if(missing(x) || !is.numeric(x) || length(x) == 0) {
    stop("x must be a numeric vector with length > 0")
  }
  lx <- length(x)
  tmp <- ""
  if(missing(dx) && lx < 2) {
    ## just a number without error
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*abs(x)) < threshold) {
      N <- N+1
    }
    tmp <- paste(format(round(x, digits=N), nsmall=N), sep="")
    if(human.readable) tmp <- convert.scientific(str=tmp)
    else tmp <- paste(format(round(x, digits=N), nsmall=N, scientific=FALSE), sep="")
  }
  else {
    ## now we need to typeset the error as well
    err <- 0.
    if(missing(dx)) err <- x[2]
    else err <- dx[1]
    if(lx > 1) x <- x[1]
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*err) < threshold) {
      N <- N+1
    }
    # if the error is large it may exceed the number of digits that one actually desires
    # also, the error may be larger or similar in size to the value itself
    # in these cases, we display it in the same format as the value, rounded to the
    # desired number of digits
    displayerr <- paste(round(10^N*err))
    if( nchar(displayerr) > digits | ceiling(log10(abs(err/x))) >= 1 ){
      displayerr <- paste(format(round(err, digits=N)))
    }

    if(human.readable){
      tmp <- convert.scientific(str = format(round(x, digits=N), nsmall=N), 
                                errstr = displayerr)
    } else {
      tmp <- paste(format(round(x, digits=N), nsmall=N, scientific=FALSE), "(", displayerr, ")", sep="")
    }
  }
  ret <- tmp
  if(with.dollar) {
    ret <- paste("$", tmp, "$", sep="")
  }
  return(ret)
}
