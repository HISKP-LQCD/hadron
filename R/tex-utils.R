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
  N <- 0
  threshold <- 10^(digits-1)

  if (is.na(x[1])) {
    base <- 0
  } else {
    base <- round(log10(x[1]))
  }
  split_base <- abs(base) > 19
  if (split_base) {
    x <- x / 10**base
    if (!missing(dx)) {
      dx <- dx / 10**base
    }
  }

  if(missing(dx) && lx < 2) {
    if( is.na(x) ) x <- 0.0
    ## just a number without error
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
    if(missing(dx)){
      if( !is.na(x[2]) ){
        err <- x[2]
      }
    } else {
      if( !is.na(dx[1]) ){
        err <- dx[1]
      }
    }
    if(lx > 1){
      x <- x[1]
      if( is.na(x) ){
        x <- 0.0
      }
    }

    if(err > 0) {
      while(round(10^N*err) < threshold) {
	N <- N+1
      }
      # if the error is large it may exceed the number of digits that one actually desires
      # also, the error may be larger or similar in size to the value itself
      # in these cases, we display it in the same format as the value, rounded to the
      # desired number of digits
      displayerr <- paste(round(10^N*err))
      if( nchar(displayerr) > digits |
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
      while(round(10^N*abs(x)) < threshold) {
	N <- N+1
      }
      displayerr <- paste(format(0))
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

escape_underscore <- function(x){
  gsub(pattern="_",
       replacement="\\_",
       x=x,
       fixed=TRUE)
}
