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

tex.catwitherror <- function(x, dx, digits=1, with.dollar=TRUE) {
  tmp <- ""
  if(missing(dx)) {
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*x) < threshold) {
      N <- N+1
    }
    tmp <- paste(format(round(x, digits=N), nsmall=N), sep="")
    tmp <- convert.scientific(str=tmp)
  }
  else {
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*dx) < threshold) {
      N <- N+1
    }
    tmp <- convert.scientific(str=format(round(x, digits=N), nsmall=N), errstr=paste(round(10^N*dx)))
  }
  ret <- tmp
  if(with.dollar) {
    ret <- paste("$", tmp, "$", sep="")
  }
  return(ret)
}
