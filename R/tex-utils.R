tex.catwitherror <- function(x, dx, digits=1, with.dollar=TRUE) {
  if(missing(dx)) {
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*x) < threshold) {
      N <- N+1
    }
    if(with.dollar){
      ret <- paste("$", format(round(x, digits=N), nsmall=N), "$", sep="")
    }
    else {
      ret <- paste(format(round(x, digits=N), nsmall=N), sep="")
    }
  }
  else {
    N <- 0
    threshold <- 10^(digits-1)
    while(round(10^N*dx) < threshold) {
      N <- N+1
    }
    if(with.dollar) {
      ret <- paste("$", format(round(x, digits=N), nsmall=N), "(", round(10^N*dx), ")$", sep="")
    }
    else {
      ret <- paste(format(round(x, digits=N), nsmall=N), "(", round(10^N*dx), ")", sep="")
    }
  }
  return(ret)
}
