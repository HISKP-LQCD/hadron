getCor <- function(T1, W, Z, type=c("cosh")) {

  ## iobs enumerating the gamma matrix combination
  ## ityp enumeratiog the smearing level
  N <- length(type)
  sign = rep(+1., times=N)
  for(i in 1:N) {
    if(type[i]=="sinh") {
      sign[i] = -1.
    }
  }
  
  for(j in 1:N) {
    for(i in 1:(T1)) {
      two <- 2.
      if(i==1 || i==(T1)) {
        ## Take care of zeros in the correlators when summing t and T-t+1
        two <- 1.
      }
      
      W[(i+(j-1)*T1),] <- (W[(i+(j-1)*T1),]
        + sign[j]*Z[(i+(j-1)*T1),])/two
    }
  }
  return(invisible(W))
}
