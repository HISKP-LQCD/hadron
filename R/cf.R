bootstrap.cf <- function(cf, boot.R=400, boot.l=2, seed=1234) {
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.cf requires an object of class cf as input! Aborting!\n")
  }
  cf$boot.samples <- TRUE
  cf$boot.R <- boot.R
  cf$boot.l <- boot.l
  cf$seed <- seed
  cf$cf0 <- apply(cf$cf, 2, mean)
  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## now we bootstrap the correlators
  cf$cf.tsboot <- tsboot(cf$cf, statistic = function(x){ return(apply(x,2,mean))},
                         R = boot.R, l=boot.l, sim="geom")
  return(invisible(cf))
}

## this is intended for instance for adding diconnected diagrams to connected ones
add.cf <- function(cf1, cf2, a=1., b=1.) {
  if(any(class(cf1) == "cf") && any(class(cf2) == "cf") &&
     all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- a*cf1$cf + b*cf2$cf
    cf$boot.samples <- FALSE
    return(cf)
  }
  else {
    stop("The two objects of class cf are not compatible\n Aborting...!\n")
  }
}

## to concatenate objects of type cf
c.cf <- function(...) {
  fcall <- match.call(expand.dots=TRUE)
  fnames <- names(fcall)
  ## first name in fnames is empty/function name
  if(length(fcall) == 2) {
    return(eval(fcall[[2]]))
  }

  cf <- eval(fcall[[2]])
  Time <- cf$Time
  cf$nrObs <- 0
  cf$sTypes <- 0
  N <- dim(cf$cf)[1]
  for(i in 2:length(fcall)) {
    if(eval(fcall[[i]])$Time != Time) {
      stop("Times must agree for different objects of type cf\n Aborting\n")
    }
    if(dim(eval(fcall[[i]])$cf)[1] != N) {
      stop("Number of measurements must agree for different objects of type cf\n Aborting\n")
    }
    cf$nrObs <- cf$nrObs + eval(fcall[[i]])$nrObs
    cf$sTypes <- cf$sTypes + eval(fcall[[i]])$sTypes
  }
  for(i in 3:length(fcall)) {
    cf$cf <- cbind(cf$cf, eval(fcall[[i]])$cf)
  }
  cf$boot.samples <- FALSE
  return(cf)
}

plot.cf <- function(cf, boot.R=400, boot.l=2, ...) {
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R, boot.l)
  }
  Err <- apply(cf$cf.tsboot$t, 2, sd)
  plotwitherror(rep(c(0:(cf$Time/2)), times=cf$nrStypes*cf$nrObs), cf$cf0, Err, ...)
}
