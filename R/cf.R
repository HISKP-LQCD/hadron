cf <- function() {
  cf <- list(cf=NULL)
  attr(cf, "class") <- c("cf", class(cf))
  return(cf)
}

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

addConfIndex2cf <- function(cf, conf.index) {
  if(is.null(cf$conf.index)) {
    cf$conf.index <- conf.index
  }
  return(cf)
}

addStat.cf <- function(cf1, cf2) {
  if(inherits(cf1, "cf") && inherits(cf2, "cf") &&
     cf1$Time == cf2$Time && dim(cf1$cf)[2] == dim(cf2$cf)[2] &&
     cf1$nrObs == cf2$nrObs && cf1$nrStypes == cf2$nrStypes
     ){
    cf <- cf1
    cf$boot.samples <- FALSE
    cf$boot.R <- NULL
    cf$boot.l <- NULL
    cf$seed <- NULL
    cf$cf <- rbind(cf1$cf, cf2$cf)
  }
  else {
    stop("addStat.cf: cf1 and cf2 not compatible. Aborting...\n")
  }
}

## averages local-smeared and smeared-local correlators in cf and adjusts
## nrStypes accordingly
## by default, assumes that LS and SL are in columns (T/2+1)+1:3*(T/2+1)
avg.ls.cf <- function(cf,cols=c(2,3)) {
  if(!any(class(cf) == "cf")) {
    stop("Input must be of class 'cf'\n")
  }
  if(cf$nrStypes < 2) {
    stop("There must be at least 2 smearing types in cf!\n")
  }
  timeslices <- cf$Time/2+1

  ind.ls <- ( (cols[1]-1)*timeslices+1 ):( cols[1]*timeslices )
  ind.sl <- ( (cols[2]-1)*timeslices+1 ):( cols[2]*timeslices )

  cf$cf[,ind.ls] <- 0.5 * ( cf$cf[,ind.ls] + cf$cf[,ind.sl] )

  cf$cf <- cf$cf[,-ind.sl]
  cf$nrStypes <- cf$nrStypes-1
  return(cf)
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

'+.cf' <- function(cf1, cf2) {
  if(all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- cf1$cf + cf2$cf
    cf$boot.samples <- FALSE
    cf$boot.R <- NULL
    cf$boot.l <- NULL
    cf$seed <- NULL

    return(cf)
  }
}

'-.cf' <- function(cf1, cf2) {
  if(all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- cf1$cf - cf2$cf
    cf$boot.samples <- FALSE
    return(cf)
  }
}

'/.cf' <- function(cf1, cf2) {
  if(all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- cf1$cf / cf2$cf
    cf$boot.samples <- FALSE
    return(cf)
  }
}

mul.cf <- function(cf, a=1.) {
  if(any(class(cf) == "cf") && is.numeric(a)) {
    cf$cf <- a*cf$cf
    cf$boot.samples=FALSE
    return(cf)
  }
  else {
    stop("Wrong classes for input objects, must be cf and numeric. Aborting...!\n")
  }
}

extractSingleCor.cf <- function(cf, id=1) {
  if(!inherits(cf, "cf")) {
    stop("extractSingleCor.cf: cf must be of class 'cf'. Aborting...\n")
  }
  cf$cf0 <- NULL
  N <- dim(cf$cf)[1]
  nCors <- N/(cf$Time/2+1)
  ii <- c(1:(cf$Time/2+1)) + (id-1)*(cf$Time/2+1)
  cf$cf <- cf$cf[,ii]
  cf$boot.samples <- FALSE
  cf$nrObs <- 1
  cf$nsStypes <- 1
  return(cf)
}

as.cf <- function(x){
  if(!inherits(x, "cf")) class(x) <- c("cf", class(x))
  x
}

is.cf <- function(x){
  inherits(x, "cf")
}

## to concatenate objects of type cf
c.cf <- function(...) {
  #fcall <- match.call(expand.dots=TRUE)
  #fnames <- names(fcall)
  ## first name in fnames is empty/function name
  fcall <- list(...)
  if(length(fcall) == 1) {
    return(eval(fcall[[1]]))
  }

  k <- -1
  for(i in 1:length(fcall)) {
    if(!is.null(fcall[[i]]$cf)) {
      k <- i
      break
    }
  }
  if(k == -1) return(eval(fcall[[1]]))
  cf <- fcall[[k]]
  Time <- cf$Time
  cf$nrObs <- 0
  cf$sTypes <- 0
  N <- dim(cf$cf)[1]
  for(i in k:length(fcall)) {
    if(!is.null(fcall[[i]]$cf)) {
      if(fcall[[i]]$Time != Time) {
        stop("Times must agree for different objects of type cf\n Aborting\n")
      }
      if(dim(fcall[[i]]$cf)[1] != N) {
        stop("Number of measurements must agree for different objects of type cf\n Aborting\n")
      }
      cf$nrObs <- cf$nrObs + fcall[[i]]$nrObs
      cf$sTypes <- cf$sTypes + fcall[[i]]$sTypes
    }
  }
  if(k < length(fcall)) {
    for(i in (k+1):length(fcall)) {
      if(!is.null(fcall[[i]]$cf)) {
        cf$cf <- cbind(cf$cf, fcall[[i]]$cf)
      }
    }
  }
  cf$boot.samples <- FALSE
  return(invisible(cf))
}

plot.cf <- function(cf, boot.R=400, boot.l=2, ...) {
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R, boot.l)
  }
  Err <- apply(cf$cf.tsboot$t, 2, sd)
  plotwitherror(rep(c(0:(cf$Time/2)), times=length(cf$cf0)/(cf$Time/2+1)), cf$cf0, Err, ...)
  return(invisible(data.frame(t=rep(c(0:(cf$Time/2)), times=length(cf$cf0)/(cf$Time/2+1)), CF=cf$cf0, Err=Err)))
}
