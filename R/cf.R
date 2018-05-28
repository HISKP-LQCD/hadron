#' Correlation function container
#'
#' This function `cf()` creates containers for correlation functions
#' of class `cf`. This class is particularly designed to deal with
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.cf} and \link{as.cf}.
#'
#' This class _must_ contain the following fields:
#'
#' - `cf`: Numeric matrix, original data for all observables and measurements.
#'
#' It _may_ also contain the following fields:
#'
#' - `icf`: Numeric matrix, imaginary part of original data. Be very careful with this as most functions just ignore the imaginary part and drop it in operations.
#' - `boot.samples`: Logical, indicating whether there are bootstrap samples available.
#' - `boot.R`: Integer, number of bootstrap samples used.
#' - `boot.l`: Integer, block length in the time-series bootstrap process.
#' - `seed`: Integer, random number generator seed used in bootstrap.
#' - `sim`: Character, `sim` argument of \link{boot::tsboot}.
#' - `cf0`: Numeric vector, mean value of original measurements.
#' - `cf.tsboot`: List, result from the \link{boot::tsboot} function.
#' - `tsboot.se`: Numeric vector, standard deviation over bootstrap samples.
#' - `jackknife.samples`: Logical, indicating whether there are jackknife samples available.
#' - `cf.jackknife`: List, containing jackknife samples:
#'   - `t`: Numeric matrix, jackknifed data sets.
#'   - `t0`: Numeric vector, copy of `cf0`.
#'   - `l`: Integer, copy of `boot.l`.
#' - `jackknife.se`: Numeric vector, standard error over jackknife samples.
#' - `nrObs`: Integer, number of different measurements contained in this correlation function. One can use \link{c.cf} to add multiple observables into one container. This is for instance needed when passing to the \link{gevp} function.
#' - `Time`: Integer, full time extent.
#' - `T`: Integer, full time extent.
#' - `nrStypes`: Integer, number of smearing types.
#' - `symmetrised`: Logical, indicating whether the correlation function has been symmetrized.
#'
#' The \link{gevp2cf} function might add the following fields:
#'
#' - `weighted`: TODO
#' - `weight.cosh`: TODO
#' - `mass1`: TODO
#' - `mass2`: TODO
#'
#' The \link{subtract.excitedstates} function adds the following fields:
#'
#' - `subtracted.values`: Numeric matrix, TODO
#' - `subtracted.ii`: Integer vector, TODO
cf <- function() {
  cf <- list(cf=NULL)
  attr(cf, "class") <- c("cf", class(cf))
  return(cf)
}

gen.block.array <- function(n, R, l, endcorr=TRUE) {
  endpt <- if (endcorr)
             n
           else n - l + 1
  nn <- ceiling(n/l)
  lens <- c(rep(l, nn - 1), 1 + (n - 1)%%l)
  st <- matrix(sample.int(endpt, nn * R, replace = TRUE),
               R)
  return(list(starts = st, lengths = lens))
}

bootstrap.cf <- function(cf, boot.R=400, boot.l=2, seed=1234, sim="geom", endcorr=TRUE) {
    
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.cf requires an object of class 'cf' as input! Aborting!\n")
  }
  boot.l <- ceiling(boot.l)
  if(boot.l < 1 || boot.l > nrow(cf$cf)) {
    stop("'boot.l' must be larger than 0 and smaller than the length of the time series! Aborting...\n")
  }
  boot.R <- floor(boot.R)
  if(boot.R < 1) {
    stop("'boot.R' must be positive!")
  }
  ## save random number generator state
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    temp <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else temp <- NULL

  cf$boot.samples <- TRUE
  cf$boot.R <- boot.R
  cf$boot.l <- boot.l
  cf$seed <- seed
  cf$sim <- sim
  cf$cf0 <- apply(cf$cf, MARGIN=2L, FUN=mean)
  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## now we bootstrap the correlators
  cf$cf.tsboot <- tsboot(cf$cf, statistic = function(x){ return(apply(x, MARGIN=2L, FUN=mean))},
                         R = boot.R, l=boot.l, sim=sim, endcorr=endcorr)
  ## the bootstrap error
  cf$tsboot.se <- apply(cf$cf.tsboot$t, MARGIN=2L, FUN=sd)
  ## restore random number generator state
  if (!is.null(temp))
    assign(".Random.seed", temp, envir = .GlobalEnv)
  else rm(.Random.seed, pos = 1)
  
  return(invisible(cf))
}

jackknife.cf <- function(cf, boot.l=2) {
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.cf requires an object of class cf as input! Aborting!\n")
  }
  if(boot.l < 1) {
    stop("boot.l must be larger than 0! Aborting...\n")
  }
  boot.l <- ceiling(boot.l)
  cf$jackknife.samples <- TRUE
  cf$boot.l <- boot.l
  ## blocking with fixed block length, but overlapping blocks
  ## number of observations
  n <- nrow(cf$cf)
  ## number of overlapping blocks
  N <- n-boot.l+1
  cf$cf0 <- apply(cf$cf, MARGIN=2L, FUN=mean)
  
  cf$cf.jackknife <- list()
  cf$cf.jackknife$t<- array(NA, dim=c(N,ncol(cf$cf)))
  cf$cf.jackknife$t0 <- cf$cf0
  cf$cf.jackknife$l <- boot.l
  for (i in 1:N) {
    ii <- c(i:(i+boot.l-1))
    ## jackknife replications of the mean
    gammai <- apply(cf$cf[-ii,], MARGIN=2L, FUN=mean)
    cf$cf.jackknife$t[i, ] <- (n*cf$cf0 - (n - boot.l)*gammai)/boot.l
  }
  ## the jackknife error
  tmp <- apply(cf$cf.jackknife$t, MARGIN=1L, FUN=function(x,y){(x-y)^2}, y=cf$cf0)
  cf$jackknife.se <- apply(tmp, MARGIN=1L,
                           FUN=function(x, l, n, N) {sqrt( l/(n-l)/N*sum( x ) ) },
                           n=n, N=N, l=boot.l)
  return(invisible(cf))
}

# Gamma method analysis on all time-slices in a 'cf' object
uwerr.cf <- function(cf, absval=FALSE){
  if(!inherits(cf, "cf")){
    stop("uwerr.cf: cf must be of class 'cf'. Aborting...\n")
  }
  uwcf <- as.data.frame( 
      t(
          apply(X=cf$cf, MARGIN=2L, 
                FUN=function(x){
                  data <- x
                  if(absval) data <- abs(x)
                  uw <- try(uwerrprimary(data=data), silent=TRUE)
                  if(any( class(uw) == 'try-error' ) ){
                    c(value=NA,
                      dvalue=NA,
                      ddvalue=NA,
                      tauint=NA,
                      dtauint=NA)
                  } else {
                    c(value=uw$value, 
                      dvalue=uw$dvalue, 
                      ddvalue=uw$ddvalue,
                      tauint=uw$tauint, 
                      dtauint=uw$dtauint)
                  }
                }
                )
      ) 
  )
  uwcf <- cbind(t=(1:ncol(cf$cf))-1,uwcf)
  return(uwcf)
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
    return(invisible(cf))
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

# "close-by-times" averaging replaces the value of the correlation function at t
# with the "hypercubic" average with the values at the neighbouring time-slices 
# with weights 0.25, 0.5 and 0.25
# it then invalidates the boundary timeslices (for all smearing types and observables)
avg.cbt.cf <- function(cf){
  if(!any(class(cf) == "cf")) {
    stop("avg.cbt.cf: Input must be of class 'cf'\n")
  }
  # copy for shifting
  cf2 <- cf
  cf <- mul.cf(cf, 0.5)
  
  # average over shifted correlation functions
  for( p in c(-1,1) ){
    cf <- cf + mul.cf(shift.cf(cf2,p),0.25)
  }
  # invalidate time slices with incorrect contributions
  for( oidx in 0:(cf$nrObs-1) ){
    for( sidx in 0:(cf$nrStypes-1) ){
      nts <- cf$Time/2+1
      if( "symmetrised" %in% names(cf) ) {
        if(!cf$symmetrised){
          nts <- cf$Time
        }
      }
      istart <- oidx*cf$nrStypes*nts + sidx*nts + 1
      iend <- istart+nts
      ii <- c(istart,iend-1)
      cf$cf[,ii] <- NA
    }
  }
  cf <- invalidate.samples.cf(cf)
  return(invisible(cf))
}

## this is intended for instance for adding diconnected diagrams to connected ones
add.cf <- function(cf1, cf2, a=1., b=1.) {
  if(any(class(cf1) == "cf") && any(class(cf2) == "cf") &&
     all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- a*cf1$cf + b*cf2$cf
    cf <- invalidate.samples.cf(cf)
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
    cf <- invalidate.samples.cf(cf)

    return(cf)
  }
}

'-.cf' <- function(cf1, cf2) {
  if(all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- cf1$cf - cf2$cf
    cf <- invalidate.samples.cf(cf)
    return(cf)
  }
}

'/.cf' <- function(cf1, cf2) {
  if(all(dim(cf1$cf) == dim(cf2$cf)) && cf1$Time == cf2$Time ) {
    cf <- cf1
    cf$cf <- cf1$cf / cf2$cf
    cf <- invalidate.samples.cf(cf)
    return(cf)
  }
}

mul.cf <- function(cf, a=1.) {
  if(any(class(cf) == "cf") && is.numeric(a)) {
    cf$cf <- a*cf$cf
    cf <- invalidate.samples.cf(cf)
    return(cf)
  }
  else {
    stop("Wrong classes for input objects, must be cf and numeric. Aborting...!\n")
  }
}

extractSingleCor.cf <- function(cf, id=c(1)) {
  if(!inherits(cf, "cf")) {
    stop("extractSingleCor.cf: cf must be of class 'cf'. Aborting...\n")
  }
  
  ii <- c()
  for(i in c(1:length(id))) {
    ii <- c(ii, c(1:(cf$Time/2+1)) + (id[i]-1)*(cf$Time/2+1))
  }

  cf$cf <- cf$cf[,ii]
  if(cf$boot.samples) {
    cf$cf0 <- cf$cf0[ii]
    cf$tsboot.se <- cf$tsboot.se[ii]
    cf$cf.tsboot$t0 <- cf$cf.tsboot$t0[ii]
    cf$cf.tsboot$t <- cf$cf.tsboot$t[,ii]
    cf$cf.tsboot$data <- cf$cf.tsboot$data[,ii]
  }
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
        cf$icf <- cbind(cf$icf, fcall[[i]]$icf)
      }
    }
  }
  cf$boot.samples <- FALSE
  return(invisible(cf))
}

plot.cf <- function(cf, boot.R=400, boot.l=2, neg.vec, rep=FALSE, ...) {
  if(is.null(cf$jackknife.samples)) {
    cf$jackknife.samples <- FALSE
  }
  if(is.null(cf$boot.samples)) {
    cf$boot.samples <- FALSE
  }
  if(!cf$boot.samples && !cf$jackknife.samples) {
    cf <- bootstrap.cf(cf, boot.R, boot.l)
  }
  Err <- numeric(0)
  if(cf$boot.samples) Err <- cf$tsboot.se
  else if(cf$jackknife.samples) Err <- cf$jackknife.se
  if(missing(neg.vec)){
    neg.vec <- rep(1,times=length(cf$cf0))
  }

  tmax <- cf$Time/2
  if( "symmetrised" %in% names(cf) ) {
    if(!cf$symmetrised){
      tmax <- cf$Time-1
    }
  }
  plotwitherror(x=rep(c(0:(tmax)), times=length(cf$cf0)/(tmax+1)), y=neg.vec*cf$cf0, dy=Err, rep=rep, ...)
  return(invisible(data.frame(t=rep(c(0:tmax), times=length(cf$cf0)/(tmax+1)), CF=cf$cf0, Err=Err)))
}

# shift a correlation function by 'places' time-slices
#   C'(t) = C(t+places)
# where places can be positive or negative as required
# this will of course mix smearings and observables
# and must be taken into account externally by
# invalidating the affected time-slices
shift.cf <- function(cf,places){
  if(!any(class(cf) == "cf")) {
    stop(".cf: Input must be of class 'cf'\n")
  }
  cf <- invalidate.samples.cf(cf)
  n <- ncol(cf$cf)
  if(places == 0){
    cf$cf <- cf$cf
  } else if ( places < 0 ){
    cf$cf <- cbind( cf$cf[, (n - abs(places) + 1):n], cf$cf[, 1:(n-abs(places))] )
  } else {
    cf$cf <- cbind( cf$cf[, (places+1):n], cf$cf[, 1:places] )
  }
  return(invisible(cf))
}

# when a correlation function is modified, any resampling should be
# invalidated
invalidate.samples.cf <- function(cf){
  cf$cf0 <- NULL
  if(cf$boot.samples) {
    cf$boot.samples <- FALSE
    cf$boot.R <- NULL
    cf$boot.l <- NULL
    cf$sim <- NULL
    cf$tsboot.se <- NULL
    cf$tsboot <- NULL
  }
  if(cf$jackknife.samples){
    cf$jackknife.samples <- FALSE
    cf$jackknife <- NULL
    cf$jackknife.se <- NULL
  }
  return(invisible(cf))
}

symmetrise.cf <- function(cf, sym.vec=c(1) ) {
  if( "symmetrised" %in% names(cf) ) {
    if(cf$symmetrised){
      message("symmetrise.cf: cf was already symmetrised\n")
      return(invisible(cf))
    }
  }
  if( cf$nrObs > 1 & length(sym.vec) == 1 ){
    sym.vec <- rep(sym.vec[1],times=cf$nrObs)
  } else if( cf$nrObs != length(sym.vec) ) {
    stop("symmetrise.cf: length of sym.vec must either be 1 or match cf$nrObs!\n")
  }

  Thalf <- cf$Time/2
  isub <- c()
  for( oidx in 0:(cf$nrObs-1) ){
    for( sidx in 0:(cf$nrStypes-1) ){
      istart <- oidx*cf$nrStypes*cf$Time + cf$Time*sidx + 1
      ihalf <- istart + Thalf
      iend <- istart + cf$Time - 1
      isub <- c(isub,(ihalf+1):iend)
      cf$cf[, (istart+1):(ihalf-1)] <- 0.5*( cf$cf[, (istart+1):(ihalf-1)] +
                                             sym.vec[oidx+1]*cf$cf[, rev((ihalf+1):iend)] )
      if( !is.null(cf$icf) ){
        cf$icf[, (istart+1):(ihalf-1)] <- 0.5*( cf$icf[, (istart+1):(ihalf-1)] + 
                                                sym.vec[oidx+1]*cf$icf[, rev((ihalf+1):iend)] )
      }
    }
  }
  # remove now unnecessary time slices 
  cf$cf <- cf$cf[, -isub]
  if( !is.null(cf$icf) ){
    cf$icf <- cf$icf[, -isub]
  }
  cf$symmetrised <- TRUE
  return(invisible(cf))
}


summary.cf <- function(cf, ...) {
  cat("T = ", cf$Time, "\n")
  cat("observations = ", dim(cf$cf)[1], "\n")
  cat("Nr Stypes = ", cf$nrStypes, "\n")
  cat("Nr Obs    = ", cf$nrObs, "\n")
  if(cf$boot.samples) {
    cat("R = ", cf$boot.R, "\n")
  }
  if(cf$boot.samples || !is.null(cf$jackknife.se)) {
    tmax <- cf$Time/2
    if( "symmetrised" %in% names(cf) ) {
      if(!cf$symmetrised){
        tmax <- cf$Time-1
      }
    }
    cat("l = ", cf$boot.l, "\n")
    out <- data.frame(t=c(0:tmax), C=cf$cf0)
  }
  if(!is.null(cf$sim)) {
    cat("sim = ", cf$sim, "\n")
  }

  if(!is.null(cf$tsboot.se)) {
    out <- cbind(out, tsboot.se=cf$tsboot.se)
  }
  if(!is.null(cf$jackknife.se)) {
    out <- cbind(out, jackknife.se=cf$jackknife.se)
  }
  if(!is.null(cf$jack.boot.se)) {
    out <- cbind(out, jab.se=cf$jack.boot.se)
  }
  if(exists("out")) {
    print(out)
  }
}

print.cf <- function(cf, ...) {
  summary(cf, ...)
}
