#' Correlation function container
#'
#' This function `cf()` creates containers for correlation functions
#' of class `cf`. This class is particularly designed to deal with
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.cf} and \link{as.cf}.
#'
#' @param nrObs Integer, number of different measurements contained in this correlation function. One can use \link{c.cf} to add multiple observables into one container. This is for instance needed when passing to the \link{gevp} function.
#' @param Time Integer, full time extent.
#' @param nrStypes Integer, number of smearing types.
#' @param symmetrised Logical, indicating whether the correlation function has been symmetrized.
#'
#' @details
#'
#' And last but not least, these are the fields that are used somewhere in the library but we have not figured out which mixin these should belong to:
#'
#' - `conf.index`: TODO
#' - `N`: Integer, number of measurements.
#' - `blockind`: TODO
#' - `jack.boot.se`: TODO
#'
#' @family cf constructors
#'
#' @export
cf <- function (nrObs = 1, Time = NA, nrStypes = 1, symmetrised = FALSE) {
  cf <- list(nrObs = nrObs, Time = Time, nrStypes = nrStypes,
             symmetrised = symmetrised)

  class(cf) <- append(class(cf), 'cf')
  return (cf)
}

#' Bootstrapped CF mixin constructor
#'
#' @param cf `cf` object to extend.
#' @param boot.R Integer, number of bootstrap samples used.
#' @param boot.l Integer, block length in the time-series bootstrap process.
#' @param seed Integer, random number generator seed used in bootstrap.
#' @param sim Character, `sim` argument of \link{boot::tsboot}.
#' @param cf.tsboot List, result from the \link{boot::tsboot} function.
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `cf0`: Numeric vector, mean value of original measurements, convenience copy of `cf.tsboot$t0`.
#' - `tsboot.se`: Numeric vector, standard deviation over bootstrap samples.
#' - `boot.samples`: Logical, indicating whether there are bootstrap samples available. This is deprecated and instead the presence of bootstrap samples should be queried with `inherits(cf, 'cf_boot')`.
#'
#' @family cf constructors
#'
#' @export
cf_boot <- function (cf, boot.R, boot.l, seed, sim, cf.tsboot) {
  cf$boot.R <- boot.R
  cf$boot.l <- boot.l
  cf$seed <- seed
  cf$sim <- sim
  cf$cf.tsboot <- cf.tsboot

  cf$cf0 <- cf.tsboot$t0
  cf$tsboot.se <- apply(cf$cf.tsboot$t, MARGIN = 2L, FUN = sd)
  cf$boot.samples <- TRUE

  class(cf) <- append(class(cf), 'cf_boot')
  return (cf)
}

#' Jackknifed CF mixin constructor
#'
#' With the `cf_jackknife` mixin, it also must have the following fields:
#'
#' @param cf `cf` object to extend.
#' @param cf0 _See above_.
#' @param boot.l _See above_.
#' @param cf.jackknife List, containing jackknife samples:
#'   - `t`: Numeric matrix, jackknifed data sets.
#'   - `t0`: Numeric vector, copy of `cf0`.
#'   - `l`: Integer, copy of `boot.l`.
#' @param jackknife.se Numeric vector, standard error over jackknife samples.
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `jackknife.samples`: Logical, indicating whether there are jackknife samples available. This is deprecated and instead the presence of jackknife samples should be queried with `inherits(cf, 'cf_jackknife')`.
#'
#' @family cf constructors
#'
#' @export
cf_jackknife <- function (cf, cf0, boot.l, cf.jackknife, jackknife.se) {
  cf$cf0 <- cf0
  cf$boot.l <- boot.l
  cf$cf.jackknife <- cf.jackknife
  cf$jackknife.se <- jackknife.se

  cf$jackknife.samples <- TRUE

  class(cf) <- append(class(cf), 'cf_jackknife')
  return (cf)
}

#' Original data CF mixin constructor
#'
#' @param .cf `cf` object to extend. Named with a leading period just to distinguish it from the member also named `cf`.
#' @param cf Numeric matrix, original data for all observables and measurements.
#' @param icf Numeric matrix, imaginary part of original data. Be very careful with this as most functions just ignore the imaginary part and drop it in operations. If it is not passed to this function, a matrix of `NA` will be created with the same dimension as `cf`.
#'
#' @family cf constructors
#'
#' @export
cf_orig <- function (.cf, cf, icf = NULL) {
  .cf$cf <- cf

  if (is.null(icf)) {
    .cf$icf <- cf
    .cf$icf[, ] <- NA
  }
  else {
    .cf$icf <- icf
  }

  class(.cf) <- append(class(.cf), 'cf_orig')
  return (.cf)
}

#' Principal correlator CF mixin constructor
#'
#' @param cf `cf` object to extend.
#' @param id Integer, number of the principal correlator from the GEVP. Ascending with eigenvalue, so `id = 1` is the lowest state.
#'
#' @family cf constructors
#'
#' @export
cf_principal_correlator <- function (cf, id) {
  cf$id <- id

  class(cf) <- append(class(cf), 'cf_principal_correlators')
  return (cf)
}

cf_shifted <- function (cf, deltat, forwardshift) {
  cf$deltat <- deltat
  cf$forwardshift <- forwardshift

  cf$shifted <- TRUE

  class(cf) <- append(class(cf), 'cf_shifted')
  return (cf)
}

#' Smeared CF mixin constructor
#'
#' @param scf Like `cf`, but with the smeared data.
#' @param iscf Like `icf`, but with the smeared data.
#' @param nrSamples TODO
#' @param obs TODO
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `smeared`: Logical, whether the correlation function has smeared data. This is deprecated and instead the presence of bootstrap samples should be queried with `inherits(cf, 'cf_smeared')`.
#'
#' @family cf constructors
#'
#' @export
cf_smeared <- function (cf, scf, iscf, nrSamples, obs) {
  cf$scf <- scf
  cf$iscf <- iscf
  cf$nrSamples <- nrSamples
  cf$obs <- obs

  cf$smeared <- TRUE

  class(cf) <- append(class(cf), 'cf_smeared')
  return (cf)
}

#' Subtracted CF mixin constructor
#'
#' @param cf `cf` object to extend.
#' @param subtracted.values Numeric matrix, TODO
#' @param subtracted.ii Integer vector, TODO
#'
#' @family cf constructors
#'
#' @export
cf_subtracted <- function (cf, subtracted.values, subtracted.ii) {
  cf$subtracted.value <- subtracted.values
  cf$subtracted.ii <- subtracted.ii

  class(cf) <- append(class(cf), 'cf_subtracted')
  return (cf)
}

#' Weighted CF mixin constructor
#'
#' @param cf `cf` object to extend.
#' @param weight.factor TODO
#' @param weight.cosh TODO
#' @param mass1 TODO
#' @param mass2 TODO
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `weighted`: Logical, indicating whether the correlation function has been weighted. This is deprecated and instead the presence of this should be queried with `inherits(cf, 'cf_weighted')`.
#'
#' @family cf constructors
#'
#' @export
cf_weighted <- function (cf, weight.factor, weight.cosh, mass1, mass2) {
  cf$weight.factor <- weight.factor
  cf$weight.cosh <- weight.cosh
  cf$mass1 <- mass1
  cf$mass2 <- mass2

  cf$weighted <- TRUE

  class(cf) <- append(class(cf), 'cf_weighted')
  return (cf)
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
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

  boot.l <- ceiling(boot.l)
  boot.R <- floor(boot.R)

  stopifnot(boot.l >= 1)
  stopifnot(boot.l <= nrow(cf$cf))
  stopifnot(boot.R >= 1)

  ## save random number generator state
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    temp <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else
    temp <- NULL

  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## now we bootstrap the correlators
  cf.tsboot <- tsboot(cf$cf, statistic = function(x){ return(apply(x, MARGIN=2L, FUN=mean))},
                         R = boot.R, l=boot.l, sim=sim, endcorr=endcorr)

  cf <- cf_boot(cf,
                boot.R = boot.R,
                boot.l = boot.l,
                seed = seed,
                sim = sim,
                cf.tsboot = cf.tsboot)

  ## restore random number generator state
  if (!is.null(temp))
    assign(".Random.seed", temp, envir = .GlobalEnv)
  else rm(.Random.seed, pos = 1)

  return(invisible(cf))
}

jackknife.cf <- function(cf, boot.l=2) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

  stopifnot(boot.l >= 1)
  boot.l <- ceiling(boot.l)

  ## blocking with fixed block length, but overlapping blocks
  ## number of observations
  n <- nrow(cf$cf)
  ## number of overlapping blocks
  N <- n-boot.l+1

  cf.jackknife <- list()
  cf.jackknife$t<- array(NA, dim=c(N,ncol(cf$cf)))
  cf.jackknife$t0 <- cf$cf0
  cf.jackknife$l <- boot.l
  for (i in 1:N) {
    ii <- c(i:(i+boot.l-1))
    ## jackknife replications of the mean
    gammai <- apply(cf$cf[-ii,], MARGIN=2L, FUN=mean)
    cf$cf.jackknife$t[i, ] <- (n*cf$cf0 - (n - boot.l)*gammai)/boot.l
  }
  ## the jackknife error
  tmp <- apply(cf$cf.jackknife$t, MARGIN=1L, FUN=function(x,y){(x-y)^2}, y=cf$cf0)
  jackknife.se <- apply(tmp, MARGIN=1L,
                           FUN=function(x, l, n, N) {sqrt( l/(n-l)/N*sum( x ) ) },
                           n=n, N=N, l=boot.l)

  cf <- cf_jackknife(cf,
                     cf0 = apply(cf$cf, MARGIN = 2, FUN = mean),
                     boot.l = boot.l,
                     cf.jackknife = cf.jackknife,
                     jackknife.se = jackknife.se)

  return (invisible(cf))
}

# Gamma method analysis on all time-slices in a 'cf' object
uwerr.cf <- function(cf, absval=FALSE){
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

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
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(cf1$Time == cf2$Time)
  stopifnot(dim(cf1$cf)[2] == dim(cf2$cf)[2])
  stopifnot(cf1$nrObs == cf2$nrObs )
  stopifnot(cf1$nrStypes == cf2$nrStypes)

  cf <- cf1

  cf$cf <- rbind(cf1$cf, cf2$cf)
  cf$icf <- rbind(cf1$icf, cf2$icf)

  cf <- invalidate.samples.cf(cf)

  return (invisible(cf))
}

## averages local-smeared and smeared-local correlators in cf and adjusts
## nrStypes accordingly
## by default, assumes that LS and SL are in columns (T/2+1)+1:3*(T/2+1)
avg.ls.cf <- function(cf, cols = c(2, 3)) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(cf$nrStypes >= 2)

  timeslices <- cf$Time/2+1

  ind.ls <- ( (cols[1]-1)*timeslices+1 ):( cols[1]*timeslices )
  ind.sl <- ( (cols[2]-1)*timeslices+1 ):( cols[2]*timeslices )

  cf$cf[,ind.ls] <- 0.5 * ( cf$cf[,ind.ls] + cf$cf[,ind.sl] )

  cf$cf <- cf$cf[,-ind.sl]
  cf$nrStypes <- cf$nrStypes-1
  return (cf)
}

# "close-by-times" averaging replaces the value of the correlation function at t
# with the "hypercubic" average with the values at the neighbouring time-slices
# with weights 0.25, 0.5 and 0.25
# it then invalidates the boundary timeslices (for all smearing types and observables)
avg.cbt.cf <- function(cf){
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

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
add.cf <- function(cf1, cf2, a=1.0, b=1.0) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- a*cf1$cf + b*cf2$cf
  cf <- invalidate.samples.cf(cf)
  return(cf)
}

'+.cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- cf1$cf + cf2$cf
  cf <- invalidate.samples.cf(cf)

  return(cf)
}

'-.cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- cf1$cf - cf2$cf
  cf <- invalidate.samples.cf(cf)
  return(cf)
}

'/.cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- cf1$cf / cf2$cf
  cf <- invalidate.samples.cf(cf)
  return (cf)
}

mul.cf <- function(cf, a=1.) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(is.numeric(a))

  cf$cf <- a*cf$cf
  cf <- invalidate.samples.cf(cf)
  return (cf)
}

extractSingleCor.cf <- function(cf, id=c(1)) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

  ii <- c()
  for(i in c(1:length(id))) {
    ii <- c(ii, c(1:(cf$Time/2+1)) + (id[i]-1)*(cf$Time/2+1))
  }

  cf$cf <- cf$cf[,ii]

  if (inherits(cf, 'cf_boot')) {
    cf$cf0 <- cf$cf0[ii]
    cf$tsboot.se <- cf$tsboot.se[ii]
    cf$cf.tsboot$t0 <- cf$cf.tsboot$t0[ii]
    cf$cf.tsboot$t <- cf$cf.tsboot$t[,ii]
    cf$cf.tsboot$data <- cf$cf.tsboot$data[,ii]
  }
  cf$nrObs <- 1
  cf$nsStypes <- 1
  return (cf)
}

is.cf <- function(x){
  inherits(x, "cf")
}

#' Concatenate correlation function objects
#'
#' @param ... One or multiple objects of type `cf_orig`.
c.cf <- function(...) {
  fcall <- list(...)

  # In case there is only one element, we do not need to do anything.
  if(length(fcall) == 1) {
    return(eval(fcall[[1]]))
  }

  # All arguments must be of type `cf` and `cf_orig` since we want to work with
  # the actual data.
  stopifnot(all(sapply(fcall, function (x) inherits(x, 'cf'))))
  stopifnot(all(sapply(fcall, function (x) inherits(x, 'cf_orig'))))

  cf <- fcall[[1]]
  Time <- cf$Time
  cf$nrObs <- 0
  cf$sTypes <- 0
  N <- dim(cf$cf)[1]
  for (i in 1:length(fcall)) {
    if (fcall[[i]]$Time != Time) {
      stop("Times must agree for different objects of type cf\n Aborting\n")
    }
    if (dim(fcall[[i]]$cf)[1] != N) {
      stop("Number of measurements must agree for different objects of type cf\n Aborting\n")
    }
    cf$nrObs <- cf$nrObs + fcall[[i]]$nrObs
    cf$sTypes <- cf$sTypes + fcall[[i]]$sTypes
  }
  if (1 < length(fcall)) {
    for (i in 2:length(fcall)) {
      cf$cf <- cbind(cf$cf, fcall[[i]]$cf)
      cf$icf <- cbind(cf$icf, fcall[[i]]$icf)
    }
  }
  cf <- invalidate.samples.cf(cf)

  return (invisible(cf))
}

plot.cf <- function(cf, neg.vec = rep(1, times = length(cf$cf0)), rep = FALSE, ...) {
  # We need to have data, otherwise we cannot plot.
  stopifnot(any(inherits(cf, c('cf_orig', 'cf_boot', 'cf_jackknife'))))

  if (inherits(cf, 'cf_boot')) {
    val <- cf$cf0
    err <- cf$tsboot.se
  } else if (inherits(cf, 'cf_jackknife')) {
    val <- cf$cf0
    err <- cf$jackknife.se
  } else {
    val <- apply(cf$cf, 2, mean)
    err <- apply(cf$cf, 2, sd) / sqrt(nrow(cf$cf))
  }

  if(!cf$symmetrised){
    tmax <- cf$Time - 1
  } else {
    tmax <- cf$Time / 2
  }

  df <- data.frame(t = rep(c(0:tmax), times = length(val)/(tmax+1)),
                   CF = val,
                   Err = err)

  plotwitherror(x = df$t, y = neg.vec * df$CF, dy = df$Err, rep = rep, ...)

  return(invisible(df))
}

# shift a correlation function by 'places' time-slices
#   C'(t) = C(t+places)
# where places can be positive or negative as required
# this will of course mix smearings and observables
# and must be taken into account externally by
# invalidating the affected time-slices
shift.cf <- function(cf, places) {
  stopifnot(all(c('cf', 'cf_orig') %in% class(cf)))

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

#' Invalidate samples
#'
#' When a correlation function is modified, any resampling should be invalidated. We could instead also choose to properly work with the samples, but most computations are done with the original data anyway.
invalidate.samples.cf <- function(cf){
  cf$boot.l <- NULL
  cf$boot.R <- NULL
  cf$boot.samples <- NULL
  cf$cf0 <- NULL
  cf$jackknife <- NULL
  cf$jackknife.samples <- NULL
  cf$jackknife.se <- NULL
  cf$seed <- NULL
  cf$sim <- NULL
  cf$tsboot <- NULL
  cf$tsboot.se <- NULL

  class(cf) <- setdiff(class(cf), c('cf_boot', 'cf_jackknife'))

  return(invisible(cf))
}

symmetrise.cf <- function(cf, sym.vec=c(1) ) {
  stopifnot(inherits(cf, 'cf'))
  stopifnot(inherits(cf, 'cf_orig'))

  if(cf$symmetrised){
    message("symmetrise.cf: cf was already symmetrised\n")
    return(invisible(cf))
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

  if (inherits(cf, 'cf_boot')) {
    cat("R = ", cf$boot.R, "\n")

    if(!cf$symmetrised){
      tmax <- cf$Time-1
    } else {
      tmax <- cf$Time/2
    }
    cat("l = ", cf$boot.l, "\n")
    out <- data.frame(t=c(0:tmax), C=cf$cf0)
    cat("sim = ", cf$sim, "\n")

    out <- cbind(out, tsboot.se=cf$tsboot.se)
  }


  if (inherits(cf, 'cf_jackknife')) {
    out <- cbind(out, jackknife.se=cf$jackknife.se)
    out <- cbind(out, jab.se=cf$jack.boot.se)
  }

  if(exists("out")) {
    print(out)
  }
}

print.cf <- function(cf, ...) {
  summary(cf, ...)
}
