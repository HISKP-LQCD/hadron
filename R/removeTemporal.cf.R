#' Remove temporal states
#'
#' Performs weighting and shifting in the rest and moving frames.
#'
#' @param cf Object of type `cf`, two-to-two particle correlation function which shall be weighted and shifted. It must be a correlation function in the frame \eqn{p_1 + p_2}.
#' @param single.cf1,single.cf2 Object of type `effectivemassfit` or `matrixfit` which contains the one particle mass in the rest frame.
#'
#' If `single.cf2` is missing, then the mass given as `single.cf1` is used as well. This is sensibly done when one scatters identical particles. But be careful: Even when `single.cf2` is missing, the `p2` is _not_ automatically copied from `p1`.
#'
#' In case `single.cf1` is missing, no weighting is performed. Instead it is assumed that the user only wants to have a simple shifting. Then this function just calls `takeTimeDiff.cf`.
#' @param p1,p2 Integer vector with three elements, containing the momenta that the one particle mass should be boosted to.
#' @param L Integer, spatial extent of the lattice.
#' @param lat.disp Logical, true when the lattice dispersion relation shall be used, otherwise continuum dispersion relation.
#' @param weight.cosh Logical, whether to use some cosh formula which might work better or something.
#'
#' @export
removeTemporal.cf <- function(cf, single.cf1, single.cf2,
                              p1=c(0,0,0), p2=c(0,0,0), L,
                              lat.disp=TRUE, weight.cosh=FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  if(missing(single.cf1)) {
    ## take the differences of C(t+1) and C(t)
    cf <- takeTimeDiff.cf(cf)
    return(invisible(cf))
  }

  Time <- cf$Time
  if(missing(L)) {
    L <- cf$Time/2
    warning("L was missing, set it to T/2\n")
  }
  mass1 <- list()
  mass2 <- list()

  if(missing(single.cf2)) {
    single.cf2 <- single.cf1
  }

  if(cf$boot.R != single.cf1$boot.R || single.cf1$boot.R != single.cf2$boot.R ||
     cf$boot.l != single.cf1$boot.l || single.cf1$boot.l != single.cf2$boot.l) {
##     cf$seed != single.cf1$seed || single.cf1$seed != single.cf2$seed) {
    stop("please provide equally bootstrapped cfs to removeTemporal.cf using the same configurations and seed for all and the same boot.R and boot.l!\n")
  }

  if(inherits(single.cf1, "effectivemassfit")) {
    mass1$t0 <- single.cf1$opt.res$par[1]
    mass1$t <- single.cf1$massfit.tsboot[,1]
  }
  else if(inherits(single.cf1, "matrixfit")) {
    mass1$t0 <- single.cf1$opt.res$par[1]
    mass1$t <- single.cf1$opt.tsboot[1,]
  }
  else {
    stop("single.cf1 must be either of class effectivemassfit or matrixfit! Aborting...\n")
  }
  if(inherits(single.cf2, "effectivemassfit")) {
    mass2$t0 <- single.cf2$opt.res$par[1]
    mass2$t <- single.cf2$massfit.tsboot[,1]
  }
  else if(inherits(single.cf1, "matrixfit")) {
    mass2$t0 <- single.cf2$opt.res$par[1]
    mass2$t <- single.cf2$opt.tsboot[1,]
  }
  else {
    stop("single.cf2 must be either of class effectivemassfit or matrixfit! Aborting...\n")
  }
  ## use momenta p1 and p2 and lattice or continuum dispersion relation to shift energies
  if(any(p1 != 0)) {
    if(lat.disp) {
      pshift <- 2*sum(sin(pi*p1/L)^2)
      mass1$t0 <- acosh( cosh(mass1$t0) + pshift )
      mass1$t <- acosh( cosh(mass1$t) + pshift )
    }
    else {
      pshift <- sum((2*pi*p1/L)^2)
      mass1$t0 <- sqrt( mass1$t0^2 + pshift )
      mass1$t <- sqrt( mass1$t^2 + pshift )
    }
  }
  if(any(p2 != 0)) {
    if(lat.disp) {
      pshift <- 2*sum(sin(pi*p2/L)^2)
      mass2$t0 <- acosh( cosh(mass2$t0) + pshift )
      mass2$t <- acosh( cosh(mass2$t) + pshift )
    }
    else {
      pshift <- sum((2*pi*p2/L)^2)
      mass2$t0 <- sqrt( mass2$t0^2 + pshift )
      mass2$t <- sqrt( mass2$t^2 + pshift )
    }
  }
  cosh.factor <- 0.
  if(weight.cosh) cosh.factor  <- 1.
  ## multiply with the exponential correction factor
  Exptt <- exp((mass2$t0-mass1$t0)*c(0:(Time/2))) + cosh.factor *exp((mass2$t0-mass1$t0)*(Time-c(0:(Time/2))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*
      (exp((mass2$t[i]-mass1$t[i])*c(0:(Time/2))) + cosh.factor *exp((mass2$t[i]-mass1$t[i])*(Time-c(0:(Time/2)))))
  }
  ## take the differences of C(t+1) and C(t)
  cf <- takeTimeDiff.cf(cf)

  ## multiply with the exponetial inverse
  Exptt <- exp(-(mass2$t0-mass1$t0)*c(-1:(Time/2-1))) + cosh.factor *exp(-(mass2$t0-mass1$t0)*(Time-c(-1:(Time/2-1))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*
      (exp(-(mass2$t[i]-mass1$t[i])*c(-1:(Time/2-1))) + cosh.factor *exp(-(mass2$t[i]-mass1$t[i])*(Time-c(-1:(Time/2-1)))) )
  }

  # We perform a clean copy of the data now to make sure that all invariants
  # hold and that no new fields have been added that we are not aware of.
  ret <- cf_meta(nrObs = cf$nrObs, Time = cf$Time, nrStypes = cf$nrStypes,
                 symmetrised = cf$symmetrised)
  ret <- cf_orig(ret,
                 cf = cf$cf)
  ret <- cf_boot(ret,
                 boot.R = cf$boot.R,
                 boot.l = cf$boot.l,
                 seed = cf$seed,
                 sim = cf$sim,
                 cf.tsboot = cf$cf.tsboot,
                 resampling_method = cf$resampling_method)
  ret <- cf_shifted(ret,
                    deltat = cf$deltat,
                    forwardshift = cf$forwardshift)
  ret <- cf_weighted(ret,
                     weight.factor = 1 / exp((mass2$t0 - mass1$t0) * 1),
                     weight.cosh = weight.cosh,
                     mass1 = mass1,
                     mass2 = mass2)

  return (invisible(ret))
}

takeTimeDiff.cf <- function (cf, deltat = 1, forwardshift = FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))

  ## number of time slices (hopefully in units of T/2+1 if the correlator has been symmetrised)
  ## and units of the time extent if it has not
  T <- cf$Time
  Nt <- dim(cf$cf)[2]

  nts <- cf$Time/2+1
    if(!cf$symmetrised){
      nts <- cf$Time
    }

  ## the time indices to be subtracted
  tt0 <- c()
  for(i in c(1:cf$nrObs)) {
    tt0 <- c(tt0, ((i-1)*(nts)+1):(i*(nts)-deltat))
  }
  tt1 <- tt0 + deltat

  ## the default is a type of backwards derivative: C'(t+deltat) = C(t) - C(t+deltat)
  ### which invalidates all time slices up to and including deltat-1
  ## alternatively, we can also do a forward derivative: C'(t) = C(t+deltat) - C(t)
  ### which will invalidate all time slices after and including tmax-(deltat-1)
  ## we do this by defining left- and right-hand indices
  # C(tlhs) = C(trhs1) - C(trhs2)
  ## note: in both cases one could in pricinple use the symmetry properties of
  ##       the full correlation function and periodic boundary condtions
  ##       to do this without any invalidation
  tlhs <- tt1
  trhs1 <- tt0
  trhs2 <- tt1
  if( forwardshift ){
    tlhs <- tt0
    trhs1 <- tt1
    trhs2 <- tt0
  }

  # take the differences, set the remaining points to NA. Apparently we don't care about the imaginary part here.
  cf$cf[,tlhs] <- cf$cf[,trhs1]-cf$cf[,trhs2]
  cf$cf[,-tlhs] <- NA

  # now the bootstrap samples
  if (inherits(cf, 'cf_boot')) {
    cf$cf.tsboot$t0[tlhs] <- cf$cf.tsboot$t0[trhs1]-cf$cf.tsboot$t0[trhs2]
    cf$cf.tsboot$t0[-tlhs] <- NA
    cf$cf.tsboot$t[,tlhs] <- cf$cf.tsboot$t[,trhs1]-cf$cf.tsboot$t[,trhs2]
    cf$cf.tsboot$t[,-tlhs] <- NA
  }

  # We perform a new construction in order to have only the fields defined that
  # we want and also to make sure that invariants are holding.
  ret <- cf_meta(nrObs = cf$nrObs, Time = cf$Time, nrStypes = cf$nrStypes,
                 symmetrised = cf$symmetrised)
  ret <- cf_orig(ret,
                 cf = cf$cf)
  if (inherits(cf, 'cf_boot')) {
    ret <- cf_boot(ret,
                   boot.R = cf$boot.R,
                   boot.l = cf$boot.l,
                   seed = cf$seed,
                   sim = cf$sim,
                   cf.tsboot = cf$cf.tsboot,
                   resampling_method = cf$resampling_method)
  }
  ret <- cf_shifted(ret,
                    deltat = deltat,
                    forwardshift = forwardshift)

  return(invisible(ret))
}
