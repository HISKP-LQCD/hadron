#' Remove temporal states
#'
#' Performs weighting and shifting in the rest and moving frames.
#'
#' @param cf Object of type `cf`, two-to-two particle correlation function which
#'   shall be weighted and shifted. It must be a correlation function in the
#'   frame \eqn{p_1 + p_2}.
#' @param single.cf1,single.cf2 Object of type `effectivemassfit` or `matrixfit`
#'   which contains the one particle mass in the rest frame.
#'
#'   If `single.cf2` is missing, then the mass given as `single.cf1` is used as
#'   well. This is sensibly done when one scatters identical particles. But be
#'   careful: Even when `single.cf2` is missing, the `p2` is _not_ automatically
#'   copied from `p1`.
#'
#'   In case `single.cf1` is missing, no weighting is performed. Instead it is
#'   assumed that the user only wants to have a simple shifting. Then this
#'   function just calls `takeTimeDiff.cf`.
#' @param p1,p2 Integer vector with three elements, containing the momenta that
#'   the one particle mass should be boosted to.
#' @param L Integer, spatial extent of the lattice.
#' @param lat.disp Logical, true when the lattice dispersion relation shall be
#'   used, otherwise continuum dispersion relation.
#' @param weight.cosh Logical, If single.cf1 is a pure cosh, the leading two
#'   thermal states also may be expressed as a cosh. If `weight.cosh` is set,
#'   they are removed simultaneously.
#' @param deltat Integer. Time shift value.
#'
#' @return
#' Returns an object of class `cf`, see \link{cf}.
#' 
#' @export
old_removeTemporal.cf <- function(cf, 
                              single.cf1, 
                              single.cf2,
                              p1=c(0,0,0), 
                              p2=c(0,0,0), 
                              L,
                              lat.disp=TRUE, 
                              weight.cosh=FALSE,
                              deltat = 1) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  warning('old_removeTemporal.cf is deprecated')

  if(missing(single.cf1)) {
    ## take the differences of C(t+1) and C(t)
    cf <- takeTimeDiff.cf(cf, deltat)
    return(invisible(cf))
  }

  Time <- cf$Time
  if(missing(L)) {
    L <- cf$Time/2
    warning("L was missing, set it to Time/2\n")
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
  if(weight.cosh) {
    cosh.factor  <- 1.
  }
  else {
    cosh.factor <- 0.
  }
  ## Multiply with the exponential correction factor
  Exptt <- exp((mass2$t0-mass1$t0)*c(0:(Time/2))) + cosh.factor *exp((mass2$t0-mass1$t0)*(Time-c(0:(Time/2))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf/t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0/Exptt
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]/
      (exp((mass2$t[i]-mass1$t[i])*c(0:(Time/2))) + cosh.factor *exp((mass2$t[i]-mass1$t[i])*(Time-c(0:(Time/2)))))
  }
  ## Take the differences C(t) - C(t + deltat)
  cf <- takeTimeDiff.cf(cf, deltat)

  ## Multiply with the exponential inverse
  ## The time has to be shifted by deltat because the first deltat timeslices are filled 
  ## up with NaN by takeTimeDiff()
  Exptt <- exp((mass2$t0-mass1$t0)*c(-deltat:(Time/2-deltat))) + cosh.factor *exp((mass2$t0-mass1$t0)*(Time-c(-deltat:(Time/2-deltat))))
  if(!is.null(cf$cf)) {
    cf$cf <- cf$cf*t(array(Exptt, dim=dim(cf$cf)[c(2,1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0*Exptt
  for(i in c(1:cf$boot.R)) {
    cf$cf.tsboot$t[i,] <- cf$cf.tsboot$t[i,]*
      (exp((mass2$t[i]-mass1$t[i])*c(-deltat:(Time/2-deltat))) + cosh.factor *exp((mass2$t[i]-mass1$t[i])*(Time-c(-deltat:(Time/2-deltat)))) )
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
                 endcorr = cf$endcorr,
                 cf.tsboot = cf$cf.tsboot,
                 resampling_method = cf$resampling_method)
  ret <- cf_shifted(ret,
                    deltat = cf$deltat,
                    forwardshift = cf$forwardshift)
  ret <- cf_weighted(ret,
                     weight.factor = 1 / exp((mass2$t0 - mass1$t0) * deltat),
                     weight.cosh = weight.cosh)

  return (invisible(ret))
}

#' Take time difference
#'
#' @description
#' Performs the calculation of the shifted correlator C_shift(t) = C(t) - C(t +/- deltat).
#'
#' @param cf Object of type `cf`, a particle correlation function which shall be shifted.
#' @param deltat integer. the time shift
#' @param forwardshift boolean. If set to `TRUE`, the forward finite
#'   difference is used instead of the backward one
#'
#' @return
#' The shifted correlator as an object of type `cf`, see \link{cf}
#'
#' @export
takeTimeDiff.cf <- function (cf, deltat = 1, forwardshift = FALSE) {
  stopifnot(inherits(cf, 'cf_meta'))

  ## number of time slices (hopefully in units of Time/2+1 if the correlator has been symmetrised)
  ## and units of the time extent if it has not
  Time <- cf$Time
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

  ## take the differences, set the remaining points to NA. Apparently we don't care about the imaginary part here.
  if( inherits(cf, 'cf_orig') ) {
    cf$cf[,tlhs] <- cf$cf[,trhs1]-cf$cf[,trhs2]
    cf$cf[,-tlhs] <- NA
  }

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
                   endcorr = cf$endcorr,
                   cf.tsboot = cf$cf.tsboot,
                   resampling_method = cf$resampling_method)
  }
  ret <- cf_shifted(ret,
                    deltat = deltat,
                    forwardshift = forwardshift)

  return(invisible(ret))
}

#' Continuum dispersion relation for CM to lattice frame
#'
#' @description
#' Converts a center of mass (CM) frame energy to the lattice frame using the
#' continuum dispersion relation.
#'
#' @param energy `double`. CM energy in lattice units, \eqn{aE}.
#' @param momentum_d `integer`. Total momentum squared of the moving frame in lattice units, \eqn{d^2}.
#' @param extent_space `integer`. Spatial extent of the lattice as a dimensionless quantity, \eqn{L/a}.
#' @param plus Boolean. Sign of a^2 artefacts.
#' @param lattice_disp Boolean. Use the lattice dispersion relation instead of the continuum one
#'
#' @return
#' `double`. Energy in the lattice frame, \eqn{aW}.
#'
#' @export
#' @family dispersion relations
dispersion_relation <- function (energy, momentum_d, extent_space, plus = TRUE, lattice_disp = FALSE) {
  sign <- if (plus) +1 else -1

  if (lattice_disp) {
    p_sq <- 2 * sum(sin(momentum_d * pi / extent_space)^2)
    energy_out <- acosh(cosh(energy) + sign * p_sq)
  } else {
    p_sq <- (2 * pi / extent_space)^2 * sum(momentum_d^2)
    energy_out <- sqrt(energy^2 + sign * p_sq)
  }

  return (energy_out)
}

#' generic function to extract a fitted mass
#'
#' @description
#' One of the main analysis tasks in \link{hadron} is the estimation
#'   of energy levels or masses from correlation functions. The
#'   corresponding analysis functions return objects, typically lists,
#'   containing the masses or energy levels. `extract_mass` is a
#'   generic function to extrac such fitted mass values.
#' 
#' @param object Object to extract the mass from.
#'
#' @return Numeric. The mass value.
#'
#' @export
extract_mass <- function (object) {
  UseMethod('extract_mass')
}

#' specialisation of \link{extract_mass} to objects of type
#' `effectivemassfit`
#' 
#' @param object Object of type `effectivemassfit` to extract the mass from.
#'
#' @return Numeric. The mass value.
#'
#' @export
extract_mass.effectivemassfit <- function (object) {
  list(t0 = object$opt.res$par[1],
       t = object$massfit.tsboot[,1])
}

#' specialisation of \link{extract_mass} to objects of type
#' `matrixfit`
#' 
#' @param object Object of type `matrixfit` to extract the mass from.
#'
#' @return Numeric. The mass value.
#'
#' @export
extract_mass.matrixfit <- function (object) {
  list(t0 = object$opt.res$par[1],
       t = object$opt.tsboot[1,])
}

make_weight_factor <- function (energy_difference, time_extent, time_start,
                                time_end, cosh_factor) {
  time_slices <- time_start:time_end
  exp(energy_difference * time_slices) +
    cosh_factor * exp(energy_difference * (time_extent - time_slices))
}

#' Weight a correlation function
#'
#' @description
#' Weights a correlation function with the given energy difference \eqn{\Delta E}{Delta E}
#' such that the function is first multiplied with
#' \eqn{\exp(\Delta E t) + c \exp(\Delta E \cdot (Time - t)}{exp(Delta E t) + c exp(Delta E(Time-t))}.
#'
#' @param cf cf_orig and possibly cf_boot object.
#' @param energy_difference_val numeric. A single energy value \eqn{\Delta E}{Delta E} for
#'   the weighting.
#' @param energy_difference_boot numeric vector. Samples for the energy
#'   difference value.
#' @param cosh_factor integer, either `+1` or `-1`. Determines the sign $c$ in
#'   the weight factor.
#' @param offset integer. Offset for the time $t$, needed for the reweighting
#'   after a shift.
#' @param inverse boolean. If `TRUE` apply inverse weight.
#'
#' @return
#' Returns an object of class `cf`, see \link{cf}.
#' 
#' @export
weight.cf <- function (cf, energy_difference_val, energy_difference_boot,
                       cosh_factor, offset = 0, inverse = FALSE) {
  Exptt <- make_weight_factor(energy_difference_val, cf$Time, offset,
                              cf$Time/2 + offset, cosh_factor)
  if (inverse) {
    Exptt <- 1 / Exptt
  }
  if (!is.null(cf$cf)) {
    cf$cf <- cf$cf * t(array(Exptt, dim = dim(cf$cf)[c(2, 1)]))
  }
  cf$cf.tsboot$t0 <- cf$cf.tsboot$t0 * Exptt
  for (i in c(1:cf$boot.R)) {
    Exptt <- make_weight_factor(energy_difference_boot[i], cf$Time, offset,
                                cf$Time/2 + offset, cosh_factor)
    if (inverse) {
      Exptt <- 1 / Exptt
    }
    cf$cf.tsboot$t[i, ] <- cf$cf.tsboot$t[i, ] * Exptt
  }

  return (cf)
}

#' Weight-shift-reweight a correlation function
#'
#' The correlation function is weighted with [`weight.cf`], then shifted, and
#' then weighted again with the inverse weighting factor.
#'
#' @inheritParams weight.cf
#'
#' @return
#' Returns an object of class `cf`, see \link{cf}.
#' 
#' @export
weight_shift_reweight.cf <- function (cf, energy_difference_val, energy_difference_boot, cosh_factor) {
  cf <- weight.cf(cf, energy_difference_val, energy_difference_boot,
                  cosh_factor, 0, TRUE)
  cf <- takeTimeDiff.cf(cf)
  cf <- weight.cf(cf, energy_difference_val, energy_difference_boot,
                  cosh_factor, -1, FALSE)

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
                 endcorr = cf$endcorr,
                 cf.tsboot = cf$cf.tsboot,
                 resampling_method = cf$resampling_method)
  ret <- cf_shifted(ret,
                    deltat = cf$deltat,
                    forwardshift = cf$forwardshift)
  ret <- cf_weighted(ret,
                     weight.factor = 1 / (exp(energy_difference_val) * 1),
                     weight.cosh = cosh_factor == +1)

  return (invisible(ret))
}

#' Remove Thermal States by Weighting and Shifting
#'
#' @param cf Object of type \link{cf}
#' @param single.cf1 Object of type \link{cf}
#' @param single.cf2 Object of type \link{cf}
#' @param p1 Numeric vector. Spatial momentum of first state
#' @param p2 Numeric vector. Spatial momentum of second state
#' @param L Integer. Spatial lattice extent.
#' @param lat.disp Boolean. Use lattice dispersion relation instead of
#'   continuum one
#' @param weight.cosh Boolean. Use cosh functional form in the
#'   weighting procedure
#'
#' @return weighted and shifted correlation function as a \link{cf} object.
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
    warning("L was missing, set it to Time/2\n")
  }

  if (missing(single.cf2)) {
    single.cf2 <- single.cf1
  }

  if(cf$boot.R != single.cf1$boot.R || single.cf1$boot.R != single.cf2$boot.R ||
     cf$boot.l != single.cf1$boot.l || single.cf1$boot.l != single.cf2$boot.l) {
##     cf$seed != single.cf1$seed || single.cf1$seed != single.cf2$seed) {
    stop("please provide equally bootstrapped cfs to removeTemporal.cf using the same configurations and seed for all and the same boot.R and boot.l!\n")
  }

  mass1 <- extract_mass(single.cf1)
  mass2 <- extract_mass(single.cf2)

  ## use momenta p1 and p2 and lattice or continuum dispersion relation to
  ## shift energies
  if (any(p1 != 0)) {
    mass1$t0 <- dispersion_relation(mass1$t0, p1, L, lattice_disp = lat.disp)
    mass1$t <- dispersion_relation(mass1$t, p1, L, lattice_disp = lat.disp)
  }

  if (any(p2 != 0)) {
    mass2$t0 <- dispersion_relation(mass2$t0, p2, L, lattice_disp = lat.disp)
    mass2$t <- dispersion_relation(mass2$t, p2, L, lattice_disp = lat.disp)
  }

  if (weight.cosh) {
    cosh_factor  <- 1.
  } else {
    cosh_factor <- 0.
  }

  weight_shift_reweight.cf(cf, mass2$t0 - mass1$t0, mass2$t - mass1$t, cosh_factor)
}
