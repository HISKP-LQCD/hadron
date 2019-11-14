#' Correlation function container
#'
#' This function `cf()` creates containers for correlation functions
#' of class `cf`. This class is particularly designed to deal with
#' correlation functions emerging in statistical and quantum field theory
#' simulations. Arithmetic operations are defined for this class in
#' several ways, as well as concatenation and \link{is.cf}.
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
cf <- function () {
  cf <- list()
  class(cf) <- append(class(cf), 'cf')
  return (cf)
}

#' CF metadata mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param nrObs Integer, number of different measurements contained in this correlation function. One can use \link{c.cf} to add multiple observables into one container. This is for instance needed when passing to the \link{gevp} function.
#' @param Time Integer, full time extent.
#' @param nrStypes Integer, number of smearing types.
#' @param symmetrised Logical, indicating whether the correlation function has been symmetrized.
#'
#' @family cf constructors
#'
#' @export
cf_meta <- function (.cf = cf(), nrObs = 1, Time = NA, nrStypes = 1, symmetrised = FALSE) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$nrObs = nrObs
  .cf$Time = Time
  .cf$nrStypes = nrStypes
  .cf$symmetrised = symmetrised

  class(.cf) <- append(class(.cf), 'cf_meta')
  return (.cf)
}

#' Bootstrapped CF mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param boot.R Integer, number of bootstrap samples used.
#' @param boot.l Integer, block length in the time-series bootstrap process.
#' @param seed Integer, random number generator seed used in bootstrap.
#' @param sim Character, `sim` argument of \link[boot]{tsboot}.
#' @param endcorr Boolean, `endcorr` argumetn of \link[boot]{tsboot}.
#' @param cf.tsboot List, result from the \link[boot]{tsboot} function for the real part.
#' @param icf.tsboot List, result from the \link[boot]{tsboot} function for the imaginay part.
#' @param resampling_method Character, either 'bootstrap' or 'jackknife'
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `cf0`: Numeric vector, mean value of original measurements, convenience copy of `cf.tsboot$t0`.
#' - `tsboot.se`: Numeric vector, standard deviation over bootstrap samples.
#' - `boot.samples`: Logical, indicating whether there are bootstrap samples available. This is deprecated and instead the presence of bootstrap samples should be queried with `inherits(cf, 'cf_boot')`.
#' - `error_fn`: Function, takes a vector of samples and computes the error. In the bootstrap case this is just the `sd` function. Use this function instead of a `sd` in order to make the code compatible with jackknife samples.
#'
#' @family cf constructors
#'
#' @export
cf_boot <- function (.cf = cf(), boot.R, boot.l, seed, sim, endcorr, cf.tsboot, icf.tsboot = NULL, resampling_method) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$boot.R <- boot.R
  .cf$boot.l <- boot.l
  .cf$seed <- seed
  .cf$sim <- sim
  .cf$endcorr <- endcorr

  .cf$cf.tsboot <- cf.tsboot
  .cf$icf.tsboot <- icf.tsboot

  if (resampling_method == 'bootstrap') {
    .cf$error_fn <- sd
    .cf$cov_fn <- cov
  }
  else if (resampling_method == 'jackknife') {
    .cf$error_fn <- function (...) jackknife_error(..., boot.l = boot.l)
    .cf$cov_fn <- jackknife_cov
  } else {
    stop('This resampling method is not implemented')
  }

  .cf$resampling_method <- resampling_method

  .cf$cf0 <- cf.tsboot$t0
  .cf$tsboot.se <- apply(.cf$cf.tsboot$t, MARGIN = 2L, FUN = .cf$error_fn)

  if( !is.null( icf.tsboot ) ){
    .cf$icf0 <- icf.tsboot$t0
    .cf$itsboot.se <- apply(.cf$icf.tsboot$t, MARGIN = 2L, FUN = .cf$error_fn)
  }
  
  .cf$boot.samples <- TRUE

  class(.cf) <- append(class(.cf), 'cf_boot')
  return (.cf)
}

#' Estimates error from jackknife samples
#'
#' Currently this uses the mean over the jackknife samples in order to compute
#' the error. It would be better in the case of a bias to use the mean over the
#' original data instead. This would require a second parameter and therefore
#' is incompatible with the previously used `sd` everywhere for the bootstrap
#' samples. As the `sd` for the bootstrap samples also does not include the
#' original data, this likely is similar in terms of bias.
#'
#' @param samples Numeric vector.
#' @param boot.l Block length for bootstrapping.
#' @param na.rm Logical. Determines whether `NA` values shall be removed, see
#' Description for details.
#'
#' @description
#' Computes the jackknife error which is just
#' \deqn{\sum_{i=0}^N (x_i - \bar x)^2 \,.}
#' Internally we use
#' \deqn{\frac{(N-1)^2}{N} \mathop{\mathrm{sd}}(X)}
#' in order to benefit from the optimized standard deviation function.
#'
#' The width of the bootstrap distribution does not change with the number of
#' elements. The jackknife distribution crucially depends on the number of
#' measurements that one started with. Therefore we cannot just drop the NA
#' values and are done with it. Instead we need to rescale with the
#' \eqn{\sqrt{N / m}} where \eqn{N} is the number of original measurements and
#' \eqn{m} is the number of non-NA values. With NA values removed we would
#' otherwise underestimate the uncertainty.
#'
#' @export
jackknife_error <- function (samples, boot.l = 1, na.rm = FALSE) {
  ## Number of jackknife samples.
  N <- length(samples)

  if (na.rm) {
    selection <- !is.na(samples)
    samples <- samples[selection]

    ## Number of non-NA samples.
    m <- sum(selection)
  } else {
    m <- N
  }

  sqrt((N - 1) * (m - 1) / (m * boot.l)) * sd(samples)
}

#' jackknife_cov
#'
#' @description
#' Computes covariance matrix for jackknife samples.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y ‘NULL’ (default) or a vector, matrix or data frame with
#'        compatible dimensions to ‘x’. The default is equivalent to 
#'        ‘y = x’ (but more efficient).
#' @param na.rm logical. The rows containing any `NA` will be deleted if this
#' option is set.
#' @param ... parameters to be forwarded to \link{cov}.
#'
#' @export
jackknife_cov <- function (x, y = NULL, na.rm = FALSE, ...) {
    factor <- 1.0
    
    if (is.null(y)) {
        N <- nrow(x)
        if (na.rm) {
            na_values <- apply(x, 1, function (row) any(is.na(row)))
            m <- sum(na_values)
            x <- x[!na_values, ]
            factor <- N / m
        }
    } else {
        N <- length(x)
        if (na.rm) {
        na_values <- is.na(x) | is.na(y)
            m <- sum(na_values)
            x <- x[!na_values]
            y <- y[!na_values]
            factor <- N / m
        }
    }
    
    (N-1)^2 / N * factor * cov(x, y, ...)
}

#' Original data CF mixin constructor
#'
#' @param .cf `cf` object to extend. Named with a leading period just to distinguish it from the member also named `cf`.
#' @param cf Numeric matrix, original data for all observables and measurements.
#' @param icf Numeric matrix, imaginary part of original data. Be very careful with this as quite a few functions just ignore the imaginary part and drop it in operations.
#'
#' @family cf constructors
#'
#' @export
cf_orig <- function (.cf = cf(), cf, icf = NULL) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$cf <- cf
  .cf$icf <- icf

  class(.cf) <- append(class(.cf), 'cf_orig')
  return (.cf)
}

#' Principal correlator CF mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param id Integer, number of the principal correlator from the GEVP.
#' Ascending with eigenvalue, so `id = 1` is the lowest state.
#' @param gevp_reference_time Integer, reference time \eqn{t_0} that has been
#' used in the GEVP.
#'
#' @family cf constructors
#'
#' @export
cf_principal_correlator <- function (.cf = cf(), id, gevp_reference_time) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$id <- id
  .cf$gevp_reference_time <- gevp_reference_time

  class(.cf) <- append(class(.cf), 'cf_principal_correlator')
  return (.cf)
}

#' Shifted CF mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param deltat TODO
#' @param forwardshift Logical, TODO
#'
#' @details
#'
#' The following fields will also be made available:
#'
#' - `shifted`: Logical, whether the correlation function has been shifted This is deprecated and instead the presence of a shift should be queried with `inherits(cf, 'cf_shifted')`.
#'
#' @family cf constructors
#'
#' @export
cf_shifted <- function (.cf = cf(), deltat, forwardshift) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$deltat <- deltat
  .cf$forwardshift <- forwardshift

  .cf$shifted <- TRUE

  class(.cf) <- append(class(.cf), 'cf_shifted')
  return (.cf)
}

#' Smeared CF mixin constructor
#'
#' @param .cf `cf` object to extend.
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
cf_smeared <- function (.cf = cf(), scf, iscf = NULL, nrSamples, obs) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$scf <- scf
  .cf$iscf <- iscf
  .cf$nrSamples <- nrSamples
  .cf$obs <- obs

  .cf$smeared <- TRUE

  class(.cf) <- append(class(.cf), 'cf_smeared')
  return (.cf)
}

#' Subtracted CF mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param subtracted.values Numeric matrix, TODO
#' @param subtracted.ii Integer vector, TODO
#'
#' @family cf constructors
#'
#' @export
cf_subtracted <- function (.cf = cf(), subtracted.values, subtracted.ii) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$subtracted.value <- subtracted.values
  .cf$subtracted.ii <- subtracted.ii

  class(.cf) <- append(class(.cf), 'cf_subtracted')
  return (.cf)
}

#' Weighted CF mixin constructor
#'
#' @param .cf `cf` object to extend.
#' @param weight.factor TODO
#' @param weight.cosh TODO
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
cf_weighted <- function (.cf = cf(), weight.factor, weight.cosh) {
  stopifnot(inherits(.cf, 'cf'))

  .cf$weight.factor <- weight.factor
  .cf$weight.cosh <- weight.cosh

  .cf$weighted <- TRUE

  class(.cf) <- append(class(.cf), 'cf_weighted')
  return (.cf)
}

#' Checks whether the cf object contains no data
#'
#' @param .cf `cf` object.
#'
#' @examples
#' # The empty cf object must be empty:
#' is_empty.cf(cf())
#'
#' # The sample cf must not be empty:
#' is_empty.cf(samplecf)
is_empty.cf <- function (.cf) {
  setequal(class(.cf), class(cf())) &&
    is.null(names(.cf))
}

#' Checks whether the resampling of two cf objects is compatible
#' 
#' @param cf1 `cf` object with `cf_boot`
#' @param cf2 `cf` object with `cf_boot`
#'
#' @return List of named booleans for each of the checked conditions
#'         with elements `boot`, `boot.R`, `boot.l`, `sim`, `endcorr`,
#'         `resampling_method`, `error_fn`, `boot_dim`, `icf` and, optionally
#'         `iboot_dim` (if both `cf1` and `cf2` contain imaginary parts).
#' @export
resampling_is_compatible(cf1, cf2){
  
  res <- list()
  res$boot <- ( inherits(cf1, 'cf_boot') & inherits(cf2, 'cf_boot') )
  res$seed <- (cf1$seed == cf2$seed)
  res$boot.R <- (cf1$boot.R == cf2$boot.R)
  res$boot.l <- (cf1$boot.l == cf2$boot.l)
  res$sim <- (cf1$sim == cf2$sim)
  res$endcorr <- (cf1$endcorr == cf2$endcorr)
  res$resampling_method <- (cf1$resampling_method == cf2$resampling_method)
  res$error_fn <- (cf1$error_fn == cf2$error_fn)
  res$boot_dim <- all(dim(cf1$cf.tsboot$t) == dim(cf2$cf.tsboot$t))
  res$icf <- (has_icf(cf1) == has_icf(cf2))
  if( has_icf(cf1) & res$icf ){
    res$iboot_dim <- all(dim(cf1$icf.tsboot$t) == dim(cf2$icf.tsboot$t))
  }

  return(res)
}

#' Checks whether the cf object contains an imaginary part
#'
#' @param .cf `cf` object
#'
#'
has_icf <- function(.cf) {
  return( !is.null(.cf$icf) )
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
  stopifnot(inherits(cf, 'cf_orig'))

  boot.l <- ceiling(boot.l)
  boot.R <- floor(boot.R)

  stopifnot(boot.l >= 1)
  stopifnot(boot.l <= nrow(cf$cf))
  stopifnot(boot.R >= 1)

  ## we set the seed for reproducibility and correlation
  old_seed <- swap_seed(seed)
  ## now we bootstrap the correlators
  cf.tsboot <- boot::tsboot(cf$cf, statistic = function(x){ return(apply(x, MARGIN=2L, FUN=mean))},
                            R = boot.R, l = boot.l, sim = sim, endcorr = endcorr)

  icf.tsboot <- NULL
  if( has_icf(cf) ){
    # no need to store the old seed again, but we definitely need to reset the RNG again!
    swap_seed(seed)
    icf.tsboot <- boot::tsboot(cf$icf, statistic = function(x){ return(apply(x, MARGIN=2L, FUN=mean)) },
                               R = boot.R, l = boot.l, sim = sim, endcorr = endcorr)
  }

  cf <- cf_boot(cf,
                boot.R = boot.R,
                boot.l = boot.l,
                seed = seed,
                sim = sim,
                endcorr = endcorr,
                cf.tsboot = cf.tsboot,
                icf.tsboot = icf.tsboot,
                resampling_method = "bootstrap")

  restore_seed(old_seed)

  return(invisible(cf))
}

jackknife.cf <- function(cf, boot.l = 1) {
  stopifnot(inherits(cf, 'cf_orig'))

  stopifnot(boot.l >= 1)
  boot.l <- ceiling(boot.l)

  ## blocking with fixed block length, but overlapping blocks
  ## number of observations
  n <- nrow(cf$cf)
  ## number of overlapping blocks
  N <- n-boot.l+1
  
  cf$cf0 <- apply(cf$cf, 2, mean)

  t <- array(NA, dim = c(N, ncol(cf$cf)))
  t0 <- cf$cf0
  for (i in 1:N) {
    ## The measurements that we are going to leave out.
    ii <- c(i:(i+boot.l-1))
    ## jackknife replications of the mean
    t[i, ] <- apply(cf$cf[-ii, ], 2L, mean)
  }

  cf.tsboot <- list(t = t,
                    t0 = t0,
                    R = N,
                    l = boot.l)

  icf.tsboot <- NULL
  if( has_icf(cf) ){
    cf$icf0 <- apply(cf$icf, 2, mean)

    t <- array(NA, dim = c(N, ncol(cf$icf)))
    t0 <- cf$icf0
    for (i in 1:N) {
      ## The measurements that we are going to leave out.
      ii <- c(i:(i+boot.l-1))
      ## jackknife replications of the mean
      t[i, ] <- apply(cf$icf[-ii, ], 2L, mean)
    }
    icf.tsboot <- list(t = t,
                       t0 = t0,
                       R = N,
                       l = boot.l)
  }

  cf <- invalidate.samples.cf(cf)
  cf <- cf_boot(cf,
                boot.R = cf.tsboot$R,
                boot.l = cf.tsboot$l,
                seed = 0,
                sim = 'geom',
                cf.tsboot = cf.tsboot,
                icf.tsboot = icf.tsboot,
                resampling_method = 'jackknife')

  return (invisible(cf))
}

#' uwerr.cf
#' @description
#' Gamma method analysis on all time-slices in a 'cf' object
#'
#' @param cf Object of type `cf` containing `cf_orig`
#'
#' @return A list with a named element `uwcf` which contains a data frame
#'         with six columns, `value`, `dvalue`, `ddvalue`, `tauint`, `dtauint`
#'         corresponding to what is returned by \link{uwerrprimary}. The sixth
#'         column, `t`, is just an index counting the columns in the original `cf$cf`.
#'         If `cf` contains an imaginary part, the return value contains another
#'         list element, `uwicf` of the same structure as `uwcf`.
#'         There are as many rows as there were columns in `cf$cf` and/or `cf$icf`.
#'         When the call to \link{uwerrprimary} fails for a particular column of `cf$cf`
#'         or `cf$icf`, the corresponding row of `uwcf` and/or `uwicf` will contain
#'         `NA`.
#' 
uwerr.cf <- function(cf){
  stopifnot(inherits(cf, 'cf_orig'))

  uw_wrapper <- function(x){
    uw_tmp <- try(uwerrprimary(data=x), silent=TRUE)
    if( any(class(uw) == "try-error") ){
      c(value=NA, dvalue=NA, ddvalue=NA, tauint=NA, dtauint=NA)
    } else {
      c(value=$uw_tmp$value, dvalue=uw_tmp$dvalue, ddvalue=uw_tmp$ddvalue,
        tauint=uw_tmp$tau_int, dtauint=uw_tmp$dtauing)
    }
  }

  res <- list()
  res[["uwcf"]] <- cbind(as.data.frame(t(apply(X=cf$cf, MARGIN=2L, FUN=uw_wrapper))),
                         t=(1:ncol(cf$cf)))
  if( has_icf(cf) ){
    res[["uwicf"]] <- cbind(as.data.frame(t(apply(X=cf$icf, MARGIN=2L, FUN=uw_wrapper))),
                            t=(1:ncol(cf$icf)))
  }
  return(res)
}

addConfIndex2cf <- function(cf, conf.index) {
  if(is.null(cf$conf.index)) {
    cf$conf.index <- conf.index
  }
  return(cf)
}

addStat.cf <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf2, 'cf'))

  if (is_empty.cf(cf1)) {
    return (invisible(cf2))
  }
  if (is_empty.cf(cf2)) {
    return (invisible(cf1))
  }

  stopifnot(inherits(cf1, 'cf_meta'))
  stopifnot(inherits(cf2, 'cf_meta'))

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

#' @title Average local-smeared and smeared-local correlators
#' @description
#' averages local-smeared and smeared-local correlators in cf and adjusts
#' nrStypes accordingly
#' by default, assumes that LS and SL are in columns (T/2+1)+1:3*(T/2+1)
#'
#' @param cf object of type \link{cf}
#' @param cols columns to be averaged over
#' 
avg.ls.cf <- function(cf, cols = c(2, 3)) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(cf$nrStypes >= 2)

  timeslices <- cf$Time/2+1

  ind.ls <- ( (cols[1]-1)*timeslices+1 ):( cols[1]*timeslices )
  ind.sl <- ( (cols[2]-1)*timeslices+1 ):( cols[2]*timeslices )

  cf$cf[,ind.ls] <- 0.5 * ( cf$cf[,ind.ls] + cf$cf[,ind.sl] )
  cf$cf <- cf$cf[,-ind.sl]

  if( has_icf(cf) ){
    cf$icf[,ind.ls] <- 0.5 * ( cf$icf[,ind.ls] + cf$icf[,ind.sl] )
    cf$icf <- cf$icf[,-ind.sl]
  }

  cf$nrStypes <- cf$nrStypes-1
  return (cf)
}

#' @title average close-by-times in a correlation function
#' @description
#' "close-by-times" averaging replaces the value of the correlation function at t
#' with the "hypercubic" average with the values at the neighbouring time-slices
#' with weights 0.25, 0.5 and 0.25
#'   C(t') = 0.25 C(t-1) + 0.5 C(t) + 0.25 C(t+1)
#' where periodic boundary conditions are assumed in shift.cf
#'
#' @param cf object of type \link{cf}
#' 
avg.cbt.cf <- function(cf){
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))

  # copy for shifting
  cf2 <- cf
  cf <- mul.cf(cf, 0.5)

  # average over shifted correlation functions
  for( p in c(-1,1) ){
    cf <- cf + mul.cf(shift.cf(cf2,p),0.25)
  }
  cf <- invalidate.samples.cf(cf)
  return(invisible(cf))
}

#' Arithmetically adds two correlation functions
#'
#' @param cf1,cf2 `cf_orig` object.
#' @param a,b Numeric. Factors that multiply the correlation function before
#' the addition.
#' 
#' Since addition is associative, this operates also on the bootstrap samples
#' and these are thus not invalidated in the process.
#'
#' @return
#' The value is
#' \deqn{a C_1 + b C_2 \,.}
#'
#' @export
add.cf <- function(cf1, cf2, a = 1.0, b = 1.0) {
  stopifnot(inherits(cf1, 'cf'))
  stopifnot(inherits(cf2, 'cf'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(inherits(cf1, 'cf_meta'))
  stopifnot(inherits(cf2, 'cf_meta'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- a*cf1$cf + b*cf2$cf

  if( has_icf(cf1) | has_icf(cf2) ){
    stopifnot( has_icf(cf1) & has_icf(cf2) )
    stopifnot( all(dim(cf1$icf) == dim(cf2$icf) ) )
    cf$icf <- a*cf1$icf + b*cf2$icf
  }

  # now reconstruct the bootstrap samples
  invalidate.samples.cf(cf)

  if( inherits(cf1, 'cf_boot') | inherits(cf2, 'cf_boot') ){
    stopifnot( all(unlist(resampling_is_compatible(cf1, cf2))) )

    cf.tsboot <- cf1$cf.tsboot
    cf.tsboot$t <- cf1$cf.tsboot$t + cf2$cf.tsboot$t
    cf.tsboot$t0 <- cf1$cf.tsboot$t0 + cf2$cf.tsboot$t0
    cf.tsboot$tseries <- cf1$cf.tsboot$tseries + cf2$cf.tsboot$tseries

    icf.tsboot <- NULL
    if( has_icf(cf1) ){
      # no further tests required as this was already tested for compatibility twice
      icf.tsboot <- cf1$icf.tsboot
      icf.tsboot$t <- a*cf1$icf.tsboot$t + b*cf2$icf.tsboot$t
      icf.tsboot$t0 <- a*cf1$icf.tsboot$t0 + b*cf2$icf.tsboot$t0
      icf.tsboot$tseries <- a*cf1$icf.tsboot$tseries + b*cf2$icf.tsboot$tseries
    }
    # use constructor to also update cf0 / icf0 and tsboot.se / itsboot.se
    cf <- cf_boot(cf = cf,
                  boot.R = cf1$boot.R,
                  boot.l = cf1$boot.l,
                  seed = cf1$seed,
                  sim = cf1$sim,
                  endcorr = cf1$endcorr,
                  cf.tsboot = cf.tsboot,
                  icf.tsboot = icf.tsboot,
                  resampling_method = cf1$resampling_method)
  }
  return(cf)
}

#' Arithmetically add correlators
#'
#' @param cf1,cf2 `cf_orig` objects.
#'
#' @export
'+.cf' <- function (cf1, cf2) {
  add.cf(cf1, cf2, a = 1.0, b = 1.0)
}

#' Arithmetically subtract correlators
#'
#' @param cf1,cf2 `cf_orig` objects.
#'
#' @export
'-.cf' <- function(cf1, cf2) {
  add.cf(cf1, cf2, a = 1.0, b = -1.0)
}

#' Divide two cf objects by each other measurement by measurement
#'
#' Note that no complex arithmetic is used, real and imaginary parts are 
#' treated as seperate and indepenent, such that the real part of one
#' is the divided by the real part of the other and similarly for the
#' imaginary parts.
#'
#' Note that this is generally only allowed on bootstrap samples and mean values,
#' although it makes sense in some exeptional circumstances. Don't use this
#' function unless you're certain that you should!
#' 
#' @param cf1,cf2 `cf_orig` objects.
#'
#' @export
'/.cf' <- function(cf1, cf2) {
  stopifnot(inherits(cf1, 'cf_meta'))
  stopifnot(inherits(cf2, 'cf_meta'))
  stopifnot(inherits(cf1, 'cf_orig'))
  stopifnot(inherits(cf2, 'cf_orig'))
  stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
  stopifnot(cf1$Time == cf2$Time)

  cf <- cf1
  cf$cf <- cf1$cf / cf2$cf

  if( has_icf(cf1) | has_icf(cf2) ){
    stopifnot(has_icf(cf1) & has_icf(cf2))
    cf$icf <- cf1$icf / cf2$icf
  }

  cf <- invalidate.samples.cf(cf)
  return (cf)
}

#' Arithmetically scale a correlator
#'
#' @param cf `cf_orig` objects.
#' @param a Numeric, scaling factor.
#'
#' @export
mul.cf <- function(cf, a=1.) {
  stopifnot(inherits(cf, 'cf_orig'))
  stopifnot(is.numeric(a))

  cf$cf <- a*cf$cf

  if( has_icf(cf) ){
    cf$icf <- a*cf$icf
  }
  if( inherits(cf, 'cf_boot') ){
    cf$cf.tsboot$t <- a*cf$cf.tsboot$t
    cf$cf.tsboot$t0 <- a*cf$cf.tsboot$t0
    cf$cf.tsboot$tseries <- a*cf$cf.tsboot$t

    if( has_icf(cf) ){
      cf$icf.tsboot$t <- a*cf$icf.tsboot$t
      cf$icf.tsboot$t0 <- a*cf$icf.tsboot$t0
      cf$icf.tsboot$tseries <- a*cf$icf.tsboot$tseries
    }
    # cf_boot will take care of cf0 / icf0 and tsboot.se / itsboot.se
    cf <- cf_boot(cf = cf,
                  boot.R = cf$boot.R,
                  boot.l = cf$boot.l,
                  seed = cf$seed,
                  sim = cf$sim,
                  endcorr = cf$endcorr,
                  cf.tsboot = cf$cf.tsboot,
                  icf.tsboot = cf$icf.tsboot,
                  resampling_method = cf$resampling_method)
  }
  return (cf)
}

extractSingleCor.cf <- function(cf, id=c(1)) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))

  ii <- c()
  for (i in c(1:length(id))) {
    num_time <- if (cf$symmetrised) cf$Time / 2 + 1 else cf$Time
    ii <- c(ii, c(1:num_time) + (id[i]-1) * num_time)
  }

  # TODO: This should be done using constructors.
  cf$cf <- cf$cf[,ii]

  if( has_icf(cf) ){
    cf$icf <- cf$icf[,ii]
  }

  if (inherits(cf, 'cf_boot')) {
    cf$cf0 <- cf$cf0[ii]
    cf$tsboot.se <- cf$tsboot.se[ii]
    cf$cf.tsboot$t0 <- cf$cf.tsboot$t0[ii]
    cf$cf.tsboot$t <- cf$cf.tsboot$t[,ii]
    cf$cf.tsboot$data <- cf$cf.tsboot$data[,ii]

    if( has_icf(cf) ){
      cf$icf0 <- cf$icf0[ii]
      cf$itsboot.se <- cf$itsboot.se[ii]
      cf$icf.tsboot$t0 <- cf$icf.tsboot$t0[ii]
      cf$icf.tsboot$t <- cf$icf.tsboot$t[,ii]
      cf$icf.tsboot$data <- cf$icf.tsboot$data[,ii]
    }
  }
  cf$nrObs <- 1
  cf$nsStypes <- 1
  return (cf)
}

#' Checks whether an object is a cf
#'
#' @param x Object, possibly of class `cf`.
#'
#' @export
is.cf <- function (x) {
  inherits(x, "cf")
}

#' Concatenate correlation function objects
#'
#' @param ... Zero or multiple objects of type `cf`.
c.cf <- function (...) {
  rval <- Reduce(concat.cf, list(...), cf())
  return (invisible(rval))
}

#' Concatenate two correlation function objects
#'
#' @param left,right `cf` objects to concatenate.
concat.cf <- function (left, right) {
  stopifnot(inherits(left, 'cf'))
  stopifnot(inherits(right, 'cf'))

  if (inherits(left, 'cf_boot') || inherits(right, 'cf_boot')) {
    warning('At least one argument of concat.cf has bootstrap/jacknife samples which cannot be concatenated. The samples will be discarded.')
  }

  # In case that one of them does not contain data, the other one is the
  # result. This satisfies the neutral element axiom of a monoid.
  if (is_empty.cf(left)) {
    return (right)
  }
  if (is_empty.cf(right)) {
    return (left)
  }

  stopifnot(inherits(left, 'cf_meta'))
  stopifnot(inherits(right, 'cf_meta'))

  # At this point both `cf` objects given here have original data, therefore we
  # need to concatenate them.

  # A few checks for compatability.
  stopifnot(left$Time == right$Time)
  stopifnot(nrow(left$cf) == nrow(right$cf))
  stopifnot(left$symmetrised == right$symmetrised)
  stopifnot(left$nrStypes == right$nrStypes)

  rval <- cf_meta(nrObs = left$nrObs + right$nrObs,
                  Time = left$Time,
                  nrStypes = left$nrStypes,
                  symmetrised = left$symmetrised)
  rval <- cf_orig(.cf = rval,
                  cf = cbind(left$cf, right$cf),
                  icf = cbind(left$icf, right$icf))
  return (invisible(rval))
}

#' Plot a correlation function
#'
#' @param x `cf_boot` object
#' @param neg.vec Numeric vector of length `cf$cf0`. This allows switching the
#' sign for certain time slices or observables such that displaying in
#' log-scale is sensible.
#' @param rep See \code{\link{plotwitherror}}.
#' @param ... Graphical parameter to be passed on to \link{plotwitherror}
#'
#' @inheritParams plotwitherror
#'
#' @export
plot.cf <- function(x, neg.vec = rep(1, times = length(cf$cf0)), rep = FALSE, ...) {
  cf <- x
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(inherits(cf, 'cf_meta'))

  val <- cf$cf0
  err <- cf$tsboot.se

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

#' shift a correlation function by 'places' time-slices
#'
#'   C'(t) = C(t+places)
#' where places can be positive or negative as required and periodic boundary conditions
#' in time are assumed
#' @param cf unsymmetrised correlation function (cf_meta and cf_orig mixins required)
#' @param places integer number of time-slices for backward (negative) or forward (positive) shifts
shift.cf <- function(cf, places) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))
  if( cf$symmetrised ){
    stop("A correlation function can only be time-shifted if it is not symmetrised!")
  }  

  if(places == 0){
    return(invisible(cf))
  }

  cf <- invalidate.samples.cf(cf)
  n <- cf$Time

  for( oidx in 0:(cf$nrObs-1) ){
    for( sidx in 0:(cf$nrStypes-1) ){
      istart <- cf$Time*cf$nrStypes*oidx + cf$Time*sidx + 1
      iend <- istart + cf$Time - 1

      if( places < 0 ){
        ishift <- c( (iend - abs(places)):iend,
                     (istart:(iend-abs(places)-1)) )
      } else {
        ishift <- c( (istart+places):iend,
                      istart:(istart+places-1) )
      }
      cf$cf[,istart:iend] <- cf$cf[,ishift, drop=FALSE]
      if( has_icf(cf) ){
        cf$icf[,istart:iend] <- cf$icf[,ishift, drop=FALSE]
      }
    }
  }
  return(invisible(cf))
}

#' Invalidate samples
#'
#' When a correlation function is modified, any resampling should be
#' invalidated. We could instead also choose to properly work with the samples,
#' but most computations are done with the original data anyway.
#'
#' @param cf `cf` object.
#'
invalidate.samples.cf <- function (cf) {
  cf$boot.l <- NULL
  cf$boot.R <- NULL
  cf$boot.samples <- NULL
  cf$seed <- NULL
  cf$sim <- NULL

  cf$cf.tsboot <- NULL
  cf$tsboot.se <- NULL
  cf$cf0 <- NULL

  if( has_icf(cf) ){
    cf$icf.tsboot <- NULL
    cf$itsboot.se <- NULL
    cf$icf0 <- NULL
  }

  class(cf) <- setdiff(class(cf), c('cf_boot', 'cf_jackknife'))

  return(invisible(cf))
}

#' symmetrise.cf
#'
#' @description
#' Symmetrises a correlation function
#' 
#' @param cf Object of type `cf`.
#' @param sym.vec Integer vector. Takes values -+1 for
#'   (anti)symmmetrisation. If longer than 1, it must be of length
#'   equal to the number of observalbes in `cf`. Otherwise, the same
#'   operation is applied to all observables in `cf`, which is the
#'   default. 
symmetrise.cf <- function(cf, sym.vec=c(1) ) {
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_orig'))

  if(cf$symmetrised){
    stop("symmetrise.cf: cf was already symmetrised")
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
      if( has_icf(cf) ){
        cf$icf[, (istart+1):(ihalf-1)] <- 0.5*( cf$icf[, (istart+1):(ihalf-1)] +
                                                sym.vec[oidx+1]*cf$icf[, rev((ihalf+1):iend)] )
      }
    }
  }
  # remove now unnecessary time slices
  cf$cf <- cf$cf[, -isub]
  if( has_icf(cf) ){
    cf$icf <- cf$icf[, -isub]
  }
  cf$symmetrised <- TRUE
  return(invisible(cf))
}

#' summary.cf
#'
#' @param object Object of type \link{cf}
#' @param ... Generic parameters to pass on.
summary.cf <- function(object, ...) {
  cf <- object
  stopifnot(inherits(cf, 'cf_meta'))

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

    if( has_icf(cf) ){
      out <- cbind(out, iC=cf$icf0, itsboot.se = cf$itsboot.se)
    }
  }

  if (inherits(cf, 'cf_jackknife')) {
    out <- cbind(out, jackknife.se=cf$jackknife.se)
    out <- cbind(out, jab.se=cf$jack.boot.se)

    if( has_icf(cf) ){
      out <- cbind(out, ijackknife.se = cf$ijackknife.se)
      out <- cbind(out, ijab.se = cf$ijack.boot.se)
    }
  }

  if(exists("out")) {
    print(out)
  }
}

#' print.cf
#'
#' @param x Object of type \link{cf}
#' @param ... Generic parameters to pass on.
print.cf <- function (x, ...) {
  summary(x, ...)
}
