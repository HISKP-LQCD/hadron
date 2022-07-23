#' Container for arbitrary timeseries
#'
#' @family tseries constructors
#'
#' @return
#' returns an object of S3 class `tseries` derived from a `list`
#'
#' @examples
#' new_tseries <- tseries()
#'
#' @export
tseries <- function() {
  tseries <- list()
  class(tseries) <- append(class(tseries), 'tseries')
  return(tseries)
}

#' Original data tseries mixin constructor
#' 
#' @param .tseries `tseries` object to extend
#' @param data Tidy data frame representing an observable or set of observables
#' of arbitrary dimensionality indexed by a "simulation time" index `md_idx`.
#' At the minimum, `data` must have two columns, one of which is named
#' `md_idx`.
#' If, for example, we are dealing with two fields, one could store their values
#' in the columns `f1` and `f2`. If these fields are not scalar, their values
#' must be indexed by further variables, see `explanatory_vars`.
#' Note that it will be internally converted to a \link{tibble}.
#'
#' @param explanatory_vars Vector of strings, names of the explanatory variables
#' in `data`. If, for example, we are dealing with a 3D field which varies in time,
#' the `explanatory_vars` might by `x`, `y` and `z`.
#'
#' @family tseries constructors
#'
#' @return
#' returns the input object of class `tseries`, extended with the elements `data`
#' and `explanatory_vars`.
#'
#' @examples
#' xyz_md <- expand.grid(x = 0:9, y = 0:9, z = 0:9, md_idx = 1:100)
#' new_tseries <- tseries_orig(data = cbind(xyz_md, val = rnorm(nrow(xyz_md))),
#'                             explanatory_vars = c("x", "y", "z") )
#'
#' @export
tseries_orig <- function(.tseries = tseries(), data, explanatory_vars = NULL) {
  stopifnot(inherits(.tseries, 'tseries'))
  stopifnot("md_idx" %in% colnames(data))
  stopifnot(ncol(data) >= 2)

  if( !is.null(explanatory_vars) ){
    stopifnot( all(explanatory_vars %in% colnames(data) ) )
  }

  .tseries$data <- tibble::as_tibble(data)
  .tseries$explanatory_vars <- explanatory_vars

  class(.tseries) <- append(class(.tseries), 'tseries_orig')
  return(.tseries)
}

#' Bootstrap tseries mixin constructor
#'
#' @details 
#' Mixin constructor to add bootstrap samples to a `tseries` object.
#'
#' @param .tseries `tseries` object
#' @param boot.R Integer, number of bootstrap samples used.
#' @param boot.l Integer, block length in the time-series bootstrap process.
#' @param seed Integer, random number generator seed used in bootstrap.
#' @param sim String, `sim` argument of \link[boot]{tsboot}.
#' @param endcorr Boolean, `endcorr` argument of \link[boot]{tsboot}.
#' @param central Numeric or complex, means of the data.
#' @param samples Numeric or complex vector, the resampling samples. 
#' @param resampling_method String, either 'bootstrap' or 'jackknife'
#' @param sample_idcs Two-dimensional array which has the values of `md_idx`
#' for each bootstrap sample as its row. The number of rows should
#' correspond to boot.R, while the number of columns should correspond
#' to the number of measurements for the original data that the
#' bootstrap samples are based on. Can be \code{NULL} when dealing
#' with parametric bootstrap or quantities which arise when combining
#' data from different statistical ensembles.
#'
#' @family tseries constructors
#' 
#' @return 
#' return the input object of class `tseries` with the bootstrap mixin added
#' which supplies the members
#' \itemize{
#'   \item{boot.R: }{see params}
#'   \item{boot.l: }{see params}
#'   \item{seed: }{see params}
#'   \item{sim: }{see params}
#'   \item{endcorr: }{see params}
#'   \item{resampling_method: }{see params}
#'   \item{samples: }{see params}
#'   \item{sample_idcs: }{see params}
#'   \item{central: }{central values of the original data or parametric central values}
#'   \item{se: }{estimates of the standard error of the means}
#'   \item{error_fn: }{function appropriate to compute the statistical error depending on \code{resampling_method}}
#'   \item{cov_fn: }{function appropriate to compute the covariance matrix depending on \code{resampling_method}}
#' }
#'
#' @export
tseries_boot <- function(.tseries = tseries(), boot.R, boot.l, seed, sim, endcorr,
                         samples, central, resampling_method,
                         explanatory_vars = NULL,
                         sample_idcs = NULL){
  stopifnot(inherits(.tseries, 'tseries'))
  stopifnot(resampling_method %in% c("bootstrap", "jackknife"))

  .tseries$boot.R <- boot.R
  .tseries$boot.l <- boot.l
  .tseries$seed <- seed
  .tseries$sim <- sim
  .tseries$endcorr <- endcorr
  .tseries$sample_idcs <- sample_idcs

  if( is.null(.tseries$explanatory_vars) ){
    .tseries$explanatory_vars <- explanatory_vars
  }

  if ( resampling_method == "bootstrap" ) {
    .tseries$error_fn <- sd
    .tseries$cov_fn <- cov
  }
  else if ( reampling_method == "jackknife" ) {
    .tseries$error_fn <- function(...) jackknife_error(..., boot.l = boot.l)
    .tseries$cov_fn <- jackknife_cov
  }

  .tseries$resampling_method <- resampling_method

  .tseries$samples <- samples
  .tseries$central <- central
  .tseries$se <- tseries_apply_reduction(data = .tseries$samples,
                                         reduction = .tseries$error_fn,
                                         explanatory_vars = .tseries$explanatory_vars)

  class(.tseries) <- append(class(.tseries), 'tseries_boot')
  return(.tseries)
}

tseries_apply_reduction <- function(data, reduction, explanatory_vars = NULL) {

  data_ref <- data
  not_across_vars <- c("md_idx", "sample_idx")
  # if the data frame represents multi-dimensional data, we group on these indices
  if( !is.null(explanatory_vars) ) { 
    data_ref <- dplyr::group_by_at(data, explanatory_vars)
    not_across_vars <- c(not_across_vars, explanatory_vars)
  }

  # and now we compute the result of the function on these groups
  return(
    dplyr::summarise(
      data_ref,
      dplyr::across(.cols = dplyr::setdiff(colnames(data_ref), not_across_vars),
                    .fns = reduction),
      .groups = "drop"
    )
  )

}


#' bootstrap a `tseries` object containing original data
#'
#' @param .tseries `tseries` object with `tseries_orig` mixin
#' @param boot.R Intenger, number of bootstrap samples to use.
#' @param boot.l Integer, block length to use in bootstrap procedure.
#' @param seed Integer, random number generator seed to be used in bootstrap.
#' @param sim String, `sim` argument passed to \link[boot]{tsboot}.
#' @param endcorr Boolean, `endcorr` argument passed to \link[boot]{tsboot}.
#' @param serial Boolean, whether to disable the use of \link[parallel]{mclapply}.
#'
#' @return `tseries` object with the additional members
#' \itemize{
#' }
#'
#' @examples
#' xyz_md <- expand.grid(x = 0:9, y = 0:9, z = 0:9, md_idx = 1:100)
#' new_tseries <- tseries_orig(data = cbind(xyz_md, val = rnorm(nrow(xyz_md))),
#'                             explanatory_vars = c("x", "y", "z") )
#'
#' bs_tseries <- bootstrap.tseries(new_tseries,
#'                                 boot.R = 500, boot.l = 5, seed = 12345,
#'                                 sim = 'geom', endcorr = TRUE, serial = FALSE)
#'
#' @export
bootstrap.tseries <- function(.tseries, 
                              boot.R, boot.l, seed, sim, endcorr,
                              serial = FALSE, progress = FALSE) {

  dplyr_version_required <- "1.0.0"
  dplyr_avail <- requireNamespace("dplyr", versionCheck = list(op = ">=", version = dplyr_version_required))
  if( !dplyr_avail ){
    stop(sprintf("The 'dplyr' package in version >= %s is required to use this function!\n",
                 dplyr_version_required))
  }

  parallel_avail <- requireNamespace("parallel")
  if( !parallel_avail & serial == FALSE ){
    stop(sprintf("The 'parallel' package is required to use this function with `serial = FALSE`!\n"))
  }

  pbmcapply_avail <- requireNamespace("pbmcapply")
  if( serial == FALSE & progress & !pbmcapply_avail ){
    message("'pbmcapply' package not found, there will be no progress information of the bootstrap.\n")
  }

  grr_avail <- requireNamespace("grr")
  if( !grr_avail ){
    stop("The 'grr' package is required to use this function!\n")
  }

  if( serial ){
    my_lapply <- lapply
  } else {
    if( pbmcapply_avail & progress ){
      my_lapply <- pbmcapply::pbmclapply
    } else {
      my_lapply <- parallel::mclapply
    }
  }
 
  stopifnot(inherits(.tseries, 'tseries_orig'))

  boot.l <- ceiling(boot.l)
  boot.R <- floor(boot.R)

  stopifnot(boot.l >= 1)
  stopifnot(boot.l <= length(unique(.tseries$data$md_idx)))
  stopifnot(boot.R >= 1)

  old_seed <- swap_seed(seed)
  # unlike in other uses of tsboot in hadron, here we bootstrap the
  # `md_idx`
  sample_idcs <- boot::tsboot(unique(.tseries$data$md_idx),
                              statistic = function(x){ x },
                              R = boot.R,
                              l = boot.l,
                              sim = sim,
                              endcorr = endcorr)$t
  restore_seed(old_seed)

  ts_samples <- my_lapply(X = 1:nrow(sample_idcs),
    FUN = function(row_idx)
      {
        idcs <- sample_idcs[row_idx,,drop=TRUE]
        # compute the row indices corresponding to the bootstrap sample indices
        # it is very surprising how slow this actually is, but considering the fact
        # that no matter how we do this, we have some 
        #   nrow(sample_idcs)*ncol(sample_idcs)*length(.tseries$data$md_idx)
        # comparisons, the slowness is not surprising
        # Seems like we've found a good reason to NOT use the long table format
        # for this kind of timeseries data...

        # BaKo: this is slower ...
        # df_sample_idcs <- unlist(lapply(idcs, function(x){ which(x == .tseries$data$md_idx) }))

        # BaKo: ... than this ...
        # df_sample_idcs <- as.vector(apply(X = array(idcs, dim=c(1,length(idcs))),
        #                                   MARGIN = 2,
        #                                   FUN = function(x){ which( x[1] == .tseries$data$md_idx ) }))

        # BaKo: ... and this seems to be the fastest solution which actually makes it somewhat bearable
        # for real datasets with hundreds of thousands of lines
        df_sample_idcs <- grr::matches(x = idcs,
                                       y = .tseries$data$md_idx,
                                       all.y = FALSE,
                                       all.x = FALSE,
                                       indexes = TRUE)$y

        # we resample the data frame, this involves copies of potentially large
        # data structures but I don't see a way around for now
        return(
          dplyr::mutate(
            tseries_apply_reduction(data = dplyr::slice(.data = .tseries$data, df_sample_idcs),
                                    reduction = mean,
                                    explanatory_vars = .tseries$explanatory_vars),
            sample_idx = row_idx
          )
        )
      }
    ) # lapply
  ts_samples <- do.call(rbind, ts_samples)

  ts_central <- tseries_apply_reduction(data = .tseries$data,
                                        reduction = mean,
                                        explanatory_vars = .tseries$explanatory_vars)

  return(
    tseries_boot(.tseries,
                 boot.R = boot.R,
                 boot.l = boot.l,
                 seed = seed,
                 sim = sim,
                 endcorr = endcorr,
                 central = ts_central,
                 samples = ts_samples,
                 resampling_method = "bootstrap",
                 sample_idcs = sample_idcs)
  )
}

#' Check whether the resampling of two tseries objects is compatible
#' 
#' @param x, y `tseries` objects with `tseries_boot`
#'
#' @details
#' Checks whether binary operations such as addition can be performed on the
#' resampling samples of `x` and `y`.
#'
#' @return
#' List of named booleans for each of the checked conditions. Since compatibility
#' may depend on the context (for example, if the samples are derived from the
#' same statistical ensembles or from different ones), it is up to the 
#' calling code to finally decide if an operation should be performed or not.
#' The list contains the elements
#' \itemize{
#'   \item{boot: }{whether both objects had resampling samples}
#'   \item{seed: }{whether the seeds were the same}
#'   \item{boot.R: }{whether the \code{boot.R} parameter was the same}
#'   \item{boot.l: }{whether the \code{boot.l} parameter was the same}
#'   \item{sim: }{whether the same type of simulation was used}
#'   \item{endcorr: }{whether the same end corrections were applied}
#'   \item{resampling_method: }{whether the same resampling method was used}
#'   \item{explanatory_vars: }{whether the names of the explanatory variables are the same (vector, one boolean per var!)}
#'   \item{samples_dims: }{whether the `samples` members have the same dimensions (vector, one boolean per dim!)}
#'   \item{sample_idcs: }{whether the `sample_idcs` members are the same (array-valued, unless either of them was \code{NULL})}
#' }
#'
#' @export
resampling_is_compatible.tseries <- function(x, y){
  res <- list()
  res$boot <- ( inherits(x, 'tseries_boot') & inherits(y, 'tseries_boot') )
  res$seed <- ( x$seed == y$seed )
  res$boot.R <- ( x$boot.R == y$boot.R )
  res$boot.l <- ( x$boot.l == y$boot.l )
  res$sim <- ( x$sim == y$sim )
  res$endcorr <- ( x$endcorr == y$endcorr )
  res$resampling_method <- ( x$resampling_method == y$resampling_method )
  res$expanatory_vars <- ( x$explanatory_vars == y$explanatory_vars )
  res$samples_dims <- ( dim(x$samples) == dim(y$samples) )
  if( is.null(x$sample_idcs) | is.null(y$sample_idcs) ){
    res$sample_idcs <- FALSE
  } else {
    res$sample_idcs <- ( x$sample_idcs == y$sample_idcs )
  }
  return(res)
}

#' Arithmetically add two timeseries
#'
#' @param x, y `tseries` objects with `tseries_boot`
#'
#' @return
#' The value is
#' \deqn{x + y \,.}
#'
#' @export
'+.tseries' <- function(x, y) {
  apply_binary_elementwise.tseries(x, y, `+`)
}

#' Arithmetically subtract two timeseries
#'
#' @param x, y `tseries` objects with `tseries_boot`
#'
#' @return
#' The value is
#' \deqn{x - y \,.}
#'
#' @export
'-.tseries' <- function(x, y) {
  apply_binary_elementwise.tseries(x, y, `-`)
}

#' Arithmetically multiply two timeseries
#'
#' @param x, y `tseries` objects with `tseries_boot`
#'
#' @return
#' The value is
#' \deqn{x * y \,.}
#'
#' @export
'*.tseries' <- function(x, y) {
  apply_binary_elementwise.tseries(x, y, `*`)
}

#' Arithmetically divide two timeseries
#'
#' @param x, y `tseries` objects with `tseries_boot`
#'
#' @return
#' The value is
#' \deqn{x / y \,.}
#'
#' @export
'/.tseries' <- function(x, y) {
  apply_binary_elementwise.tseries(x, y, `/`)
}

#' apply a unary function element-wise to a `tseries` object
#'
#' @param ts `tseries` object with the `tseries_boot` mixin.
#' @param fn Function to be applied.
#'
#' @return
#' The `tseries` object with the function applied to `samples` and `central`.
#' `se` is updated and if the original \code{ts} had the `tseries_orig` mixin,
#' it is removed and the `data` field set to \code{NULL}.
#'
#' @export
apply_unary_elementwise.tseries <- function(ts, fn){

  dplyr_version_required <- "1.0.0"
  dplyr_avail <- requireNamespace("dplyr", versionCheck = list(op = ">=", version = dplyr_version_required))
  if( !dplyr_avail ){
    stop(sprintf("The 'dplyr' package in version >= %s is required to use this function!\n",
                 dplyr_version_required))
  }
  
  stopifnot(inherits(ts, 'tseries_boot'))
  stopifnot('function' %in% class(fn))
 
  obs_vars <- dplyr::setdiff(colnames(ts$samples), c(ts$explanatory_vars, "sample_idx"))

  samples <- dplyr::mutate(ts$samples,
                           dplyr::across(.cols = obs_vars, .fns = fn))
  central <- dplyr::mutate(ts$central,
                           dplyr::across(.cols = obs_vars, .fns = fn))

  return(
    tseries_boot(samples = samples,
                 central = central,
                 boot.R = ts$boot.R,
                 boot.l = ts$boot.l,
                 seed = ts$seed,
                 sim = ts$sim,
                 endcorr = ts$endcorr,
                 resampling_method = ts$resampling_method,
                 sample_idcs = ts$sample_idcs,
                 explanatory_vars = ts$explanatory_vars)
  )
}

#' @export
apply_binary_elementwise.tseries <- function(x, y, `%op%`) {
  stopifnot(inherits(x, 'tseries_boot'))
  stopifnot(inherits(y, 'tseries_boot'))

  res_compat <- resampling_is_compatible(x, y)
  stopifnot( res_compat$boot )
  stopifnot( res_compat$resampling_method )
  stopifnot( all(res_compat$explanatory_vars) )
  stopifnot( all(res_compat$samples_dims) )

  # In most applications, the two timeseries that we're combining come
  # from the same statistical ensemble.
  # However, we may also be in the situation where we want to combine
  # timeseries from different statistical ensembles.
  # In those cases, it is only the dimensions of the resampling
  # samples which are relevant.
  # Other details such as `boot.l`, `sim`, `endcorr` get lost.
  samples <- x$samples

  obs_vars <- dplyr::setdiff(colnames(samples), c(x$explanatory_vars, "sample_idx"))

  samples[,obs_vars] <- samples[,obs_vars] %op% y$samples[,obs_vars]

  central <- x$central
  central[,obs_vars] <- central[,obs_vars] %op% y$central[,obs_vars]

  # we need this style because `ifelse` can't deal with NULL in either
  # of its branches
  if( all(res_compat$sample_idcs) ){
    res_sample_idcs <- x$sample_idcs
  } else {
    res_sample_idcs <- NULL
  }

  return(
    tseries_boot(boot.R = x$boot.R,
                 boot.l = ifelse(res_compat$boot.l, x$boot.l, NA),
                 sim = ifelse(res_compat$sim, x$sim, NA),
                 endcorr = ifelse(res_compat$endcorr, x$endcorr, NA),
                 seed = ifelse(res_compat$seed, x$seed, NA),
                 samples = samples,
                 central = central,
                 explanatory_vars = x$explanatory_vars,
                 sample_idcs = res_sample_idcs,
                 resampling_method = x$resampling_method)
  )
}

#' @export
apply_reduce_plan.tseries <- function(ts, reduce_vars, plan){ 
  stopifnot( inherits(ts, "tseries_boot") | inherits(ts, "tseries_orig") )
  stopifnot( all( reduce_vars %in% ts$explanatory_vars ) )

  remaining_explanatory_vars <- dplyr::setdiff(ts$explanatory_vars, reduce_vars)
  
  # this function can work either on the original data or the samples
  # and central values
  # on the original data, the reduction is applied *per measurement*
  # on the samples, the reduction is applied *per sample*
  # of course, it is the caller who has to decide if a *per measurement*
  # evaluation is a valid thing to do
  dat <- list()
  if( inherits(ts, "tseries_boot") ){
    sample_vars <- c("sample_idx", remaining_explanatory_vars)
    # group_by_at doesn't mind NULL or empty groupings
    dat[["central"]] <- dplyr::group_by_at(ts$central, remaining_explanatory_vars)
    dat[["samples"]] <- dplyr::group_by_at(ts$samples, sample_vars)
  }
  if( inherits(ts, "tseries_orig") ){
    md_vars <- c("md_idx", remaining_explanatory_vars)
    dat[["orig"]] <- dplyr::group_by_at(ts$data, md_vars)
  }
  
  out <- lapply(
    X = dat,
    FUN = function(x){
      # for each element of the reduction plan, we create a new summary variable
      # accordining to the name of the element of the plan list
      # and the expression for this element
      res_lst <- lapply(
        X = names(plan),
        FUN = function(out_var){
          dplyr::summarise(
            x,
            !!as.character(out_var) := eval(plan[[out_var]]),
            .groups = "drop"
          )
        }
      )
      Reduce(auto_full_join,res_lst)
    }
  )
  names(out) <- names(dat)

  # special treatment here because ifelse cannot deal with `NULL` in any of its branches
  if( length(remaining_explanatory_vars) == 0 ){
    res_explanatory_vars <- NULL
  } else {
    res_explanatory_vars <- remaining_explanatory_vars
  }

  if( !inherits(ts, "tseries_boot") ){
    return( tseries_orig(data = out$orig,
                         explanatory_vars = res_explanatory_vars) )
  } else {
    if( inherits(ts, "tseries_orig") ){
      ts_base <- tseries_orig(data = out$orig,
                              explanatory_vars = res_explanatory_vars)
    } else {
      ts_base = tseries()
    }
    return(
      tseries_boot(.tseries = ts_base,
                   boot.R = ts$boot.R,
                   boot.l = ts$boot.l,
                   seed = ts$seed,
                   sim = ts$sim,
                   endcorr = ts$endcorr,
                   samples = out$samples,
                   central = out$central,
                   resampling_method = ts$resampling_method,
                   explanatory_vars = res_explanatory_vars,
                   sample_idcs = ts$sample_idcs)
    )
  }
}

#' @export
apply_transmute_plan.tseries <- function(ts, plan){ 
  stopifnot( inherits(ts, "tseries_boot") )

  dat <- list()
  dat[["central"]] <- ts$central
  dat[["samples"]] <- dplyr::group_by_at(ts$samples, "sample_idx")

  # we save some typing here by using a list for the central values and samples
  # for each element of this list, we apply the mutation plan
  out <- lapply(
    X = dat,
    FUN = function(x){
      # for each element of the reduction plan, we create a new variable
      # accordining to the name of the element of the plan list
      # and the expression for this element
      res_lst <- lapply(
        X = names(plan),
        FUN = function(out_var){
          dplyr::ungroup(
            dplyr::transmute(
              x,
              !!as.character(out_var) := eval(plan[[out_var]])
            )
          )
        }
      )
      Reduce(auto_full_join,res_lst)
    }
  )
  names(out) <- names(dat)

  return(
    tseries_boot(boot.R = ts$boot.R,
                 boot.l = ts$boot.l,
                 seed = ts$seed,
                 sim = ts$sim,
                 endcorr = ts$endcorr,
                 samples = out$samples,
                 central = out$central,
                 resampling_method = ts$resampling_method,
                 explanatory_vars = ts$explanatory_vars,
                 sample_idcs = ts$sample_idcs)
  )
}

apply_binary_mutate_plan.tseries <- function(x, y, `%op%`, plan){
  stopifnot(inherits(x, 'tseries_boot'))
  stopifnot(inherits(y, 'tseries_boot'))
  stopifnot( all(c("vars_out", "x_vars", "y_vars") %in% names(plan)) )

  scale1 <- rep(1.0, nrow(plan))
  if( "x_scale" %in% names(plan) ){
    scale1 <- plan$x_scale
  }
  scale2 <- rep(1.0, nrow(plan))
  if( "y_scale" %in% names(plan) ){
    scale2 <- plan$y_scale
  }
   
  res_compat <- resampling_is_compatible(x, y)
  stopifnot( res_compat$boot )
  stopifnot( res_compat$resampling_method )
  stopifnot( all(res_compat$explanatory_vars) )

  stopifnot( all(plan$x_vars %in% names(x$samples)) )
  stopifnot( all(plan$y_vars %in% names(y$samples)) )
  
  stopifnot( all(plan$x_vars %in% names(x$central)) )
  stopifnot( all(plan$y_vars %in% names(y$central)) )

  sample_vars <- c("sample_idx", x$explanatory_vars)
  samples <- x$samples[,sample_vars]
  central <- x$central[,x$explanatory_vars]

  for( idx in 1:nrow(plan) ){
    samples <- 
      dplyr::mutate(
        samples, 
        !!as.character(plan$vars_out[idx]) := (scale1[idx] * x$samples[[ plan$x_vars[idx] ]]) %op% 
                                              (scale2[idx] * y$samples[[ plan$y_vars[idx] ]])
      )

    central <- 
      dplyr::mutate(
        central,
        !!as.character(plan$vars_out[idx]) := (scale1[idx] * x$central[[ plan$x_vars[idx] ]]) %op% 
                                              (scale2[idx] * y$central[[ plan$y_vars[idx] ]])
      )
  }
  
  # we need this style because `ifelse` can't deal with NULL in either
  # of its branches
  if( all(res_compat$sample_idcs) ){
    res_sample_idcs <- x$sample_idcs
  } else {
    res_sample_idcs <- NULL
  }
  return(
    tseries_boot(boot.R = x$boot.R,
                 boot.l = ifelse(res_compat$boot.l, x$boot.l, NA),
                 sim = ifelse(res_compat$sim, x$sim, NA),
                 endcorr = ifelse(res_compat$endcorr, x$endcorr, NA),
                 seed = ifelse(res_compat$seed, x$seed, NA),
                 samples = samples,
                 central = central,
                 explanatory_vars = x$explanatory_vars,
                 sample_idcs = res_sample_idcs,
                 resampling_method = x$resampling_method)
  )
}

#' @export
add.tseries <- function(x, y, plan){
  apply_binary_mutate_plan.tseries(x, y, `+`, plan)
}

#' @export
subtract.tseries <- function(x, y, plan){
  apply_binary_mutate_plan.tseries(x, y, `-`, plan)
}

#' @export
mul.tseries <- function(x, y, plan){
  apply_binary_mutate_plan.tseries(x, y, `*`, plan)
}

#' @export
div.tseries <- function(x, y, plan){
  apply_binary_mutate_plan.tseries(x, y, `/`, plan)
}

#' @export
subset.tseries <- function(x, subset) {
  stopifnot( inherits(x, "tseries") )
  stopifnot( inherits(x, "tseries_orig") | inherits(x, "tseries_boot") )

  x_new <- x

  if( inherits(x, "tseries_orig") ){
    cols_keep <- c("md_idx", x$explanatory_vars, subset)
    stopifnot( all(subset %in% colnames(x$data) ) )
    x_new$data <- x$data[,cols_keep] 
  }

  if( inherits(x, "tseries_boot") ){
    stopifnot( all(subset %in% colnames(x$samples) ) )
    stopifnot( all(subset %in% colnames(x$central) ) )
    stopifnot( all(subset %in% colnames(x$se) ) )
    
    cols_keep <- c(x$explanatory_vars, subset)
    x_new$central <- x$central[,cols_keep]
    x_new$se <- x$se[,cols_keep]

    cols_keep <- c("sample_idx", cols_keep)
    x_new$samples <- x$samples[,cols_keep]
  }

  return(x_new)
}

#' plot_timeseries
#' 
#' @description
#' function to plot timeseries data, a corresponding histogram
#' and an error shading for an error analysis via uwerr
#'
#' @param dat Timeseries to analyse.
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param plotsize Width and Height of plot.
#' @param titletext Text in the plot title.
#' @param hist.by Numeric. Stepping to compute the histogram breaks.
#' @param stat_range Optional integer vector of length 2. Start and end indices
#'        of the subset of `dat` to be plotted. If left empty, all of `dat` will be
#'        plotted.
#' @param pdf.filename String. PDF filename.
#' @param name String. Timeseries name.
#' @param hist.probs Optional numeric vector of length 2. Probability extrema to limit the width
#'        of the histogram or smoothed density plots. By default all data is used. Note: this
#'        has not effect on the analysis as a whole or other plots.
#' @param smooth_density Boolean. Instead of plotting a histogram, use a smoothed density.
#' @param errorband_color String. Colour of the error band.
#' @param type String. Plot type, see \link{plot} for details.
#' @param uwerr.S Numeric. `S` of the \link{uwerr} method to be used.
#' @param time_factor Numeric. Factor by which any auto-correlation times shoud be
#'                             multiplied. Used, for example, when running non-unit-length 
#'                             trajectories or employing strided measurements.
#' @param periodogram Boolean. Whether to show a periodogram.
#' @param debug Boolean. Generate debug output.
#' @param uw.summary Boolean. Generate an \link{uwerr} summary plot.
#' @param ... Generic graphical parameters to be passed on.
#'
#' @return
#' Returns a \link{data.frame} with named columns `val`, `dval`, `tauint`, `dtauint`, `Wopt`
#' and `stringsAsFactors`, see \link{uwerr}.
#' 
#' @export
plot_timeseries <- function(dat, 
                            ylab, plotsize, titletext, hist.by,
                            stat_range = c(1.0, length(dat$y)),
                            pdf.filename,
                            name="", xlab="$t_\\mathrm{MD}$", 
                            hist.probs=c(0.0,1.0), errorband_color=rgb(0.6,0.0,0.0,0.6),
                            type='l',
                            uwerr.S=2,
                            time_factor=1.0,
                            smooth_density=FALSE,
                            periodogram=FALSE,debug=FALSE,uw.summary=TRUE,...) {

  stopifnot(length(stat_range) == 2)
  stopifnot(length(hist.probs) == 2)

  yrange <- range(dat$y)

  stat_y <- dat$y[ seq(stat_range[1], stat_range[2]) ]
 
  dat$t <- dat$t*time_factor

  uw.data <- uwerrprimary(stat_y, S=uwerr.S)
  if(debug) {
    print(paste("uw.",name,sep=""))
    print(summary(uw.data))
  }
  
  tikzfiles <- NULL
  if( !missing(pdf.filename) ){
    tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  }

  op <- par(family="Palatino",cex.main=0.8,font.main=1)
  on.exit(par(op))
  par(mgp=c(2,1.0,0))

  # plot the timeseries
  plot(x=dat$t,
       xlim=range(dat$t),
       y=dat$y,
       ylab=ylab,
       type=type,
       xlab=xlab,
       main=titletext, ...)

  rect(xleft=dat$t[ stat_range[1] ],
       xright=dat$t[ stat_range[2] ],
       ytop=uw.data$value+uw.data$dvalue,
       ybottom=uw.data$value-uw.data$dvalue,
       border=FALSE, col=errorband_color)
  abline(h=uw.data$value,col="black",lwd=2)

  legend(x="topright",
         legend=sprintf("%s $=%s$",
                         ylab,
                         tex.catwitherror(x=uw.data$value,
                                          dx=uw.data$dvalue,
                                          digits=2,
                                          with.dollar=FALSE,
                                          with.cdot=FALSE)
                         ),
         lty=1,
         pch=NA,
         col="red",
         bty='n',
         cex = 0.8)
  
  # plot the corresponding histogram
  hist.data <- NULL
  
  hist.breaks <- floor( ( max(stat_y)-min(stat_y) ) / uw.data$dvalue )
  if(!missing(hist.by)){
    hist.breaks <- floor( ( max(stat_y)-min(stat_y) ) / hist.by )
  } else {
    if(hist.breaks < 10 || hist.breaks > 150){
      hist.breaks <- 70
    }
  }
 
  if( smooth_density ){
    # determine an appropriate bandwidth for the density estimate
    # by finding the closest power of 2 smaller than half the
    # number of measurements
    n <- ifelse(length(stat_y)/2 > 512, 
                2^floor(log2(length(stat_y)/2)),
                512)
    d <- density(stat_y, bw = "SJ", n = n)
    plot(d, 
         xlim=quantile(stat_y,probs=hist.probs),
         main=titletext,
         xlab=ylab, 
         breaks=hist.breaks)
    ytop <- max(d$y)
    ybottom <- 0.0
  } else { 
    hist.data <- hist(stat_y,
                      xlim=quantile(stat_y,probs=hist.probs),
                      main=titletext,
                      xlab=ylab, 
                      breaks=hist.breaks)
    ytop <- max(hist.data$counts)
    ybottom <- 0.0
  }
  rect(ytop=ytop,
       ybottom=0,
       xright=uw.data$value+uw.data$dvalue,
       xleft=uw.data$value-uw.data$dvalue,
       border=FALSE, col=errorband_color)
  abline(v=uw.data$value, col="black", lwd=2) 

  # and a periodogram
  if(periodogram)
  {
    spec.pgram(x=stat_y,
               main=paste(ylab,paste("raw periodogram",titletext)))
  }

  # and the uwerr plots
  if(uw.summary){
    plot(uw.data, 
         main=paste(ylab,paste("UWErr analysis",titletext)),
         plot.hist=FALSE)
  }
   
  if(!missing(pdf.filename)){
    tikz.finalize(tikzfiles)
  }

  return(t(data.frame(val=uw.data$value, dval=uw.data$dvalue, tauint=uw.data$tauint*time_factor, 
                      dtauint=uw.data$dtauint*time_factor, Wopt=uw.data$Wopt*time_factor, stringsAsFactors=FALSE)))
}

#' plot_eigenvalue_timeseries
#' 
#' @description
#' function to plot timeseries of eigenvlues, including minimum and maximum eigenvalue bands 
#'  as found in the monomial_0x.data files produced by tmLQCD
#' 
#' @param dat Timeseries to analyse.
#' @param ylab Y-axis label.
#' @param plotsize Width and Height of plot.
#' @param filelabel String. Label of the file.
#' @param titletext Text in the plot title.
#' @param stat_range range of statistics to use.
#' @param time_factor Numeric. Factor by which any auto-correlation times shoud be
#'                             multiplied. Used, for example, when running non-unit-length 
#'                             trajectories or employing strided measurements.
#' @param pdf.filename String. PDF filename.
#' @param errorband_color String. Colour of the error band.
#' @param debug Boolean. Generate debug output.
#'
#' @return
#' Returns a list with two named elements `mineval` and `maxeval` for the minimal
#' and the maximal eigenvalue, see \link{plot_timeseries}.
#' 
#' @export
plot_eigenvalue_timeseries <- function(dat,
                                       stat_range,
                                       time_factor = 1.0,
                                       ylab, plotsize, filelabel,titletext,
                                       pdf.filename,
                                       errorband_color=rgb(0.6,0.0,0.0,0.6),
                                       debug=FALSE) {
  if( missing(stat_range) ) { stat_range <- c(1,nrow(dat)) }
  yrange <- range(dat[,2:5])

  dat$t <- dat$t*time_factor

  stat_min_ev <- dat[ seq(stat_range[1], stat_range[2]), "min_ev" ]
  stat_max_ev <- dat[ seq(stat_range[1], stat_range[2]), "max_ev" ]

  uw.min_ev <- uwerrprimary( stat_min_ev )
  uw.max_ev <- uwerrprimary( stat_max_ev )

  if(debug){
    print("uw.eval.min_ev")
    print(summary(uw.min_ev))
    print("uw.eval.max_ev")
    print(summary(uw.max_ev))
  }

  tikzfiles <- NULL
  if(!missing(pdf.filename)){
    tikzfiles <- tikz.init(basename=pdf.filename,width=plotsize,height=plotsize)
  }
  par_save <- par(mgp=c(2,1,0))
  on.exit(par(par_save))

  # plot the timeseries
  plot(x=dat$traj, xlim=range(dat$traj), 
       y=dat$min_ev, ylim=yrange, 
       t='l', 
       ylab=ylab, 
       xlab=expression(t[MD]), 
       main=titletext, log='y', tcl=0.02)

  lines(x=dat$traj, y=dat$max_ev)
 
  ## add the approximation interval
  lines(x=dat$traj, y=dat$ev_range_min, lty=2, col="darkgreen")
  lines(x=dat$traj, y=dat$ev_range_max, lty=2, col="darkgreen")
  
  # plot the corresponding histograms with error bands and mean values
  hist.min_ev <- hist(stat_min_ev, main=paste("min. eval",titletext),xlab="min. eval",tcl=0.02)
  rect(ytop=max(hist.min_ev$counts),
       ybottom=0,
       xright=uw.min_ev$value+uw.min_ev$dvalue,
       xleft=uw.min_ev$value-uw.min_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.min_ev$value,col="black")                                                                                                   

  hist.max_ev <- hist(stat_max_ev, main=paste("max. eval",titletext), xlab="max. eval")
  rect(ytop=max(hist.max_ev$counts),
       ybottom=0,
       xright=uw.max_ev$value+uw.max_ev$dvalue,
       xleft=uw.max_ev$value-uw.max_ev$dvalue,border=FALSE,col=errorband_color)
  abline(v=uw.max_ev$value,col="black")                                                                                                   

  if(!missing(pdf.filename)){
    tikz.finalize(tikzfiles)
  }

  return(list(mineval=t(data.frame(val=uw.min_ev$value, dval=uw.min_ev$dvalue,
                                   tauint=uw.min_ev$tauint*time_factor, 
                                   dtauint=uw.min_ev$dtauint*time_factor, 
                                   Wopt=uw.min_ev$Wopt*time_factor, stringsAsFactors=FALSE)),
              maxeval=t(data.frame(val=uw.max_ev$value, dval=uw.max_ev$dvalue,
                                   tauint=uw.max_ev$tauint*time_factor, 
                                   dtauint=uw.max_ev$dtauint*time_factor, 
                                   Wopt=uw.max_ev$Wopt*time_factor, stringsAsFactors=FALSE)) ) )
              
}

