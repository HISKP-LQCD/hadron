#' @title disconnected contribution to current insertion three-point function
#' @description Computes the quark-line disconnected contribution to a three-point function
#'              of the form
#'                \deqn{ C_3(t, \Delta t = t_{snk} - t_{src}) = C_2(t_{snk}, t_{src}) * L(t) }
#'              \eqn{\forall t} considering only the case t_{snk} > t_{src}.
#' @param cf_2pt 'raw_cf' container holding two-point part of three-point function
#'               in lattice-absolute coordinates (not relative to source!)         
#' @param loop 'raw_cf' container holding loop contribution, suitably spin-projected
#'             and averaged over stochastic samples.
#' @param src_ts Integer vector, the source time slices that were used for the computation
#'               of the two-point function in lattice-absolute coordinates.
#'               Must be of the same length as the number of measurements in \code{cf_2pt}.
#' @param dt Integer, the source-sink separation that should be computed.
#' @param reim_loop String, one of 'real', 'imag' or 'both'. Specifies whether
#'                  just the real or imaginary part should be considered when
#'                  constructing the correlation with the two-point function.
#' @param reim_2pt String, same as \code{reim_loop} but for the two-point
#'                 contribution to the three-point function.
#' @param vev_subtract Boolean, whether the loop contains a vev which should
#'                     be subtracted.
#' @return `raw_cf` container with the product of loop and 2pt function, shifted
#'         in time to be relative to source using the info from \code{src_ts}
#'
#' @export
disc_3pt <- function(cf_2pt, 
                     loop,
                     src_ts,
                     dt,
                     reim_loop = 'both',
                     reim_2pt = 'both',
                     vev_subtract = FALSE)
{
  stopifnot( inherits(cf_2pt, 'raw_cf_data') )
  stopifnot( inherits(cf_2pt, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( cf_2pt$nts == loop$nts )
  stopifnot( cf_2pt$nrStypes == 1 )
  stopifnot( cf_2pt$nrObs == 1 )

  if( dim(cf_2pt$data)[1] %% dim(loop$data)[1] != 0 ){
    stop(paste("The number of measurements of the 2pt function must be",
               "divisible by the number of measurements of the loop contribution."))
  }
  
  if( length(dim(loop$data)) == 5 ){
    stop("The loop data should already be averaged over stochastic samples")
  }
  if( length(loop$dim) == 2 ){
    if( any( loop$dim != 1 ) ){
      stop("The loop data should be spin-projected")
    }
  }
  if( length(src_ts) != dim(cf_2pt$data)[1] ){
    stop(paste("The number of source times 'src_ts' provided should be",
               "the same as the number of measurements for the 2pt function"))
  }
  if( length(dt) != 1 ){
    stop("One single source-sink separation 'dt' should be provided!")
  }
  if( cf_2pt$Time != cf_2pt$nts ){
    stop("For this calculation, 'cf_2pt$Time == cf_2pt$nts' is required!")
  }
  
  if( !(reim_loop %in% c('real','imag','both')) | 
      !(reim_2pt %in% c('real', 'imag', 'both')) ){
    stop("'reim_loop' and 'reim_2pt' can only be one of 'real', 'imag' or 'both'")
  }
  
  if( reim_loop == 'real' ){
    reim_loop_fn <- function(x){ Re(x) }
  } else if ( reim_loop == 'imag' ) {
    reim_loop_fn <- function(x){ Im(x) }
  } else {
    reim_loop_fn <- function(x){ x }
  }
  
  if( reim_2pt == 'real' ){
    reim_2pt_fn <- function(x){ Re(x) }
  } else if ( reim_2pt == 'imag' ) {
    reim_2pt_fn <- function(x){ Im(x) }
  } else {
    reim_2pt_fn <- function(x){ x }
  }
  
  if( vev_subtract ){
    loop$data <- loop$data - mean(loop$data)
  }

  # in general, multiple source locations have been used per configuration
  # to compute the 2pt function
  stat_rep <- dim(cf_2pt$data)[1] / dim(loop$data)[1]

  # we need to replicate the loop measurements as many times as we have
  # measurements of the 2pt function per gauge configuration
  data <- loop$data[unlist(lapply(1:dim(loop$data)[1],
                                  function(x){ rep(x, stat_rep) })), , , ,drop=FALSE]

  # and now we correlate the loop with the 2pt function for this particular source
  # sink separation on each configuration
  # note that we don't shift the data to be relative to source, so this still has
  # to be done later
  for( i in 1:length(src_ts) ){
    data[i, , , ] <- reim_loop_fn(data[i, , , ,drop=FALSE]) * 
                       reim_2pt_fn( cf_2pt$data[i, ((src_ts[i]+cf_2pt$Time+dt) %% cf_2pt$Time) + 1, 1, 1, drop=TRUE] )
  }
  
  cf_3pt <- raw_cf_meta(Time=cf_2pt$Time, dim=c(1,1))
  cf_3pt <- raw_cf_data(cf_3pt,
                        data = data
                        )
  cf_3pt <- shift.raw_cf(cf_3pt, places = src_ts)

  return(cf_3pt)
}

#' @title compute two-point correlation function between quark loops
#' @param loop_snk 'raw_cf' container with spin-projected quark loop
#'                 at sink (see \link{loop_spin_project})
#' @param loop_src 'raw_cf' container with spin-projected quark loop 
#'                 at source (see \link{loop_spin_project})
#' @param random_vectors_outer_product Boolean. 
#'               For testing purposes, the average over all
#'               random sample combinations can be carried out explicitly
#'               as \eqn{ \sum_{i \neq j} Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_j ] }
#'               instead of the (much faster) equivalent
#'               \eqn{ (\sum_i Tr[ \Gamma_snk M_i ])*(\sum_j Tr[ \Gamma_src M_j ]) -
#'                     \sum_i ( Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_i ] ) }.
#' @param nstoch_to_avg String or integer, how many of the available stochastic
#'                      samples should be averaged over. See \link{loop_stochav} for
#'                      details.
#' @return 'raw_cf' container with two-point function of these two quark loops.
#'         In the calculation, both averaging over all source locations and
#'         the average over all stochastic sample combinations are performed.
#' @export
loop_2pt <- function(loop_snk,
                     loop_src,
                     random_vectors_outer_product = FALSE,
                     nstoch_to_avg = 'all'){
  stopifnot( inherits(loop_snk, 'raw_cf_meta') )
  stopifnot( inherits(loop_snk, 'raw_cf_data') )
  stopifnot( inherits(loop_src, 'raw_cf_meta') )
  stopifnot( inherits(loop_src, 'raw_cf_data') )
  stopifnot( all(dim(loop_snk$data) == dim(loop_src$data)) )
  if( length(dim(loop_snk$data)) != 5 ){
    stop("'loop_snk' should have 5 dimensions")
  }
  if( length(dim(loop_src$data)) != 5 ){
    stop("'loop_src' should have 5 dimensions")
  }

  Time <- loop_snk$Time

  if( nstoch_to_avg == "all" ){
    nstoch_to_avg <- loop_src$dim[1]
  }
  if( nstoch_to_avg >= 1 && nstoch_to_avg <= loop_src$dim[1] ){
    loop_src$data <- loop_src$data[,,c(1:nstoch_to_avg),,,drop=FALSE]
    loop_src$dim[1] <- nstoch_to_avg
    loop_snk$data <- loop_snk$data[,,c(1:nstoch_to_avg),,,drop=FALSE]
    loop_snk$dim[1] <- nstoch_to_avg
  } else {
    stop("The only valid values for 'nstoch_to_avg' are 'all' or an integer smaller or equal to 'loop_src$dim[1]'")
  }

  loop_src_avg <- loop_stochav(loop_src)

  cf_2pt <- raw_cf_meta(Time = Time,
                        dim = c(1,1))
  cf_2pt <- raw_cf_data(cf_2pt,
                        data = array(as.complex(0.0), dim = dim(loop_src_avg$data)))

  dims <- loop_snk$dim

  random_combinations <- t(combn(x=1:dims[1], m=2))
  if( random_vectors_outer_product ){
    random_idcs <- data.frame(r1=random_combinations[,1], 
                              r2=random_combinations[,2])
  }

  for( t_sep in 0:(Time-1) ){
    loop_snk_shift <- shift.raw_cf(loop_snk, t_sep)
    # for testing, we can compute the 2pt function explicitly using i != j
    # random combinations
    # of course, this is extremely slow
    if( random_vectors_outer_product ){
      cf_2pt$data[,(t_sep+1),,] <- apply(loop_src$data[,,random_idcs$r1,,,drop=FALSE] * 
                                         loop_snk_shift$data[,,random_idcs$r2,,,drop=FALSE], 
                                         c(1,4,5), 
                                         mean)
    } else {
      loop_snk_avg_shift <- loop_stochav(loop_snk_shift)
      
      # Compute 1/N_r^2 (\sum_i Tr[ \Gamma_snk M_i ])*(\sum_j Tr[ \Gamma_src M_j ])
      c_sep_full <- loop_src_avg * loop_snk_avg_shift

      # Compute 1/N_r \sum_i ( Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_i ] )
      c_sep <- loop_src * loop_snk_shift
      c_sep_diag <- loop_stochav(c_sep)

      # mean over all sources for this source-sink separation
      # at the same time, use the identity
      # \sum_{i != j} Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_j ] =
      # (\sum_i Tr[ \Gamma_snk M_i ])*(\sum_j Tr[ \Gamma_src M_j ]) -
      #  \sum_i ( Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_i ] )
      # note that since we work with averages (rather than sums), we need to account
      # for the relative normalisation of the full and diagonal random vector products
      # as well as the relative weight of the full N_r^2 combinations
      # and the unique C(N_r,2) combinations, such that the mean that we take
      # in the end is over those
      cf_2pt$data[,(t_sep+1),,] <- 0.5/nrow(random_combinations)*
                                      apply(dims[1]^2*c_sep_full$data - dims[1]*c_sep_diag$data, 
                                            c(1,3,4),
                                            mean) 
    }
  }
  return(cf_2pt)
}

#' @title average over stochastic samples of loop
#' @description Perform mean over the third dimension of the loop data.
#' @param loop 'raw_cf' container with loop data
#' @param nstoch_to_avg String or integer, number of stochastic samples 
#'                      to average over. Only possible string is 'all'.
#'                      If an integer is supplied, it must be at least '1'
#'                      and at most consistent with the number of stochastic
#'                      samples in \code{loop}.
#' @return
#' Returns the input `loop` object with named elemens `data` and
#' `dim` added.
#' 
#' @export
loop_stochav <- function(loop, nstoch_to_avg = 'all'){
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( length(dim(loop$data)) == 5 )
  dims <- dim(loop$data)

  if(is.character(nstoch_to_avg)){
    if(nstoch_to_avg != 'all'){
      stop("The only string input supported fro 'nstoch_to_avg' is 'all'")
    }
  }
  if(nstoch_to_avg == 'all'){
    nstoch_to_avg <- dims[3]
  } else if( nstoch_to_avg < 1 | nstoch_to_avg > dims[3] ){
    message <- sprintf(paste("invalid value for 'nstoch_to_avg': %d. Must be at least '1'",
                             "and at most the number of stochastic samples in 'loop': %d"),
                       nstoch_to_avg, dims[3])
    stop(message)
  }
  loop$data <- apply( loop$data[,,1:nstoch_to_avg,,,drop=FALSE], c(1,2,4,5), mean )
  loop$dim <- dims[c(4,5)]
  return(loop)
}

#' @title subtract vev from loop data
#' @description Convenience function to subtract any possible vacuum-expectation
#'              value from a loop matrix. The expectation value of each component
#'              of the internal dimensions is subtracted individually.
#'              Averaging over stochstic samples can be restricted to a subset, see
#'              \code{nstoch_to_avg} input parameter.
#' @param loop 'raw_cf' container with loop data
#' @param nstoch_to_avg String or integer, number of stochastic samples 
#'                      to average over. Only possible string is 'all'. If an
#'                      integer is provided it must be at least '1' and at most
#'                      consistent with the number of stochastic samples in \code{loop}.
#'
#' @return
#' Returns the input `loop` object with added data.
#' 
#' @export
loop_vev_subtract <- function(loop, nstoch_to_avg = 'all'){
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( length(dim(loop$data)) == 5 )

  dims <- dim(loop$data)
  if(is.character(nstoch_to_avg)){
    if(nstoch_to_avg != 'all'){
      stop("The only string input supported fro 'nstoch_to_avg' is 'all'")
    }
  }
  if(nstoch_to_avg == 'all'){
    nstoch_to_avg <- dims[3]
  } else if( nstoch_to_avg < 1 | nstoch_to_avg > dims[3] ){
    message <- sprintf(paste("invalid value for 'nstoch_to_avg': %d. Must be at least '1'",
                             "and at most the number of stochastic samples in 'loop': %d"),
                       nstoch_to_avg, dims[3])
    stop(message)
  }

  vevs <- apply(loop$data[ , ,1:nstoch_to_avg,,,drop=FALSE], c(4,5), mean )
  for( d1 in 1:dims[4] ){
    for( d2 in 1:dims[5] ){
      loop$data[1:dims[1], 1:dims[2], 1:dims[3], d1, d2] <- loop$data[1:dims[1], 1:dims[2], 1:dims[3], d1, d2] - vevs[d1,d2]
    }
  }
  return(loop)
}
  

#' @title spin projection of quark loop data
#' @description Implements the operation
#'   \deqn{ L = a*( \Gamma_{ik} M_{ki} ) }
#' to give the trace of a quark loop \eqn{M} multiplied by a gamma structure \eqn{\Gamma}
#' and scaled by a complex factor \eqn{a}.
#' @param loop 'raw_cf' container with loop data
#' @param gamma 4x4 complex matrix
#' @param reim String, one of 'real', 'imag' or 'both'. After the spin projection and trace,
#'             the result can be restricted to just the real or imaginary part, if desired.
#'             Useful for the cases in which it is clear that only one or the other contains
#'             any signal.
#' @param stochav Boolean, specifies whether the average over stochastic samples should be
#'                performed. This makes the projection much faster
#'                but of course prevents the projected loop data to be
#'                used for the construction of diagrams with multiple quark loops.
#' @param scale_factor Complex scaling factor to be applied.
#' @param herm_conj Boolean, optionally the loop matrix \eqn{M} can be hermitian
#'                  conjugated before the spin projection is performed.
#'
#' @return
#' Returns an object of class \link{raw_cf}.
#' 
#' @export
loop_spin_project <- function(loop,
                              gamma,
                              reim = 'both',
                              stochav = FALSE, 
                              scale_factor = as.complex(1.0),
                              herm_conj = FALSE){
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( all(dim(gamma) == 4) & length(dim(gamma)) == 2 )
  stopifnot( length(factor) == 1 )
  
  if( stochav ){
    stopifnot( length(dim(loop$data)) == 5 )
  } else {
    stopifnot( length(dim(loop$data)) == 4 | length(dim(loop$data)) == 5 )
  }

  if( !(reim %in% c('real','imag','both') ) ){
    stop("'reim' can only be one of 'real', 'imag' or 'both'")
  }

  ## in some cases, we want to extract just the real or imaginary
  ## part of the estimate of the loop trace
  if( reim == 'real' ){
    reim_fn <- function(x){ Re(x) }
  } else if ( reim == 'imag' ) {
    reim_fn <- function(x){ Im(x) }
  } else {
    reim_fn <- function(x){ x }
  }
  
  if( stochav ){
    loop <- loop_stochav(loop)
  }

  dims <- dim(loop$data)
  if( any( dims[ c(length(dims)-1, length(dims)) ] != 4 ) ){
    stop("The two right-most dimensions of 'loop$data' should both be '4'.")
  }

  # perform gamma projection
  # using apply along the dimensions that we want to preserve (i.e., anything
  # but the spin degrees of freedom) 
  proj_loop <- array(apply(X = loop$data,
                           MARGIN = c(1:(length(dims)-2)),
                           FUN = function(x){
                                   # the as.complex is required to retain the same data type going in
                                   # and coming out. If that is not done, proj_loop will
                                   # be coerced to something that is not an array.
                                   if( herm_conj ){
                                     return( as.complex(reim_fn(scale_factor*sum(diag( gamma %*% Conj(t(x)))) ) ) )
                                   }else{
                                     return( as.complex(reim_fn(scale_factor*sum(diag( gamma %*% x )) ) ) )
                                   }
                                 }),
                           dim = c(dims[1:(length(dims)-2)],1,1)
                           )

  proj_dims <- dim(proj_loop)

  rval <- raw_cf_meta(Time = loop$Time,
                      dim = proj_dims[3:length(proj_dims)] )
  rval <- raw_cf_data(rval, data = proj_loop)
  return(rval)
}

