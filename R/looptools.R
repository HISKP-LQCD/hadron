#' @title compute two-point correlation function between quark loops
#' @param loop_snk 'raw_cf' container with spin-projected quark loop
#'                 at sink (see \link{loop_spin_project})
#' @param loop_src 'raw_cf' container with spin-projected quark loop 
#'                 at source (see \link{loop_spin_project})
#' @return 'raw_cf' container with two-point function of these two quark loops.
#'         In the calculation, both averaging over all source locations and
#'         the average over all stochastic sample combinations are performed.
#' @export
loop_2pt <- function(loop_snk, loop_src){
  stopifnot( inherits(loop_snk, 'raw_cf_meta') )
  stopifnot( inherits(loop_snk, 'raw_cf_data') )
  stopifnot( inherits(loop_src, 'raw_cf_meta') )
  stopifnot( inherits(loop_src, 'raw_cf_data') )
  stopifnot( all(dim(loop_snk$data) == dim(loop_src$data)) )

  Time <- loop_snk$Time

  loop_src_avg <- loop_stochav(loop_src)

  cf_2pt <- raw_cf_meta(Time = Time,
                        dim = c(1,1))
  cf_2pt <- raw_cf_data(cf_2pt,
                        data = array(as.complex(0.0), dim = dim(loop_src_avg$data)))

  for( t_sep in 0:(Time-1) ){
    loop_snk_shift <- shift.raw_cf(loop_snk, t_sep)
    loop_snk_avg_shift <- loop_stochav(loop_snk_shift)
    
    # Compute (\sum_i Tr[ \Gamma_snk M_i ])*(\sum_j Tr[ \Gamma_src M_j ])
    c_sep_cross <- loop_src_avg * loop_snk_avg_shift
    
    # Compute \sum_i Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_i ]
    c_sep <- loop_src * loop_snk_shift
    c_sep_diag <- loop_stochav(c_sep)

    # mean over all sources for this source-sink separation
    # at the same time, use the identity
    # \sum_{i != j} Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_j ] =
    # (\sum_i Tr[ \Gamma_snk M_i ])*(\sum_j Tr[ \Gamma_src M_j ]) -
    #  \sum_i Tr[ \Gamma_snk M_i ] Tr[ \Gamma_src M_i ]
    cf_2pt$data[,(t_sep+1),,] <- (1.0/Time)*( apply(c_sep_cross$data, c(1,3,4), mean) - 
                                              apply(c_sep_diag$data, c(1,3,4), mean) )
  }
  return(cf_2pt)
}

#' @title average over stochastic samples of loop
#' @description Perform mean over the third dimension of the loop data.
#' @param loop 'raw_cf' container with loop data
#' @export
loop_stochav <- function(loop){
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( length(dim(loop$data)) == 5 )
  loop$data <- apply( loop$data, c(1,2,4,5), mean )
  loop$dim <- c(4,4)
  return(loop)
}

#' @title spin projection of quark loop data
#' @description Implements the operation
#'   \deqn{ L = a*( \Gamma_{ik} M_{ki} ) }
#' to give the trace of a quark loop \eqn{M} multiplied by a gamma structure \eqn{\Gamma}
#' and scaled by a complex factor \eqn{a}.
#' @param loop 'raw_cf' container with loop data
#' @param gamma 4x4 complex matrix
#' @param stochav Boolean, specified whether the average over stochastic samples should be
#'                   performed. This makes the projection much faster
#'                   but of course prevents the projected loop data to be
#'                   used for the construction of diagrams with multiple quark loops.
#' @param factor Complex scaling factor to be applied. 
#' @export
loop_spin_project <- function(loop, gamma, stochav = FALSE, factor = as.complex(1.0) ){
  stopifnot( inherits(loop, 'raw_cf_meta') )
  stopifnot( inherits(loop, 'raw_cf_data') )
  stopifnot( all(dim(gamma) == 4) & length(dim(gamma)) == 2 )
  stopifnot( length(dim(loop$data)) == 5 )
  stopifnot( length(factor) == 1 )
  
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
  proj_loop <- factor*array(apply(X = loop$data,
                                  MARGIN = c(1:(length(dims)-2)),
                                  FUN = function(x){
                                          sum(diag( gamma %*% x ))
                                        }),
                            dim = c(dims[1:(length(dims)-2)],1,1)
                           )
  proj_dims <- dim(proj_loop)

  rval <- raw_cf_meta(Time = loop$Time,
                      dim = proj_dims[3:length(proj_dims)] )
  rval <- raw_cf_data(rval, data = proj_loop)
  return(rval)
}
