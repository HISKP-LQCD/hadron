#' @title Generate HDF5 key for a momentum and spin-projected CVC loop 
#' @param loop_type String, loop type.
#' @param i_stoch Integer, stochastic sample id number.
#' @param gamma Integer, CVC convention gamma matrix identifier.
#' @param p Integer vector of length 3, (x,y,z) components of the momentum
#'          vector in lattice units.
cvc_local_loop_key <- function(loop_type, i_stoch, gamma, p)
{
  stopifnot( length(p) == 3 )

  sprintf("/%s/istoch_%04d/g%d/px%dpy%dpz%d",
          loop_type,
          i_stoch,
          gamma,
          p[1], p[2], p[3])
}

#' @title Generate key to identify a momentum and spin-projected loop 
#' @param loop_type String, loop type.
#' @param gamma Integer, CVC convention gamma matrix identifier.
#' @param p Integer vector of length 3, (x,y,z) components of the momentum
#'          vector in lattice units.
cvc_local_loop_key <- function(loop_type, gamma, p)
{
  stopifnot( length(p) == 3 )

  sprintf("/%s/g%d/px%dpy%dpz%d",
          loop_type,
          gamma,
          p[1], p[2], p[3])
}

#' @title Generate HDF5 key for CVC 'correlators' meson 2pt function
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param src_ts Integer, source time slice.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
correlators_key_meson_2pt <- function(fwd_flav, bwd_flav, src_ts, snk_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )

  sprintf("/%s+-g-%s-g/t%d/gf%d/pfx%dpfy%dpfz%d/gi%d/pix%dpiy%dpiz%d",
          bwd_flav,
          fwd_flav,
          src_ts,
          snk_gamma,
          snk_p[1],snk_p[2],snk_p[3],
          src_gamma,
          src_p[1],src_p[2],src_p[3])
}

#' @title Generate key string to identify a meson 2pt function
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
cf_key_meson_2pt <- function(fwd_flav, bwd_flav, snk_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  sprintf("/%s+-g-%s-g/gf%d/pfx%dpfy%dpfz%d/gi%d/pix%dpiy%dpiz%d",
          bwd_flav,
          fwd_flav,
          snk_gamma,
          snk_p[1], snk_p[2], snk_p[3],
          src_gamma,
          src_p[1],src_p[2],src_p[3])
}

#' @title Generate HDF5 key for CVC 'correlators' meson 3pt function with a local insertion
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param seq_flav String, "sequential" quark flavour identifier.
#' @param src_ts Integer, source time slice.
#' @param dt Integer, source-sink separation.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param cur_gamma Integer, CVC convention gamma matrix identified at the insertion.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
correlators_key_meson_3pt_local <- function(fwd_flav, bwd_flav, seq_flav, src_ts, dt, snk_gamma, cur_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  sprintf("/s%s%s+-g-%s-g/t%d/dt%d/gf%d/pfx%dpfy%dpfz%d/gc%d/gi%d/pix%dpiy%dpiz%d",
          bwd_flav,
          seq_flav,
          fwd_flav,
          src_ts,
          dt,
          snk_gamma,
          snk_p[1],snk_p[2],snk_p[3],
          cur_gamma,
          src_gamma,
          src_p[1],src_p[2],src_p[3])
}

#' @title Generate key to identify a meson 3pt funtion with a local insertion
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param seq_flav String, "sequential" quark flavour identifier.
#' @param dt Integer, source-sink separation.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param cur_gamma Integer, CVC convention gamma matrix identified at the insertion.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
cf_key_meson_3pt_local <- function(fwd_flav, bwd_flav, seq_flav, dt, snk_gamma, cur_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  sprintf("/s%s%s+-g-%s-g/dt%d/gf%d/pfx%dpfy%dpfz%d/gc%d/gi%d/pix%dpiy%dpiz%d",
          bwd_flav,
          seq_flav,
          fwd_flav,
          dt,
          snk_gamma,
          snk_p[1],snk_p[2],snk_p[3],
          cur_gamma,
          src_gamma,
          src_p[1],src_p[2],src_p[3])
}

#' @title Convert correlation funtion in CVC HDF5 format to 'raw_cf'
#' @description Given a numeric vector of alternating real and imaginary parts of a
#'              correaltion function, creates an object of class 'raw_cf' with
#'              a single measurement, inferring \code{Time} and internal dimensions
#'              from the passed numeric vector.
#' @param cf_dat Numeric vector of alternating real and imaginary parts of a
#'               correlation function.
cvc_to_raw_cf <- function(cf_dat)
{
  # idcs for real and imaginary parts
  ridcs <- seq(1,length(cf_dat),2)
  iidcs <- seq(2,length(cf_dat),2)

  # turn it into a complex 4D array of 'internal' size c(1,1)
  # the dimensions correspond to measurement, time, internal
  cf_dat <- array(complex(real = cf_dat[ridcs],
                          imaginary = cf_dat[iidcs]),
                          dim = c(1,length(ridcs),c(1,1)))

  cf <- raw_cf_data(cf = raw_cf_meta(Time=length(ridcs), dim=c(1,1)),
                    data = cf_dat)
  return(cf)
}


