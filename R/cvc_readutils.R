#' @title Generate HDF5 key for a momentum and spin-projected CVC loop 
#' @param loop_type String, loop type.
#' @param i_stoch Integer, stochastic sample id number.
#' @param gamma Integer, CVC convention gamma matrix identifier.
#' @param p Integer vector of length 3, (x,y,z) components of the momentum
#'          vector in lattice units.
cvc_local_loop_key <- function(loop_type, istoch, gamma, p)
{
  stopifnot( length(p) == 3 )

  sprintf("/%s/istoch_%04d/g%d/px%dpy%dpz%d",
          loop_type,
          istoch,
          gamma,
          p[1], p[2], p[3])
}

#' @title Generate key to identify a momentum and spin-projected loop 
#' @param loop_type String, loop type.
#' @param istoch Integer, index of the stochastic sample.
#' @param gamma Integer, CVC convention gamma matrix identifier.
#' @param p Integer vector of length 3, (x,y,z) components of the momentum
#'          vector in lattice units.
cvc_local_loop_key <- function(loop_type, istoch, gamma, p)
{
  stopifnot( length(p) == 3 )

  sprintf("/%s/istoch_%04d/g%d/px%dpy%dpz%d",
          loop_type,
          istoch,
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

#' @title read HDF5 loop files in the CVC loop format
#' @description The CVC naive_loops code produces HDF5 files which contain
#'              a matrix of momenta and the data for the loops (without
#'              spin projection) organised by stochastic sample. Currently, the
#'              reading code assumes that there is a single configuration stored per
#'              file and the "trajectory" parameter in CalcLoops is assumed
#'              to take its default value of '4'.
#' @param selections Named list with names from the list 'Naive', 'Scalar', 'dOp', 'Loops'
#'                   'LpsDw', 'LpsDwCv', 'LoopsCv' specifying the requesetd loop types. 
#'                   The elements of this list are in turn expected
#'                   be data frames of the form
#'                     \tabular{rrr}{
#'                       \strong{qx} \tab \strong{qy} \tab \strong{qz} \cr
#'                       0           \tab 0           \tab 1           \cr
#'                       -2          \tab 1           \tab -3          \cr
#'                       ...         \tab ...         \tab ...
#'                     }
#'                   specifying the momentum combinations to be extracted for each
#'                   loop type.
#' @param files Vector of strings, list of HDF5 files to be processed.
#' @param Time Integer, time extent of the lattice.
#' @param nstoch Integer, number of stochastic samples to be expected in file.
#' @param verbose Boolean, output I/O time per file. Requires 'tictoc' package.
#' @return Named nested list of the same length as \code{selections} containg the loop data
#'         in the \link{raw_cf} format. Each named element corresponds to one loop
#'         type and each element of the underlying numbered list corresponds to one momentum
#'         combination as specified via \code{selections} for this loop type in the same order.
#'         
#' @export
cvc_read_loops <- function(selections, files, Time, nstoch, verbose = FALSE){
  rhdf5_avail <- requireNamespace("rhdf5")
  dplyr_avail <- requireNamespace("dplyr")
  if( !rhdf5_avail | !dplyr_avail ){
    stop("The 'dplyr' and 'rhdf5' packages are required to use this function!\n")
  }
  if( verbose ){
    tictoc_avail <- requireNamespace("tictoc")
    if( !tictoc_avail ){
      stop("Time reporting requires the 'tictoc' package!")
    }
  }

  rval <- list()
  selected_loop_types <- names(selections)

  if( any( !(selected_loop_types %in% c("naive") ) ) ){
    stop("Currently only the 'naive' loop types is supported!")
  }

  for( ifile in 1:length(files) ){
    f <- files[ifile]
    if(verbose) tictoc::tic(sprintf("Reading %s",f))
    
    # The file names are of the form 
    # path/loops.XXXX.h5
    # and we want to recover XXXX
    cid_in_filename <- as.integer(strsplit(basename(f), split = ".", fixed = TRUE)[[1]][2])

    h5f <- rhdf5::H5Fopen(f, flags = "H5F_ACC_RDONLY")
    
    group_names <- h5ls(h5f)$name
    
    avail_loop_types <- unlist( lapply( selected_loop_types, function(x){ x %in% group_names } ) )
    if( any( !avail_loop_types ) ){
      msg <- sprintf("Some selected loop types could not be found in %s :\n %s",
                     f,
                     do.call(paste, as.list( selected_loop_types[!avail_loop_types] ) )
                     )
      stop(msg)
    }
    
    # how many stochastic samples are available and does it match out expectation?
    stoch_avail <- sort(as.numeric(unlist(lapply(X = strsplit(unique(group_names[ grepl("istoch", group_names) ]),
                                       "_"),
                                 FUN = function(x){ x[2] })
                          )))
    if( length(stoch_avail) != nstoch ){
      stop(sprintf("Number of stochastic samples in file %s :\n%d, expected %d!",
                   f, length(stoch_avail), nstoch))
    }

    for( loop_type in selected_loop_types ){
      if( !(loop_type %in% names(rval) ) ){
        rval[[ loop_type ]] <- list()
      }
      
      selected_moms <- selections[[ loop_type ]][, c("qx","qy","qz")]
      for( igamma in 1:16 ){
        gamma <- igamma-1
        if( length(rval[[ loop_type ]]) < igamma ){
          rval[[loop_type]][[igamma]] <- list()
        }
        for( imom in 1:nrow(selected_moms) ){
          if( length(rval[[loop_type]][[igamma]]) < imom){
            rval[[loop_type]][[igamma]][[imom]] <- array(as.complex(NA),
                                                         dim = c(length(files),
                                                                 Time,
                                                                 nstoch,
                                                                 1, 1))
          }
          for( istoch in 1:nstoch ){
            key <- cvc_local_loop_key(loop_type = loop_type,
                                      istoch = istoch-1,
                                      gamma = gamma,
                                      p = as.integer(selected_moms[imom,]))

            single_meas <- cvc_to_raw_cf(h5_get_dataset(h5f, key))

            rval[[loop_type]][[igamma]][[imom]][ifile, 1:Time, istoch, 1, 1] <-
              single_meas$data 
          } # istoch
        } # imom
      } # igamma
    } # loop_type
    H5Fclose(h5f)
    if(verbose) tictoc::toc()
  } # ifile
  
  # finally make complete `raw_cf` objects
  for( itype in 1:length(rval) ){
    for( igamma in 1:length(rval[[itype]]) ){
      for( imom in 1:length(rval[[itype]][[igamma]]) ){
        rval[[itype]][[igamma]][[imom]] <- 
              raw_cf_data(raw_cf_meta(Time = Time,
                                      nrObs = 1,
                                      nrStypes = 1,
                                      dim = c(nstoch, 1, 1),
                                      nts = Time),
                          data = rval[[loop_type]][[igamma]][[imom]])
      } # imom
    } # igamma
  } # itype
  return(rval)
}
