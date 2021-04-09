#' @title Generate HDF5 key for a momentum and spin-projected CVC loop 
#' @param loop_type String, loop type.
#' @param istoch Integer, stochastic sample id number.
#' @param gamma Integer, CVC convention gamma matrix identifier.
#' @param p Integer vector of length 3, (x,y,z) components of the momentum
#'          vector in lattice units.
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
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
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
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
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
correlators_key_meson_2pt <- function(fwd_flav, bwd_flav, src_ts, snk_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  stopifnot( is.character(fwd_flav) )
  stopifnot( is.character(bwd_flav) )
  stopifnot( is.integer(src_ts) )
  stopifnot( is.integer(snk_gamma) )
  stopifnot( is.integer(src_gamma) )

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
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
cf_key_meson_2pt <- function(fwd_flav, bwd_flav, snk_gamma, src_gamma, src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  stopifnot( is.character(fwd_flav) )
  stopifnot( is.character(bwd_flav) )
  stopifnot( is.integer(snk_gamma) )
  stopifnot( is.integer(src_gamma) )
  
  sprintf("/%s+-g-%s-g/gf%d/pfx%dpfy%dpfz%d/gi%d/pix%dpiy%dpiz%d",
          bwd_flav,
          fwd_flav,
          snk_gamma,
          snk_p[1], snk_p[2], snk_p[3],
          src_gamma,
          src_p[1],src_p[2],src_p[3])
}

#' @title Generate HDF5 key for CVC 'correlators' meson 3pt function with a local or derivative insertion
#' 
#' The key for a meson three-point function has the form:
#' 
#'  /sud+-g-u-g/t10/dt12/gf5/pfx0pfy0pfz0/gc0/Ddim0_dir0/Ddim1_dir1/D[...]/gi5/pix0piy0piz0
#' 
#' where, from left to right:
#' \itemize{
#'  \item 'u' is the flavour of the "backward" propagator
#'  \item 'd' is the flavour of the "sequential" propagator
#'  \item '+' indicates that 'sud' is daggered
#'  \item 'g' indicates a gamma insertion 
#'  \item 'u' is the flavour of the foward propagator
#'  \item 'g' indicates a Dirac structure at the source
#'  \item 'tXX' is the source time slice
#'  \item 'dtYY' is the source-sink separation
#'  \item 'gfN' gamma structure at the sink in CVC indexing
#'  \item 'pfxXpfyYpfzZ' is the sink momentum in CVC convention (sink and source phases are both e^{ipx})
#'  \item 'gcN' gamma structure at the current insertion point in CVC indexing
#'  \item 'DdimJ_dirK' covariant displacement applied in dimension 'J', direction 'K'
#'         where it should be noted that this is. in operator notation, i.e., the right-most
#'         displacement is the one applied first.
#'  \item [...]
#'  \item 'giN' gamma structure at the souce in CVC indexing
#'  \item 'pixXpiyYpizZ' at the source in CVC convention
#' }
#'
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param seq_flav String, "sequential" quark flavour identifier.
#' @param src_ts Integer, source time slice.
#' @param dt Integer, source-sink separation.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param cur_gamma Integer, CVC convention gamma matrix identified at the insertion.
#' @param cur_displ_dim Integer vector of dimensions (0,1,2,3 <-> t,x,y,z) in which
#'                      covariant displacements have been applied. This vector will be
#'                      parsed in reverse order, such that the first element here
#'                      is the first displacement applied to the spinor in the calculation
#'                      and the right-most element in the key.
#'                      Length must be matched to 'cur_displ_dir'.
#'                      Defaults to 'NA' for no displacements.
#' @param cur_displ_dir Integer vector of directions (forward, backward) <-> (0,1) in
#'                      which the covariant displacements have been applied. Parsing
#'                      as for 'cur_displ_dim'. Length must be matched to 'cur_displ_dim'.
#'                      Defaults to 'NA' for no displacements.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
correlators_key_meson_3pt <- function(fwd_flav, bwd_flav, seq_flav,
                                      src_ts, dt, 
                                      snk_gamma,
                                      cur_gamma, cur_displ_dim = NA, cur_displ_dir = NA,
                                      src_gamma, 
                                      src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  stopifnot( is.character(fwd_flav) )
  stopifnot( is.character(bwd_flav) )
  stopifnot( is.integer(src_ts) )
  stopifnot( is.integer(dt) )
  stopifnot( is.integer(snk_gamma) )
  stopifnot( is.integer(src_gamma) )
  
  if( !is.na(cur_displ_dim) | !is.na(cur_displ_dir) ){
    stopifnot( all(cur_displ_dim %in% c(0:3)) )
    stopifnot( all(cur_displ_dir %in% c(0,1)) )
    stopifnot( length(cur_displ_dim) == length(cur_displ_dir) )
  }

  key <- sprintf("/s%s%s+-g-%s-g/t%d/dt%d/gf%d/pfx%dpfy%dpfz%d/gc%d",
                 bwd_flav,
                 seq_flav,
                 fwd_flav,
                 src_ts,
                 dt,
                 snk_gamma,
                 snk_p[1],snk_p[2],snk_p[3],
                 cur_gamma)

  # parse the displacements in reverse order such as to obtain operator notation
  # in the key, with the right-most key value representing the first operator
  # applied 
  if( !is.na(cur_displ_dim) ){
    for( i_displ in c(length(cur_displ_dim):1) ){
      key <- sprintf("%s/Ddim%d_dir%d", 
                     key, cur_displ_dim[i_displ], cur_displ_dir[i_displ])
    }
  }
  key <- sprintf("%s/gi%d/pix%dpiy%dpiz%d",
                 key,
                 src_gamma,
                 src_p[1],src_p[2],src_p[3])

  return(key)
}

#' @title Generate HDF5 key for CVC 'correlators' meson 3pt function with a local or derivative insertion
#' 
#' The key for a meson three-point function has the form:
#' 
#'  /sud+-g-u-g/t10/dt12/gf5/pfx0pfy0pfz0/gc0/Ddim0_dir0/Ddim1_dir1/D[...]/gi5/pix0piy0piz0
#' 
#' where, from left to right:
#' \itemize{
#'  \item 'u' is the flavour of the "backward" propagator
#'  \item 'd' is the flavour of the "sequential" propagator
#'  \item '+' indicates that 'sud' is daggered
#'  \item 'g' indicates a gamma insertion 
#'  \item 'u' is the flavour of the foward propagator
#'  \item 'g' indicates a Dirac structure at the source
#'  \item 'tXX' is the source time slice
#'  \item 'dtYY' is the source-sink separation
#'  \item 'gfN' gamma structure at the sink in CVC indexing
#'  \item 'pfxXpfyYpfzZ' is the sink momentum in CVC convention (sink and source phases are both e^{ipx})
#'  \item 'gcN' gamma structure at the current insertion point in CVC indexing
#'  \item 'DdimJ_dirK' covariant displacement applied in dimension 'J', direction 'K'
#'         where it should be noted that this is. in operator notation, i.e., the right-most
#'         displacement is the one applied first.
#'  \item [...]
#'  \item 'giN' gamma structure at the souce in CVC indexing
#'  \item 'pixXpiyYpizZ' at the source in CVC convention
#' }
#'
#' @param fwd_flav String, "forward" quark flavour identifier.
#' @param bwd_flav String, "backward" quark flavour identifier.
#' @param seq_flav String, "sequential" quark flavour identifier.
#' @param dt Integer, source-sink separation.
#' @param snk_gamma Integer, CVC convention gamma matrix identifier at the source.
#' @param cur_gamma Integer, CVC convention gamma matrix identified at the insertion.
#' @param cur_displ_dim Integer vector of dimensions (0,1,2,3 <-> t,x,y,z) in which
#'                      covariant displacements have been applied. This vector will be
#'                      parsed in reverse order, such that the first element here
#'                      is the first displacement applied to the spinor in the calculation
#'                      and the right-most element in the key.
#'                      Length must be matched to 'cur_displ_dir'.
#'                      Defaults to 'NA' for no displacements.
#' @param cur_displ_dir Integer vector of directions (forward, backward) <-> (0,1) in
#'                      which the covariant displacements have been applied. Parsing
#'                      as for 'cur_displ_dim'. Length must be matched to 'cur_displ_dim'.
#'                      Defaults to 'NA' for no displacements.
#' @param src_gamma Integer, CVC convention gamma matrix identified at the sink.
#' @param src_p Integer vector of length 3. (x,y,z) components of the source momentum
#'              vector in lattice units.
#' @param snk_p Integer vector of length 3. (x,y,z) components of the sink momentum
#'              vector in lattice units.
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
cf_key_meson_3pt <- function(fwd_flav, bwd_flav, seq_flav,
                             dt, 
                             snk_gamma,
                             cur_gamma, cur_displ_dim = NA, cur_displ_dir = NA,
                             src_gamma, 
                             src_p, snk_p)
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  stopifnot( is.character(fwd_flav) )
  stopifnot( is.character(bwd_flav) )
  stopifnot( is.integer(dt) )
  stopifnot( is.integer(snk_gamma) )
  stopifnot( is.integer(src_gamma) )
  if( !is.na(cur_displ_dim) | !is.na(cur_displ_dir) ){
    stopifnot( all(cur_displ_dim %in% c(0:3)) )
    stopifnot( all(cur_displ_dir %in% c(0,1)) )
    stopifnot( length(cur_displ_dim) == length(cur_displ_dir) )
  }

  key <- sprintf("/s%s%s+-g-%s-g/dt%d/gf%d/pfx%dpfy%dpfz%d/gc%d",
                 bwd_flav,
                 seq_flav,
                 fwd_flav,
                 dt,
                 snk_gamma,
                 snk_p[1],snk_p[2],snk_p[3],
                 cur_gamma)

  # parse the displacements in reverse order such as to obtain operator notation
  # in the key, with the right-most key value representing the first operator
  # applied 
  if( !is.na(cur_displ_dim) ){
    for( i_displ in c(length(cur_displ_dim):1) ){
      key <- sprintf("%s/Ddim%d_dir%d", 
                     key, cur_displ_dim[i_displ], cur_displ_dir[i_displ])
    }
  }
  key <- sprintf("%s/gi%d/pix%dpiy%dpiz%d",
                 key,
                 src_gamma,
                 src_p[1],src_p[2],src_p[3])

  return(key)
}

#' @title Convert correlation function read from CVC HDF5 or AFF format to 'raw_cf'
#' @description Given a numeric vector of alternating real and imaginary parts of a
#'              correlation function, creates an object of class 'raw_cf' with
#'              a single measurement, inferring \code{Time} from the passed numeric vector
#'              while the shape of the internal dimensions has to be specified explicitly
#'              if larger than `one by one` (\code{c(1,1)}).
#' @param cf_dat Numeric vector of alternating real and imaginary parts of a
#'               correlation function. Ordering of the input should be complex, internal dimensions,
#'               time (left to right, fastest to slowest).
#' @param dims Integer vector with the sizes of the internal dimensions. For example,
#'             \code{c(4,4)} for spin correlators.
#' @return `raw_cf` object with a \code{data} member which contains the data (as complex numbers)
#'         in the shape \code{c(1,nts,dims)}, where `nts` is the number of time slices
#'         inferred from the length of \code{cfdat} and the product of the internal dimensions \code{dims}.
#'                 
#' @export
cvc_to_raw_cf <- function(cf_dat, dims = c(1,1))
{
  number_of_internal_dims <- prod(dims)

  # idcs for real and imaginary parts
  ridcs <- seq(1,length(cf_dat),2)
  iidcs <- seq(2,length(cf_dat),2)

  nts <- length(ridcs)/number_of_internal_dims

  # turn it into a complex 4D array of with ordering of
  # internal dimensions and time in the 'wrong' order
  cf_dat <- array(complex(real = cf_dat[ridcs],
                          imaginary = cf_dat[iidcs]),
                          dim = c(1,dims,nts))
  cf_dims <- dim(cf_dat)
  # reshape the array, ordering 'measurement' (a single one), 'time',
  # internal dimensions
  cf_dat <- aperm(cf_dat, c(1,length(cf_dims),2:(length(cf_dims)-1)))

  cf <- raw_cf_data(cf = raw_cf_meta(Time=nts, dim=dims),
                    data = cf_dat)
  return(cf)
}



#### FIXME FIXME FIXME 
#### THIS NOT NOT UP TO DATE
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
#' @param verbose Boolean, output I/O time per file. Requires 'tictoc' package. Default FALSE.
#' @param check_group_names Boolean, check if the group names that we're about to read actually
#'                          exist in the file. This is quite slow because it uses \code{rhdf5::h5ls}.
#'                          Default FALSE.
#' @return Named nested list of the same length as \code{selections} containg the loop data
#'         in the \link{raw_cf} format. Each named element corresponds to one loop
#'         type and each element of the underlying numbered list corresponds to one momentum
#'         combination as specified via \code{selections} for this loop type in the same order.
#'         
#' @export
cvc_read_loops <- function(selections, files, Time, nstoch, verbose = FALSE, check_group_names = FALSE){
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
   
    if( check_group_names ){ 
      group_names <- rhdf5::h5ls(h5f)$name
      
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
    rhdf5::H5Fclose(h5f)
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
