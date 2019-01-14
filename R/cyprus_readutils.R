#' @title HDF5 key for Cyprus CalcLoops scalar-type loops
#' @description Generates an HDF5 key (full path) for the scalar
#'              type loops from the Cyprus CalcLoops application.
#' @param istoch Integer, index of the stochastic sample that the key should
#'               be generated for.
#' @param loop_type String, name of loop type. Allowed values:
#'                  'Scalar', 'dOp'
#' @param cid Integer, configuration number, internally produced by the CalcLoops
#'            tool via the "trajectory" input flag. The default is '4' as this is
#'            often not used in practice.
#' @export
cyprus_make_key_scalar <- function(istoch, loop_type, cid = 4){
  if( any( !(loop_type %in% c("Scalar","dOp")) ) ){
    stop("The only scalar loop types are 'Scalar' and 'dOp'")
  }
  return(sprintf("/conf_%04d/Nstoch_%04d/%s/loop",
                 cid, istoch, loop_type)
        )
}

#' @title HDF5 key for Cyprus CalcLoops derivative-type loops
#' @description Generates an HDF5 key (full path) for the derivative
#'              type loops from the Cyprus CalcLoops application.
#' @param istoch Integer, index of the stochastic sample that the key should
#'               be generated for.
#' @param loop_type String, name of loop type. Allowed values:
#'                  'Loops', 'LpsDw', 'LpsDwCv', 'LoopsCv'
#' @param dir Integer, lattice direction of the derivative. Allowed values:
#'            \code{0 == x}, \code{1 == y}, \code{2 == z}, \code{3 == t}.
#' @param cid Integer, configuration number, internally produced by the CalcLoops
#'            tool via the "trajectory" input flag. The default is '4' as this is
#'            often not used in practice.
#' @export
cyprus_make_key_deriv <- function(istoch, loop_type, dir, cid = 4){
  deriv_loop_types <- c("LpsDw", "Loops", "LpsDwCv", "LoopsCv")
  if( any( !(loop_type %in% deriv_loop_types ) ) ) {
    stop(sprintf("The only derivative loop types are %s",
                 do.call(paste, list(deriv_loop_types))
                 )
        )
  }
  stopifnot( dir %in% c(0:3) )
  return(sprintf("/conf_%04d/Nstoch_%04d/%s/dir_%02d/loop",
                 cid, istoch, loop_type, dir))
}


#' @title read HDF5 loop files in the Cyprus CalcLoops format
#' @description The CalcLoops code produces HDF5 files which contain
#'              a matrix of momenta and the data for the loops (without
#'              spin projection) organised by stochastic sample. Currently, the
#'              reading code assumes that there is a single configuration stored per
#'              file and the "trajectory" parameter in CalcLoops is assumed
#'              to take its default value of '4'.
#' @param selections Named list with names from the list 'Scalar', 'dOp', 'Loops'
#'                   'LpsDw', 'LpsDwCv', 'LoopsCv' specifying the requesetd loop types. 
#'                   The elements of this list are in turn expected
#'                   be data frames of the form
#'
#'                     \tabular{rrr}{
#'                       \strong{qx} \tab \strong{qy} \tab \strong{qz} \cr
#'                       0           \tab 0           \tab 1           \cr
#'                       "..."       \tab "..."       \tab "..."
#'                     }
#'
#'                   specifying the momentum combinations to be extracted for each
#'                   loop type.
#' @param files Vector of strings, list of HDF5 files to be processed.
#' @param accumulated Boolean, specifies whether the loops, as organised by stochastic sample,
#'                    are accumulated, such that, say, element \code{n} corresponds to the 
#'                    sum over the first \code{n} stochastic samples. If specified as \code{TRUE},
#'                    the data is post-processed to recover the measurements for the particular
#'                    samples.
#' @return Named list of the same length as \code{selections} containg the loop data
#'         in the \link{raw_cf} format.
#' @export
cyprus_read_loops <- function(selections, files, accumulated = TRUE){
  rhdf5_avail <- requireNamespace("rhdf5")
  dplyr_avail <- requireNamespace("dplyr")
  if( !rhdf5_avail | !dplyr_avail ){
    stop("The 'dplyr' and 'rhdf5' packages are required to use this function!\n")
  }

  rval <- list()
  selected_loop_types <- names(selections)

  for( ifile in 1:length(files) ){
    f <- files[ifile]
    h5f <- rhdf5::H5Fopen(f)
    avail_loop_types <- h5_names_exist(h5f, selected_loop_types)
    if( any( !avail_loop_types ) ){
      msg <- sprintf("Some selected loop types could not be found in %s:\n %s",
                     f,
                     do.call(paste, as.list( selected_loop_types[!avail_loop_types] ) )
                     )
      stop(msg)
    }
    if( !H5Lexists(h5f, "Momenta_list_xyz") ){
      stop(sprintf("'Momenta_list_xyz' could not be found in %s!", f))
    }
    # we transpose this to get the momenta as the rows of a matrix
    momenta_avail <- as.data.frame(t(h5_get_dataset(h5f, "Momenta_list_xyz")))
    colnames(momenta_avail) <- c("qx","qy","qz")
    # index the momentum combinations
    momenta_avail <- cbind(momenta_avail, idx = 1:nrow(momenta_avail))
    
    group_names <- h5ls(h5f)$name
    # how many stochastic sources are available?
    stoch_avail <- sort(as.numeric(unlist(lapply(X = strsplit(unique(group_names[ grepl("Nstoch", group_names) ]),
                                       "_"),
                                 FUN = function(x){ x[2] })
                          )))

    for( loop_type in selected_loop_types ){
      # check if all the momenta that we want are in the file
      # we do this per loop_type as we could have different selections
      # for different loop types
      missing_momenta <- dplyr::anti_join(x = selections[[ loop_type ]],
                                          y = momenta_avail,
                                          by = c("qx","qy","qz"))
      if( nrow(missing_momenta) > 0 ){
        msg <- sprintf("\nMomenta\n%s\ncould not be found in %s!",
                      do.call(paste,
                              list(...=apply(X = missing_momenta,
                                             MARGIN = 1,
                                             FUN = function(x){
                                                 sprintf("(%d,%d,%d)", x[1], x[2], x[3])
                                               }
                                           ),
                                    sep = '\n'
                                  )
                              ),
                      f
                      )
        stop(msg)
      }
      # select which elements we need to read
      selected_momenta <- dplyr::inner_join(x = selections[[ loop_type ]],
                                            y = momenta_avail,
                                            by = c("qx","qy","qz"))

      for( istoch in stoch_avail ){
        key <- cyprus_make_key_scalar(istoch = istoch,
                                      loop_type = loop_type)

        # read the data, which comes in the ordering
        #   complex, gamma, mom_idx, time
        # we permute it to
        #   time, gamma, complex, mom_idx
        # this is quite expensive, but it makes filling the target
        # array much easier below
        data <- aperm(h5_get_dataset(h5f, key),
                      perm = c(4,2,1,3))
        dims <- dim(data)
        nts <- dims[1]

        if( !(loop_type %in% names(rval) ) ){
          rval[[ loop_type ]] <- list()
        }
        for( mom_idx in 1:nrow(selected_momenta) ){
          if( length(rval[[ loop_type ]]) < mom_idx ){
            rval[[loop_type]][[mom_idx]] <- array(as.complex(NA), dim=c(length(files),
                                                                        nts,
                                                                        length(stoch_avail),
                                                                        4,
                                                                        4)
                                                  )
          }
          rval[[loop_type]][[mom_idx]][ifile, 1:nts, istoch, 1:4, 1:4] <-
            complex(real = data[1:nts, 1:16, 1, selected_momenta$idx[mom_idx] ],
                    imaginary = data[1:nts, 1:16, 2, selected_momenta$idx[mom_idx] ])
        }
      }
    }
  }
  return(rval)
}

