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
#' @param accumulated Boolean, depending on whether the loop data was accumulated
#'                    over the stochastic source d.o.f. or not, the keys are
#'                    different. Default: FALSE
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
cyprus_make_key_scalar <- function(istoch, loop_type, cid = 4, accumulated = FALSE){
  if( any( !(loop_type %in% c("Scalar","dOp","Naive")) ) ){
    stop("The only scalar loop types are 'Scalar', 'Naive' and 'dOp'")
  }

  nchar <- "n"
  if( accumulated ){ nchar <- "N" }

  return(sprintf("/conf_%04d/%sstoch_%04d/%s/loop",
                 cid, nchar, istoch, loop_type)
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
#' @param accumulated Boolean, depending on whether the loop data was accumulated
#'                    over the stochastic source d.o.f. or not, the keys are
#'                    different. Default: FALSE
#' @return
#' A character vector with the HDF5 key.
#' 
#' @export
cyprus_make_key_vector <- function(istoch, loop_type, dir, cid = 4, accumulated = FALSE){
  vector_loop_types <- c("LpsDw", "Loops", "LpsDwCv", "LoopsCv")
  if( any( !(loop_type %in% vector_loop_types ) ) ) {
    stop(sprintf("The only derivative loop types are %s",
                 do.call(paste, list(vector_loop_types))
                 )
        )
  }
  stopifnot( dir %in% c(0:3) )

  nchar <- "n"
  if( accumulated ){ nchar <- "N" }

  return(sprintf("/conf_%04d/%sstoch_%04d/%s/dir_%02d/loop",
                 cid, nchar, istoch, loop_type, dir))
}


#' @title read HDF5 loop files in the Cyprus CalcLoops format
#' @description The CalcLoops code produces HDF5 files which contain
#'              a matrix of momenta and the data for the loops (without
#'              spin projection) organised by stochastic sample. Currently, the
#'              reading code assumes that there is a single configuration stored per
#'              file.
#' @param selections Named list with names from the list 'Naive', 'Scalar', 'dOp', 'Loops'
#'                   'LpsDw', 'LpsDwCv', 'LoopsCv' specifying the requesetd loop types. 
#'                   The elements of this list are in turn expected
#'                   be data frames of the form
#'                     \tabular{rrr}{
#'                       \strong{px} \tab \strong{py} \tab \strong{pz} \cr
#'                       0           \tab 0           \tab 1           \cr
#'                       -2          \tab 1           \tab -3          \cr
#'                       ...         \tab ...         \tab ...
#'                     }
#'                   specifying the momentum combinations to be extracted for each
#'                   loop type.
#' @param files Vector of strings, list of HDF5 files to be processed.
#' @param Time Integer, time extent of the lattice.
#' @param nstoch Integer, number of stochastic samples to be expected in file.
#' @param accumulated Boolean or vector of boolean, specifies whether the loops, 
#'                    as organised by stochastic sample,
#'                    are accumulated, such that, say, element \code{n} corresponds to the 
#'                    sum over the first \code{n} stochastic samples. If specified as \code{TRUE},
#'                    the data is post-processed to recover the measurements for the particular
#'                    samples. In case this is specified as a vector, it must be of the same
#'                    length as \code{files}. Default: TRUE.
#' @param legacy_traj Boolean. The root group for the loop data is 'conf_xxxx', where 'xxxx'
#'                    corresponds to what is passed via the 'traj' flag to CalcLoops. When
#'                    left empty, this defaults to '0004'. If this was left emtpy when
#'                    the loop files were generated, set this to \code{TRUE} and the paths
#'                    will be constructed with 'conf_0004' as their root group. 
#'                    When specified as a vector, it must be of length \code{length(files)}
#'                    giving the integer configuration indices, such as \code{c(0, 2, 4, 6)}
#'                    Default: TRUE.
#' @param verbose Boolean, output I/O time per file. Requires 'tictoc' package. Default FALSE.
#' @param check_group_names Boolean, employ \code{rhdf5::h5ls} to check if all the group names
#'                          that we want to read are actually in the file. This can be slow
#'                          for large files. Default: FALSE.
#' @param spin_project Boolean, whether the loops should be spin projected after being read.
#'                     Must be provided to together with \code{project_gamma}! Default: FALSE
#' @param project_gamma Named list of the same length as \code{selections} containing, for each
#'                      selected loop type a 4x4 complex-valued projection matrix. For vector
#'                      loop types, one matrix must be provided per direction (so \code{project_gamma$loop_type}
#'                      is a numbered list with indices \code{c(1,2,3,4)}. Default: NULL
#' @param use_parallel Boolean, whether to parallelise over the files using the `mclapply` from the `parallel`
#'                     package.
#' @return Named nested list of the same length as \code{selections} containg the loop data
#'         in the \link{raw_cf} format. Each named element corresponds to one loop
#'         type.
#'         For scalar-valued loop types, each element of the underlying numbered list
#'         corresponds to one momentum combination as specified via \code{selections} 
#'         for this loop type in the same order.
#'         For the vector-valued loop types, the first level of the underlying numbered
#'         list has four elements corresponding to the four derivative directions 
#'         in the order t,x,y,z.
#'         At the next level, the underlying numbered list corresponds to the momentum
#'         combinations for this loop type and derivative direction, just as for the
#'         scalar type.
#'         
#' @export
cyprus_read_loops <- function(selections, files, Time, nstoch,
                              accumulated = TRUE,
                              legacy_traj = TRUE, 
                              verbose = FALSE,
                              check_group_names = FALSE,
                              spin_project = FALSE,
                              project_gamma = NULL,
                              use_parallel = TRUE){
  rhdf5_avail <- requireNamespace("rhdf5")
  dplyr_avail <- requireNamespace("dplyr")
  if( !rhdf5_avail | !dplyr_avail ){
    stop("The 'dplyr' and 'rhdf5' packages are required to use this function!\n")
  }

  if( verbose ){
    tictoc_avail <- requireNamespace("tictoc")
    if( !tictoc_avail ){
      stop("Time reporting (verbose = TRUE) requires the 'tictoc' package!")
    }
  }
  if( use_parallel ){
    parallel_avail <- requireNamespace("parallel")
    if( !parallel_avail ){
      stop("Running in parallel requires the 'parallel' package!")
    }
    lapply_fn <- parallel::mclapply
  } else {
    lapply_fn <- lapply
  }

  files_not_exist <- unlist(lapply(X = files,
                                   FUN = 
                                     function(x){
                                       !file.exists(x)
                                     }
                                   )
                           )
  if( any( files_not_exist ) ){
    warning("The files: ", files[files_not_exist], " do not exist! \n")
    stop("There were missing files!")
  }
 
  # we require 'accumulated' and 'legacy_traj' to be vectors below
  # so if only one value was passed, replicate these here
  if( length(accumulated) == 1 ){
    accumulated <- rep(x = accumulated,
                       times = length(files))
  }
  if( length(legacy_traj) == 1 ){
    legacy_traj <- rep(x = legacy_traj,
                       times = length(files))
  }

  scalar_loop_types <- c("Scalar", "dOp", "Naive")
  vector_loop_types <- c("Loops", "LpsDw", "LpsDwCv", "LoopsCv")
  supported_loop_types <- c(scalar_loop_types, vector_loop_types)

  rval <- list()
  selected_loop_types <- names(selections)

  if( any( !(selected_loop_types %in% supported_loop_types ) ) ){
    invalid_loop_types <- 
      selected_loop_types[ !(selected_loop_types %in% supported_loop_types)]
    args <- lapply(X=1:length(invalid_loop_types), FUN=function(x){ invalid_loop_types[x] })
    args$sep <- ", "
    stop(sprintf("Loop types '%s' are not supported!",
                do.call(paste, args)))
  }

  # if we want to spin project, projection gamma matrices should be provided
  # for all selected loop types
  if( spin_project ){
    stopifnot( !is.null(project_gamma) )
    stopifnot( all( names(selections) == names(project_gamma) ) )
    selected_vec_loop_types <- which( names(selections) %in% vector_loop_types )

    incorrectly_sized_projectors <- NULL
    for( loop_type in selected_vec_loop_types ){
      if( length(project_gamma[[loop_type]] != 4) ){
        incorrectly_sized_projectors <- c(incorrectly_sized_projectors, loop_type)
      }
    }
    if( length(incorrectly_sized_projectors) > 0 ){
      stop(sprintf("'project_gamma' for loop types '%s' were not of the right length!",
                   do.call(paste, incorrectly_sized_projectors)
                   )
          )
    }
    
    # when spin-projecting, we pre-allocate an internal dimension of 1x1
    internal_dim <- 1
  } else {
    # otherwise, we pre-allocate 4x4 matrices internally
    internal_dim <- 4
  }

  loops_per_file <- lapply_fn(
    X = 1:length(files),
    FUN = function(ifile){
      f <- files[ifile]
      if(verbose){
        message(sprintf("Opening %s\n",f))
        tictoc::tic("Loop extraction")
      }

      # our output for this file
      this_file_loops <- list()
      
      # The file names are of the form 
      # path/MG_loop_FLAVORquark_conf_conf.XXXX_runtype_probD8_part1_stoch_NeV0_NsYYYY_step0001_QsqZZ.h5
      # and we want to recover XXXX
      tokens <- unlist(strsplit(basename(f), split = ".", fixed = TRUE))
      cid_in_filename <- as.integer(strsplit(tokens, split = "_", fixed = TRUE)[[2]][1])

      # In the generation of loop files via CalcLoops in the "flavour singlet" project running
      # on PizDaint from Oct. 2018 - Oct. 2019, there was an issue with the naming of
      # the "trajectory" group name in the HDF5 files. If files were affected,
      # the stored "trajectory" id was always "4" instead of the actual configuration number.
      if( legacy_traj[ifile] ){
        cid_to_read <- 4
      } else {
        cid_to_read <- cid_in_filename
      }

      h5f <- rhdf5::H5Fopen(f, flags = "H5F_ACC_RDONLY")
   
      if( check_group_names ){ 
        group_names <- rhdf5::h5ls(h5f)$name

        # are all requested loop types really available?
        avail_loop_types <- unlist( lapply( selected_loop_types, function(x){ x %in% group_names } ) )
        if( any( !avail_loop_types ) ){
          msg <- sprintf("Some selected loop types could not be found in %s:\n %s",
                         f,
                         do.call(paste, as.list( selected_loop_types[!avail_loop_types] ) )
                         )
          stop(msg)
        }
        # how many stochastic samples are available and does it match our expectation?
        if( accumulated[ifile] ){
          stoch_text <- "Nstoch"
        } else {
          stoch_text <- "nstoch"
        }
        stoch_avail <- sort(as.numeric(
                         unlist(
                            lapply(X = strsplit(unique(group_names[ grepl(stoch_text, group_names) ]), "_"),
                                   FUN = function(x){ x[2] })
                          )
                       ))
        if( length(stoch_avail) != nstoch ){
          stop(sprintf("Number of stochastic samples in file %s :\n%d, expected %d!",
                       f, length(stoch_avail), nstoch))
        }

        # check if there are multiple instances of 'conf_xxxx' group names
        # in the file
        if( length( unique( group_names[ grepl("conf", group_names) ] ) ) > 1 ){
          warning(sprintf(paste("The file\n%s\ncontains more than one 'conf_xxxx'",
                                "group names.\nThis is currently not really supported,",
                                "but we will continue and attempt to read 'conf_%04d' (and no others!)"),
                          f, cid_to_read),
                  immediate. = TRUE)
        }
      } # if(check_group_names)

      if( !rhdf5::H5Lexists(h5f, "Momenta_list_xyz") ){
        stop(sprintf("'Momenta_list_xyz' could not be found in %s!", f))
      }
      # we transpose this to get the momenta as the rows of a matrix
      momenta_avail <- as.data.frame(t(h5_get_dataset(h5f, "Momenta_list_xyz")))
      colnames(momenta_avail) <- c("px","py","pz")
      # index the momentum combinations
      momenta_avail <- cbind(momenta_avail, idx = 1:nrow(momenta_avail))

      for( loop_type in selected_loop_types ){
        # check if all the momenta that we want are in the file
        # we do this per loop_type as we could have different selections
        # for different loop types
        missing_momenta <- dplyr::anti_join(x = selections[[ loop_type ]],
                                            y = momenta_avail,
                                            by = c("px","py","pz"))
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
        # select which elements we need to read from the available momenta
        # this fixes the set of momenta_avail$idx which are going to be processed
        selected_momenta <- dplyr::inner_join(x = selections[[ loop_type ]],
                                              y = momenta_avail,
                                              by = c("px","py","pz"))

        for( istoch in 1:nstoch ){
          if( loop_type %in% scalar_loop_types ){
            key <- cyprus_make_key_scalar(istoch = istoch,
                                          loop_type = loop_type,
                                          cid = cid_to_read,
                                          accumulated = accumulated[ifile])

            # read the data, which comes in the ordering
            #   complex, gamma, mom_idx, time
            # we permute it to
            #   time, gamma, complex, mom_idx
            # this is quite expensive, but it makes filling the target
            # array much easier below
            # Note that 'gamma' is of length 16
            data <- h5_get_dataset(h5f, key)
            # first we select the momenta that we actually want
            data <- data[,,selected_momenta$idx,,drop=FALSE]
            # and then perform the reshaping on this (possibly) reduced array
            data <- aperm(data, perm = c(4,2,1,3))
            dims <- dim(data)
            nts <- dims[1]

            if( nts != Time ){
              stop(sprintf("For %s, the number of time slices %d does not match the expected time extent %d!",
                           key, nts, Time))
            }

            if( !(loop_type %in% names(this_file_loops) ) ){
              this_file_loops[[ loop_type ]] <- list()
            }
            for( mom_idx in 1:nrow(selected_momenta) ){
              # if the list element corresponding to 'mom_idx' is empty,
              # we pre-allocate some space for an array there
              # note that the first dimension is of size 1, we will concatenate over this dimension
              # after all files have been read
              if( length(this_file_loops[[ loop_type ]]) < mom_idx ){
                this_file_loops[[loop_type]][[mom_idx]] <- 
                  array(as.complex(NA), dim=c(1,
                                              Time,
                                              nstoch,
                                              internal_dim,
                                              internal_dim)
                                             )
              }
              # in order to make spin projection easier, we create an intermediate raw_cf object
              # with a single measurement
              temp_raw_cf <- raw_cf_data(raw_cf_meta(Time = Time,
                                                     nrObs = 1,
                                                     nrStypes = 1,
                                                     dim = c(1, 4, 4),
                                                     nts = Time),
                                         data = array(complex(real = data[1:Time, 1:16, 1, mom_idx ],
                                                              imaginary = data[1:Time, 1:16, 2, mom_idx ]),
                                                      dim = c(1, Time, 1, 4, 4)))
              if(spin_project){
                temp_raw_cf <- loop_spin_project(loop = temp_raw_cf,
                                                 gamma = project_gamma[[loop_type]])
              }
              # and now we extract the data into our pre-allocated tensor again
              # further below, the entire tensor will be turned back into raw_cf
              # This is wasteful, but allows us to reuse loop_spin_project here
              this_file_loops[[loop_type]][[mom_idx]][1, 1:Time, istoch, 1:internal_dim, 1:internal_dim] <-
                temp_raw_cf$data
            }
          } else if( loop_type %in% vector_loop_types ){
            if( !(loop_type %in% names(this_file_loops) ) ){
              this_file_loops[[loop_type]] <- list()
            }
            if( length(this_file_loops[[loop_type]]) == 0 ){
              for( dir_idx in 1:4 ){
                this_file_loops[[loop_type]][[dir_idx]] <- list()
              }
            }
            # we are going to interchange t <-> x as direction 0 below
            deriv_dirs_out <- c(1,2,3,0)
            for( deriv_dir_in in c(0,1,2,3) ){
              key <- cyprus_make_key_vector(istoch = istoch,
                                            loop_type = loop_type,
                                            dir = deriv_dir_in,
                                            cid = cid_to_read,
                                            accumulated = accumulated[ifile])

              # read the data, which comes in the ordering
              #   complex, gamma, mom_idx, time
              # we permute it to
              #   time, gamma, complex, mom_idx
              # this is quite expensive, but it makes filling the target
              # array much easier below
              # Note that 'gamma' is of length 16
              data <- h5_get_dataset(h5f, key)
              # first we select the momenta that we actually want
              data <- data[,,selected_momenta$idx,,drop=FALSE]
              # and then perform the reshaping on this (possibly) reduced array
              data <- aperm(data, perm = c(4,2,1,3))
              dims <- dim(data)
              nts <- dims[1]
            
              if( nts != Time ){
                stop(sprintf("For %s, the number of time slices %d does not match the expected time extent %d!",
                             key, nts, Time))
              }

              # we interchange t <-> x as direction 0
              # and assign an idx in 1:4 -> (t,1),(x,2),(y,3),(z,4) 
              deriv_dir_out_idx <- deriv_dirs_out[ deriv_dir_in+1 ]+1
              for( mom_idx in 1:nrow(selected_momenta) ){
                if( length(this_file_loops[[loop_type]][[deriv_dir_out_idx]]) < mom_idx ){
                  this_file_loops[[loop_type]][[deriv_dir_out_idx]][[mom_idx]] <- 
                    array(as.complex(NA), dim=c(1,
                                                Time,
                                                nstoch,
                                                internal_dim,
                                                internal_dim)
                          )
                }
                # in order to make spin projection easier, we create an intermediate raw_cf object
                # with a single measurement
                temp_raw_cf <- raw_cf_data(raw_cf_meta(Time = Time,
                                                       nrObs = 1,
                                                       nrStypes = 1,
                                                       dim = c(1, 4, 4),
                                                       nts = Time),
                                           data = array(complex(real = data[1:Time, 1:16, 1, mom_idx ],
                                                                imaginary = data[1:Time, 1:16, 2, mom_idx ]),
                                                        dim = c(1,Time,1,4,4)))
                if(spin_project){
                  temp_raw_cf <- loop_spin_project(loop = temp_raw_cf,
                                                   gamma = project_gamma[[loop_type]][[deriv_dir_out_idx]])
                }
                # And now we extract the data into our pre-allocated tensor again
                # further below, the entire tensor will be turned back into raw_cf
                # This is wasteful, but allows us to reuse loop_spin_project here
                this_file_loops[[loop_type]][[deriv_dir_out_idx]][[mom_idx]][1, 1:Time, istoch, 1:internal_dim, 1:internal_dim] <-
                  temp_raw_cf$data
              } # for(mom_idx)
            } # for(deriv_dir_in)
          } # if(loop_type)
        } # for(istoch)
      } # for(loop_type)
      rhdf5::H5Fclose(h5f)

      # if the current file was produced in the "accumulated" convention -> de-accumulate it
      if( accumulated[ifile] ){
        for( loop_type in selected_loop_types ){
          # recover measurements from individual stochastic samples
            for( mom_idx in 1:nrow(selected_momenta) ){
              if( loop_type %in% scalar_loop_types ){
                temp <- this_file_loops[[loop_type]][[mom_idx]]
                for( istoch in 2:nstoch ){
                  this_file_loops[[loop_type]][[mom_idx]][1,,istoch,,] <- temp[1,,istoch,,] - temp[1,,(istoch-1),,]
                }
              } else if( loop_type %in% vector_loop_types ){
                # for the vector-valued loop types we need to duplicate some code
                for( dir_idx in 1:4 ){
                  temp <- this_file_loops[[loop_type]][[dir_idx]][[mom_idx]]
                  for( istoch in 2:nstoch ){
                    this_file_loops[[loop_type]][[dir_idx]][[mom_idx]][1,,istoch,,] <- temp[1,,istoch,,] - temp[1,,(istoch-1),,]
                  }
                } # for(dir_idx)
              } # if(loop_type)
            } # for(mom_idx) 
          } # for(selected_loop_types)
      } # if(accumulated)

      if(verbose) tictoc::toc()
      return(this_file_loops)
    } # FUN(ifile)
  ) # lapply

  # this is our final return value for which we will now collect the loop data
  # from all files, converting it to raw_cf in the process
  rval <- list()
  for( loop_type in selected_loop_types ){
    rval[[loop_type]] <- list()
    if( loop_type %in% scalar_loop_types ){
      for( mom_idx in 1:nrow(selections[[loop_type]]) ){
        temp <- array(as.complex(NA), dim=c(length(files),
                                            Time,
                                            nstoch,
                                            internal_dim,
                                            internal_dim)
                      )
        for( ifile in 1:length(files) ){
          temp[ifile, , , , ] <- 
            loops_per_file[[ifile]][[loop_type]][[mom_idx]][1, , , , ]
        }
        rval[[ loop_type ]][[ mom_idx ]] <-
          raw_cf_data(raw_cf_meta(Time = Time,
                                  nrObs = 1,
                                  nrStypes = 1,
                                  dim = c(nstoch, internal_dim, internal_dim),
                                  nts = Time),
                      data = temp)
        
      }
      # and the same for the vector loop types
    } else if( loop_type %in% vector_loop_types ) { 
      for( dir_idx in 1:4 ){
        rval[[loop_type]][[dir_idx]] <- list()
        for( mom_idx in 1:nrow(selections[[loop_type]][[dir_idx]]) ){
          temp <- array(as.complex(NA), dim=c(length(files),
                                              Time,
                                              nstoch,
                                              internal_dim,
                                              internal_dim)
                        )
          for( ifile in 1:length(files) ){
            temp[ifile, , , , ] <- 
              loops_per_file[[ifile]][[loop_type]][[dir_idx]][[mom_idx]][1, , , , ]
          }
          rval[[loop_type]][[dir_idx]][[mom_idx]] <- 
            raw_cf_data(raw_cf_meta(Time = Time,
                                    nrObs = 1,
                                    nrStypes = 1,
                                    dim = c(nstoch, internal_dim, internal_dim),
                                    nts = Time),
                        data = temp)
        } # for(mom_idx)
      } # for(dir_idx)
    } # if( loop_type )
  } # for( loop_type )
  return(rval)
}

