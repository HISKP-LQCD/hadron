#' @title get dataset from HDF5 file
#' @param h5f HDF5 file opened with \code{rhdf5::H5Fopen}
#' @param key String, full path to dataset.
#' @param check_exists Boolean, check if key actually exists (keep in mind overhead).
#'
#' @return
#' Returns the requested dataset, if successfully read from file.
#' 
#' @export
h5_get_dataset <- function(h5f, key, check_exists = TRUE)
{
  rhdf5_avail <- requireNamespace("rhdf5")
  stopifnot( rhdf5_avail )
  exists <- ifelse(check_exists, rhdf5::H5Lexists(h5f, key), TRUE)
  if( exists ){
    h5d <- rhdf5::H5Dopen(h5f, key)
    rval <- rhdf5::H5Dread(h5d)
    rhdf5::H5Dclose(h5d)
  } else {
    stop(sprintf("Dataset %s could not be found!", key))
  }
  return(rval)
}

#' @title check if group names exist in HDF5 file
#' @description The group names in an HDF5 file are stored as full paths
#'              as well as a flat vector. It is thus possible to check
#'              if a particular set of group names exist in the file
#'              by parsing the \code{name} member of the output
#'              of \code{rhdf5::h5ls}. This function does just that.
#' @param h5f HDF5 file handle openend with \code{rhdf5::H5Fopen}
#' @param nms_to_find Vector of strings, group names (not full paths) which
#'                    are to be located in the file.
#' @return Vector of booleans of the same length as \code{nms_to_find}
#'         indicating whether the name at the same index position
#'         was located in the file.
h5_names_exist <- function(h5f, nms_to_find){
  rhdf5_avail <- requireNamespace("rhdf5")
  stopifnot( rhdf5_avail )
  nms <- rhdf5::h5ls(h5f)$name
  unlist( lapply( nms_to_find, function(x){ x %in% nms } ) )
}

