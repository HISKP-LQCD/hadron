#' @title get dataset from HDF5 file
#' @param h5f HDF5 file opened with \code{H5Fopen}
#' @param key String, full path to dataset.
#' @export
h5_get_dataset <- function(h5f,key)
{
  requireNamespace("rhdf5")
  exists <- rhdf5::H5Lexists(h5f, key)
  if( exists ){
    h5d <- rhdf5::H5Dopen(h5f, key)
    rval <- rhdf5::H5Dread(h5d)
    rhdf5::H5Dclose(h5d)
  } else {
    stop(sprintf("Dataset %s could not be found!", key))
  }
  return(rval)
}
