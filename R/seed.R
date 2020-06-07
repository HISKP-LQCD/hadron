#' Set seed and store a seed which can be used to
#' reset the random number generator 
#'
#' @param new_seed integer. The new seed that is to be set. In case this is
#' parameter is missing, no changes are made and the function just returns
#' `NULL`. This is useful because a function can just pass on its own `seed`
#' argument and therefore control whether the seed shall be fixed or left
#' as-is.
#'
#' @return
#' The generated seed is returned if it exists. Otherwise `NULL`. In case that
#' `new_seed` was missing, `NULL` is returned.
swap_seed <- function (new_seed) {
  seed_range <- 2^31-1
  if (missing(new_seed)) {
    return (NULL)
  }

  temp <- sample.int(size=1, n=seed_range)
  set.seed(new_seed)

  return (temp)
}

#' Restore random number generator state
#'
#' @return
#' No return value, but the random seed is reset to
#' `old_seed`.
#' 
#' @param old_seed integer. Previous seed that should be restored globally.
restore_seed <- function (old_seed) {
  if (!is.null(old_seed))
    set.seed(old_seed)
}
