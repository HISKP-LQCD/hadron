#' Save random number generator state
#'
#' @param new_seed integer. The new seed that is to be set. In case this is
#' parameter is missing, no changes are made and the function just returns
#' `NULL`. This is useful because a function can just pass on its own `seed`
#' argument and therefore control whether the seed shall be fixed or left
#' as-is.
#'
#' @return
#' The old random seed is returned if it exists. Otherwise `NULL`. In case that
#' `new_seed` was missing, `NULL` is returned.
swap_seed <- function (new_seed) {
  if (missing(new_seed)) {
    return (NULL)
  }

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    temp <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else
    temp <- NULL

  set.seed(new_seed)

  return (temp)
}

#' Restore random number generator state
#'
#' @param old_seed integer. Previous seed that should be restored globally.
restore_seed <- function (old_seed) {
  if (!is.null(old_seed))
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  else rm(.Random.seed, pos = 1)
}
