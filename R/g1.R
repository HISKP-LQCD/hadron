#' g1
#'
#' @description
#' Implementation of the Gasser-Leutwyler function g_1 for
#' computing finite volume effects.
#' 
#' @param x Numeric. x-value
#'
#' @export
g1 <- function(x) {

  weights <- c(6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.)
  ex <- c(1:20)
  res <- x
  for( i in 1:length(x)) {
    sex <- x[i]*sqrt(ex)
    res[i] <- sum(4*weights*besselK(sex, 1)/(sex))
  }
  return(res)
}
