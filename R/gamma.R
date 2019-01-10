#' Accessor function for \code{\link{gm}}
#'
#' Retrieve the entries of the \code{\link[hadron]{gm}} three-index array containing various
#' gamma structures in a natural indexing.
#'
#' @param mu Number or string denoting
#' \itemize{
#'   \item Lorentz index (0,1,2,3) for \eqn{\gamma^\mu}
#'   \item 5 for \eqn{\gamma^5}
#'   \item "Pp" or "Pm" for the positive and negative parity projectors respectively
#' }
#' 
#' @concept Dirac
#' @concept gamma
#' @export
gm_mu <- function(mu){
  if( is.numeric(mu) ){
    if( mu < 0 | mu > 5 ){
      stop(sprintf("mu=%d is not a valid input", mu))
    } else {
      if( mu == 5 ){
        return(gm[mu,,])
      } else {
        return(gm[mu+1,,])
      }
    }
  } else if( is.character(mu) ){
    if( mu == "Pp" ){
      return(gm[6,,])
    } else if ( mu == "Pm" ){
      return(gm[7,,])
    } else {
      stop(sprintf("mu=%s is not a valid input", mu))
    }
  } else {
    # whatever was given, let's transform it into a string and attempt to print it
    stop(sprintf("mu=%s is not a valid input", mu))
  }
}

#' Array of gamma structures
#'
#' Array of 4x4 complex gamma matrices in the tmLQCD chiral gamma basis,
#' where \eqn{\gamma^5 = \gamma^0 \gamma^1 \gamma^2 \gamma^3 = } diag(c(1,1,-1,-1))
#'
#' The index mapping is as follows
#' \itemize{
#'   \item \code{gm[1,,]} \eqn{\gamma^0}
#'   \item \code{gm[2,,]} \eqn{\gamma^1}
#'   \item \code{gm[3,,]} \eqn{\gamma^2}
#'   \item \code{gm[4,,]} \eqn{\gamma^3}
#'   \item \code{gm[5,,]} \eqn{\gamma^5}
#'   \item \code{gm[6,,]} positive parity projector \eqn{ \frac{1}{2} (1 + \gamma^0) }
#'   \item \code{gm[7,,]} negative parity projector \eqn{ \frac{1}{2} (1 - \gamma^0) }
#' }
#'
#' The function \code{\link{gm_mu}} can be used to access its elements using a more
#' "natural" indexing.
#.
#' @concept Dirac
#' @concept gamma
#' @export
gm <- array(complex(16), dim=c(7,4,4))

## gamma0
gm[1,3,1] <- -1.
gm[1,4,2] <- -1.
gm[1,1,3] <- -1.
gm[1,2,4] <- -1.

## gamma1
gm[2,4,1] <- -1i
gm[2,3,2] <- -1i
gm[2,2,3] <- +1i
gm[2,1,4] <- +1i

## gamma2
gm[3,4,1] <- -1.
gm[3,3,2] <- +1.
gm[3,2,3] <- +1.
gm[3,1,4] <- -1.

## gamma3
gm[4,3,1] <- -1i
gm[4,4,2] <- +1i
gm[4,1,3] <- +1i
gm[4,2,4] <- -1i

# gamma5
gm[5,1,1] <- 1
gm[5,2,2] <- 1
gm[5,3,3] <- -1
gm[5,4,4] <- -1

# parity projectors
gm[6,,] <- 0.5*( diag(1,4,4) + gm_mu(0) )
gm[7,,] <- 0.5*( diag(1,4,4) - gm_mu(0) )

sigmamunu <- array(complex(16), dim=c(4,4,4,4))

tmatrix <- array(0., dim=c(4,4))
tmatrix[1,1] <- 1
tmatrix[1,2] <- 10
tmatrix[2,1] <- 100
tmatrix[2,2] <- 1000

tmatrix[3,3] <- 10000
tmatrix[3,4] <- 100000
tmatrix[4,3] <- 1000000
tmatrix[4,4] <- 10000000

for(mu in 1:4) {
  for(nu in 1:4) {
    sigmamunu[mu, nu,,] <- gm[mu,,] %*% gm[nu,,]
  }
}

Tr <- function(x) {
  return(sum(diag(x)))
}
