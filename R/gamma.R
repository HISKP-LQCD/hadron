#' @title Accessor function for \code{\link{gm}}
#'
#' @description
#' Retrieve the entries of the \code{\link[hadron]{gm}} list of 
#' three-index arrays containing various
#' gamma structures in a natural indexing.
#'
#' @param mu Number or string denoting
#' \itemize{
#'   \item Lorentz index (0,1,2,3,4) for \eqn{\gamma^\mu}
#'   \item 5 for \eqn{\gamma^5}
#'   \item "Pp" or "Pm" for the positive and negative parity projectors respectively
#' }
#' @param basis String, gamma basis to use. Possible values
#'              \describe{
#'                \item{'ukqcd':}{UKQCD gamma basis with \eqn{\gamma^i, i \in [1,2,3,4]} 
#'                                and \eqn{\gamma^5 = \gamma^1 \gamma^2 \gamma^3 \gamma^4},
#'                                such that \eqn{1 = x}, \eqn{4 = t}.}
#'                \item{'chiral_tmlqcd':}{Chiral gamma basis used by tmLQCD with 
#'                                        \eqn{\gamma^\mu, \mu \in [0,1,2,3]} and 
#'                                        \eqn{\gamma^5 = \gamma^0 \gamma^1 \gamma^2 \gamma^3},
#'                                        such that \eqn{0 = t}, \eqn{3 = z}.}
#'              }
#'
#' @concept Dirac
#' @concept gamma
#' 
#' @return
#' Returns the requested \eqn{\gamma}{gamma} matrix as a 4x4 complex valued
#' array, see \link{gm}.
#' 
#' @export
gm_mu <- function(mu, basis = 'chiral_tmlqcd'){
  if( !(basis %in% c('ukqcd', 'chiral_tmlqcd')) ){
    stop("Only 'ukqcd' and 'chiral_tmlqcd' gamma bases are currently supported!")
  }
  if( is.numeric(mu) ){
    if( basis == 'chiral_tmlqcd' ){
      if( mu < 0 | mu > 5 ){
        stop(sprintf("mu=%d is not a valid input", mu))
      } else {
        if( mu == 5 ){
          return(gm[[basis]][mu,,])
        } else {
          return(gm[[basis]][mu+1,,])
        }
      }
    } else if ( basis == 'ukqcd' ){
      if( mu < 1 | mu > 5 ){
        stop(sprintf("mu=%d is not a valid input", mu))
      } else {
        return(gm[[basis]][mu,,])
      }
    }
  } else if( is.character(mu) ){
    if( mu == "Pp" ){
      return(gm[[basis]][6,,])
    } else if ( mu == "Pm" ){
      return(gm[[basis]][7,,])
    } else {
      stop(sprintf("mu=%s is not a valid input", mu))
    }
  } else {
    # whatever was given, let's transform it into a string and attempt to print it
    stop(sprintf("mu=%s is not a valid input", mu))
  }
}

#' @name gm
#' 
#' @title List of arrays of gamma structures
#'
#' @description
#' List of arrays of 4x4 complex gamma matrices in the tmLQCD chiral gamma basis,
#' where \eqn{\gamma^5 = \gamma^0 \gamma^1 \gamma^2 \gamma^3 = } \code{diag(c(1,1,-1,-1))}
#' and the UKQCD gamma basis, where \eqn{\gamma^5 = \gamma^0 \gamma^1 \gamma^2 \gamma^3}.
#'
#' The index mappings are as follows
#' \itemize{
#'   \item \code{gm[['chiral_tmlqcd']][1,,]} \eqn{\gamma^0}
#'   \item \code{gm[['chiral_tmlqcd']][2,,]} \eqn{\gamma^1}
#'   \item \code{gm[['chiral_tmlqcd']][3,,]} \eqn{\gamma^2}
#'   \item \code{gm[['chiral_tmlqcd']][4,,]} \eqn{\gamma^3}
#'   \item \code{gm[['chiral_tmlqcd']][5,,]} \eqn{\gamma^5}
#'   \item \code{gm[['chiral_tmlqcd']][6,,]} positive parity projector \eqn{ \frac{1}{2} (1 + \gamma^0) }
#'   \item \code{gm[['chiral_tmlqcd']][7,,]} negative parity projector \eqn{ \frac{1}{2} (1 - \gamma^0) }
#' }
#' \itemize{
#'   \item \code{gm[['ukqcd']][1,,]} \eqn{\gamma^1}
#'   \item \code{gm[['ukqcd']][2,,]} \eqn{\gamma^2}
#'   \item \code{gm[['ukqcd']][3,,]} \eqn{\gamma^3}
#'   \item \code{gm[['ukqcd']][4,,]} \eqn{\gamma^4}
#'   \item \code{gm[['ukqcd']][5,,]} \eqn{\gamma^5}
#'   \item \code{gm[['ukqcd']][6,,]} positive parity projector \eqn{ \frac{1}{2} (1 + \gamma^4) }
#'   \item \code{gm[['ukqcd']][7,,]} negative parity projector \eqn{ \frac{1}{2} (1 - \gamma^4) }
#' }
#'
#' The function \code{\link{gm_mu}} can be used to access its elements using a more
#' "natural" indexing.
#'
#' @docType data
#' @concept Dirac
#' @concept gamma
NULL

gm <- list()
  
gm[['chiral_tmlqcd']] <- array(complex(16), dim=c(7,4,4))
## gamma0
gm[['chiral_tmlqcd']][1,3,1] <- -1.
gm[['chiral_tmlqcd']][1,4,2] <- -1.
gm[['chiral_tmlqcd']][1,1,3] <- -1.
gm[['chiral_tmlqcd']][1,2,4] <- -1.
## gamma1
gm[['chiral_tmlqcd']][2,4,1] <- -1i
gm[['chiral_tmlqcd']][2,3,2] <- -1i
gm[['chiral_tmlqcd']][2,2,3] <- +1i
gm[['chiral_tmlqcd']][2,1,4] <- +1i
## gamma2
gm[['chiral_tmlqcd']][3,4,1] <- -1.
gm[['chiral_tmlqcd']][3,3,2] <- +1.
gm[['chiral_tmlqcd']][3,2,3] <- +1.
gm[['chiral_tmlqcd']][3,1,4] <- -1.
## gamma3
gm[['chiral_tmlqcd']][4,3,1] <- -1i
gm[['chiral_tmlqcd']][4,4,2] <- +1i
gm[['chiral_tmlqcd']][4,1,3] <- +1i
gm[['chiral_tmlqcd']][4,2,4] <- -1i
# gamma5
gm[['chiral_tmlqcd']][5,1,1] <- 1
gm[['chiral_tmlqcd']][5,2,2] <- 1
gm[['chiral_tmlqcd']][5,3,3] <- -1
gm[['chiral_tmlqcd']][5,4,4] <- -1
# parity projectors
gm[['chiral_tmlqcd']][6,,] <- 0.5*( diag(1,4,4) + gm_mu(0,'chiral_tmlqcd') )
gm[['chiral_tmlqcd']][7,,] <- 0.5*( diag(1,4,4) - gm_mu(0,'chiral_tmlqcd') )

gm[['ukqcd']] <- array(complex(16), dim=c(7,4,4))
## gamma1
gm[['ukqcd']][1,,] <- -gm_mu(1, 'chiral_tmlqcd')
## gamma2
gm[['ukqcd']][2,,] <- -gm_mu(2, 'chiral_tmlqcd')
## gamma3
gm[['ukqcd']][3,,] <- -gm_mu(3, 'chiral_tmlqcd')
## gamma4
gm[['ukqcd']][4,1,1] <- +1
gm[['ukqcd']][4,2,2] <- +1
gm[['ukqcd']][4,3,3] <- -1
gm[['ukqcd']][4,4,4] <- -1
# gamma5
gm[['ukqcd']][5,,] <- gm_mu(1, 'ukqcd') %*% gm_mu(2, 'ukqcd') %*%
                      gm_mu(3, 'ukqcd') %*% gm_mu(4, 'ukqcd')
# parity projectors
gm[['ukqcd']][6,,] <- 0.5*( diag(1,4,4) + gm_mu(4,'ukqcd') )
gm[['ukqcd']][7,,] <- 0.5*( diag(1,4,4) - gm_mu(4,'ukqcd') )

sigmamunu <- list()
  
sigmamunu[['chiral_tmlqcd']] <- array(complex(16), dim=c(4,4,4,4))
sigmamunu[['ukqcd']] <- array(complex(16), dim=c(4,4,4,4))

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
    sigmamunu[['chiral_tmlqcd']][mu, nu,,] <- gm[['chiral_tmlqcd']][mu,,] %*% gm[['chiral_tmlqcd']][nu,,]
    sigmamunu[['ukqcd']][mu, nu,,] <- gm[['ukqcd']][mu,,] %*% gm[['ukqcd']][nu,,]
  }
}

Tr <- function(x) {
  return(sum(diag(x)))
}
