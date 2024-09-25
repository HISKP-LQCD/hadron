#' @title Lanczos method for LQCD correlators
#' 
#' @description
#'   blub ...
#'
#' @param cf object of type \link{cf}
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Lanczos analysis
#' @return
#'   tbw
#' 
#' @family lanczos
#' @export
bootstrap.lanczos <- function(cf, N = (cf$Time/2+1)) {
  ## wrapper function, not yet bootstrapping...
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  res <- lanczos.solve(cf=cf$cf0, N=N)
  return(res)
}

#' @title Lanczos solver
#' 
#' @description
#'   blub ...
#'
#' @param cf Numeric vector (this will generally be a correlation function or a bootstrap sample thereof).
#' @param N Integer. Maximal time index in correlation function to be used in
#'                   Lanczos analysis
#' @return
#'   tbw
#' 
#' @family lanczos
lanczos.solve <- function(cf, N) {
  ## container for the eigenvalues per m
  evs <- rep(NA, times=N/2)
  for(m in c(1:(N/2))) {
    Aj <- cf[2:N]/cf[1]
    Bj <- rep(0, times=length(Aj))
    Bjp1 <- Bj
    Gj <- Bj
    Gjp1 <- Gj
    Ajm1 <- Bj
    Ajp1 <- Bj
    
    alpha <- rep(NA, times=m)
    beta <- alpha
    gamma <- alpha
    alpha[1] <- Aj[1]
    beta[1] <- 0.
    gamma[1] <- 0.

    if(m > 1) {
      for(j in c(1:(m-1))) {
        ## eq.(53) suppl. mat.
        sr <- Aj[2]-alpha[j]^2-beta[j]*gamma[j]
        ## eq.(41) suppl. mat.
        rhojp1 <- sqrt(abs(sr))
        taujp1 <- sr/rhojp1
        ## below eq.(28) suppl. mat.
        kmax <- 2*(m-j) + 1
        ## B_{j+1}^k and G_{j+1}^k
        for(k in c(1:kmax)) {
          ## eqs. (44) and (45) suppl. mat.
          Gjp1[k] <- Aj[k+1] - alpha[j]*Aj[k] - gamma[j]*Bj[k]
          Bjp1[k] <- Aj[k+1] - alpha[j]*Aj[k] - beta[j]*Gj[k]
        }
        Gjp1 <- Gjp1/taujp1
        Bjp1 <- Bjp1/rhojp1
        ## A_{j+1}^k eq. (46) suppl. mat.
        for(k in c(1:kmax)) {
          Ajp1[k] <-  Aj[k+2] - 2*alpha[j]*Aj[k+1] + alpha[j]^2*Aj[k] +
            alpha[j]*(beta[j]*Gj[k] + gamma[j]*Bj[k]) -
            (beta[j]*Gj[k+1] + gamma[j]*Bj[k+1]) +
            gamma[j]*beta[j]*Ajm1[k]
        }
        Ajp1 <- Ajp1/rhojp1/taujp1
        ## eq. (38) suppl. mat.
        alpha[j + 1] <- Ajp1[1]
        beta[j + 1] <- Bjp1[1]
        gamma[j + 1] <- Gjp1[1]
        ## don't do the copying if not needed
        if(j == m-1) break
        Ajm1 <- Aj
        Aj <- Ajp1
        Bjm1 <- Bj
        Bj <- Bjp1
        Gjm1 <- Gj
        Gj <- Gjp1
      }
    }
    if(m == 1) {
      ## 1x1 case
      evs[m] <- alpha[1]
    }
    else{
      ##cat("m= ", m, "\n\n")
      ## construct mxm tridiagonal matrix
      M <- diag(alpha[1:m])
      ii <- c(1:(m-1))
      M[row(M) - col(M) == -1] <- beta[2:m]
      M[row(M) - col(M) == +1] <- gamma[2:m]
      ##print(M)
      
      ## eigensolve and extract the lowest eigenvalue
      eigvalues <- eigen(M, symmetric=FALSE, only.values = TRUE, EISPACK = FALSE)$values
      ##cat("\n")
      eigvalues <- Re(eigvalues[Im(eigvalues) == 0])
      eigvalues <- eigvalues[eigvalues > 0 & eigvalues < 1]
      ##print(eigvalues)
      evs[m] <- sort(eigvalues, decreasing=TRUE)[1]
    }
  }
  return(evs)
}
