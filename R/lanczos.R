bootstrap.lanczos <- function(cf, m = 2, N = (cf$Time/2+1)) {
  ## wrapper function, not yet bootstrapping...
  stopifnot(inherits(cf, 'cf_meta'))
  stopifnot(inherits(cf, 'cf_boot'))

  res <- lanczos.solve(cf=cf0, N=N)
  return(res)
}

lanczos.solve <- function(cf, N) {
  N <- length(cf)
  Aj <- cf[2:N]/cf[1]
  ## container for the eigenvalues per m
  evs <- rep(NA, times=Time/2)
  for(m in c(1:N)) {
    B <- rep(0, times=length(Aj))
    Bjp1 <- B
    G <- B
    Gjp1 <- G
    Ajm1 <- B
    Ajp1 <- B
    
    alpha <- rep(NA, times=m)
    beta <- alpha
    gamma <- alpha
    alpha[1] <- Aj[1]
    beta[1] <- 0.
    gamma[1] <- 0.

    if(m > 1) {
      for(jp1 in c(2:m)) {
        ## eq.(53) suppl. mat.
        sr <- A[2]-alpha[jp1-1]-beta[jp1-1]*gamma[jp1-1]
        ## eq.(41) suppl. mat.
        rhojp1 <- sqrt(abs(sr))
        taujp1 <- sr/rhojp1
        ## below eq.(28) suppl. mat.
        kmax <- 2*(m-jp1) + 1
        ## B_{j+1}^k and G_{j+1}^k
        for(k in c(2:kmax)) {
          ## eqs. (44) and (45) suppl. mat.
          Gjp1[k] <- Aj[k+1] - alpha[jp1-1]*Aj[k] - gamma[jp1-1]*Bj[k]
          Bjp1[k] <- Aj[k+1] - alpha[jp1-1]*Aj[k] - beta[jp1-1]*Gj[k]
        }
        Gjp1 <- Gjp1/taujp1
        Bjp1 <- Bjp1/rhojp1
        ## A_{j+1}^k eq. (46) suppl. mat.
        for(k in c(2:kmax)) {
          Ajp1[k] <-  Aj[k+2] - 2*alpha[jp1-1]*Aj[k+1] + alpha[jp1-1]^2*Aj[k] +
            alpha[jp1-1]*(beta[jp1-1]*Gj[k] + gamma[jp1-1]*Bj[k]) -
            (beta[jp1-1]*Gj[k+1] + gamma[jp1-1]*Bj[k+1]) +
            gamma[jp1-1]*beta[jp1-1]*Ajm1[k]
        }
        Ajp1 <- Ajp1/rhojp1/taujp1
        ## eq. (38) suppl. mat.
        alpha[jp1] <- Ajp1[1]
        beta[jp1] <- Bjp1[1]
        gamma[jp1] <- Gjp1[1]
        ## don't do the copying if not needed
        if(jp1 == m) break
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
      ## construct mxm tridiagonal matrix
      M <- diag(alpha[1:m])
      ii <- c(1:(m-1))
      M[row(M) - col(M) == -1] <- beta[2:m]
      M[row(M) - col(M) == +1] <- gamma[2:m]
      ## eigensolve and extract the lowest eigenvalue
      eigvalues <- eigen(M, symmetric=FALSE, only.values = TRUE, EISPACK = FALSE)$values
      evs[m] <- sort(eigvalues[eigvalues > 0])[1]
    }
  }
  return(evs)
}
