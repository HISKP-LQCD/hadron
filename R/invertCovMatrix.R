invertCovMatrix <- function(cf, boot.l=1, boot.samples=FALSE) {
  ## compute compute the correctly normalised inverse of a noisy covariance matrix
  ## see C. Michael hep-lat/9412087

  ## block data first
  ncf <- cf
  if(boot.l > 1) {
    ncf <- block.ts(cf, l=boot.l)
  }
  ## compute covariance matrix and invert
  CovMatrix <- cov(ncf)
  ## we have a real, symmetric square matrix with dimension n=length(ncf[1,])
  n <- length(ncf[1,])
  ## the number of observations
  N <- length(ncf[,1])
  M <- matrix()
  if(floor(sqrt(n)) < 0*N) {
    ## use singular value decomposition
    cov.svd <- svd(CovMatrix)
    ## replace smallest singular values by their mean, if needed
    cov.svd$d[floor(sqrt(n)):n] <-
      mean(cov.svd$d[floor(sqrt(n)):n])
    ## construct the inverse
    D <- diag(1./cov.svd$d)
    M <- cov.svd$v %*% D %*% t(cov.svd$u)
  }
  else {
    ## use cholesky decomposition for real symmetric matrix
    M <- chol2inv(chol(CovMatrix))
  }
  if(!boot.samples) {
    M <- N*M
  }
  return(invisible(M))
}
