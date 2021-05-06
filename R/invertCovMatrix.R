#' Inverts the covariance matrix for noisy data
#' 
#' The covariance matrix of noisy data is inverted. Special care is taken in
#' treating spurious small modes of the matrix, which are likely to arise from
#' too much noise in the data.
#' 
#' The inverse covariance matrix is estimated. If the number of observations is
#' too small the procedure described in the reference is used to remove
#' spuriously small eigenvalues of the covariance matrix.
#' 
#' We always keep the \eqn{\sqrt{R}}{sqrt(R)} largest eigenvalues exactly and
#' replace the remaining smallest ones by their mean.
#' 
#' @param cf The data for which the covariance matrix is to be computed. It is
#' expected to be an array or matrix with dimension RxN, where R is the number
#' of observations and N the number of observables.
#' 
#' \code{cf} can be either real data or bootstrap data. In the latter case
#' \code{boot.samples=TRUE} must be set for proper normalisation of the inverse
#' matrix.
#' @param boot.l If set to a value larger than 1 the data will be blocked with
#' blocklength \code{boot.l} before the covariance matrix is computed.
#' @param boot.samples If set to \code{TRUE} the data is treated a pseudo data
#' from a bootstrap procedure.
#' @param cov_fn Function that computes the covariance matrix from the given
#' samples.
#' @return Returns the inverse covariance matrix as an object of class
#' \code{\link{matrix}}.
#' @author Carsten Urbach, \email{curbach@@gmx.de}
#' @seealso \code{\link{cov}}, \code{\link{matrix}}
#' @references C.Michael, A.McKerrell, Phys.Rev. D51 (1995) 3745-3750,
#' hep-lat/9412087
#' @keywords covariance matrix correlated chisqr
#' @examples
#' X <- array(rnorm(4000), dim=c(1000, 4))
#' invertCovMatrix(cf=X, boot.samples=TRUE)
#' M <- invertCovMatrix(cf=X, boot.samples=TRUE) 
#' M
#' 
#' @export invertCovMatrix
invertCovMatrix <- function(cf, boot.l=1, boot.samples=FALSE, cov_fn = cov) {
  ## compute compute the correctly normalised inverse of a noisy covariance matrix
  ## see C. Michael hep-lat/9412087

  ## block data first, this should only be done if we're not dealing with boostrap samples
  ## because these already stem from an appropriate block sampling procedure
  ncf <- cf
  if(boot.l > 1 && boot.samples == FALSE) {
    ncf <- block.ts(cf, l=boot.l)
  }
  ## compute covariance matrix and invert
  CovMatrix <- cov_fn(ncf)
  ## we have a real, symmetric square matrix with dimension n=length(ncf[1,])
  n <- length(ncf[1,])
  ## the number of observations
  N <- length(ncf[,1])
  M <- matrix()

  if(n > floor(sqrt(N))) {
    ## use singular value decomposition
    cov.svd <- svd(CovMatrix)
    ## replace smallest singular values by their mean, if needed
    ## we keep floor(sqrt(N)) exact eigenvalues and replace all smaller once
    ## by their average value
    cov.svd$d[floor(sqrt(N)):n] <-
      mean(cov.svd$d[floor(sqrt(N)):n])
    ## construct the inverse
    D <- diag(1./cov.svd$d)
    M <- cov.svd$v %*% D %*% t(cov.svd$u)
  }
  else {
    ## use cholesky decomposition for real symmetric matrix
    M <- chol2inv(chol(CovMatrix))
  }
  ## for bootstrap samples the error is equal to sd
  ## otherwise to sd/sqrt(N)
  if(!boot.samples) {
    M <- N*M
  }
  return(invisible(M))
}
