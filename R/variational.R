#' variational analysis of a correlator matrix
#' 
#' variational analysis of a correlator matrix
#' 
#' 
#' @param Cor correlator matrix
#' @param ta first time value for eigenvector determination
#' @param tb second time value for eigenvector determination
#' @param tmax maximal time value to be considered in the analysis
#' @param N matrix size. must be <= to the size of the correlator matrix
#' @param T1 maximal time value in correlator matrix
#' @param matrix.size matrix size to be used for start value determination
#' (e.g. for pionfit)
#' @param no.masses number of mass values for start value determination
#' @return returns a list with following entries \item{t}{ the list of time
#' values } \item{res.values}{ the mass values in a \code{tmax-ta} times
#' \code{N} array for all analysed t-values } \item{par}{ startvalues for a fit
#' } \item{variational.masses}{ the list of mass values as determined from the
#' variational analysis }
#' @author Carsten Urbach, \email{carsten.urbach@@physik.hu-berlin.de}
#' @seealso pion
#' @keywords variational
variational <- function(Cor, ta, tb, tmax, N, T1, matrix.size, no.masses) {

  ta=ta+1
  tb=tb+1
  if(ta > tb) {
    tmp <- ta
    ta <- tb
    tb <- tmp
  }
  res.values <- array(0., dim=c(tmax-ta,N))
  Ca <- getNxNmatrix(Cor=Cor, t=ta, T1=T1, N=N)
  Cb <- getNxNmatrix(Cor=Cor, t=tb, T1=T1, N=N)
  C3 <- getNxNmatrix(Cor=Cor, t=tb, T1=T1, N=N)

                                        # first index is rows, second columns
  Ca.inv <- solve(Ca)
  C3 <- Ca.inv %*% Cb
  variational.solve <- eigen(C3, symmetric=FALSE, only.values = FALSE, EISPACK=FALSE)
  # check and sort eigenvalues
  for(i in 1:N) {
    if(abs(variational.solve$values[i]) > 0.95/(tb-ta)) {
      variational.solve$values[i] <- 0.0001
    }
  }
  sortindex <- order(-log(abs(variational.solve$values)*(tb-ta)))
  # get the left eigenvectors, the eigenvectors have unit length
  # this does not quite work for the pion ?
  left.vectors <- crossprod(Ca, variational.solve$vectors)

  X <- crossprod(left.vectors, variational.solve$vectors)
  for(i in 1:N) {
#    left.vectors[,i] <- left.vectors[,i]/X[i,i]
    variational.solve$vectors[,i] <- variational.solve$vectors[,i]/X[i,i]
  }
  variational.masses <-  -log(abs(variational.solve$values[sortindex]))/(tb-ta)
  for(t in tb:(tmax)) {
    values <- numeric(N)
    ca <- getNxNmatrix(Cor=Cor, t=t-1, T1=T1, N=N)
    cb <- getNxNmatrix(Cor=Cor, t=t, T1=T1, N=N)
                                        # first index is rows, second columns
    ca.inv <- solve(ca)
    c3 <- ca.inv %*% cb
    X <- crossprod(left.vectors, c3 %*% variational.solve$vectors)
    for(j in 1:N) {
      values[j] <- X[j,j]
    }
    rm(X)
                                        #  cat(ta, tb, variational.masses, variational.solve$values, "\n")
    res.values[t-ta,] <- -log(abs(values[sortindex]))
  }
      
  par <- c(2*left.vectors[(1:matrix.size),sortindex[1]],
           -log(abs(variational.solve$values[sortindex[1]]))/(tb-ta))
  if(no.masses > 1) {
    for(i in 2:(no.masses)) {
      par <- c(par,
               2*left.vectors[(1:matrix.size),sortindex[i]],
               -log(abs(variational.solve$values[sortindex[i]]))/(tb-ta)) 
    }
  }

  return(invisible(list(t=c(tb:tmax), res.values=res.values, par=par,
                        variational.masses = variational.masses)))
  
}

#variational.analysis <- function(cmicor,ta=3, tb=4, tmax, N) {


#}
