## routine to block a vector or a two dimensional array in blocks of length l
## discards anything that doesn't fit into exact multiples of the block length
## this is much faster than the blocking
## of tsboot
block.ts <- function(data, l=2) {
  if(l == 1) {
    return(invisible(data))
  }
  if(is.vector(data)) {
    N <- floor(length(data)/l)*l
    return( apply(array(data, dim=c(l, N/l)), 2, mean))
  }
  if(length(dim(data))!=2) {
    stop("block.ts currently only implemented for vectors of 2-dim arrays\n")
  }
  N <- floor(length(data[,1])/l)*l
  ncf <- array(0, dim=c(N/l,length(data[1,])))
  j <- 1
  for ( i in seq(1,N,l)) {
    ncf[j,] <- apply(data[i:(i+l-1),], 2, mean)
    j <- j+1
  }
  return(invisible(ncf))
}
