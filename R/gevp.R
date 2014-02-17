## t0 start counting at 0!
##                           ( 1 , 2 )
## order c(1,2,3,4) goes to 
##                           ( 3 , 4 )

gevp <- function(cf, Time, t0=1, matrix.size=2, element.order=c(1,2,3,4), for.tsboot=TRUE) {
  
  Thalf <- Time/2
  if(length(dim(cf)) == 2) {
    Cor <- apply(cf, 2, mean)
  }
  else {
    Cor <- cf
  }

  ## need to check consistency of cf here!
  ## can only operate on a square matrix
  if(length(Cor)/(Thalf+1) != matrix.size^2) {
    stop("gevp can only operate on square matrices! Aborting!\n")
  }

  ## symmetrise Cor
  tt <- c(1:(Thalf+1))
  for(i in 1:matrix.size) {
    for(j in i:matrix.size) {
      Cor[((i-1)*matrix.size + (j-1))*(Thalf+1)+tt] <- 0.5*(Cor[((i-1)*matrix.size + (j-1))*(Thalf+1)+tt]+
                                                            Cor[((j-1)*matrix.size + (i-1))*(Thalf+1)+tt])
      Cor[((j-1)*matrix.size + (i-1))*(Thalf+1)+tt] <- Cor[((i-1)*matrix.size + (j-1))*(Thalf+1)+tt]
    }
  }

  ii <- c()
  for(i in 1:matrix.size^2) {
    ii <- c(ii, (i-1)*(Thalf+1)+1)
  }
  ## re-order as to match the input order
  ii <- ii[element.order]
  
  evalues <-  array(NA, dim=c(Thalf+1, matrix.size))
  evectors <- array(NA, dim=c(Thalf+1, matrix.size, matrix.size))
  ## matrix at t=t0
  cM <- matrix(Cor[ii+t0], nrow=matrix.size, ncol=matrix.size)
  L <- chol(cM)
  invL <- solve(L)
  ## now the time dependence
  ## we need to multiply from the left with t(invL) and from the right with invL
  for(t in (t0+1):(Thalf)) {
    variational.solve <- eigen(t(invL) %*% matrix(Cor[ii+t], nrow=matrix.size, ncol=matrix.size) %*% invL,
                               symmetric=TRUE, only.values = FALSE, EISPACK=FALSE)
    sortindex <- order(variational.solve$values)
    evalues[t,] <- variational.solve$values[sortindex]
    evectors[t,,] <- variational.solve$vectors[, sortindex]
  }
  if(for.tsboot) {
    return(c(as.vector(evalues), as.vector(evectors)))
  }
  
  return(list(evalues=evalues, evectors=evectors))
}


bootstrap.gevp <- function(cf, t0, boot.R=400, boot.l=2, matrix.size=2, element.order=c(1,2,3,4), seed=1234) {
  ## number of measurements
  N <- length(cf$cf[,1])
  res <- gevp(cf$cf, Time=cf$Time, t0, matrix.size, element.order, for.tsboot=FALSE)
  ## we set the seed for reproducability and correlation
  set.seed(seed)
  ## and bootstrap the GEVP
  gevp.tsboot <- tsboot(tseries=cf$cf, statistic=gevp, R=boot.R, l=boot.l, sim="geom",
                        Time=cf$Time, t0=t0, matrix.size=matrix.size, element.order=element.order,
                        for.tsboot=TRUE)
  

}
