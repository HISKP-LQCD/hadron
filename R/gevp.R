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
    sortindex <- order(variational.solve$values, decreasing=TRUE)
    evalues[t,] <- variational.solve$values[sortindex]
    evectors[t,,] <- variational.solve$vectors[, sortindex]
  }
  ## in case of bootstrapping everything (eigenvalues and eigenvectors)
  ## is concatenated into a single vector
  if(for.tsboot) {
    return(c(as.vector(evalues), as.vector(evectors)))
  }
  
  return(list(evalues=evalues, evectors=evectors))
}


bootstrap.gevp <- function(cf, t0=1, boot.R=400, boot.l=2, matrix.size=2, element.order=c(1,2,3,4), seed=1234) {
  ## number of measurements
  if(!any(class(cf) == "cf")) {
    stop("bootstrap.gevp requires an object of class cf as input! Aborting!\n")
  }
  N <- length(cf$cf[,1])
  if(!cf$boot.samples) {
    cf <- bootstrap.cf(cf, boot.R=boot.R, boot.l=boot.l, seed=seed)
  }
  res <- gevp(cf$cf0, Time=cf$Time, t0, matrix.size, element.order, for.tsboot=FALSE)


  gevp.tsboot <- t(apply(cf$cf.tsboot$t, 1, gevp, Time=cf$Time, t0=t0,
                         matrix.size=matrix.size, element.order=element.order,
                         for.tsboot=TRUE))

  ## gevp.tsboot contains first the N*(Thalf+1) eigenvalues
  ## and the the N*N*(Thalf+1) eigenvectors
  
  ret <- list(cf=cf, res.gevp=res, gevp.tsboot=gevp.tsboot)
  class(ret) <- c("gevp", class(ret))
  return(invisible(ret))
}

gevp2cf <- function(gevp, id=1) {
  cf <- list()
  cf$cf0 <- gevp$res.gevp$evalues[,id]
  cf$boot.samples <- TRUE
  cf$nrStypes <- 1
  cf$nrObs <- 1
  cf$Time <- gevp$cf$Time
  tt <- (id-1)*(cf$Time/2+1)+seq(1, cf$Time/2+1)
  cf$cf.tsboot <- list()
  cf$cf.tsboot$t <- gevp$gevp.tsboot[,tt]
  cf$cf.tsboot$t0 <- gevp$res.gevp$evalues[,id]
  attr(cf, "class") <- c("cf", class(cf))
  return(invisible(cf))
}
