## creates a square Hankel matrix of size nxn
hankel.matrix <- function(n, seq){
  outer(1:n,1:n,
        function(x,y) seq[x+y-1])
}

gevp.hankel.old <- function(cf, t0, deltat = 1, n, N, id=c(1),
                            range=c(0.,1.), eps=0.001, debug=FALSE, cor.id=1) {
  if(t0+2*n+deltat+1 > N) {
    stop("t0+n+deltat > N\n")
  }
  ## build the two Hankel matrices
  cM1 <- hankel.matrix(n, cf[((cor.id-1)*N+t0):((cor.id)*N)])
  cM2 <- hankel.matrix(n, cf[(((cor.id-1)*N+t0)+deltat):((cor.id)*N)])

  qr.cM1 <- qr(cM1)
  ##M <- try(qr.coef(qr.cM1, cM2), TRUE)
  M <- qr.coef(qr.cM1, cM2)
  if(inherits(M, "try-error")) {
    ##cM1.svd <- svd(cM1)
    ##D <- diag(1./cM1.svd$d)
    ##M <- cM1.svd$v %*% D %*% t(cM1.svd$u)
    
    return(rep(NA, times=length(id)))
  }

  M.eigen <- eigen(M, symmetric=FALSE, only.values=TRUE)

  ## reduce to real eigenvalues with real part larger than 1
  ii <- which(abs(Im(M.eigen$values)) < eps*abs(Re(M.eigen$values)) & Re(M.eigen$values) > range[1] & Re(M.eigen$values) < range[2] )
  if(debug) print(M.eigen$values[which(abs(Im(M.eigen$values)) < eps*abs(Re(M.eigen$values)))])
  return(sort(Re(M.eigen$values[ii]), decreasing=TRUE)[id])
  
}

gevp.hankel <- function(cf, t0=1, deltat=1, n, N, eps=0.0001, range=c(0,1),
                        submatrix.size=1, element.order=c(1,2,3,4), id=c(1),
                        debug=FALSE) {
  if(t0+2*n+deltat > N) {
    stop("t0+n+deltat > N\n")
  }
  t0 <- t0+1
  cM1 <- array(NA, dim=c(submatrix.size*n, submatrix.size*n))
  cM2 <- cM1
  ii <- seq(from=1, to=submatrix.size*n, by=submatrix.size)

  for(i in c(1:submatrix.size)) {
    for(j in c(1:submatrix.size)) {
      cor.id <- element.order[(i-1)*submatrix.size+j]
      cM1[ii+i-1,ii+j-1] <- hankel.matrix(n, cf[((cor.id-1)*N+t0):((cor.id)*N)])
      cM2[ii+i-1,ii+j-1] <- hankel.matrix(n, cf[(((cor.id-1)*N+t0)+deltat):((cor.id)*N)])
    }
  }
  ## symmetrise matrices
  cM1 <- 0.5*(cM1 + t(cM1))
  cM2 <- 0.5*(cM2 + t(cM2))

  qr.cM1 <- qr(cM1)

  M <- try(qr.coef(qr.cM1, cM2), TRUE)
  if(inherits(M, "try-error")) {
    ##cM1.svd <- svd(cM1)
    ##D <- diag(1./cM1.svd$d)
    ##M <- cM1.svd$v %*% D %*% t(cM1.svd$u)
    
    return(c(NA))
  }
  M[abs(M) < 1.e-12] <- 0
  M.eigen <- try(eigen(M, symmetric=FALSE, only.values=TRUE), TRUE)
  if(inherits(M.eigen, "try-error")) {
    return(NA)
  }
  ii <- which(abs(Im(M.eigen$values)) < eps & Re(M.eigen$values) > range[1]
              & Re(M.eigen$values) < range[2])
  if(id == "all") {
    return(invisible(Re(M.eigen$values[ii])))
  }
  return(invisible(Re(M.eigen$values[ii[id]])))
}

bootstrap.hankel <- function(cf, t0, deltat=1, n, N, id=c(1), range=c(0,1), eps=0.001,
                             element.order=c(1,2,3,4), submatrix.size=1) {
  if(t0+2*n+deltat >= N) {
    stop("t0+2*n+1+deltat must be smaller than N\n")
  }

  if(cf$boot.samples == FALSE) {
    stop("Need to call bootstrap.cf first. Aborting...\n")
  }

  evs <- gevp.hankel(cf$cf0, t0=t0, n=n, N=N, id=id, range=range, eps=eps, debug=FALSE, deltat=deltat,
                     submatrix.size=submatrix.size, element.order=element.order)
  l <- length(evs)

  evs.tsboot <- apply(cf$cf.tsboot$t, 1, gevp.hankel, t0=t0, n=n, N=N, id=id, range=range, eps=eps, deltat=deltat,
                      submatrix.size=submatrix.size, element.order=element.order)
  res <- list(evs=evs, evs.tsboot=evs.tsboot)
}


analyse.hankel <- function(cf, t0=1, range=c(1./exp(1),1.), n=5,
                           element.order=c(1,2,3,4), submatrix.size=1,
                           deltat, id="all") {

  ##def.par <- par(no.readonly = TRUE) # save default, for resetting...
  N=cf$Time/2+1
  if(missing(deltat)) {
    ##par(mfrow=c(length( c(1:(N-2*n-1-t0)) ), 1))
    for(deltat in c(1:(N-2*n-1-t0))) {
      bla <- bootstrap.hankel(cf, t0=t0, n=n, N=N, range=c(1./exp(2*deltat),1), deltat=deltat,
                              submatrix.size=submatrix.size, element.order=element.order, id=id)
      ##cat(bla$evs, "\n")
      ##cat("n=", n, -log(bla$evs), apply(-log(bla$evs.tsboot), 1, sd, na.rm=TRUE), "\n")
      new_window_if_appropriate()
      hist(-log(unlist(bla$evs.tsboot))/deltat, breaks=200, xlim=c(0,2), xlabel=c(""), ylabel=c(""))
    }
  }
  else {
    bla <- bootstrap.hankel(cf, t0=t0, n=n, N=N, range=c(1./exp(2*deltat),1), deltat=deltat,
                            submatrix.size=submatrix.size, element.order=element.order, id=id)
    ##cat(bla$evs, "\n")
    ##cat("n=", n, -log(bla$evs), apply(-log(bla$evs.tsboot), 1, sd, na.rm=TRUE), "\n")
    new_window_if_appropriate()
    hist(-log(unlist(bla$evs.tsboot))/deltat, breaks=200, xlim=c(0,2))
  }
  ##par(def.par)  #- reset to default
}
