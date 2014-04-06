## creates a square Hankel matrix of size nxn
hankel.matrix <- function(n, seq){
  outer(1:n,1:n,
        function(x,y) seq[x+y-1])
}

ev.hankel <- function(cf, t0, n, N, id=1, cutoff=1, eps=0.001) {
  if(t0+n > N) {
    stop("t0+n > N\n")
  }
  ## build the two Hankel matrices
  cM2 <- hankel.matrix(n, cf[t0:N])
  cM1 <- hankel.matrix(n, cf[(t0+1):N])

  qr.cM1 <- qr(cM1)
  M <- try(qr.solve(qr.cM1, cM2), TRUE)
  if(inherits(M, "try-error")) {
    return(NA)
  }

  M.eigen <- eigen(M, symmetric=FALSE, only.values=TRUE)

  ## reduce to real eigenvalues with real part larger than 1
  ii <- which(abs(Im(M.eigen$values)) < eps*abs(Re(M.eigen$values)) & Re(M.eigen$values) > cutoff )

  return(sort(Re(M.eigen$values[ii]), decreasing=FALSE)[id])
  
}

bootstrap.hankel <- function(cf, t0, n, N, id=1, cutoff=1, eps=0.001) {
  if(t0+n > N) {
    stop("t0+n > N\n")
  }

  if(cf$boot.samples == FALSE) {
    stop("Need to call bootstrap.cf first. Aborting...\n")
  }

  evs <- ev.hankel(cf$cf0, t0, n, N, id=id)
  l <- length(evs)

  evs.tsboot <- t(apply(cf$cf.tsboot$t, 1, ev.hankel, t0=t0, n=n, N=N, id=id, cutoff=cutoff, eps=eps))
  res <- list(evs=evs, evs.tsboot=evs.tsboot)
}
