## interface function for the C implementation of Lueschers Zeta function

LuescherZeta <- function(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1, lambda=1, tol=0.000001, verbose=FALSE) {
  if(!is.numeric(qsq) || !is.numeric(l) || !is.numeric(m) || !is.numeric(dvec) || !is.numeric(gamma) || !is.numeric(lambda) || !is.numeric(tol))
    stop("argument qsq to LuescherZeta must be numeric!\n")
  n <- length(qsq)
  if(length(gamma) == 1) {
    gamma <- rep(gamma, times=n)
  }
  if(n != length(gamma)) {
    stop("gamma must be either a scalar or a vector of length identical with the length of qsq\n")
  }
  verb <- 0
  if(verbose) verb <- 1
  return(.Call("LuscherZetaArray", as.double(qsq), as.integer(n), as.integer(l[1]), as.integer(m[1]), as.double(dvec[1:3]), as.double(gamma), as.double(lambda[1]), as.double(tol), as.integer(verb)))
}
