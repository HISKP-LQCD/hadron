## interface function for the C implementation of Lueschers Zeta function

LuescherZeta <- function(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1, lambda=1, N=8, N2=3) {
  if(!is.numeric(qsq) || !is.numeric(l) || !is.numeric(m) || !is.numeric(dvec) || !is.numeric(gamma) || !is.numeric(lambda) || !is.numeric(N) || !is.numeric(N2))
    stop("argument qsq to LuescherZeta must be numeric!\n")
  n <- length(qsq)
  if(n==1) return(.Call("LuscherZeta", as.double(qsq), as.integer(l[1]), as.integer(m[1]), as.double(dvec[1:3]), as.double(gamma[1]), as.double(lambda[1]), as.integer(N[1]), as.integer(N2[1])))

  return(.Call("LuscherZetaArray", as.double(qsq), as.integer(n), as.integer(l[1]), as.integer(m[1]), as.double(dvec[1:3]), as.double(gamma[1]), as.double(lambda[1]), as.integer(N[1]), as.integer(N2[1])))
}
