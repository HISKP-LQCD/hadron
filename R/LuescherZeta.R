## interface function for the C implementation of Luescher Zeta function

LuescherZeta <- function(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1, lambda=1, N=10) {
  if(length(qsq)==1) return(.Call("LuscherZeta", qsq, l, m, dvec, gamma, lambda, N))
  n <- length(qsq)
  return(.Call("LuscherZetaArray", qsq, n, l, m, dvec, gamma, lambda, N))
}
