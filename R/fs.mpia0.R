
fs.mpia0 <- function(mps, fps, L) {

  fn <- function(n, cn, mpsL) {
    cn*exp(-n*mpsL)/sqrt(n*mpsL)*(1-17/(8*n*mpsL))
  }
  cn <- c(6, 12, 8, 6, 24, 24, 0, 12)
  n <- c(1, sqrt(2), sqrt(3), 2, sqrt(5), sqrt(6), sqrt(7), sqrt(8))
  S <- sum(fn(n, cn, mpsL=mps*L))

  return(mps^4/fps^4/2^(13/2)/pi^(5/2)*S)
}
