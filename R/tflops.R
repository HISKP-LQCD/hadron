tflops <- function(L, T, N, tau, nconf=1000) {
  res <- L^3*T*nconf*tau*N*(1356+168)/(1000^4*365*24*60*60)
  return(res)
}
