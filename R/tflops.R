tflops <- function(L, Time, N, tau, nconf=1000) {
  res <- L^3*Time*nconf*tau*N*(1356+168)/(1000^4*365*24*60*60)
  return(res)
}
