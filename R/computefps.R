computefps <- function(PP, mass, mu1, mu2, Kappa) {
  return(4*Kappa*(mu1+mu2)/2/sqrt(2)*abs(PP)/sqrt(mass^3))
}
