computefps <- function(mfit, PP, mass, mu1, mu2, Kappa, normalisation="cmi", reduce.latarts=FALSE) {
  if(missing(mu1)) {
    stop("computefps: mu1  must be specified! Aborting...\n")
  }
  if(missing(Kappa)) {
    Kappa <- NA
  }
  if(missing(mu2)) mu2 <- mu1
  if(missing(mfit)) {
    denominator <- sqrt(mass^3)
    if(reduce.latarts) denominator <- sinh(mass)*sqrt(mass)
    
    if(normalisation != "cmi") return((mu1+mu2)*abs(PP)/denominator)
    else return(sqrt(2)*Kappa*(mu1+mu2)*abs(PP)/denominator)
  }
  else if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa
    
    denominator <- sqrt(mfit$opt.res$par[1]^3)
    if(reduce.latarts) denominator <- sinh(mfit$opt.res$par[1])*sqrt(mfit$opt.res$par[1])
    mfit$fps <- 2*k*(mu1+mu2)/sqrt(2)*abs(mfit$opt.res$par[2])/denominator

    denominator <- sqrt(mfit$opt.tsboot[1,]^3)
    if(reduce.latarts) denominator <- sinh(mfit$opt.tsboot[1,])*sqrt(mfit$opt.tsboot[1,])
    mfit$fps.tsboot <- sqrt(2.)*k*(mu1+mu2)*abs(mfit$opt.tsboot[2,])/denominator

    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$reduce.latarts <- reduce.latarts
    return(invisible(mfit))
  }
  else if(inherits(mfit, "gevp.amplitude")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa

    denominator <- sqrt(mfit$m0^3)
    if(reduce.latarts) denominator <- sinh(mfit$m0)*sqrt(mfit$m0)
    mfit$fps <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude)/denominator

    denominator <- sqrt(mfit$m0.tsboot^3)
    if(reduce.latarts) denominator <- sinh(mfit$m0.tsboot)*sqrt(mfit$m0.tsboot)
    mfit$fps.tsboot <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude.tsboot[,1])/denominator
    
    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$reduce.latarts <- reduce.latarts
    return(invisible(mfit))
  }
  else{
    stop("computefps: expecting object of class matrixfit, gevp.amplitude or at least an amplitude\nAborting...!\n")
  }
}
