computefps <- function(mfit, PP, mass, mu1, mu2, Kappa, normalisation="cmi") {
  if(missing(mu2)) mu2 <- mu1
  if(missing(mfit)) {
    if(normalisation != "cmi") return((mu1+mu2)/sqrt(2)*abs(PP)/sqrt(mass^3))
    else return(2*Kappa*(mu1+mu2)/sqrt(2)*abs(PP)/sqrt(mass^3))
  }
  else if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- 0.5
    else  mfit$kappa <- Kappa
    mfit$fps <- 2*k*(mu1+mu2)/sqrt(2)*abs(mfit$opt.res$par[2])/sqrt(mfit$opt.res$par[1]^3)
    mfit$fps.tsboot <- 2*k*(mu1+mu2)/sqrt(2)*abs(mfit$opt.tsboot[2,])/sqrt(mfit$opt.tsboot[1,]^3)
    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    return(mfit)
  }
  else{
    stop("computefps expecting object of clas matrixfit\nAborting...!\n")
  }
}
