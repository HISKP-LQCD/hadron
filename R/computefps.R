computefps <- function(mfit, PP, mass, mu1, mu2, Kappa, normalisation="cmi", disprel="continuum", boot.fit=TRUE) {
  if(missing(mu1)) {
    stop("computefps: mu1  must be specified! Aborting...\n")
  }
  if(missing(Kappa)) {
    Kappa <- NA
  }
  if(missing(mu2)) mu2 <- mu1
  if(missing(mfit)) {
    denominator <- sqrt(mass^3)
    if(disprel == "lattice") denominator <- sinh(mass)*sqrt(mass)
    
    if(normalisation != "cmi") return((mu1+mu2)*abs(PP)/denominator)
    else return(sqrt(2)*Kappa*(mu1+mu2)*abs(PP)/denominator)
  }
  else if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa
    
    denominator <- sqrt(mfit$opt.res$par[1]^3)
    if(disprel == "lattice") denominator <- sinh(mfit$opt.res$par[1])*sqrt(mfit$opt.res$par[1])
    mfit$fps <- 2*k*(mu1+mu2)/sqrt(2)*abs(mfit$opt.res$par[2])/denominator

    if(boot.fit) {
      denominator <- sqrt(mfit$opt.tsboot[1,]^3)
      if(disprel == "lattice") denominator <- sinh(mfit$opt.tsboot[1,])*sqrt(mfit$opt.tsboot[1,])
      mfit$fps.tsboot <- sqrt(2.)*k*(mu1+mu2)*abs(mfit$opt.tsboot[2,])/denominator
    } else {
      mfit$fps.tsboot <- NA
    }

    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$disprel <- disprel
    return(invisible(mfit))
  }
  else if(inherits(mfit, "gevp.amplitude")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa

    denominator <- sqrt(mfit$m0^3)
    if(disprel == "lattice") denominator <- sinh(mfit$m0)*sqrt(mfit$m0)
    mfit$fps <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude)/denominator

    if(boot.fit) {
      denominator <- sqrt(mfit$m0.tsboot^3)
      if(disprel == "lattice") denominator <- sinh(mfit$m0.tsboot)*sqrt(mfit$m0.tsboot)
      mfit$fps.tsboot <- sqrt(2)*k*(mu1+mu2)*abs(mfit$meanAmplitude.tsboot[,1])/denominator
    } else {
      mfit$fps.tsboot <- NA
    }
    
    mfit$mu1 <- mu1
    mfit$mu2 <- mu2
    mfit$normalisation <- normalisation
    mfit$disprel <- disprel
    return(invisible(mfit))
  }
  else{
    stop("computefps: expecting object of class matrixfit, gevp.amplitude or at least an amplitude\nAborting...!\n")
  }
}

computefpsOS <- function(mfit, Kappa=sqrt(0.5), normalisation="cmi", boot.fit=TRUE, ZA=1, ZAboot, dZA) {
  if(any(class(mfit) == "matrixfit")) {
    k <- Kappa
    if(normalisation != "cmi") k <- sqrt(0.5)
    else  mfit$kappa <- Kappa

    mfit$ZA <- ZA
    zab <- rep(1, times=mfit$boot.R)
    if(!missing(ZAboot) || !missing(dZA)) {
      if(missing(ZAboot)){
        zab <- rnorm(n=mfit$boot.R, mean=ZA, sd=dZA)
      }
      else {
        zab <- ZAboot
        if(length(ZAboot) != mfit$boot.R) {
          dZA <- sd(ZAboot)
          zab <- rnorm(n=mfit$boot.R, mean=ZA, sd=dZA)
        }
      }
    }
    mfit$ZAboot <- zab
    
    mfit$fpsOS <- ZA*sqrt(2)*k*mfit$opt.res$par[3]*sqrt(mfit$opt.res$par[1])/mfit$opt.res$par[1]    
    if(boot.fit) {
      mfit$fpsOS.tsboot <- zab*sqrt(2)*k*mfit$opt.tsboot[3,]*sqrt(mfit$opt.tsboot[1,])/mfit$opt.tsboot[1,]
    }
    else {
      mfit$fps.tsboot <- NA
    }
    mfit$normalisationOS <- normalisation
    return(invisible(mfit))
  }
  else {
    stop("computefps: expecting object of class matrixfit\nAborting...!\n")
  }
}
