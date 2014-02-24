computefps <- function(mfit, PP, mass, mu1, mu2, Kappa) {
  if(missing(mfit)) {
    return(4*Kappa*(mu1+mu2)/2/sqrt(2)*abs(PP)/sqrt(mass^3))
  }
  else{
    return(c(4*Kappa*(mu1+mu2)/2/sqrt(2)*abs(mfit$opt.res$par[2])/sqrt(mfit$opt.res$par[1]^3),
             sd(4*Kappa*(mu1+mu2)/2/sqrt(2)*abs(mfit$opt.tsboot[,2])/sqrt(mfit$opt.tsboot[,1]^3)))
           )
  }
}
