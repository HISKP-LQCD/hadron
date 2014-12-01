cat("...determining phase shift from energy levels for", pc, "\n")
pionfilelist <- Sys.glob("pion.p0.effectivemass.*.Rdata")
pipifilelist <- Sys.glob(paste(pc, ".effectivemass.*.Rdata", sep=""))

res <- array(0., dim=c(boot.R+1, length(pionfilelist), length(pipifilelist), 6))

## interpolate Zeta function based on first available fits for pi and pipi
if(interpolate) {
  cat("...computing spline interpolation for Zeta function\n")
  load(pionfilelist[1])
  load(pipifilelist[1])
  
  qsq <- pc.effectivemass$massfit.tsboot[,1]^2/4. - pion.effectivemass$massfit.tsboot[,1]^2
  qtsq <- L^2/4/pi^2*qsq
  rg <- range(qtsq)
  d <- rg[2]-rg[1]
  rg[1] <- rg[1]-3*d/2.
  rg[2] <- rg[2]+3*d/2.
  cat("range chosen to be from", rg[1], "to", rg[2], "with spacing", d/200, "\n")
  x <- seq(rg[1], rg[2], d/200)
  y <- Re(LuescherZeta(x, l=l, m=m, gamma=gammaboost, dvec=dvec))
  SplineReZ <- splinefun(x,y)
  cat("...done\n")
}

for(i in c(1:length(pionfilelist))) {
  load(pionfilelist[i])
  for(j in c(1:length(pipifilelist))) {
    cat(i, j, "\n")
    if(debug) cat(i, j, pionfilelist[i], pipifilelist[j], "\n")
    load(pipifilelist[j])
    if(pc.effectivemass$boot.R != pion.effectivemass$boot.R) stop("number of boostrap samples does not match. Aborting ...\n")
    ## q^2
    qsq <- pc.effectivemass$opt.res$par[1]^2/4. - pion.effectivemass$opt.res$par[1]^2
    ## \tilde q^2
    qtsq <- L^2/4/pi^2*qsq
    
    Z <- Re(LuescherZeta(qtsq, l=l, m=m, gamma=gammaboost, dvec=dvec))
    qcotdelta <- 2.*Z/(gammaboost*L*sqrt(pi))
    delta <- atan2(gammaboost*pi^(3/2)*sqrt(qtsq), Z)*180/pi
    res[1, i, j, ] <- c(qsq, qtsq, qcotdelta, delta, pion.effectivemass$opt.res$par[1], pc.effectivemass$opt.res$par[1])
    
    qsq <- pc.effectivemass$massfit.tsboot[,1]^2/4. - pion.effectivemass$massfit.tsboot[,1]^2
    qtsq <- L^2/4/pi^2*qsq
    if(interpolate) {
      if(rg[1] <= min(qtsq) && rg[2] >= max(qtsq))
        Z <- SplineReZ(qtsq)
      else {
        warning(paste("approximation range of spline exceeded, using Zeta function itself", i, j, min(qtsq), max(qtsq), "\n"))
        Z <- Re(LuescherZeta(qtsq, l=l, m=m, gamma=gammaboost, dvec=dvec))
      }
    }
    else Z <- Re(LuescherZeta(qtsq, l=l, m=m, gamma=gammaboost, dvec=dvec))
    qcotdelta <- 2.*Z/(gammaboost*L*sqrt(pi))
    delta <- atan2(gammaboost*pi^(3/2)*sqrt(qtsq), Z)*180/pi
    res[c(2:(boot.R+1)), i, j, 1] <- qsq
    res[c(2:(boot.R+1)), i, j, 2] <- qtsq
    res[c(2:(boot.R+1)), i, j, 3] <- qcotdelta
    res[c(2:(boot.R+1)), i, j, 4] <- delta
    res[c(2:(boot.R+1)), i, j, 5] <- pion.effectivemass$massfit.tsboot[,1]
    res[c(2:(boot.R+1)), i, j, 6] <- pc.effectivemass$massfit.tsboot[,1]
    
    if(debug) cat(qcotdelta, "and", delta, "\n")
  }
  save(res, file=paste("res.",pc,".CM.Rdata", sep=""))
}
