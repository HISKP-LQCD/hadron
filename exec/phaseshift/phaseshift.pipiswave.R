## version with continuums dispersion relation
compute.gamma.free2 <- function(mpisq, p1=c(0,0,0), p2=c(0,0,0), L, dvec) {
  E <- sqrt(mpisq + sum((p1*2*pi/L)^2)) + sqrt(mpisq+sum((p2*2*pi/L)^2))
  ##cat("E non-interacting", sqrt(mpisq + sum((p1*2*pi/L)^2)) + sqrt(mpisq+sum((p2*2*pi/L)^2)), "\n")
  return(E/sqrt(E^2 - sum((2*pi*dvec/L)^2)))
}

## version with lattice dispersion relation
compute.gamma.free <- function(mpisq, p1, p2, L, dvec) {
  cmpi <- cosh(sqrt(mpisq))
  E <- acosh(cmpi + 2*sum(sin(pi*p1/L)^2)) + acosh(cmpi + 2*sum(sin(pi*p2/L)^2))
  Ecm <- acosh(cosh(E) - 2*sum(sin(pi*dvec/L)^2))
  return(E/Ecm)
}

phaseshift.pipi.swave <- function(PC="pc1", tp="TP0",
                                  p1=c(0,0,0), p2=c(0,0,0), dr=c(-1, 1),
                                  boot.R=400, boot.l=1, L=32, T=64, dvec=c(0,0,0), debug=FALSE
                                  ) {
  l <- 0
  m <- 0
  ## dvec needs to be defined as the direction vector for the total momentum
  cat("...determining phase shift from energy levels for", PC, "\n")
  cat("momentum vector:", dvec, "\n")
  pionfilelist <- Sys.glob("pion.p0.effectivemass.*.Rdata")
  pipifilelist <- Sys.glob(paste(PC, ".", tp, ".effectivemass.*.Rdata", sep=""))
  
  res <- array(0., dim=c(boot.R+1, length(pionfilelist), length(pipifilelist), 8))
  
  cat("\n...determining shifts for", length(pionfilelist), "x", length(pipifilelist), "combinations\n")
  for(i in c(1:length(pionfilelist))) {
    load(pionfilelist[i])
    for(j in c(1:length(pipifilelist))) {
      cat(i, j, pionfilelist[i], pipifilelist[j], "\n")

      load(pipifilelist[j])
      if(pc.effectivemass$boot.R != pion.effectivemass$boot.R) stop("number of boostrap samples does not match. Aborting ...\n")

      qtsq <- compute.qtildesq(pc.effectivemass$opt.res$par[1], dvec=dvec, L=L, mpi=pion.effectivemass$opt.res$par[1])

      Z <- try(LuescherZeta(qtsq$qtsq, l=l, m=m, gamma=qtsq$gammaboost, dvec=dvec))
      if(inherits(Z, "try-error")) Z <- SplineReZ(qtsq$qtsq)
      Z <- Re(Z)
      qcotdelta <- 2.*Z/(qtsq$gammaboost*L*sqrt(pi))
      delta <- atan(qtsq$q/qcotdelta)*180/pi
      res[1, i, j, ] <- c(qtsq$q^2, qtsq$qtsq, qcotdelta, delta, pion.effectivemass$opt.res$par[1], pc.effectivemass$opt.res$par[1],
                          pion.effectivemass$Qval, pc.effectivemass$Qval)
      
      ## now we bootstrap
      qtsq <- compute.qtildesq(pc.effectivemass$massfit.tsboot[,1], dvec=dvec, L=L, mpi=pion.effectivemass$massfit.tsboot[,1])

      Z <- Re(LuescherZeta(qtsq$qtsq, l=l, m=m, gamma=qtsq$gammaboost, dvec=dvec))
      res[c(2:(boot.R+1)), i, j, 1] <- qtsq$q^2
      res[c(2:(boot.R+1)), i, j, 2] <- qtsq$qtsq
      ## q cot(delta)
      res[c(2:(boot.R+1)), i, j, 3] <- 2.*Z/(qtsq$gammaboost*L*sqrt(pi))
      ## delta
      res[c(2:(boot.R+1)), i, j, 4] <- atan(qtsq$q/res[c(2:(boot.R+1)), i, j, 3])*180/pi
      ## Epi
      res[c(2:(boot.R+1)), i, j, 5] <- pion.effectivemass$massfit.tsboot[,1]
      ## Epipi
      res[c(2:(boot.R+1)), i, j, 6] <- pc.effectivemass$massfit.tsboot[,1]
      ## p-value pion effective mass fit
      res[c(2:(boot.R+1)), i, j, 7] <- 1-pchisq(pion.effectivemass$massfit.tsboot[,2], pion.effectivemass$dof)
      ## p-value pipi effective mass fit
      res[c(2:(boot.R+1)), i, j, 8] <- 1-pchisq(pc.effectivemass$massfit.tsboot[,2], pc.effectivemass$dof)
      if(debug) cat("qcotdelta = ", qcotdelta, "; delta = ", delta, "\n")
    }
    save(res, file=paste("res.", PC, ".", tp, ".Rdata", sep=""))
  }
  return(invisible(res))
}
