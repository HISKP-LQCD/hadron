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

## version with continuums dispersion relation
compute.qtildesq2 <- function(E, dvec, mpi, L) {
  Ecmsq <- E^2 - sum((2*pi*dvec/L)^2)
  qsq <- Ecmsq/4.-mpi^2
  return(data.frame(gammaboost=E/sqrt(Ecmsq), qtsq=qsq*(L/2./pi)^2, qsq=qsq))
}

## version with lattice dispersion relation
compute.qtildesq <- function(E, dvec, mpi, L) {
  ## cosh(Ecm) = cosh(E) - 2 sum sin^2(Pi/2)
  Ecm <- acosh(cosh(E) - 2*sum(sin(pi*dvec/L)^2))
  ## cosh(Ecm/2) = 2 sin^2(q*/2) + cosh(mpi)
  qs = 2*asin(sqrt( (cosh(Ecm/2) - cosh(mpi))/2. ))
  ## qs = 2 pi q / L
  return(data.frame(gammaboost=E/Ecm, qtsq=(L*qs/2./pi)^2, q=qs))
}

phaseshift.pipi.swave <- function(PC="pc1", tp="TP0", interpolate=TRUE, interpolation.n=100,
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
  
  ## interpolate Zeta function based on first available fits for pi and pipi
  if(interpolate) {
    cat("...computing spline interpolation for Zeta function\n")
    load(pionfilelist[1])
    load(pipifilelist[1])

    ## \tilde q^2    
    qtsq <- compute.qtildesq(E=pc.effectivemass$massfit.tsboot[,1], dvec=dvec, mpi=pion.effectivemass$massfit.tsboot[,1], L=L)
    cat("E:", pc.effectivemass$opt.res$par[1], "mpi:", pion.effectivemass$opt.res$par[1], "gamma:", qtsq$gammaboost[1], "\n")
    
    rg <- range(qtsq$qtsq)
    d <- rg[2]-rg[1]
    rg[1] <- rg[1] + dr[1]*d
    rg[2] <- rg[2] + dr[2]*d
    if(rg[1] < 0) rg[1] <- 0.0001
    ## Zeta function has poles for q^2=0,1,2,3,4...
    ## we don't want to interpolate over a pole
    if(rg[1] > 0 && ceiling(rg[1]) < rg[2]) {
      if(ceiling(rg[1])-rg[1] > rg[2] - ceiling(rg[1])) {
        rg[2] <- ceiling(rg[1])-0.0001
      }
      else {
        rg[1] <- ceiling(rg[1])+0.0001
      }
    }
    else if(rg[1] < 0 && rg[2] > 0) {
      if(abs(rg[1]) > abs(rg[2])) {
        rg[2] <- -0.0001
      }
      else {
        rg[1] <- 0.0001
      }
    }
    cat("range chosen to be from", rg[1], "to", rg[2], "with spacing", d/interpolation.n, "\n")
    x <- seq(rg[1], rg[2], d/interpolation.n)
    y <- Re(LuescherZeta(x, l=l, m=m, gamma=qtsq$gammaboost, dvec=dvec))
    SplineReZ <- splinefun(x,y)
    cat("...done\n")
  }
  cat("\n...determining shifts for", length(pionfilelist), "x", length(pipifilelist), "combinations\n")
  for(i in c(1:length(pionfilelist))) {
    load(pionfilelist[i])
    for(j in c(1:length(pipifilelist))) {
      cat(i, j, pionfilelist[i], pipifilelist[j], "\n")

      load(pipifilelist[j])
      if(pc.effectivemass$boot.R != pion.effectivemass$boot.R) stop("number of boostrap samples does not match. Aborting ...\n")

      qtsq <- compute.qtildesq(pc.effectivemass$opt.res$par[1], dvec=dvec, L=L, mpi=pion.effectivemass$opt.res$par[1])
      qsq <- qtsq$qsq

      Z <- try(LuescherZeta(qtsq$qtsq, l=l, m=m, gamma=qtsq$gammaboost, dvec=dvec))
      if(inherits(Z, "try-error")) Z <- SplineReZ(qtsq$qtsq)
      Z <- Re(Z)
      qcotdelta <- 2.*Z/(qtsq$gammaboost*L*sqrt(pi))
      delta <- atan(qtsq$q/qcotdelta)*180/pi
      res[1, i, j, ] <- c(qtsq$q^2, qtsq$qtsq, qcotdelta, delta, pion.effectivemass$opt.res$par[1], pc.effectivemass$opt.res$par[1],
                          pion.effectivemass$Qval, pc.effectivemass$Qval)
      
      ## now we bootstrap
      qtsq <- compute.qtildesq(pc.effectivemass$massfit.tsboot[,1], dvec=dvec, L=L, mpi=pion.effectivemass$massfit.tsboot[,1])

      if(interpolate) {
        Z <- SplineReZ(qtsq$qtsq)
        if(rg[1] > min(qtsq$qtsq, na.rm=TRUE) && rg[2] < max(qtsq$qtsq, na.rm=TRUE)) {
          ii <- which(qtsq$qtsq < rg[1])
          jj <- which(qtsq$qtsq > rg[2])
          warning(paste("approximation range of spline exceeded", i, j, min(qtsq$qtsq, na.rm=TRUE), max(qtsq$qtsq, na.rm=TRUE), ", using Zeta function itself for", length(c(ii,jj)), "points \n"))
          Z[c(ii,jj)] <- Re(LuescherZeta(qtsq$qtsq[c(ii,jj)], l=l, m=m, gamma=qtsq$gammaboost[c(ii,jj)], dvec=dvec))
        }
      }
      else Z <- Re(LuescherZeta(qtsq$qtsq, l=l, m=m, gamma=qtsq$gammaboost, dvec=dvec))
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
