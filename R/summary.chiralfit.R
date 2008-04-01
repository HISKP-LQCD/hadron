plot.chiralfit <- function(fit, ...) {
  N <- length(fit$data)
  npar <- length(fit$par)
  par <- fit$par
  X11()
  fplot <- dev.cur()
  X11()
  mpsplot <- dev.cur()
  X11()
  mpsmuplot <- dev.cur()
  X11()
  mNplot <- dev.cur()
  
  for(i in 1:N) {
    dev.set(fplot)
    ij <- ii[[i]]
    if(fit$fit.l12) {
      aLamb1=par[8+2*N]/par[4+i]
      aLamb2=par[9+2*N]/par[4+i]
    }
    else {
                                        #        aLamb1=sqrt(exp(-0.4+log((0.1396*a_fm/0.1973)^2)))
      aLamb1=sqrt(exp(-0.4)*(0.1396*fit$result$a[i]/0.1973)^2)
                                        #        aLamb2=sqrt(exp(4.3+log((0.1396*a_fm/0.1973)^2)))
      aLamb2=sqrt(exp(4.3)*(0.1396*fit$result$a[i]/0.1973)^2)
    }
    
    res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
               aLamb4=par[2]/par[4+i], ampiV=fit$data[[i]]$mps[ij], afpiV=fit$data[[i]]$fps[ij],
               aF0=fit$data[[i]]$fps[ij], a_fm=fit$result$a[i], L=fit$data[[i]]$L[ij], rev=-1, printit=FALSE)

    mpsV <- res$mpiFV

    fpsV <- res$fpiFV
    xfit <- seq(from=0., to=1.05*max(fit$data[[i]]$mu*fit$par[4+i]/fit$par[4+N+i], na.rm=TRUE),
                length.out=500)

    r0TwoB <- par[4]
    r0sqTwoBmu <- r0TwoB*xfit
    msq <- getmpssq(r0sqTwoBmu, fit$par, N, fit$fit.l12)
    f <- getfps(r0sqTwoBmu, fit$par, N, fit$fit.l12)
    mN <- getmN(r0sqTwoBmu, fit$par, N)
    color <- c("red", "blue", "black")
    if(i == 1) {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    fit$r0data$r0[i]*fpsV, fit$r0data$r0[i]*fit$data[[i]]$dfps[ij],
                    ylim=c(0.85*min(fit$r0data$r0[i]*fpsV, na.rm=TRUE), 1.1*max(fit$r0data$r0[i]*fpsV, na.rm=TRUE)),
                    xlim=c(0.,1.1*max(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu, na.rm=TRUE)), col=color[i],
                    pch=i, ylab="r0 fps", xlab="r0 mu / ZP")
      lines(spline(xfit, f), lty=2)
    }
    else{
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    fit$r0data$r0[i]*fpsV, fit$r0data$r0[i]*fit$data[[i]]$dfps[ij], rep=TRUE,
                    col=color[i], pch=i)
    }
    dev.set(mpsplot)
    if(i == 1) {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    (fit$r0data$r0[i]*mpsV)^2, 2*fit$r0data$r0[i]*fit$r0data$r0[i]*mpsV*fit$data[[i]]$dmps[ij],
                    ylim=c(0., 1.1*max((fit$r0data$r0[i]*mpsV)^2, na.rm=TRUE)),
                    xlim=c(0.,1.1*max(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij], na.rm=TRUE)), col=color[i],
                    pch=i, ylab="(r0 mps)^2", xlab="r0 mu / ZP")
      lines(spline(xfit, msq), lty=2)
    }
    else {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    (fit$r0data$r0[i]*mpsV)^2, 2*fit$r0data$r0[i]*fit$r0data$r0[i]*mpsV*fit$data[[i]]$dmps[ij],
                    rep=TRUE, col=color[i], pch=i)      
    }
    dev.set(mpsmuplot)
    if(i == 1) {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    (fit$r0data$r0[i]*mpsV)^2/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                    2*fit$r0data$r0[i]*fit$r0data$r0[i]*mpsV*fit$data[[i]]$dmps[ij]/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                    ylim=c(0.85*min((fit$r0data$r0[i]*mpsV)^2/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                      na.rm=TRUE),
                      1.1*max((fit$r0data$r0[i]*mpsV)^2/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                              na.rm=TRUE)),
                    xlim=c(0.,1.1*max(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij], na.rm=TRUE)), col=color[i],
                    pch=i, ylab="(r0 mps)^2/mu", xlab="r0 mu / ZP")
      lines(spline(xfit, msq/xfit), lty=2)
    }
    else {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    (fit$r0data$r0[i]*mpsV)^2/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                    2*fit$r0data$r0[i]*fit$r0data$r0[i]*mpsV*fit$data[[i]]$dmps[ij]/(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij]),
                    rep=TRUE, col=color[i], pch=i)      
    }
    dev.set(mNplot)
    if(i==1) {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    fit$r0data$r0[i]*fit$data[[i]]$mN[ij], fit$r0data$r0[i]*fit$data[[i]]$dmN[ij],
                    ylim=c(0.85*min(fit$r0data$r0[i]*fit$data[[i]]$mN[ij], na.rm=TRUE),
                      1.1*max(fit$r0data$r0[i]*fit$data[[i]]$mN[ij], na.rm=TRUE)),
                    xlim=c(0.,1.1*max(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij], na.rm=TRUE)), col=color[i],
                    pch=i, ylab="r0 mN", xlab="r0 mu / ZP")
      lines(spline(xfit, mN), lty=2)
    }
    else {
      plotwitherror(fit$r0data$r0[i]/fit$ZPdata$ZP[i]*fit$data[[i]]$mu[ij],
                    fit$r0data$r0[i]*fit$data[[i]]$mN[ij], fit$r0data$r0[i]*fit$data[[i]]$dmN[ij],
                    rep=TRUE, col=color[i], pch=i)
    }
  }
}

print.chiralfit <- function(fit, ...) {
  summary(fit, ...)
}

summary.chiralfit <- function(fit, show.input=FALSE, show.chis=FALSE) {
  N <- length(fit$data)
  npar <- length(fit$par)
  if(!is.null(fit$boot.result)) {
    cat("Errors computed using", fit$boot.R, "bootsamples \n")
  }
  cat("chisqr       = ", fit$result$chisqr, "\n")
  cat("dof          = ", fit$result$dof, "\n")
  cat("red. chisqr  = ", fit$result$chisqr/fit$result$dof, "\n")
  if(!is.null(fit$boot.result)) {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "\n")
      cat("l2           = ", fit$result$l2, "+-", sd(fit$boots[,(9+2*N)], na.rm=TRUE), "\n")
    }
    cat("l3           = ", fit$result$l3, "+-", sd(fit$boots[,(1+2*N)], na.rm=TRUE), "\n")
    cat("l4           = ", fit$result$l4, "+-", sd(fit$boots[,(2+2*N)], na.rm=TRUE), "\n")
    cat("F0           = ", fit$result$F, "+-", sd(fit$boots[,(5+2*N)], na.rm=TRUE), "MeV \n")
    cat("B0           = ", fit$result$B0, "+-", sd(fit$boots[,(6+2*N)], na.rm=TRUE), "MeV \n")
    cat("mN0          = ", fit$result$mN0, "+-", sd(fit$boots[,(7+2*N)], na.rm=TRUE), "MeV \n")
    cat("mN           = ", fit$result$mN, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "MeV \n")
    cat("c1           = ", fit$result$c1, "+-", sd(fit$boots[,(8+2*N+6+2*N)], na.rm=TRUE), "\n")
    cat("gA           = ", fit$result$gA, "+-", sd(fit$boots[,(8+2*N+7+2*N)], na.rm=TRUE), "\n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "+-", sd(fit$boots[,1], na.rm=TRUE), "MeV \n\n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$r0data$r0[i], ": a = ", fit$result$a[i], "+-",
          sd(fit$boots[,(N+i)], na.rm=TRUE),"fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw r0 data used:\n")
        cat(fit$r0data$r0[i], fit$r0data$dr0[i], "\n")
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
  }
  else {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "\n")
      cat("l2           = ", fit$result$l2, "\n")
    }
    cat("l3           = ", fit$result$l3, "\n")
    cat("l4           = ", fit$result$l4, "\n")
    cat("F0           = ", fit$result$F, "MeV \n")
    cat("B0           = ", fit$result$B0, "MeV \n")
    cat("mN0          = ", fit$result$mN0, "MeV \n")
    cat("mN           = ", fit$result$mN, "MeV \n")
    cat("c1           = ", fit$result$c1, "\n")
    cat("gA           = ", fit$result$gA, "\n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "MeV \n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$r0data$r0[i], ": a = ", fit$result$a[i], "fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw r0 data used:\n")
        cat(fit$r0data$r0[i], fit$r0data$dr0[i], "\n")
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
  }
  if(show.chis) {
    getchi.Na.withr0ZP(par=fit$par, data=fit$data, ii=fit$ii, r0data=fit$r0data, ZPdata=fit$ZPdata,
                       fsmethod="cdh", a=fit$result$a, fit.l12=fit$fit.l12)
  }
}
