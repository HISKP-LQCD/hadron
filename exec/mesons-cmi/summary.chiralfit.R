print.chiralfit <- function(fit, ...) {
  summary(fit, ...)
}

summary.chiralfit <- function(fit, show.input=FALSE, show.chis=FALSE) {
  N <- length(fit$data)
  npar <- length(fit$par)
  Hinv <- try(solve( fit$fit$hessian ))
  if(inherits(Hinv, "try-error")) {
    Hinv <- diag(nrow=npar, ncol=npar)
  }
  dpar <- numeric(npar)
  for(i in 1:npar) {
    dpar[i] <- sqrt(2*Hinv[i,i])
  }
  if(fit$method != "no") {
    cat("Errors computed using", fit$boot.R, "bootsamples \n")
  }
  cat("FS           = ", fit$fsmethod, "\n")
  cat("chisqr       = ", fit$result$chisqr, "\n")
  cat("dof          = ", fit$result$dof, "\n")
  cat("red. chisqr  = ", fit$result$chisqr/fit$result$dof, "\n")
  cat("Qval         = ", 1-pgamma(fit$result$chisqr/2, fit$result$dof/2), "\n")
  if(fit$method != "no") {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "+-", sd(fit$boots[,(3+3*N)], na.rm=TRUE), "\n")
      cat("l2           = ", fit$result$l2, "+-", sd(fit$boots[,(4+3*N)], na.rm=TRUE), "\n")
    }
    cat("l3           = ", fit$result$l3, "+-", sd(fit$boots[,(1+2*N)], na.rm=TRUE), "\n")
    cat("l4           = ", fit$result$l4, "+-", sd(fit$boots[,(2+2*N)], na.rm=TRUE), "\n")
    cat("F0           = ", fit$result$F, "+-", sd(fit$boots[,(5+2*N)], na.rm=TRUE), "GeV \n")
    cat("Fpi/F0       = ", 1./fit$result$F*0.1307, "+-", sd(1./fit$boots[,(5+2*N)], na.rm=TRUE)*0.1307, "\n")
    cat("B0           = ", fit$result$B0, "+-", sd(fit$boots[,(6+2*N)], na.rm=TRUE), "GeV \n")
    cat("Sigma^(1/3)  = ", fit$result$Sigma, "+-", sd(fit$boots[,(10+2*N)], na.rm=TRUE), "GeV \n")
    cat("<r^2>_s      = ", fit$result$rssq, "+-", sd(fit$boots[,(11+2*N)], na.rm=TRUE), "fm^2 \n")
    if(fit$fit.mN) {
      cat("mN0          = ", fit$result$mN0, "+-", sd(fit$boots[,(7+2*N)], na.rm=TRUE), "GeV \n")
      cat("mN           = ", fit$result$mN, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "GeV \n")
      cat("mN/mN0       = ", fit$result$mN/fit$result$mN0, "+-", sd(fit$boots[,(8+2*N)]/fit$boots[,(7+2*N)], na.rm=TRUE),
          "\n")
      cat("c1           = ", fit$result$c1, "+-", sd(fit$boots[,(9+2*N)], na.rm=TRUE), "GeV^(-1)\n")
      cat("gA           = ", fit$result$gA, "+-", sd(fit$boots[,(13+2*N+7+2*N)], na.rm=TRUE), "\n")
      cat("sigma(0)     = ", fit$result$s0, "+-", sd(fit$boots[,(13+2*N)], na.rm=TRUE), "\n")
    }
    cat("r0           = ", fit$result$r0, "+-", sd(fit$boots[,(12+2*N)], na.rm=TRUE), "fm \n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "+-", sd(fit$boots[,1], na.rm=TRUE), "GeV \n\n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$par[4+i], ": a = ", fit$result$a[i], "+-",
          sd(fit$boots[,(N+i)], na.rm=TRUE),"fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "+-", sd(fit$boots[,(13+2*N+4+i)], na.rm=TRUE), "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "+-", sd(fit$boots[,(13+2*N+4+N+i)], na.rm=TRUE), "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
    cat("Fit Parameter:\n")
    print(data.frame(Par=fit$par, Err=fit$boot.result[(2*N+14):(2*N+13+npar),2],
                     Bias=fit$par- fit$boot.result[(2*N+14):(2*N+13+npar),1], ErrHessian=dpar))
  }
  else {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "\n")
      cat("l2           = ", fit$result$l2, "\n")
    }
    cat("l3           = ", fit$result$l3, "\n")
    cat("l4           = ", fit$result$l4, "\n")
    cat("F0           = ", fit$result$F, "MeV \n")
    cat("Fpi/F0       = ", 1./fit$result$F*0.1307, "\n")
    cat("B0           = ", fit$result$B0, "MeV \n")
    cat("Sigma^(1/3)  = ", fit$result$Sigma, "GeV \n")
    cat("<r^2>_s      = ", fit$result$rssq, "fm^2 \n")
    if(fit$fit.mN) {
      cat("mN0          = ", fit$result$mN0, "MeV \n")
      cat("mN           = ", fit$result$mN, "MeV \n")
      cat("mN/mN0       = ", fit$result$mN/fit$result$mN0, "\n")
      cat("c1           = ", fit$result$c1, "\n")
      cat("gA           = ", fit$result$gA, "\n")
      cat("sigma(0)     = ", fit$result$s0, "\n")
    }
    cat("r0           = ", fit$result$r0, "fm \n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "MeV \n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$par[4+i], ": a = ", fit$result$a[i], "fm \n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
    cat("Fit Parameter:\n")
    print(data.frame(Par=fit$par, ErrHessian=dpar))
  }
  if(show.chis) {
    getchi.Na.withr0ZP(par=fit$par, data=fit$data, ii=fit$ii, r0data=fit$r0data, ZPdata=fit$ZPdata,
                       fsmethod="cdh", a=fit$result$a, fit.l12=fit$fit.l12)
  }
}

