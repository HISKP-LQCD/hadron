print.pionChiPTfit <- function(fit, ...) {
  summary(fit, ...)
}

summary.pionChiPTfit <- function(fit, show.input=FALSE, show.chis=FALSE) {
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
    cat("l3           = ", fit$result$l3, "+-", sd(fit$boots[,(1+3*N)], na.rm=TRUE), "\n")
    cat("l4           = ", fit$result$l4, "+-", sd(fit$boots[,(2+3*N)], na.rm=TRUE), "\n")
    cat("F0           = ", fit$result$F, "+-", sd(fit$boots[,(5+3*N)], na.rm=TRUE), "GeV \n")
    cat("Fpi/F0       = ", 1./fit$result$F*0.1307, "+-", sd(1./fit$boots[,(5+3*N)], na.rm=TRUE)*0.1307, "\n")
    cat("B0           = ", fit$result$B0, "+-", sd(fit$boots[,(6+3*N)], na.rm=TRUE), "GeV \n")
    cat("Sigma^(1/3)  = ", fit$result$Sigma, "+-", sd( ((fit$boots[,(6+3*N)]*fit$boots[,(5+3*N)]^2)/2)^(1/3), na.rm=TRUE), "GeV \n")
    cat("<r^2>_s      = ", fit$result$rssq, "+-", sd(fit$boots[,(9+3*N)], na.rm=TRUE), "fm^2 \n")
    if(fit$fit.mN) {
      cat("mN0          = ", fit$result$mN0, "+-", sd(fit$boots[,(7+2*N)], na.rm=TRUE), "GeV \n")
      cat("mN           = ", fit$result$mN, "+-", sd(fit$boots[,(8+2*N)], na.rm=TRUE), "GeV \n")
      cat("mN/mN0       = ", fit$result$mN/fit$result$mN0, "+-", sd(fit$boots[,(8+2*N)]/fit$boots[,(7+2*N)], na.rm=TRUE),
          "\n")
      cat("c1           = ", fit$result$c1, "+-", sd(fit$boots[,(9+2*N)], na.rm=TRUE), "GeV^(-1)\n")
      cat("gA           = ", fit$result$gA, "+-", sd(fit$boots[,(13+2*N+7+2*N)], na.rm=TRUE), "\n")
      cat("sigma(0)     = ", fit$result$s0, "+-", sd(fit$boots[,(13+2*N)], na.rm=TRUE), "\n")
    }
    cat("r0           = ", fit$result$r0, "+-", sd(fit$boots[,(7+3*N)], na.rm=TRUE), "fm \n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "+-", sd(fit$boots[,1], na.rm=TRUE), "GeV \n\n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$par[4+i], ": a = ", fit$result$a[i], "+-",
          sd(fit$boots[,(N+i)], na.rm=TRUE),"fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "+-", sd(fit$boots[,(9+3*N+4+i)], na.rm=TRUE), "\n")
      cat("            fitted ZP   = ", fit$par[6+N+i], "+-", sd(fit$boots[,(9+3*N+6+N+i)], na.rm=TRUE), "\n")
      if(show.input) {
        cat("Raw data used:\n")
        print(fit$data[[i]])
        cat("Raw ZP data used:\n")
        cat(fit$ZPdata$ZP[i], fit$ZPdata$dZP[i], "\n")
      }
      cat("\n")
    }
    cat("Fit Parameter:\n")
    print(data.frame(Par=fit$par, Err=fit$boot.result[(3*N+10):(3*N+9+npar),2], Bias=fit$par- fit$boot.result[(3*N+10):(3*N+9+npar),1], ErrHessian=dpar))
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
      cat("            fitted ZP   = ", fit$par[6+N+i], "\n")
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

printtab.old <- function(x, dx, digits=2, c=1.) {
  if(!is.null(dx)) {
    x <- c*x
    dx <- c*dx
    n <- 1.
    while(n*dx < 10^(digits-1)) {
      n=n*10.
    }
    op <- options(scipen=6)
    cat("$", round(n*x)/n, "(", signif(n*dx, digits=digits), ")$\n", sep="")
    options(op)
  }
  else {
    cat("$", c*x,"()$\n", sep="")
  }
}

printtab <- function(x, dx, digits=2, c=1., endl=TRUE) {
  if(!is.null(dx)) {

    x <- c*x
    dx <- c*dx
    n <- 1.
    i <- 0
    while(n*dx < 10^(digits-1)) {
      i <- i+1
      n=n*10.
    }
    y <- round(n*abs(x))/n
    z <- abs(round((abs(x)-floor(abs(x))), digits=i))
    if(x < 0.) cat("$-", floor(y), sep="")
    else cat("$", floor(y), sep="")
    if(n > 1.) cat(".", sep="")
    while(i > 0) {
      z <- z*10.000001
      cat(floor(z))
      z <- z-floor(z)
      i <- i-1
    }
    op <- options(scipen=6)
    if(floor(dx) > 0 && n > 1.) cat("(", signif(dx, digits=digits), ")$", sep="")
    else cat("(", signif(n*dx, digits=digits), ")$", sep="")
    options(op)
  }
  else {
    cat("$", c*x,"()$", sep="")
  }
  if(endl) cat(" \\\\\n", sep="")
}


tab <- function(fitlist) {
  nl <- length(fitlist)
  
  par <- fitlist[[1]]$par
  br <- fitlist[[1]]$boot.result
  npar <- length(par)
  N <- length(fitlist[[1]]$data)
  #cat("FS    ", sep="")
  #for(k in 1:nl) {
  #  cat(" & ", sep="")
  #  cat(fitlist[[1]]$fsmethod)
  #}
  #cat(" \\\\\n", sep="")
  cat("$2r_0B_0$   ", sep="")
  nm <- 3*N+9

  ## 2 r0 B0
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$par[4], fitlist[[k]]$boot.result[nm+4,2], endl=FALSE)
  }
  cat(" \\\\\n$r_0f_0$   ", sep="")
  ## r0 f0
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$par[3], fitlist[[k]]$boot.result[nm+3,2], endl=FALSE)
  }
  cat(" \\\\\n", sep="")
  
  if(fitlist[[1]]$fit.mN) {
    ## r0 mN
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$par[2*N+5], fitlist[[k]]$boot.result[nm+2*N+5,2], endl=FALSE)
    }
    cat(" \\\\\n", sep="")
  }
  # r0
  for(i in 1:N) {
    cat("$r_0/a$   ", sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$par[4+i], fitlist[[k]]$boot.result[nm+4+i,2], endl=(k==nl))
    }
  }
  # sr0 (slope)
  for(i in 1:2) {
    cat("$sr_1}$   ", sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$par[4+N+i], fitlist[[k]]$boot.result[nm+4+N+i,2], endl=(k==nl))
    }
  }
  # ZP
  for(i in 1:N) {
    cat("$Z_\\mathrm{P}$   ", sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$par[6+N+i], fitlist[[k]]$boot.result[nm+6+N+i,2], endl=FALSE)
    }
    cat(" \\\\\n", sep="")
  }
  # l1 l2
  for(i in 1:2) {
    cat(c("$r_0\\Lambda_1$   ","$r_0\\Lambda_2$   ")[i], sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      if(fitlist[[k]]$fit.l12) printtab(fitlist[[k]]$par[6+2*N+i], fitlist[[k]]$boot.result[nm+6+2*N+i,2], endl=FALSE)
      else cat("$-$")
    }
    if(k==nl)cat(" \\\\\n", sep="")
  }
  # l3 l4
  for(i in 1:2) {
    cat(c("$r_0\\Lambda_3$   ","$r_0\\Lambda_4$   ")[i], sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(abs(fitlist[[k]]$par[i]), fitlist[[k]]$boot.result[nm+i,2], endl=(k==nl))
    }
  }
  # kM, kF
  for(i in 1:2) {
    cat(c("$k_M$   ","$k_F$   ")[i], sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      if(fitlist[[k]]$fit.kmf) printtab(fitlist[[k]]$par[7+2*N+i], fitlist[[k]]$boot.result[nm+7+2*N+i,2], endl=FALSE)
      else cat("$-$")
      if(k==nl) cat(" \\\\\n", sep="")
    }
  }
#  for(i in 1:2) {
#    if(fitlist[[k]]$fit.mN) printtab(fitlist[[k]]$par[5+2*N+i], fitlist[[k]]$boot.result[nm+5+2*N+i,2])
#  }
  # Dm, Df
  for(i in 1:0) {
    cat(c("$D_f$   ","$D_m$   ")[i+1], sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      if(fitlist[[k]]$fit.asq) {
        printtab(fitlist[[k]]$par[npar-i], fitlist[[k]]$boot.result[nm+npar-i,2], endl=(k==nl))
      }
      else {
        cat("$-$")
      }
    }
  }
  cat("$\\chi^2/\\mathrm{dof}$   ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    cat("$",round(fitlist[[k]]$result$chisqr, digits=2),"/",fitlist[[k]]$result$dof,"$", sep="")
  }
  cat(" \\\\\nCL   ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    cat("$",round(1-pchisq(fitlist[[k]]$result$chisqr, fitlist[[k]]$result$dof), digits=4), "$", sep="")
  }
  cat(" \\\\\n", sep="")
  cat("\\hline\n$m_{u,d}\\ [\\mathrm{MeV}]$   ")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$mu.phys, fitlist[[k]]$boot.result[1,2], c=1000., endl=(k==nl))
  }
  for(i in 1:N) {
    cat("$a(xxx)\\ [\\mathrm{fm}]$   ", sep="")
    for(k in 1:nl) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$result$a[i], fitlist[[k]]$boot.result[N+i,2], endl=(k==nl))
    }
  }
  cat("$r_0\\ [\\mathrm{fm}]$   ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$r0, fitlist[[k]]$boot.result[7+3*N,2], endl=(k==nl))
  }
  cat("$\\bar{\\ell}_1$   ", sep="")
  for(k in 1:nl) {
    if(fitlist[[k]]$fit.l12) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$result$l1, fitlist[[k]]$boot.result[3+3*N, 2], endl=(k==nl))
    }
    else cat("$-$", sep="")
  }
  cat("$\\bar{\\ell}_2$   ", sep="")
  for(k in 1:nl) {
    if(fitlist[[k]]$fit.l12) {
      cat(" & ", sep="")
      printtab(fitlist[[k]]$result$l2, fitlist[[k]]$boot.result[4+3*N, 2], endl=(k==nl))
    }
    else cat("$-$", sep="")
  }
  cat("$\\bar{\\ell}_3$   ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$l3, fitlist[[k]]$boot.result[1+3*N, 2], endl=(k==nl))
  }
  cat("$\\bar{\\ell}_4$   ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$l4, fitlist[[k]]$boot.result[2+3*N, 2], endl=(k==nl))
  }
  cat("$f_0\\ [\\mathrm{MeV}]$    ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$F, fitlist[[k]]$boot.result[5+3*N, 2], c=1000., endl=(k==nl))
  }
                                        #  try(stddev=sd(1./fitlist[[k]]$boots[,(5+3*N)])*0.1307)
                                        #  if(inherits(stddev, "try-error")) stddev=NA
                                        #  printtab(1./fitlist[[k]]$result$F*0.1307, 0.001)
  cat("$f_\\pi/f_0$    ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(1./fitlist[[k]]$result$F*0.1307, sd(1./fitlist[[k]]$boots[,(5+3*N)]*0.1307, na.rm=TRUE), endl=(k==nl))
  }
                                        #  cat("-----\n")
  cat("$B_0\\ [\\mathrm{MeV}]$    ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$B0, fitlist[[k]]$boot.result[6+3*N, 2], c=1000., endl=(k==nl))
  }
  cat("$\\Sigma^{1/3}\\ [\\mathrm{MeV}]$    ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$Sigma, sd( ((fitlist[[k]]$boots[,(6+3*N)]*fitlist[[k]]$boots[,(5+3*N)]^2)/2)^(1/3), na.rm=TRUE), c=1000., endl=(k==nl))
  }
  cat("$\\langle r^2\\rangle_s\\ [\\mathrm{fm}^2]$    ", sep="")
  for(k in 1:nl) {
    cat(" & ", sep="")
    printtab(fitlist[[k]]$result$rssq, fitlist[[k]]$boot.result[9+3*N, 2], endl=(k==nl))
  }
  if(fitlist[[k]]$fit.mN) {
    printtab(fitlist[[k]]$result$mN0, fitlist[[k]]$boot.result[7+2*N, 2], c=1000.)
    printtab(fitlist[[k]]$result$c1, fitlist[[k]]$boot.result[9+2*N, 2])
    printtab(fitlist[[k]]$result$s0, fitlist[[k]]$boot.result[13+2*N, 2], c=1000)
    printtab(fitlist[[k]]$result$mN, fitlist[[k]]$boot.result[8+2*N, 2], c=1000)
    printtab(fitlist[[k]]$result$mN/fitlist[[k]]$result$mN0, sd(fitlist[[k]]$boots[,(8+2*N)]/fitlist[[k]]$boots[,(7+2*N)], na.rm=TRUE))
  }
}
