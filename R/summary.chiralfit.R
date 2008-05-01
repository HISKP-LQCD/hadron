plot.chiralfit <- function(fit, write.data=FALSE, plot.file=FALSE, plot.all=FALSE,...) {
  N <- length(fit$data)
  npar <- length(fit$par)
  fit.a <- -1.
  par <- fit$par
  X11()
  par(list(las=1, tck=0.01, mgp=c(3,1.,0)))
  fplot <- dev.cur()
  X11()
  par(list(las=1, tck=0.01, mgp=c(3,1.,0)))
  mpsplot <- dev.cur()
  X11()
  par(list(las=1, tck=0.01, mgp=c(3,1.,0)))
  mpsmuplot <- dev.cur()
  if(fit$fit.mN) {
    X11()
    par(list(las=1, tck=0.01, mgp=c(3,1.,0)))
    mNplot <- dev.cur()
  }

  
  for(i in 1:N) {
    if(fit$fit.asq) {
      fit.a <- i
    }
    dev.set(fplot)
    ij <- fit$ii[[i]]
    if(plot.all) {
      ij <- c(1:length(fit$data[[i]]$mu))
    }
    if(fit$fsmethod == "cdh") {
      if(fit$fit.l12) {
        aLamb1=par[8+2*N]/par[4+i]
        aLamb2=par[9+2*N]/par[4+i]
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*fit$result$a[i]/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*fit$result$a[i]/0.1973)^2)
      }
      
      
      res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                 aLamb4=par[2]/par[4+i], ampiV=fit$data[[i]]$mps, afpiV=fit$data[[i]]$fps,
                 aF0=fit$data[[i]]$fps, a_fm=fit$result$a[i], L=fit$data[[i]]$L, rev=-1, printit=FALSE)
      
      mpsV <- res$mpiFV
      
      fpsV <- res$fpiFV
    }
    else {
      mpsv <- fit$data[[i]]$L[ij]*fit$data[[i]]$mps/par[4+i]
      r <-  - fit$data[[i]]$mps^2/(4.0*pi*par[3])^2*g1(mpsv)/par[4+i]
      
      mpsV <- fit$data[[i]]$mps*(1.+0.5*r)/par[4+i]
      fpsV <- fit$data[[i]]$fps*(1.0-2.0*r)/par[4+i]
    }
    
    xfit <- seq(from=0., to=1.05*max(fit$data[[i]]$mu[ij]*fit$par[4+i]/fit$par[4+N+i], na.rm=TRUE),
                length.out=500)
    
    r0TwoB <- par[4]
    r0sqTwoBmu <- r0TwoB*xfit
    msq <- getmpssq(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=fit.a)
    f <- getfps(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=fit.a)
    mN <- getmN(r0sqTwoBmu, fit$par, N, fit.asq=fit.a)
    color <- c("red", "blue", "darkseagreen", "darksalmon", "darkviolet")
    xmu <- fit$par[4+i]/fit$par[4+N+i]*fit$data[[i]]$mu
                                        # continuum curves
    if(fit$fit.asq && i == 1) {
      f2 <- getfps(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=-1.)
      mN2 <- getmN(r0sqTwoBmu, fit$par, N, fit.asq=-1.)
      msq2 <- getmpssq(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=-1.)
    }
    if(write.data) {
      file <- paste("curve_",i,".dat", sep="")
      write.table(data.frame(xfit, msq, f, mN),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
      file <- paste("data_cor_", fit$fsmethod, i,".dat", sep="")
      write.table(data.frame(fit$data[[i]]$mu, fit$data[[i]]$L,
                             fit$par[4+i], fit$par[4+N+i],
                             xmu, fpsV, fit$data[[i]]$dfps,
                             mpsV, fit$data[[i]]$dmps),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
      file <- paste("data_", i,".dat", sep="")
      write.table(data.frame(fit$data[[i]]$mu, fit$data[[i]]$L,
                             fit$par[4+i], fit$par[4+N+i],
                             fit$data[[i]]$fps, fit$data[[i]]$dfps,
                             fit$data[[i]]$mps, fit$data[[i]]$dmps),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
    }
    if(i == 1) {
      plotwitherror(xmu[ij],
                    fit$par[4+i]*fpsV[ij], fit$par[4+i]*fit$data[[i]]$dfps[ij],
                    ylim=c(0.9*par[3], 1.1*max(fit$par[4+i]*fpsV[ij], na.rm=TRUE)),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20+i, ylab=expression(r[0]*f[PS]), xlab=expression(r[0]*mu[R]))
      lines(xfit, f, lty=2)
      points(fit$result$r0*fit$result$mu.phys/0.1973, fit$result$r0*0.1307/0.1973, pch=25, bg="turquoise")
      if(!is.null(fit$boot.result)) {
        arrows((fit$result$r0*fit$result$mu.phys-sd(fit$boots[,1], na.rm=TRUE))/0.1973,
               fit$result$r0*0.1307/0.1973,
               (fit$result$r0*fit$result$mu.phys+sd(fit$boots[,1], na.rm=TRUE))/0.1973,
               fit$result$r0*0.1307/0.1973,
               length=0.01,angle=0,code=3)
      }
    }
    else {
      plotwitherror(xmu[ij],
                    fit$par[4+i]*fpsV[ij], fit$par[4+i]*fit$data[[i]]$dfps[ij], rep=TRUE,
                    col=color[i], pch=20+i, bg=color[i])
      if(fit$fit.asq) {
        lines(xfit, f, lty=2)
      }
    }
    if(fit$fit.asq && i==1) {
      lines(xfit, f2, lty=1)
      if(write.data) {
        file <- paste("cont_curve_",i,".dat", sep="")
        write.table(data.frame(xfit, msq2, f2, mN2),
                    row.names=FALSE, col.names=FALSE,
                    file=file)
      }
    }
    dev.set(mpsplot)
    if(i == 1) {
      plotwitherror(xmu[ij],
                    (fit$par[4+i]*mpsV[ij])^2, 2*fit$par[4+i]^2*mpsV[ij]*fit$data[[i]]$dmps[ij],
                    ylim=c(0., 1.1*max((fit$par[4+i]*mpsV[ij])^2, na.rm=TRUE)),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20+i, ylab=expression((r[0]*m[PS])^2), xlab=expression(r[0]*mu[R]))
      lines(xfit, msq, lty=2)
    }
    else {
      plotwitherror(xmu[ij],
                    (fit$par[4+i]*mpsV[ij])^2, 2*fit$par[4+i]^2*mpsV[ij]*fit$data[[i]]$dmps[ij],
                    rep=TRUE, col=color[i], pch=20+i, bg=color[i])
      if(fit$fit.asq) {
        lines(xfit, msq, lty=2)
      }
    }
    if(fit$fit.asq && i==1) {
      lines(xfit, msq2, lty=1)
    }
    dev.set(mpsmuplot)
    if(i == 1) {
      plotwitherror(xmu[ij],
                    (fit$par[4+i]*mpsV[ij])^2/(xmu[ij]),
                    2*fit$par[4+i]^2*mpsV[ij]*fit$data[[i]]$dmps[ij]/(xmu[ij]),
                    ylim=c(0.95*min((fit$par[4+i]*mpsV[ij])^2/(xmu[ij]),
                      na.rm=TRUE),
                      1.1*fit$par[4]),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20+i, ylab=expression((r[0]*m[PS])^2/(r[0]*mu[R])), xlab=expression(r[0]*mu[R]))
      lines(xfit, msq/xfit, lty=2)
      points(0., fit$par[4], pch=25, bg="sandybrown")
      if(!is.null(fit$boot.result)) {
        arrows(0,
               fit$par[4]-fit$boot.result[(2*N+13+4),2],
               0,
               fit$par[4]+fit$boot.result[(2*N+13+4),2],
               length=0.01,angle=90,code=3)
      }
    }
    else {
      plotwitherror(xmu[ij],
                    (fit$par[4+i]*mpsV[ij])^2/(xmu[ij]),
                    2*fit$par[4+i]^2*mpsV[ij]*fit$data[[i]]$dmps[ij]/(xmu[ij]),
                    rep=TRUE, col=color[i], pch=20+i, bg=color[i])
      if(fit$fit.asq) {
        lines(xfit, msq/xfit, lty=2)
      }
    }
    if(fit$fit.asq && i==1) {
      lines(xfit, msq2/xfit, lty=1)
    }
    if(fit$fit.mN) {
      dev.set(mNplot)
      if(i==1) {
        plotwitherror(xmu[ij],
                      fit$par[4+i]*fit$data[[i]]$mN[ij], fit$par[4+i]*fit$data[[i]]$dmN[ij],
                      ylim=c(0.90*par[5+2*N],
                        1.1*max(fit$par[4+i]*fit$data[[i]]$mN[ij], na.rm=TRUE)),
                      xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                      pch=20+i, ylab=expression(r[0]*m[N]), xlab=expression(r[0]*mu[R]))
        lines(xfit, mN, lty=2)
        points(fit$result$r0*fit$result$mu.phys/0.1973, fit$result$r0*fit$result$mN/0.1973,
               pch=25, bg="sandybrown")
        if(!is.null(fit$boot.result)) {
          arrows(fit$result$r0*fit$result$mu.phys/0.1973,
                 fit$result$r0*fit$result$mN/0.1973-fit$boot.result[(2*N+13+5+2*N),2],
                 fit$result$r0*fit$result$mu.phys/0.1973,
                 fit$result$r0*fit$result$mN/0.1973+fit$boot.result[(2*N+13+5+2*N),2],
                 length=0.01,angle=90,code=3)
        }
      }
      else {
        plotwitherror(xmu[ij],
                      fit$par[4+i]*fit$data[[i]]$mN[ij], fit$par[4+i]*fit$data[[i]]$dmN[ij],
                      rep=TRUE, col=color[i], pch=20+i, bg=color[i])
        if(fit$fit.asq) {
          lines(xfit, mN, lty=2)
        }
      }
      if(fit$fit.asq && i==1) {
        lines(xfit, mN2, lty=1)
      }
    }
  }
  if(plot.file) {
    dev.set(fplot)
    dev.copy2eps(file="fpsvmu.eps")
    dev.set(mpsplot)
    dev.copy2eps(file="mpssqvmu.eps")
    dev.set(mpsmuplot)
    dev.copy2eps(file="mpssqovmuvmu.eps")
    if(fit$fit.mN) {
      dev.set(mNplot)
      dev.copy2eps(file="mNvmu.eps")
    }
  }
}

print.chiralfit <- function(fit, ...) {
  summary(fit, ...)
}

summary.chiralfit <- function(fit, show.input=FALSE, show.chis=FALSE) {
  N <- length(fit$data)
  npar <- length(fit$par)
  if(fit$fit.mN) {
    Hinv <- try(solve( fit$fit$hessian ))
  }
  else Hinv <- diag(1., nrow=npar, ncol=npar)
  dpar <- numeric()
  for(i in 1:npar) {
    dpar[i] <- sqrt(2*Hinv[i,i])
  }
  if(!is.null(fit$boot.result)) {
    cat("Errors computed using", fit$boot.R, "bootsamples \n")
  }
  cat("FS           = ", fit$fsmethod, "\n")
  cat("chisqr       = ", fit$result$chisqr, "\n")
  cat("dof          = ", fit$result$dof, "\n")
  cat("red. chisqr  = ", fit$result$chisqr/fit$result$dof, "\n")
  cat("Qval         = ", 1-pgamma(fit$result$chisqr/2, fit$result$dof/2), "\n")
  if(!is.null(fit$boot.result)) {
    if(fit$fit.l12) {
      cat("l1           = ", fit$result$l1, "+-", sd(fit$boots[,(3+2*N)], na.rm=TRUE), "\n")
      cat("l2           = ", fit$result$l2, "+-", sd(fit$boots[,(4+2*N)], na.rm=TRUE), "\n")
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
    cat("r0           = ", fit$result$r0, "+-", sd(fit$boots[,(12+2*N)], na.rm=TRUE), "\n")
    cat("mu_phys      = ", fit$result$mu.phys[1], "+-", sd(fit$boots[,1], na.rm=TRUE), "GeV \n\n")
    for(i in 1:N) {
      cat("lattice spacing", i, ":\n")
      cat("lattice spacing at r0/a = ",fit$r0data$r0[i], ": a = ", fit$result$a[i], "+-",
          sd(fit$boots[,(N+i)], na.rm=TRUE),"fm \n")
      cat("            fitted r0/a = ", fit$par[4+i], "+-", sd(fit$boots[,(13+2*N+4+i)], na.rm=TRUE), "\n")
      cat("            fitted ZP   = ", fit$par[4+N+i], "+-", sd(fit$boots[,(13+2*N+4+N+i)], na.rm=TRUE), "\n")
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
    cat("Fit Parameter:\n")
    print(data.frame(Par=fit$par, Err=fit$boot.result[(2*N+14):(2*N+13+npar),2], Bias=fit$par- fit$boot.result[(2*N+14):(2*N+13+npar),1], ErrHessian=dpar))
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
    cat("r0           = ", fit$result$r0, "\n")
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
    cat("Fit Parameter:\n")
    print(data.frame(Par=fit$par, ErrHessian=dpar))
  }
  if(show.chis) {
    getchi.Na.withr0ZP(par=fit$par, data=fit$data, ii=fit$ii, r0data=fit$r0data, ZPdata=fit$ZPdata,
                       fsmethod="cdh", a=fit$result$a, fit.l12=fit$fit.l12)
  }
}

printtab <- function(x, dx, digits=2, c=1.) {
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

tab <- function(fit) {
  par <- fit$par
  br <- fit$boot.result
  npar <- length(par)
  N <- length(fit$data)
  cat(fit$fsmethod, "\n")
  nm <- 2*N+13

  printtab(par[4], br[nm+4,2])
  printtab(par[3], br[nm+3,2])
  if(fit$fit.mN) printtab(par[2*N+5], br[nm+2*N+5,2])
  for(i in 1:N) {
    printtab(par[4+i], br[nm+4+i,2])
  }
  for(i in 1:N) {
    printtab(par[4+N+i], br[nm+4+N+i,2])
  }
  for(i in 1:2) {
    if(fit$fit.l12) printtab(par[7+2*N+i], br[nm+7+2*N+i,2])
    else cat("$-$\n")
  }
  for(i in 1:2) {
    printtab(par[i], br[nm+i,2])
  }
  for(i in 3:4) {
    if(fit$fit.kmf) printtab(par[7+2*N+i], br[nm+7+2*N+i,2])
    else cat("$-$\n")
  }
  for(i in 1:2) {
    if(fit$fit.mN) printtab(par[5+2*N+i], br[nm+5+2*N+i,2])
  }
  for(i in 2:0) {
    if(fit$fit.asq) {
      printtab(par[npar-i], br[nm+npar-i,2])
    }
  }
  cat("$",fit$result$chisqr,"/",fit$result$dof,"$\n", sep="")
  cat("---\n")
  printtab(fit$result$mu.phys*1000, br[1,2], c=1000.)
  for(i in 1:2) {
    printtab(fit$result$a[i], br[N+i,2])
  }
  printtab(fit$result$r0, br[12+2*N,2])
  if(fit$fit.l12) {
    printtab(fit$result$l1, br[3+2*N, 2])
    printtab(fit$result$l2, br[4+2*N, 2])
  }
  else cat("$-$\n$-$\n")
  printtab(fit$result$l3, br[1+2*N, 2])
  printtab(fit$result$l4, br[2+2*N, 2])
  printtab(fit$result$F, br[5+2*N, 2], c=1000.)
  printtab(1./fit$result$F*0.1307, sd(1./fit$boots[,(5+2*N)], na.rm=TRUE)*0.1307)
  printtab(fit$result$B0, br[6+2*N, 2], c=1000.)
  printtab(fit$result$Sigma, br[10+2*N, 2], c=1000.)
  printtab(fit$result$rssq, br[11+2*N, 2])
  if(fit$fit.mN) {
    printtab(fit$result$mN0, br[7+2*N, 2], c=1000.)
    printtab(fit$result$c1, br[9+2*N, 2])
    printtab(fit$result$s0, br[13+2*N, 2], c=1000)
    printtab(fit$result$mN, br[8+2*N, 2], c=1000)
    printtab(fit$result$mN/fit$result$mN0, sd(fit$boots[,(8+2*N)]/fit$boots[,(7+2*N)], na.rm=TRUE))
  }
}
