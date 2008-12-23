plot.pionChiPTfit <- function(fit, write.data=FALSE, plot.file=FALSE, plot.all=FALSE,
                              filesuffix="", rawpoints=TRUE,
                              leg.text=c(expression(beta==3.9), expression(beta==4.05)), ...) {
  N <- length(fit$data)
  npar <- length(fit$par)
  dpar <- numeric(npar)
  Hinv <- try(solve( fit$fit$hessian ))
  if(inherits(Hinv, "try-error")) {
    Hinv <- diag(nrow=npar, ncol=npar)
    dpar <- rep(0.0001, times=npar)
  }
  else {
    for(i in 1:npar) {
      dpar[i] <- sqrt(2*Hinv[i,i])
    }
  }
  fit.a <- -1.
  par <- fit$par
  if(!plot.file) {
    X11()
    par(list(las=1, tck=0.015, mgp=c(3.,.5,0)))
    fplot <- dev.cur()
    X11()
    par(list(las=1, tck=0.015, mgp=c(3.,.5,0)))
    mpsplot <- dev.cur()
    X11()
    par(list(las=1, tck=0.015, mgp=c(3.,.5,0)))
    mpsmuplot <- dev.cur()
    X11()
    par(list(las=1, tck=0.015, mgp=c(3.,.5,0)))
    r0plot <- dev.cur()
    if(fit$fit.mN) {
      X11()
      par(list(las=1, tck=0.015, mgp=c(3.,.5,0)))
      mNplot <- dev.cur()
    }
  }
  else {
    postfix <- paste(filesuffix, ".eps", sep="")
    if(fit$fit.asq) postfix <- paste(".asq", postfix, sep="")
    if(fit$fit.nnlo) postfix <- paste(".nnlo", postfix, sep="")
    if(fit$fit.kmf) postfix <- paste(".kmf", postfix, sep="")
    postfix <- paste(".",fit$fsmethod, postfix, sep="")
    postscript(file = paste("fpsvmu", postfix, sep=""), family="NimbusRom", paper="special",
               horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
    par(tck=0.02, mgp=c(2.5,.5,0), las=1, lwd=0.5, cex.lab=1.2)
    fplot <- dev.cur()
    postscript(file = paste("mpssqvmu", postfix, sep=""), family="NimbusRom", paper="special",
               horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
    par(tck=0.02, mgp=c(2.5,.5,0), las=1, lwd=0.5, cex.lab=1.2)
    mpsplot <- dev.cur()
    postscript(file = paste("mpssqovmuvmu", postfix, sep=""), family="NimbusRom", paper="special",
               horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
    par(tck=0.02, mgp=c(2.5,.5,0), las=1, lwd=0.5, cex.lab=1.2)
    mpsmuplot <- dev.cur()
    postscript(file = paste("r0vmu", postfix, sep=""), family="NimbusRom", paper="special",
               horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
    par(tck=0.02, mgp=c(2.5,.5,0), las=1, lwd=0.5, cex.lab=1.2)
    r0plot <- dev.cur()
    if(fit$fit.mN) {
      postscript(file = paste("mNvmu", postfix, sep=""), family="NimbusRom", paper="special",
                 horizontal=FALSE, onefile=FALSE, width=6.996766, height=6.996766)
      par(tck=0.02, mgp=c(2.5,.5,0), las=1, lwd=0.5, cex.lab=1.2)
      mNplot <- dev.cur()
    }
  }

  
  for(i in 1:N) {
    ur0 <- fit$par[4+i]
    dr0 <- dpar[4+i]
    mZP <- fit$ZPdata$ZP[i]
    fZP <- fit$par[4+2*N+i]
    dZP <- fit$ZPdata$dZP[i]
    if(rawpoints) {
      uZP <- mZP
    }
    else {
      uZP <- fZP
    }
    if(fit$fit.asq) {
      fit.a <- i
    }
    ij <- fit$ii[[i]]
    if(plot.all) {
      ij <- c(1:length(fit$data[[i]]$mu))
    }
    if(fit$fsmethod == "cdh" || fit$fsmethod == "cdhnew") {
      if(fit$fit.l12) {
        aLamb1=par[5+3*N]/ur0
        aLamb2=par[6+3*N]/ur0
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*fit$result$a[i]/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*fit$result$a[i]/0.1973)^2)
      }
      
      if(fit$fsmethod == "cdhnew") {
        res <- cdhnew(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/ur0,
                      aLamb4=par[2]/ur0, ampiV=fit$data[[i]]$mps,
                      afpiV=fit$data[[i]]$fps, aF0=par[3]/ur0,
                      a2B0mu=par[4]*fit$data[[i]]$mu/ur0/mZP, L=fit$data[[i]]$L, rev=-1,
                      printit=FALSE)
      }
      else {
        res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/ur0,
                   aLamb4=par[2]/ur0, ampiV=fit$data[[i]]$mps, afpiV=fit$data[[i]]$fps,
                   aF0=fit$data[[i]]$fps, a_fm=fit$result$a[i], L=fit$data[[i]]$L, rev=-1,
                   printit=FALSE)
      }
      
      mpsV <- res$mpiFV
      
      fpsV <- res$fpiFV
    }
    else {
      mpsv <- fit$data[[i]]$L*fit$data[[i]]$mps
      r <-  - fit$data[[i]]$mps^2/(4.0*pi*par[3])^2*g1(mpsv)
      
      mpsV <- fit$data[[i]]$mps*(1.+0.5*r)
      fpsV <- fit$data[[i]]$fps*(1.0-2.0*r)
    }
    
#    xfit <- seq(from=0., to=1.05*max(fit$data[[i]]$mu[ij]*ur0/fZP, na.rm=TRUE),
#                length.out=500)
    xfit <- seq(from=0., to=1.05*max(fit$data[[i]]$mu*ur0/fZP, na.rm=TRUE),
                length.out=500)

    
    r0TwoB <- par[4]
    r0sqTwoBmu <- r0TwoB*xfit
    msq <- getmpssq.pion(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=fit.a)
    f <- getfps.pion(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=fit.a)
    mN <- getmN(r0sqTwoBmu, fit$par, N, fit.asq=fit.a)
    color <- c("red", "blue", "darkseagreen", "darksalmon", "darkviolet")
    xmu <- ur0/uZP*fit$data[[i]]$mu
    fitr0 <- par[4+i] + par[4+N+i]*xfit^2/ur0^2*fZP^2
    pxmu <- split(xmu[ij], fit$data[[i]]$L[ij])
    pxmuall <- split(xmu, fit$data[[i]]$L)
    pdfps <- split(fit$data[[i]]$dfps[ij], fit$data[[i]]$L[ij])
    pfpsV <- split(fpsV[ij], fit$data[[i]]$L[ij])
    pmpsV <- split(mpsV[ij], fit$data[[i]]$L[ij])
    pdmps <- split(fit$data[[i]]$dmps[ij], fit$data[[i]]$L[ij])
    pmN <- split(fit$data[[i]]$mN[ij], fit$data[[i]]$L[ij])
    pdmN <- split(fit$data[[i]]$dmN[ij], fit$data[[i]]$L[ij])
    r0a <- split(fit$data[[i]]$r0a, fit$data[[i]]$L)
    dr0a <- split(fit$data[[i]]$dr0a, fit$data[[i]]$L)
    k <- length(pfpsV)
                                        # continuum curves
    if(fit$fit.asq && i == 1) {
      f2 <- getfps.pion(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=-1.)
      mN2 <- getmN(r0sqTwoBmu, fit$par, N, fit.asq=-1.)
      msq2 <- getmpssq.pion(r0sqTwoBmu, fit$par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=-1.)
    }
    if(write.data) {
      file <- paste("curve_",i,".dat", sep="")
      #r0*mu*ZP, (r0mps)^2, r0fps, r0mN
      write.table(data.frame(xfit, msq, f, mN),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
      file <- paste("data_cor_", fit$fsmethod, i,".dat", sep="")
      write.table(data.frame(fit$data[[i]]$mu, fit$data[[i]]$L,
                             ur0, fZP,
                             xmu, fpsV, fit$data[[i]]$dfps,
                             mpsV, fit$data[[i]]$dmps,
                             fit$data[[i]]$mN, fit$data[[i]]$dmN,
                             ur0, dr0,
                             mZP, dZP),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
      file <- paste("data_", i,".dat", sep="")
      write.table(data.frame(fit$data[[i]]$mu, fit$data[[i]]$L,
                             ur0, fZP,
                             fit$data[[i]]$fps, fit$data[[i]]$dfps,
                             fit$data[[i]]$mps, fit$data[[i]]$dmps,
                             fit$data[[i]]$mN, fit$data[[i]]$dmN,
                             ur0, dr0,
                             mZP, dZP),
                  row.names=FALSE, col.names=FALSE,
                  file=file)
    }
    # fps as a function of mu
    dev.set(fplot)
    if(i == 1) {
      j <- 1
      plotwitherror(pxmu[[j]],
                    ur0*pfpsV[[j]],
                    sqrt((ur0*pdfps[[j]])^2 + (dr0*pfpsV[[j]])^2),
                    ylim=c(0.9*par[3], 1.1*max(ur0*fpsV[ij], na.rm=TRUE)),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20, ylab=expression(r[0]*f[PS]), xlab=expression(r[0]*mu[R]),
                    axes=F, rep=rep)
      for(j in 2:k) {
        plotwitherror(pxmu[[j]],
                      ur0*pfpsV[[j]],
                      sqrt((ur0*pdfps[[j]])^2 + (dr0*pfpsV[[j]])^2),
                      col=color[i], bg=color[i],
                      pch=20+(i-1)*2+(j-1), rep=TRUE)
      }
      rm(j)
      axis(1, lwd=0.5)
      axis(2, lwd=0.5)
      axis(3, lwd=0.5, labels=F)
      axis(4, lwd=0.5, labels=F)
      legend(x="bottomright", legend=c(leg.text, expression(r[0]*f[pi]), "Fit"),
             pch=c(21:(21+N-1), 25, -1), col=c(color[1:N],"turquoise", "black"),
             pt.bg=c(color[1:N],"turquoise", "black"),
             inset=.05, lty=c(rep(0,times=N+1), 1))
      box()
      if(fit$fit.asq) lines(xfit, f, lty="solid", col=color[i])
      else lines(xfit, f, lty="solid")
      points(fit$result$r0*fit$result$mu.phys/0.1973, fit$result$r0*0.1307/0.1973, pch=25,
             bg="turquoise", col="turquoise")
      if(fit$method != "no") {
        arrows((fit$result$r0*fit$result$mu.phys-sd(fit$boots[,1], na.rm=TRUE))/0.1973,
               fit$result$r0*0.1307/0.1973,
               (fit$result$r0*fit$result$mu.phys+sd(fit$boots[,1], na.rm=TRUE))/0.1973,
               fit$result$r0*0.1307/0.1973,
               length=0.01,angle=0,code=3)
      }
    }
    else {
      for(j in 1:length(pfpsV)) {
        plotwitherror(pxmu[[j]],
                      ur0*pfpsV[[j]],
                      sqrt((ur0*pdfps[[j]])^2 + (dr0*pfpsV[[j]])^2),
                      rep=TRUE, col=color[i], pch=20+(i-1)*2+(j-1), bg=color[i])
      }
      if(fit$fit.asq) {
        lines(xfit, f, lty="solid", col=color[i])
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
      j <- 1
      plotwitherror(pxmu[[j]],
                    (ur0*pmpsV[[j]])^2,
                    sqrt((2*ur0^2*pmpsV[[j]]*pdmps[[j]])^2 + (2*dr0*ur0*pmpsV[[j]]^2)^2),
                    ylim=c(0., 1.1*max((ur0*mpsV[ij])^2, na.rm=TRUE)),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20, ylab=expression((r[0]*m[PS])^2), xlab=expression(r[0]*mu[R]),
                    axes=F)
      for(j in 2:k) {
        plotwitherror(pxmu[[j]],
                      (ur0*pmpsV[[j]])^2,
                      sqrt((2*ur0^2*pmpsV[[j]]*pdmps[[j]])^2 + (2*dr0*ur0*pmpsV[[j]]^2)^2),
                      col=color[i], bg=color[i],
                      pch=20+(j-1), rep=TRUE)
      }
      axis(1, lwd=0.5)
      axis(2, lwd=0.5)
      axis(3, lwd=0.5, labels=F)
      axis(4, lwd=0.5, labels=F)
      legend(x="bottomright", legend=c(leg.text, "Fit"),
             pch=c(21:(21+N-1), -1), col=c(color[1:N], "black"),
             pt.bg=c(color[1:N], "black"),
             inset=.05, lty=c(rep(0,times=N), 1))
      box()
      if(fit$fit.asq) lines(xfit, msq, lty="solid", col=color[i])
      else lines(xfit, msq, lty="solid")
    }
    else {
      for(j in 1:k) {
        plotwitherror(pxmu[[j]],
                      (ur0*pmpsV[[j]])^2,
                      sqrt((2*ur0^2*pmpsV[[j]]*pdmps[[j]])^2 + (2*dr0*ur0*pmpsV[[j]]^2)^2),
                      col=color[i], bg=color[i],
                      pch=20+(i-1)*2+(j-1), rep=TRUE)
      }
      if(fit$fit.asq) {
        lines(xfit, msq, lty="solid", col=color[i])
      }
    }
    if(fit$fit.asq && i==1) {
      lines(xfit, msq2, lty=1)
    }
    dev.set(mpsmuplot)
    if(i == 1) {
      j <- 1
      plotwitherror(pxmu[[j]],
                    (ur0*pmpsV[[j]])^2/(pxmu[[j]]),
                    sqrt((2*ur0*pmpsV[[j]]*pdmps[[j]])^2 + (dr0*pmpsV[[j]]^2)^2
                         + (ur0*pmpsV[[j]]^2*dZP/uZP)^2 )*ur0/(pxmu[[j]]),
                    ylim=c(0.95*min((ur0*mpsV[ij])^2/(xmu[ij]),
                      na.rm=TRUE),
                      1.1*fit$par[4]),
                    xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20, ylab=expression((r[0]*m[PS])^2/(r[0]*mu[R])), xlab=expression(r[0]*mu[R]),
                    axes=F)
      for(j in 2:k) {
        plotwitherror(pxmu[[j]],
                      (ur0*pmpsV[[j]])^2/(pxmu[[j]]),
                      sqrt((2*ur0*pmpsV[[j]]*pdmps[[j]])^2 + (dr0*pmpsV[[j]]^2)^2
                           + (ur0*pmpsV[[j]]^2*dZP/uZP)^2 )*ur0/(pxmu[[j]]),
                      col=color[i], bg=color[i],
                      pch=20+(j-1), rep=TRUE)
      }
      axis(1, lwd=0.5)
      axis(2, lwd=0.5)
      axis(3, lwd=0.5, labels=F)
      axis(4, lwd=0.5, labels=F)
      legend(x="topright", legend=c(leg.text, expression(2*r[0]*B[0]), "Fit"),
             pch=c(21:(21+N-1), 25, -1), col=c(color[1:N],"sandybrown", "black"),
             pt.bg=c(color[1:N],"sandybrown", "black"),
             inset=.05, lty=c(rep(0,times=N+1), 1))      
      box()
      if(fit$fit.asq) lines(xfit, msq/xfit, lty="solid", col=color[i])
      else lines(xfit, msq/xfit, lty="solid")
      points(0., fit$par[4], pch=25, col="sandybrown", bg="sandybrown")
      if(fit$method != "no") {
        # 2*r0*B0
        arrows(0,
               fit$par[4]-fit$boot.result[(3*N+9+4),2],
               0,
               fit$par[4]+fit$boot.result[(3*N+9+4),2],
               length=0.01,angle=90,code=3)
      }
    }
    else {
      for(j in 1:k) {
        plotwitherror(pxmu[[j]],
                      (ur0*pmpsV[[j]])^2/(pxmu[[j]]),
                      sqrt((2*ur0*pmpsV[[j]]*pdmps[[j]])^2 + (dr0*pmpsV[[j]]^2)^2
                           + (ur0*pmpsV[[j]]^2*dZP/uZP)^2 )*ur0/(pxmu[[j]]),
                      col=color[i], bg=color[i],
                      pch=20+(i-1)*2+(j-1), rep=TRUE)
      }
      if(fit$fit.asq) {
        lines(xfit, msq/xfit, lty="solid", col=color[i])
      }
    }
    if(fit$fit.asq && i==1) {
      lines(xfit, msq2/xfit, lty=1)
    }
    dev.set(r0plot)
    if(i==1) {
      j <- 1
      plotwitherror(pxmuall[[j]]^2, r0a[[j]]/fit$par[4+i], dr0a[[j]]/fit$par[4+i],
                    ylim=c(0.85*min(fit$data[[i]]$r0a/fit$par[4+i], na.rm=TRUE),
                      1.05*max(fit$data[[i]]$r0a/fit$par[4+i], na.rm=TRUE)),
                    xlim=c(0.,1.1*max(xmu^2, na.rm=TRUE)), col=color[i], bg=color[i],
                    pch=20, ylab=expression(r[0]/r[0](mu==0)), xlab=expression((r[0]*mu[R])^2),
                    axes=F)
      for(j in 2:length(pxmuall)) {
        plotwitherror(pxmuall[[j]]^2, r0a[[j]]/fit$par[4+1], dr0a[[j]]/fit$par[4+1],
                      col=color[i], bg=color[i],
                      pch=20+(j-1), rep=TRUE)
      }
      axis(1, lwd=0.5)
      axis(2, lwd=0.5)
      axis(3, lwd=0.5, labels=F)
      axis(4, lwd=0.5, labels=F)
      legend(x="bottomright", legend=c(leg.text, expression(r[0]/r[0](mu==0)), "Fit"),
             pch=c(21:(21+N-1), 25, -1), col=c(color[1:N],"sandybrown", "black"),
             pt.bg=c(color[1:N],"sandybrown", "black"),
             inset=.05, lty=c(rep(0,times=N+1), 1))      
      box()
      lines(xfit^2, fitr0/fit$par[4+i], lty="solid", col=c(color[i]))
      points((fit$result$r0*fit$result$mu.phys/0.1973)^2, fit$par[4+i]/fit$par[4+i],
             pch=25, bg="sandybrown", col="sandybrown")
      if(fit$method != "no") {

      }
    }
    else {
      for(j in 1:length(pxmuall)) {
        plotwitherror(pxmuall[[j]]^2, r0a[[j]]/fit$par[4+i], dr0a[[j]]/fit$par[4+i],
                      col=color[i], bg=color[i],
                      pch=20+(i-1)*2+(j-1), rep=TRUE)
      }
      lines(xfit^2, fitr0/fit$par[4+i], lty="solid", col=c(color[i]))
    }
    if(fit$fit.mN) {
      dev.set(mNplot)
      if(i==1) {
        j <- 1
        plotwitherror(pxmu[[j]],
                      ur0*pmN[[j]],
                      sqrt((ur0*pdmN[[j]])^2 + (dr0*pmN[[j]])^2),
                      ylim=c(0.90*par[5+2*N],
                        1.1*max(ur0*fit$data[[i]]$mN[ij], na.rm=TRUE)),
                      xlim=c(0.,1.1*max(xmu, na.rm=TRUE)), col=color[i], bg=color[i],
                      pch=20, ylab=expression(r[0]*m[N]), xlab=expression(r[0]*mu[R]),
                      axes=F)
        for(j in 2:k) {
          plotwitherror(pxmu[[j]],
                        ur0*pmN[[j]],
                        sqrt((ur0*pdmN[[j]])^2 + (dr0*pmN[[j]])^2),
                        col=color[i], bg=color[i],
                        pch=20+(j-1), rep=TRUE)
        }
        axis(1, lwd=0.5)
        axis(2, lwd=0.5)
        axis(3, lwd=0.5, labels=F)
        axis(4, lwd=0.5, labels=F)
        legend(x="bottomright", legend=c(leg.text, expression(r[0]*m[N]), "Fit"),
               pch=c(21:(21+N-1), 25, -1), col=c(color[1:N],"sandybrown", "black"),
               pt.bg=c(color[1:N],"sandybrown", "black"),
               inset=.05, lty=c(rep(0,times=N+1), 1))      
        box()
        if(fit$fit.asq) lines(xfit, mN, lty="solid", col=color[i])
        else lines(xfit, mN, lty="solid")
        points(fit$result$r0*fit$result$mu.phys/0.1973, fit$result$r0*fit$result$mN/0.1973,
               pch=25, bg="sandybrown", col="sandybrown")
        if(fit$method != "no") {
          arrows(fit$result$r0*fit$result$mu.phys/0.1973,
                 fit$result$r0*fit$result$mN/0.1973-fit$boot.result[(3*N+9+5+2*N),2],
                 fit$result$r0*fit$result$mu.phys/0.1973,
                 fit$result$r0*fit$result$mN/0.1973+fit$boot.result[(3*N+9+5+2*N),2],
                 length=0.01,angle=90,code=3)
        }
      }
      else {
        for(j in 1:k) {
          plotwitherror(pxmu[[j]],
                        ur0*pmN[[j]],
                        sqrt((ur0*pdmN[[j]])^2 + (dr0*pmN[[j]])^2),
                        col=color[i], bg=color[i],
                        pch=20 + (i-1)*2+(j-1), rep=TRUE)
        }
        plotwitherror(xmu[ij],
                      ur0*fit$data[[i]]$mN[ij],
                      sqrt((ur0*fit$data[[i]]$dmN[ij])^2 + (dr0*fit$data[[i]]$mN[ij])^2),
                      rep=TRUE, col=color[i], pch=20+i, bg=color[i])
        if(fit$fit.asq) {
          lines(xfit, mN, lty="solid", col=color[i])
        }
      }
      if(fit$fit.asq && i==1) {
        lines(xfit, mN2, lty=1)
      }
    }
    rm(dr0, mZP, fZP, dZP, ur0, uZP)
  }
  if(plot.file) {
    dev.set(fplot)
    dev.off()
    dev.set(mpsplot)
    dev.off()
    dev.set(mpsmuplot)
    dev.off()
    dev.set(r0plot)
    dev.off()
    if(fit$fit.mN) {
      dev.set(mNplot)
      dev.off()
    }
  }
}
