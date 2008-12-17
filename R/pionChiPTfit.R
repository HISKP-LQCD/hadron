fovermps.pion <- function(x, par, N, fit.nnlo=FALSE, fit.kmf=fit.kmf, fit.asq=-1.) {
  r0sqTwoBmu <- par[4]*x
  mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  if(any(mpssq <= 0) || any(is.nan(mpssq))) return(NaN)
  fpsV <- getfps.pion(r0sqTwoBmu, par=par, N=N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  mpsV <- sqrt(mpssq)
  return(fpsV/mpsV - (130.7/135))
}


# order of fit parameters (N no of lattice spacings)
# 1:              r0*L3
# 2:              r0*L4
# 3:              r0*f0
# 4:              2*r0*B0
# 5 to 4+N:       r0
# 5+N to 4+2*N    sr0 (slope of r0 wrt mu^2 )
# 5+2*N to 4+3*N: ZP
# 4+3*N+1:        r0*L1
# 4+3*N+2:        r0*L2
# 4+3*N+3:        k_M
# 4+3*N+4:        k_F
# npar-1:         D_m^0
# npar:           D_f^0

pionChiPTfit <- function(data, startvalues, bootsamples, fsmethod="gl", a.guess,
                         ZPdata, ZPbootsamples, method="no",
                         fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                         fit.corr=FALSE, ii, boot.R=100, debug=FALSE) {
  if(missing(data) || missing(startvalues) || missing(ZPdata)) {
    stop("data and startvalues must be provided!\n")
  }
  if(missing(bootsamples)) {
    cat("bootstrap samples missing, will not compute error estimates!\n")
    method <- "no"
  }
  else if(missing(boot.R)) {
    boot.R <- length(bootsamples[[1]][,1,1])
  }
  if(fit.nnlo) {
    fit.l12 <- TRUE
  }
  if(!fit.nnlo) {
    fit.kmf=FALSE
  }
  N <- length(data)
  if(missing(ii)) {
    ii <- list(c(1:length(data[[1]]$mu)))
    for (i in 2:N) {
      ii <- c(ii, list(c(1:length(data[[i]]$mu))))
    }
  }
  else {
    N <- length(ii)
  }
  np <- 4+3*N
  if(fit.l12) np=np+2
  if(fit.asq) np=np+2
  if(fit.nnlo && fit.kmf) np=np+2
  if(length(startvalues) != np) stop("length of startvalues", length(startvalues),
             "must match number of fit parameters ", np, "!\n")
  if(fsmethod=="cdh" || fsmethod=="cdhnew") {
    if(missing(a.guess)) {
      cat("a.guess was missing! Guessing myself!\n")
      a.guess <- rep(0.1, times=N) 
    }
  }

  corrmatrix <- NULL
  if(!missing(bootsamples) && fit.corr==TRUE) {
    corrmatrix <- list()
    for(i in 1:N) {
      corrmatrix[[i]] <- array(0., dim=c(2,2,length(data[[i]]$mu)))
      for(j in 1:length(data[[i]]$mu)) {
        corr <- cor(bootsamples[[i]][,1:2,j])
        cat(i, " ", j, "\n")
        print(corr)
        corr[1,1] <- corr[1,1]*data[[i]]$dmps[j]^2
        corr[2,2] <- corr[2,2]*data[[i]]$dfps[j]^2
        corr[2,1] <- corr[2,1]*data[[i]]$dmps[j]*data[[i]]$dfps[j]
        corr[1,2] <- corr[1,2]*data[[i]]$dmps[j]*data[[i]]$dfps[j]
        corrmatrix[[i]][,,j] <- solve(corr)
      }
    }
  }
  else {
    warning("no bootstrapsamples to compute the correlation matrix\nPerform uncorrelated fit\n")
    fit.corr = FALSE
  }
  
  if(debug) {
    chisqr.piononly(par=startvalues, data=data, ii=ii,fsmethod=fsmethod,
                    a.guess=a.guess, ZPdata=ZPdata,
                    fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, cm=corrmatrix,
                    printit=debug)
  }
  parscale <- startvalues
  for(i in 1:length(startvalues)) {
    if(parscale[i] < 0.001) parscale[i] <- 1.
  }
  mini <- optim(par=startvalues, fn=chisqr.piononly, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, parscale=parscale),
                data=data, ii=ii,
                fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                fit.mN=fit.mN, cm=corrmatrix)
  mini <- optim(par=mini$par, fn=chisqr.piononly, method="BFGS", hessian=TRUE,
                control=list(maxit=500, trace=debug, parscale=mini$par),
                data=data, ii=ii,
                fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                fit.mN=fit.mN, cm=corrmatrix)
  if(mini$convergence != 0) {
    warning("Attention: optim did not converge in initial run!\n Please adjust initial parameter values!")
  }
  if(debug) {
    chisqr.piononly(par=mini$par, data=data, ii=ii,fsmethod=fsmethod,
                    a.guess=a.guess, ZPdata=ZPdata,
                    fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, cm=corrmatrix,
                    printit=debug)
  }
  dof <- 0

  for(i in 1:length(data)) {
    dof <- dof + length( na.omit(data[[i]]$mps[ii[[i]]])) +
      length( na.omit( data[[i]]$fps[ii[[i]]] ) ) +
        length( na.omit( data[[i]]$r0a ) )
  }
  dof <- dof + length(ZPdata$ZP) - length(startvalues)
  if(fit.l12) {
    dof <- dof + 2
  }
  if(fit.kmf) {
    dof <- dof + 2
  }
  chisqr <- mini$value
  par <- mini$par
  
                                        # compute observables
  mu.phys <- numeric()
  a <- numeric()
  r0 <- numeric()
  l3 <- numeric()
  l4 <- numeric()
  l1 <- numeric()
  l2 <- numeric()
  F <- numeric()
  B0 <- numeric()
  R0 <- numeric()
  Zp <- numeric()
  sr0 <- numeric()
  Sigma <- numeric()
  rssq <- numeric()
  
  r0mu.phys <- uniroot(fovermps.pion, c(0.00001, 0.050),
                     tol=1.e-12, par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1., N=N)$root
  r0TwoB <- par[4]
  r0sqTwoBmu <- r0TwoB*r0mu.phys
  r0F <- par[3]
  mpisq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
  fpi <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
  lmpssq <-  log( 0.1396^2 )
  l3 <- log((par[1]/fpi*0.1307)^2) - lmpssq
  l4 <- log((par[2]/fpi*0.1307)^2) - lmpssq
  if(fit.l12) {
    l1 <- log((par[5+3*N]/fpi*0.1307)^2) - lmpssq
    l2 <- log((par[6+3*N]/fpi*0.1307)^2) - lmpssq
  }
  F <- r0F/fpi*0.1307
  B0 <- r0TwoB/fpi*0.1307/2.
  r0 <- fpi/0.1307*0.1973
  mu.phys <-  r0mu.phys/fpi*0.1307
  Sigma <- (B0*F^2/2)^(1/3)
  rssq <- 12./(4*pi*F)^2*(l4-12./13.)*0.1973^2
  for(i in 1:N) {
    a[i] <- fpi/0.1307*0.1973/par[4+i]
    sr0[i] <- par[4+i+N]
    Zp[i] <- par[4+i+2*N]
  }
  
  boot.result <- NULL
  boots <- NULL
  if(method != "no") {

    if(missing(ZPbootsamples)) {
      ZPbootsamples <- array(0., dim=c(boot.R,N))
      for(i in 1:N) {
        ZPbootsamples[,i] <- rnorm(boot.R, mean=ZPdata$ZP[i], sd=ZPdata$dZP[i])
      }
    }

    boots <- array(0., dim=c(boot.R, (3*N+length(startvalues)+9)))
    for(s in 1:boot.R) {
      df <- list(data.frame(mu=data[[1]]$mu, mps=bootsamples[[1]][s, 1,],
                            dmps=data[[1]]$dmps,
                            fps=bootsamples[[1]][s, 2,],
                            dfps=data[[1]]$dfps, L=data[[1]]$L,
                            r0a=rnorm(length(data[[1]]$r0a), mean=data[[1]]$r0a, sd=data[[1]]$dr0a),
                            dr0a=data[[1]]$dr0a))
      ZPdf <- data.frame(ZP=ZPbootsamples[s,], dZP=ZPdata$dZP[1:N])
      for(i in 2:N) {
        df[[i]] <- data.frame(mu=data[[i]]$mu, mps=bootsamples[[i]][s,1,],
                              dmps=data[[i]]$dmps,
                              fps=bootsamples[[i]][s,2,],
                              dfps=data[[i]]$dfps, L=data[[i]]$L,
                              r0a=rnorm(length(data[[i]]$r0a), mean=data[[i]]$r0a, sd=data[[i]]$dr0a),
                              dr0a=data[[i]]$dr0a)
      }

      mini.boot <- optim(par=par, fn=chisqr.piononly, method="BFGS", hessian=FALSE,
                         control=list(maxit=150, parscale=par, trace=TRUE),
                         data=df, ii=ii,
                         fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdf, cm=corrmatrix,
                         fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf)
      if(mini.boot$convergence != 0) {
        warning("optim did not converge for one bootsample")
        par.boot <- rep(NA, length(mini.boot$par))
        boots[s,] <- NA
      }
      else {
        par.boot <- mini.boot$par
        r0F <- numeric()
        mpssqr <- numeric()
        boots[s,1] <- uniroot(fovermps.pion, c(0.0001, 0.030),
                              tol=1.e-12, par=par.boot, fit.nnlo=fit.nnlo,
                              fit.kmf=fit.kmf, fit.asq=-1., N=N)$root
        
        r0sqTwoBmu <- par.boot[4]*boots[s,1]
        r0F <- par.boot[3]
        mpssqr <- getmpssq.pion(r0sqTwoBmu, par.boot, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
        fpi <- getfps.pion(r0sqTwoBmu, par.boot, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
        #mu.phys
        boots[s,1] <- boots[s,1]/fpi*0.1307
        for(i in 1:N) {
          boots[s,(i+N)] <- fpi/0.1307*0.1973/par.boot[4+i]
          boots[s,(i+2*N)] <- par.boot[4+N+i]
        }
        lmpssqr <- log( 0.1396^2 )
        # l3 and l4
        boots[s,(1+3*N)] <- log((par.boot[1]/fpi*0.1307)^2) - lmpssqr
        boots[s,(2+3*N)] <- log((par.boot[2]/fpi*0.1307)^2) - lmpssqr
        if(fit.l12) {
          # l1 and l2
          boots[s,(3+3*N)] <- log((par.boot[5+3*N]/fpi*0.1307)^2) - lmpssqr
          boots[s,(4+3*N)] <- log((par.boot[6+3*N]/fpi*0.1307)^2) - lmpssqr
        }
        # r0*f0 and 2*r0*B0
        boots[s,(5+3*N)] <- par.boot[3]/fpi*0.1307
        boots[s,(6+3*N)] <- par.boot[4]/fpi*0.1307/2.
        # r0
        boots[s,(7+3*N)] <- fpi/0.1307*0.1973
        # <r^2>_s
        boots[s,(8+3*N)] <- 12./(4*pi*boots[s,(5+3*N)])^2*(boots[s,(2+3*N)]-12./13.)*0.1973^2
        # sigma
        boots[s,(9+3*N)] <- (boots[s,(5+3*N)]*(boots[s,(6+3*N)])^2/2)^(1/3)
        for(i in 1:length(par.boot)) {
          boots[s,(3*N+9+i)] <- par.boot[i]
        }
      }
    }
    boot.result <- array(0., dim=c(length(boots[1,]), 2))
    for(i in 1:length(boots[1,])) {
      boot.result[i,1] <-  mean(boots[,i], na.rm=TRUE)
      boot.result[i,2] <-  sd(boots[,i], na.rm=TRUE)
    }
  }

  if(missing(bootsamples)) {
    bootsamples <- NULL

  }
  if(method == "no") {
    ZPbootsamples <- NULL
  }
  
  result <- list(par=par, result=list(mu.phys=mu.phys, F=F, a=a, l1=l1, l2=l2, l3=l3, l4=l4, r0=r0,
                            B0=B0, ZP=Zp, Sigma=Sigma, rssq=rssq, sr0=sr0,
                            chisqr=chisqr, dof=dof), fit=mini,
                 data=data, boot.result=boot.result, boots=boots, ZPdata=ZPdata,
                 bootsamples=bootsamples, ZPbootsamples=ZPbootsamples, method=method,
                 ii=ii, fit.l12=fit.l12, boot.R=boot.R, fsmethod=fsmethod, fit.asq=fit.asq,
                 fit.kmf=fit.kmf, fit.nnlo=fit.nnlo, fit.mN=FALSE, fit.corr=fit.corr)
  attr(result, "class") <- c("chiralfit", "list")  
  return(invisible(result))
}

chisqr.piononly <- function(par, data, ii, ZPdata, fsmethod="gl", a.guess,
                            fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                            printit=FALSE, fit.mN=TRUE, cm=NULL) {

  fit.a <- -1.
  N <- length(ii)
  chisum <- 0.
  r0TwoB <- par[4]
  r0F <- par[3]
  fpi <- numeric()

  r0mu.phys <- try(uniroot(fovermps.pion,
                           c(0.0001, 0.030), tol=1.e-12,
                           par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1.,
                           fit.kmf=fit.kmf, N=N)$root, silent=T)
  for( i in 1:N) {
    if(fit.asq) {
      fit.a <- i
    }
    ij <- ii[[i]]
    if(any(r0TwoB <= 0)) {
      return(invisible(NaN))
    }
# fit r0/a as a function of (a mu)^2 first
    if(inherits(r0mu.phys, "try-error") || is.nan(r0mu.phys)) {
      r0 <- par[4+i]  + (data[[i]]$mu^2)*par[4+N+i]
    }
    else {
      r0 <- par[4+i]  + (data[[i]]$mu^2 -r0mu.phys/par[4+i])*par[4+N+i]
    }
    chisum <- chisum + sum(((data[[i]]$r0a-r0)/data[[i]]$dr0)^2, na.rm=TRUE)
    
# now the rest
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+2*N+i]
    mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
    fpsV <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)

    if(any(is.nan(mpssq)) || any(mpssq <= 0)) {
      return(NaN)
    }
    if(fsmethod=="cdh" || fsmethod=="cdhnew") {
      if(inherits(r0mu.phys, "try-error") || is.nan(r0mu.phys)) a_fm <- a.guess[i]
      else {
        a_fm <- getfps.pion(r0sqTwoBmu=r0TwoB*r0mu.phys, par, fit.nnlo=fit.nnlo,
                            fit.kmf=fit.kmf, fit.asq=FALSE)/0.1307*0.1973/par[4+i]
      }
      if(fit.l12) {
        aLamb1=par[5+3*N]/par[4+i]
        aLamb2=par[6+3*N]/par[4+i]
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*a_fm/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*a_fm/0.1973)^2)
      }

      if(fsmethod=="cdhnew") {
        res <- cdhnew(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                      aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i],
                      afpiV=fpsV/par[4+i], aF0=par[3]/par[4+i],
                      a2B0mu=r0sqTwoBmu/par[4+i]^2, L=data[[i]]$L[ij], rev=1,
                      printit=FALSE)
      }
      else {
        res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                   aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i], afpiV=fpsV/par[4+i],
                   aF0=fpsV/par[4+i], a_fm=a_fm, L=data[[i]]$L[ij], rev=1, printit=FALSE)
      }
      mpsV <- res$mpiFV
      fpsV <- res$fpiFV
    }
    else {
      mpsv <- data[[i]]$L[ij]*sqrt(mpssq)/par[4+i]
      r <-  mpssq/(4.0*pi*par[3])^2*g1(mpsv)/par[4+i]

      mpsV <- sqrt(mpssq)*(1.+0.5*r)/par[4+i]
      fpsV <- fpsV*(1.0-2.0*r)/par[4+i]
    }
    if(is.null(cm)) {
      chisum <- chisum + (sum(((data[[i]]$mps[ij]-mpsV)/(data[[i]]$dmps[ij]))^2) +
                          sum(((data[[i]]$fps[ij]-fpsV)/(data[[i]]$dfps[ij]))^2))
    }
    else {
      chisum <- chisum + (sum(((data[[i]]$mps[ij]-mpsV))^2*cm[[i]][1,1,ij]) +
                          sum(((data[[i]]$fps[ij]-fpsV))^2*cm[[i]][2,2,ij]) +
                          2*sum((data[[i]]$mps[ij]-mpsV)*(data[[i]]$fps[ij]-fpsV)*cm[[i]][2,1,ij]))
    }
#    cat((sum(((data[[i]]$mps[ij]-mpsV)/(data[[i]]$dmps[ij]))^2) + sum(((data[[i]]$fps[ij]-fpsV)/(data[[i]]$dfps[ij]))^2)), "\n")
#    cat(sum((data[[i]]$mps[ij]-mpsV)^2*cm[[i]][1,1,ij] + (data[[i]]$fps[ij]-fpsV)^2*cm[[i]][2,2,ij] + (data[[i]]$mps[ij]-mpsV)*(data[[i]]$fps[ij]-fpsV)*cm[[i]][2,1,ij]), "\n",cm[[i]][1,1,ij], "\n\n")
    if(fit.l12 && i==1) {
      chisum <- chisum + ((-0.4-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/0.6)^2 +
        ((4.3-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/0.1)^2
    }
    if(fit.kmf && i==1) {
      chisum <- chisum + (par[7+2*N]/5.)^2 + (par[8+2*N]/5.)^2
    }
    if(printit) {
      cat("r0chiral", par[4+i], "\n")
      cat("datar0  ", data[[i]]$r0a[ij], "\n")
      cat("modelr0 ", r0, "\n")
      cat("chir0   ", (data[[i]]$r0a[ij]-r0)/data[[i]]$dr0[ij], "\n")
      cat("modelZP ", par[4+2*N+i], "\n")
      cat("ZPdata  ", ZPdata$ZP[i], "\n")
      cat("chiZP   ", (ZPdata$ZP[i]-par[4+2*N+i])/ZPdata$dZP[i], "\n")
      cat("modelm  ", mpsV, "\n")
      cat("datam   ", data[[i]]$mps[ij], "\n")
      cat("errm    ", data[[i]]$dmps[ij], "\n")
      cat("chim    ", ((data[[i]]$mps[ij]-mpsV)/data[[i]]$dmps[ij]), "\n")
      cat("modelf  ", fpsV, "\n")
      cat("dataf   ", data[[i]]$fps[ij], "\n")
      cat("errf    ", data[[i]]$dfps[ij], "\n")
      cat("chif    ", ((data[[i]]$fps[ij]-fpsV)/data[[i]]$dfps[ij]), "\n")
    }

  }
  # ZP
  chisum <- chisum + sum(((par[(5+2*N):(4+3*N)]-ZPdata$ZP)/ZPdata$dZP)^2)
  if(printit) {
    cat("chisqr ", chisum, "\n\n")
  }
  return(invisible(chisum))
}

getchi.piononly <- function(par, data, ii, ZPdata, fsmethod="gl", a, fit.l12=FALSE)
{
  N <- length(ii)
  chisum <- 0.
  for( i in 1:N) {
    a_fm = a[i]
    ij <- ii[[i]]
    r0TwoB <- par[4]
    if(any(r0TwoB < 0)) {
      return(invisible(NaN))
    }
    r0F <- par[3]
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+N+i]
    mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.l12)
    fpsV <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.l12)
    if(any(is.nan(mpssq)) || any(mpssq < 0)) {
      return(NaN)
    }
    if(fsmethod=="cdh" || fsmethod=="cdhnew") {
      if(fit.l12) {
        aLamb1=par[8+2*N]/par[4+i]
        aLamb2=par[9+2*N]/par[4+i]
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*a_fm/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*a_fm/0.1973)^2)
      }
      if(fsmethod=="cdhnew") {
        res <- cdhnew(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                      aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i],
                      afpiV=fpsV/par[4+i], aF0=par[3]/par[4+i],
                      a2B0mu=r0sqTwoBmu/par[4+i]^2, L=data[[i]]$L[ij], rev=1,
                      printit=FALSE)
      }
      else {
        res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                   aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i], afpiV=fpsV/par[4+i],
                   aF0=fpsV/par[4+i], a_fm=a_fm, L=data[[i]]$L[ij], rev=1, printit=FALSE)
      }
      mpsV <- res$mpiFV
      fpsV <- res$fpiFV
    }
    else {
      mpsv <- L[ij]*sqrt(mpssq)/par[4+i]
      r <-  mpssq/(4.0*pi*par[3])^2*g1(mpsv)/par[4+i]

      mpsV <- sqrt(mpssq)*(1.+0.5*r)/par[4+i]
      fpsV <- fpsV*(1.0-2.0*r)/par[4+i]
    }
    cat("\n*** lattice spacing", i, "a=", a_fm,"fm ***\n")
    cat("mps:\n")
    print(data.frame(data=data[[i]]$mps[ij], Err=data[[i]]$dmps[ij], Fit=mpsV, Chi=(data[[i]]$mps[ij]-mpsV)/(data[[i]]$dmps[ij])))
    cat("fps:\n")
    print(data.frame(data=data[[i]]$fps[ij], Err=data[[i]]$dfps[ij], Fit=fpsV, Chi=(data[[i]]$fps[ij]-fpsV)/(data[[i]]$dfps[ij])))
    if(fit.l12) {
      cat("l12:\n")
      print(data.frame(data=c(-0.4,4.3), Err=c(0.6,0.1),
                       Fit=c((log((aLamb1)^2/(0.1396*a_fm/0.1973)^2)), (log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))),
                       Chi=c((-0.4-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/0.6, (4.3-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/0.1)
                       ))
    }
#    cat("r0:\n")
#    print(data.frame(data=r0data$r0[i],Err=r0data$dr0[i], Fit=par[4+i], Chi=(r0data$r0[i]-par[4+i])/r0data$dr0[i]))
    cat("ZP:\n")
    print(data.frame(data=ZPdata$ZP[i],Err=ZPdata$dZP[i], Fit=par[4+N+i], Chi=(ZPdata$ZP[i]-par[4+N+i])/ZPdata$dZP[i]))
    
# the nucleon
    mN <- getmN(r0sqTwoBmu=r0sqTwoBmu, par, N)/par[4+i]
    cat("mN:\n")
    print(data.frame(data=data[[i]]$mN, Err=data[[i]]$dmN, Fit=mN, Chi=(data[[i]]$mN-mN)/data[[i]]$dmN))
  }
}

getmpssq.pion <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1) {

  npar <- length(par)
  xi <- r0sqTwoBmu/(4.0*pi*par[3])^2
  rln3 <- log(r0sqTwoBmu/par[1]^2)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
  }
  if(fit.nnlo) {
    fit.k <- 0.
    if(fit.kmf) {
      fit.k <- 1.
    }
    rln1 <- log(r0sqTwoBmu/par[3*N+5]^2)
    rln2 <- log(r0sqTwoBmu/par[3*N+6]^2)
    # this line makes it work even if we do not fit k_M,F
    # however, fit.k=0. then...
    np <- ifelse((3*N+7)<=npar, 3*N+7, npar)
    return(r0sqTwoBmu*(1. + xi*rln3 + asq/par[4+fit.asq]^2*par[npar-1] +
                       17./2.*xi^2*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)^2
                                    +fit.k*8./17.*par[np]
                                    )
                       )
           )
  }
  return(r0sqTwoBmu*(1.0 + asq/par[4+fit.asq]^2*par[npar-1] + r0sqTwoBmu*(rln3/(4.0*pi*par[3])^2 )))
}

getfps.pion <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1.) {

  npar <- length(par)
  xi <- r0sqTwoBmu/(4.0*pi*par[3])^2
  rln3 <- log(r0sqTwoBmu/par[1]^2)
  rln4 <- log(r0sqTwoBmu/par[2]^2)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
  }
  if(fit.nnlo) {
    fit.k <- 0.
    if(fit.kmf) {
      fit.k <- 1.
    }
    rln1 <- log(r0sqTwoBmu/par[3*N+5]^2)
    rln2 <- log(r0sqTwoBmu/par[3*N+6]^2)
    np <- ifelse((3*N+8)<=npar, 3*N+8, npar)
    return(par[3]*(1. - 2.*xi*rln4 +  asq/par[4+fit.asq]^2*par[npar] -
                   5*xi^2*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)^2
                           + fit.k*4./5.*par[np]
                           )
                   )
           )
  }
  return(par[3]*(1.0 + asq/par[4+fit.asq]^2*par[npar]  - 2.0*r0sqTwoBmu*(rln4)/(4.0*pi*par[3])^2 ))
}
