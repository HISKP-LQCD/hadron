fovermps.Na.withr0 <- function(x, par, N, fit.nnlo=FALSE, fit.kmf=fit.kmf, fit.asq=-1.) {
  r0sqTwoBmu <- par[4]*x
  mpssq <- getmpssq(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  if(any(mpssq <= 0) || any(is.nan(mpssq))) return(NaN)
  fpsV <- getfps(r0sqTwoBmu, par=par, N=N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  mpsV <- sqrt(mpssq)
  return(fpsV/mpsV - (130.7/135))
}


# order of fit parameters (N no of lattice spacings)
# 1:              r0*L3
# 2:              r0*L4
# 3:              r0*f0
# 4:              2*r0*B0
# 5 to 4+N:       r0
# 5+N to 4+2*N:   ZP
# 5+2*N:          r0*m_N
# 6+2*N:          r0c_1
# 7+2*N:          g_A
# 7+2*N+1:        r0*L1
# 7+2*N+2:        r0*L2
# 7+2*N+3:        k_M
# 7+2*N+4:        k_F
# npar-2:         D_m^0
# npar-1:         D_f^0
# npar:           D_N^0

chiralfit <- function(data, startvalues, bootsamples, fsmethod="gl", a.guess,
                      r0bootsamples, r0data, ZPdata, ZPbootsamples, method="no",
                      fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                      ii, boot.R=100, debug=FALSE, fit.mN=TRUE, fit.corr=FALSE) {
  if(missing(data) || missing(startvalues)) {
    stop("data and startvalues must be provided!\n")
  }
  if(missing(bootsamples) && fit.corr) {
    warning("bootstrap samples missing, will not compute correlation matrix!\n")
    fit.corr=FALSE
  }
  if(!missing(bootsamples) && missing(boot.R)) {
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
  np <- 7+2*N
  if(fit.l12) np=np+2
  if(fit.asq) np=np+3
  if(fit.nnlo && fit.kmf) np=np+2
  if(length(startvalues) != np) stop("length of startvalues", length(startvalues),
             "must match number of fit parameters ", np, "!\n")
  if(fsmethod=="cdh") {
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
    chisqr.Na.withr0ZP(par=startvalues, data=data, ii=ii,fsmethod=fsmethod,
                       a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                       fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                       printit=debug)
  }
  parscale <- startvalues
  for(i in 1:length(startvalues)) {
    if(parscale[i] < 0.001) parscale[i] <- 1.
  }
  mini <- optim(par=startvalues, fn=chisqr.Na.withr0ZP, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, parscale=parscale, REPORT=50),
                data=data, ii=ii,
                fsmethod=fsmethod, a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, fit.mN=fit.mN)
  mini <- optim(par=mini$par, fn=chisqr.Na.withr0ZP, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii,
                fsmethod=fsmethod, a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, fit.mN=fit.mN)

  if(mini$convergence != 0) {
    warning("Attention: optim did not converge in initial run!\n Please adjust start values!")
  }
  if(debug) {
    chisqr.Na.withr0ZP(par=mini$par, data=data, ii=ii,fsmethod=fsmethod,
                       a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                       fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                       printit=debug)
    
  }
  dof <- 0
  for(i in 1:length(data)) {
    dof = dof + 2*length( ii[[i]] ) 
    if(fit.mN){
      dof = dof + length(na.omit(data[[i]]$mN[ii[[i]]]))
    }
  }
  dof <- dof + length(r0data$r0) + length(ZPdata$ZP) - length(startvalues)
  if(fit.l12) {
    dof <- dof + 2
  }
  if(fit.kmf) {
    dof <- dof + 2
  }
  if(!fit.mN) {
    dof <- dof + 3
  }
  chisqr <- mini$value
  par <- mini$par

  # compute observables
  mu.phys <- numeric()
  a <- numeric()
  r0 <- numeric()
  l3 <- numeric()
  l4 <- numeric()
  l1 <- 0.
  l2 <- 0.
  F <- numeric()
  B0 <- numeric()
  R0 <- numeric()
  Zp <- numeric()
  mN <- numeric()
  c1 <- numeric()
  gA <- numeric()
  Sigma <- numeric()
  rssq <- numeric()
  s0 <- numeric()
  
  r0mu.phys <- uniroot(fovermps.Na.withr0, c(0.00001, 0.050),
                     tol=1.e-12, par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1., N=N)$root
  r0TwoB <- par[4]
  r0sqTwoBmu <- r0TwoB*r0mu.phys
  r0F <- par[3]
  r0mN <- par[5+2*N]
  mpisq <- getmpssq(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
  fpi <- getfps(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
  mN <- getmN(r0sqTwoBmu, par, N)/fpi*0.1307
  lmpssq <-  log( 0.1396^2 )
  l3 <- log((par[1]/fpi*0.1307)^2) - lmpssq
  l4 <- log((par[2]/fpi*0.1307)^2) - lmpssq
  if(fit.l12) {
    l1 <- log((par[8+2*N]/fpi*0.1307)^2) - lmpssq
    l2 <- log((par[9+2*N]/fpi*0.1307)^2) - lmpssq
  }
  F <- r0F/fpi*0.1307
  mN0 <- r0mN/fpi*0.1307
  c1 <- par[6+2*N]*fpi/0.1307
  gA <- par[7+2*N]
  B0 <- r0TwoB/fpi*0.1307/2.
  r0 <- fpi/0.1307*0.1973
  mu.phys <-  r0mu.phys/fpi*0.1307
  Sigma <- (B0*F^2/2)^(1/3)
  rssq <- 12./(4*pi*F)^2*(l4-12./13.)*0.1973^2
  s0 <- -4.*c1*0.1396^2 - 9*gA^2/32/pi/F^2*0.1396^3
  for(i in 1:N) {
    a[i] <- fpi/0.1307*0.1973/par[4+i]
    Zp[i] <- par[4+i+N]
  }
  
  boot.result <- NULL
  boots <- NULL
  if(method != "no") {
    if(missing(r0bootsamples)) {
      r0bootsamples <- array(0., dim=c(boot.R,N))
      for(i in 1:N) {
        r0bootsamples[,i] <- rnorm(boot.R, mean=r0data$r0[i], sd=r0data$dr0[i])
      }
    }

    if(missing(ZPbootsamples)) {
      ZPbootsamples <- array(0., dim=c(boot.R,N))
      for(i in 1:N) {
        ZPbootsamples[,i] <- rnorm(boot.R, mean=ZPdata$ZP[i], sd=ZPdata$dZP[i])
      }
    }

    boots <- array(0., dim=c(boot.R, (2*N+length(startvalues)+15)))
    for(s in 1:boot.R) {
      if(missing(bootsamples)) {
        bmps=rnorm(length(data[[1]]$mps), mean=data[[1]]$mps, sd=data[[1]]$dmps)
        bfps=rnorm(length(data[[1]]$fps), mean=data[[1]]$fps, sd=data[[1]]$dfps)
      }
      else {
        bmps=bootsamples[[1]][s, 1,]
        bfps=bootsamples[[1]][s, 2,]
      }
      bmN=rnorm(length(data[[1]]$mN), mean=data[[1]]$mN, sd=data[[1]]$dmN)
      df <- list(data.frame(mu=data[[1]]$mu, mps=bmps, dmps=data[[1]]$dmps,
                            fps=bfps, dfps=data[[1]]$dfps, L=data[[1]]$L,
                            mN=bmN, dmN=data[[1]]$dmN))
      r0df <- data.frame(r0=r0bootsamples[s,], dr0=r0data$dr0[1:N])
      ZPdf <- data.frame(ZP=ZPbootsamples[s,], dZP=ZPdata$dZP[1:N])
      i <- 2
      while(i <= N) {
        if(missing(bootsamples)) {
          bmps=rnorm(length(data[[i]]$mps), mean=data[[i]]$mps, sd=data[[i]]$dmps)
          bfps=rnorm(length(data[[i]]$fps), mean=data[[i]]$fps, sd=data[[i]]$dfps)
        }
        else {
          bmps=bootsamples[[i]][s, 1,]
          bfps=bootsamples[[i]][s, 2,]
        }
        bmN=rnorm(length(data[[i]]$mN), mean=data[[i]]$mN, sd=data[[i]]$dmN)
        df[[i]] <- data.frame(mu=data[[i]]$mu, mps=bmps, dmps=data[[i]]$dmps,
                              fps=bfps, dfps=data[[i]]$dfps, L=data[[i]]$L,
                              mN=bmN, dmN=data[[i]]$dmN)
        i <- i+1
      }

      mini.boot <- optim(par=par, fn=chisqr.Na.withr0ZP, method="BFGS", hessian=FALSE,
                         control=list(maxit=500, trace=1, parscale=par, REPORT=50),
                         data=df, ii=ii,
                         fsmethod=fsmethod, a.guess=a.guess, r0data=r0df, ZPdata=ZPdf,
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
        boots[s,1] <- uniroot(fovermps.Na.withr0, c(0.0001, 0.030),
                              tol=1.e-12, par=par.boot, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1., N=N)$root
        
        r0sqTwoBmu <- par.boot[4]*boots[s,1]
        r0F <- par.boot[3]
        mpssqr <- getmpssq(r0sqTwoBmu, par.boot, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
        fpi <- getfps(r0sqTwoBmu, par.boot, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
        boots[s,1] <- boots[s,1]/fpi*0.1307
        for(i in 1:N) {
          boots[s,(i+N)] <- fpi/0.1307*0.1973/par.boot[4+i]
        }
        lmpssqr <- log( 0.1396^2 )
        #l3 l4
        boots[s,(1+2*N)] <- log((par.boot[1]/fpi*0.1307)^2) - lmpssqr
        boots[s,(2+2*N)] <- log((par.boot[2]/fpi*0.1307)^2) - lmpssqr
        #l1 l2
        if(fit.l12) {
          boots[s,(3+2*N)] <- log((par.boot[8+2*N]/fpi*0.1307)^2) - lmpssqr
          boots[s,(4+2*N)] <- log((par.boot[9+2*N]/fpi*0.1307)^2) - lmpssqr
        }
        #f0
        boots[s,(5+2*N)] <- par.boot[3]/fpi*0.1307
        #2B0
        boots[s,(6+2*N)] <- par.boot[4]/fpi*0.1307/2.
        #mN0
        boots[s,(7+2*N)] <- par.boot[5+2*N]/fpi*0.1307
        #mN
        boots[s,(8+2*N)] <- getmN(r0sqTwoBmu, par.boot, N)/fpi*0.1307
        #c1
        boots[s,(9+2*N)] <- par.boot[6+2*N]*fpi/0.1307
        #sigma
        boots[s,(10+2*N)] <- (boots[s,(6+2*N)]*boots[s,(5+2*N)]^2/2)^(1/3)
        #<r^2>
        boots[s,(11+2*N)] <- 12./(4*pi*boots[s,(5+2*N)])^2*(boots[s,(2+2*N)]-12./13.)*0.1973^2
        #r0
        boots[s,(12+2*N)] <- fpi/0.1307*0.1973
        # sigma(0)
        boots[s,(13+2*N)] <- -4.*boots[s,(9+2*N)]*0.1396^2 - 9*par[7+2*N]^2/32/pi/boots[s,(5+2*N)]^2*0.1396^3
        for(i in 1:length(par.boot)) {
          boots[s,(2*N+13+i)] <- par.boot[i]
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
    r0bootsamples <- NULL
    ZPbootsamples <- NULL
  }
  
  result <- list(par=par, result=list(mu.phys=mu.phys, F=F, a=a, l1=l1, l2=l2, l3=l3, l4=l4, r0=r0,
                            B0=B0, mN0=mN0, mN=mN, c1=c1, gA=gA, ZP=Zp, Sigma=Sigma, rssq=rssq, s0=s0,
                            chisqr=chisqr, dof=dof), fit=mini,
                 data=data, boot.result=boot.result, boots=boots, r0data=r0data, ZPdata=ZPdata,
                 bootsamples=bootsamples, r0bootsamples=r0bootsamples, ZPbootsamples=ZPbootsamples,
                 ii=ii, fit.l12=fit.l12, boot.R=boot.R, fsmethod=fsmethod, fit.asq=fit.asq,
                 fit.kmf=fit.kmf, fit.nnlo=fit.nnlo, fit.mN=fit.mN, method=method, fit.corr=fit.corr)
  attr(result, "class") <- c("chiralfit", "list")  
  return(invisible(result))
}

chisqr.Na.withr0ZP <- function(par, data, ii, r0data, ZPdata, fsmethod="gl", a.guess,
                               fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                               printit=FALSE, fit.mN=TRUE) {

  fit.a <- -1.
  N <- length(ii)
  chisum <- 0.
  r0TwoB <- par[4]
  r0F <- par[3]
  fpi <- numeric()
  if(fsmethod=="cdh") {
    r0mu.phys <- try(uniroot(fovermps.Na.withr0,
                             c(0.0001, 0.030), tol=1.e-12,
                             par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1.,
                             fit.kmf=fit.kmf, N=N)$root, silent=T)
#    fpi <- getfps(r0sqTwoBmu=r0TwoB*r0mu.phys, par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=FALSE)
  }
  for( i in 1:N) {
    if(fit.asq) {
      fit.a <- i
    }
    ij <- ii[[i]]
    if(any(r0TwoB <= 0)) {
      return(invisible(NaN))
    }
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+N+i]
    mpssq <- getmpssq(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
    fpsV <- getfps(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)

    if(any(is.nan(mpssq)) || any(mpssq <= 0)) {
      return(NaN)
    }
    if(fsmethod=="cdh") {
      if(inherits(r0mu.phys, "try-error") || is.nan(r0mu.phys)) a_fm <- a.guess[i]
      else {
        a_fm <- getfps(r0sqTwoBmu=r0TwoB*r0mu.phys, par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=FALSE)/0.1307*0.1973/par[4+i]
      }
      if(fit.l12) {
        aLamb1=par[8+2*N]/par[4+i]
        aLamb2=par[9+2*N]/par[4+i]
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*a_fm/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*a_fm/0.1973)^2)
      }

      res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                 aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i], afpiV=fpsV/par[4+i],
                 aF0=fpsV/par[4+i], a_fm=a_fm, L=data[[i]]$L[ij], rev=1, printit=FALSE)
      mpsV <- res$mpiFV
      fpsV <- res$fpiFV
    }
    else {
      mpsv <- data[[i]]$L[ij]*sqrt(mpssq)/par[4+i]
      r <-  mpssq/(4.0*pi*par[3])^2*g1(mpsv)/par[4+i]

      mpsV <- sqrt(mpssq)*(1.+0.5*r)/par[4+i]
      fpsV <- fpsV*(1.0-2.0*r)/par[4+i]
    }
    chisum <- chisum + (sum(((data[[i]]$mps[ij]-mpsV)/(data[[i]]$dmps[ij]))^2) +
                        sum(((data[[i]]$fps[ij]-fpsV)/(data[[i]]$dfps[ij]))^2))
    if(fit.l12 && i==1) {
      chisum <- chisum + ((-0.4-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/0.6)^2 +
        ((4.3-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/0.1)^2
    }
    if(fit.kmf && i==1) {
      chisum <- chisum + (par[10+2*N])^2 + (par[11+2*N])^2
    }
# the nucleon
    if(fit.mN) {
      mN <- getmN(r0sqTwoBmu=r0sqTwoBmu, par, N, fit.asq=fit.a)/par[4+i]
      chisum <- chisum + sum(((data[[i]]$mN[ij]-mN)/data[[i]]$dmN[ij])^2, na.rm=TRUE)
    }

    if(printit) {
      cat("mu-val  ", data[[i]]$mu[ij], "\n")
      cat("r0model ", par[4+i], "\n")
      cat("r0data  ", r0data$r0[i], "\n")
      cat("chir0   ", (r0data$r0[i]-par[4+i])/r0data$dr0[i], "\n\n")
      cat("ZPmodel ", par[4+N+i], "\n")
      cat("ZPdata  ", ZPdata$ZP[i], "\n")
      cat("chiZP   ", (ZPdata$ZP[i]-par[4+N+i])/ZPdata$dZP[i], "\n\n")
      cat("mpsmodel", mpsV, "\n")
      cat("mpsdata ", data[[i]]$mps[ij], "\n")
      cat("chimps  ", ((data[[i]]$mps[ij]-mpsV)/data[[i]]$dmps[ij]), "\n\n")
      cat("fpsmodel", fpsV, "\n")
      cat("fpsdata ", data[[i]]$fps[ij], "\n")
      cat("chif    ", ((data[[i]]$fps[ij]-fpsV)/data[[i]]$dfps), "\n\n")
      if(fit.mN) {
        cat("mNmodel ", mN, "\n")
        cat("mNdata  ", data[[i]]$mN[ij], "\n")
        cat("chimN   ", ((data[[i]]$mN[ij]-mN)/data[[i]]$dmN), "\n\n")
      }
    }

  }
  chisum <- chisum + sum(((par[(5):(4+N)]-r0data$r0)/r0data$dr0)^2) +
    sum(((par[(5+N):(4+2*N)]-ZPdata$ZP)/ZPdata$dZP)^2)
  if(printit) {
    cat("Total Chisqr ", chisum, "\n\n")
  }
  return(invisible(chisum))
}

getchi.Na.withr0ZP <- function(par, data, ii, r0data, ZPdata, fsmethod="gl", a, fit.l12=FALSE) {


  N <- length(ii)
  chisum <- 0.
  for( i in 1:N) {
    a_fm = a[i]
    ij <- ii[[i]]
    r0TwoB <- par[4]
    if(any(r0TwoB < 0)) {
      return(invisible(NaN))
    }
#    r0 <- par[4+i]
#    ZP <- par[4+N+i]
    r0F <- par[3]
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+N+i]
    mpssq <- getmpssq(r0sqTwoBmu, par, N, fit.nnlo=fit.l12)
    fpsV <- getfps(r0sqTwoBmu, par, N, fit.nnlo=fit.l12)
    if(any(is.nan(mpssq)) || any(mpssq < 0)) {
      return(NaN)
    }
    if(fsmethod=="cdh") {
      if(fit.l12) {
        aLamb1=par[8+2*N]/par[4+i]
        aLamb2=par[9+2*N]/par[4+i]
      }
      else {
        aLamb1=sqrt(exp(-0.4)*(0.1396*a_fm/0.1973)^2)
        aLamb2=sqrt(exp(4.3)*(0.1396*a_fm/0.1973)^2)
      }

      res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                 aLamb4=par[2]/par[4+i], ampiV=sqrt(mpssq)/par[4+i], afpiV=fpsV/par[4+i],
                 aF0=fpsV/par[4+i], a_fm=a_fm, L=data[[i]]$L[ij], rev=1, printit=FALSE)
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
    cat("r0:\n")
    print(data.frame(data=r0data$r0[i],Err=r0data$dr0[i], Fit=par[4+i], Chi=(r0data$r0[i]-par[4+i])/r0data$dr0[i]))
    cat("ZP:\n")
    print(data.frame(data=ZPdata$ZP[i],Err=ZPdata$dZP[i], Fit=par[4+N+i], Chi=(ZPdata$ZP[i]-par[4+N+i])/ZPdata$dZP[i]))
    
# the nucleon
    mN <- getmN(r0sqTwoBmu=r0sqTwoBmu, par, N)/par[4+i]
    cat("mN:\n")
    print(data.frame(data=data[[i]]$mN, Err=data[[i]]$dmN, Fit=mN, Chi=(data[[i]]$mN-mN)/data[[i]]$dmN))
  }
}

getmpssq <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1) {

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
    rln1 <- log(r0sqTwoBmu/par[2*N+8]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+9]^2)
    np <- ifelse((2*N+10)<=npar, 2*N+10, npar)
    return(r0sqTwoBmu*(1. + xi*rln3 + asq/par[4+fit.asq]^2*par[npar-2] +
                       17./2.*xi^2*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)^2
                                    +fit.k*8./17.*par[np]
                                    )
                       )
           )
  }
  return(r0sqTwoBmu*(1.0 + asq/par[4+fit.asq]^2*par[npar-2] + r0sqTwoBmu*(rln3/(4.0*pi*par[3])^2 )))
}

getfps <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1.) {

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
    rln1 <- log(r0sqTwoBmu/par[2*N+8]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+9]^2)
    np <- ifelse((2*N+11)<=npar, 2*N+11, npar)
    return(par[3]*(1. - 2.*xi*rln4 +  asq/par[4+fit.asq]^2*par[npar-1] -
                   5*xi^2*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)^2
                           + fit.k*4./5.*par[np]
                           )
                   )
           )
  }
  return(par[3]*(1.0 + asq/par[4+fit.asq]^2*par[npar-1]  - 2.0*r0sqTwoBmu*(rln4)/(4.0*pi*par[3])^2 ))
}

getmN <- function(r0sqTwoBmu, par, N, fit.asq=-1.) {
  npar <- length(par)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
  }
  return(par[5+2*N]*(1 + asq/par[4+fit.asq]^2*par[npar]) - 4.*par[6+2*N]*r0sqTwoBmu - 6.*par[7+2*N]^2/(32*pi*par[3]^2)*r0sqTwoBmu^(3/2))
}
