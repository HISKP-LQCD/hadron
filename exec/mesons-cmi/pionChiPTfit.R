fovermps.pion <- function(x, par, N, fit.nnlo=FALSE, fit.kmf=fit.kmf, fit.asq=-1.) {
  ##r0sqTwoBmu <- par[4]*x
  ##mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  ##fpsV <- getfps.pion(r0sqTwoBmu, par=par, N=N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  ##mpsV <- sqrt(mpssq)
  ##if(any(mpssq <= 0) || any(is.nan(mpssq))) return(NaN)
  lres <- getfpsmps(par[4]*x, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.asq)
  return(lres$fps/lres$mps - (130.7/135))
  ##return(fpsV/mpsV - (130.7/135))
}


# order of fit parameters (N no of lattice spacings)
# 1:              r0*L3
# 2:              r0*L4
# 3:              r0*f0
# 4:              2*r0*B0
# 5 to 4+N:       r0/a chiral
# 5+N to 4+2*N    sr0 (slope of r0 wrt mu^2 )
# 5+2*N to 4+3*N: ZP
# 4+3*N+1:        r0*L1
# 4+3*N+2:        r0*L2
# 4+3*N+3:        k_M
# 4+3*N+4:        k_F
# npar-1:         D_m^0
# npar:           D_f^0

pionChiPTfit <- function(data, startvalues_, bootsamples, fsmethod="gl", a.guess,
                         ZPdata, method="no", seed=123456, r0exp=2,
                         priors = list(l1=-0.4, dl1=0.6, l2=4.3, dl2=0.1, kM=0., dkM=10., kF=0., dkF=10.),
                         fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                         fit.corr=FALSE, ii, boot.R=100, debug=FALSE) {
  
  if(missing(data) || missing(startvalues_) || missing(ZPdata)) {
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
    if(N > 1) for (i in 2:N) {
      ii <- c(ii, list(c(1:length(data[[i]]$mu))))
    }
  }
  else {
    N <- length(ii)
  }
  np <- 6+2*N
  if(fit.l12) np=np+2
  if(fit.asq) np=np+2
  if(fit.nnlo && fit.kmf) np=np+2
  if(length(startvalues_) < np) stop("length of startvalues", length(startvalues_),
             "must be larger than number of fit parameters ", np, "!\n")
  startvalues <- startvalues_[1:np]
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
    chisqr.piononly(par=startvalues, data=data, ii=ii, fsmethod=fsmethod,
                    a.guess=a.guess, ZPdata=ZPdata, priors=priors, r0exp=r0exp,
                    fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, cm=corrmatrix,
                    printit=debug)
  }
  parscale <- startvalues
  for(i in 1:length(startvalues)) {
    if(parscale[i] < 0.001) parscale[i] <- 1.
  }
  mini <- optim(par=startvalues, fn=chisqr.piononly, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, parscale=parscale, REPORT=50),
                data=data, ii=ii, priors=priors,  r0exp=r0exp,
                fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                fit.mN=fit.mN, cm=corrmatrix)
  mini <- optim(par=mini$par, fn=chisqr.piononly, method="BFGS", hessian=TRUE,
                control=list(maxit=500, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, priors=priors,  r0exp=r0exp,
                fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                fit.mN=fit.mN, cm=corrmatrix)
  mini <- optim(par=mini$par, fn=chisqr.piononly, method="BFGS", hessian=TRUE,
                control=list(maxit=500, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, priors=priors,  r0exp=r0exp,
                fsmethod=fsmethod, a.guess=a.guess, ZPdata=ZPdata,
                fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf,
                fit.mN=fit.mN, cm=corrmatrix)
  if(mini$convergence != 0) {
    warning("Attention: optim did not converge in initial run!\n Please adjust initial parameter values!")
  }
  if(debug) {
    chisqr.piononly(par=mini$par, data=data, ii=ii,fsmethod=fsmethod,
                    a.guess=a.guess, ZPdata=ZPdata, priors=priors, r0exp=r0exp,
                    fit.nnlo=fit.nnlo, fit.l12=fit.l12, fit.asq=fit.asq, fit.kmf=fit.kmf, cm=corrmatrix,
                    printit=debug)
  }
  dof <- 0

  for(i in 1:length(data)) {
    dof <- dof + length( na.omit(data[[i]]$mps[ii[[i]]])) +
      length( na.omit( data[[i]]$fps[ii[[i]]] ) ) +
        length( na.omit( data[[i]]$r0a ) )
  }
  dof <- dof + N - length(startvalues)
  if(fit.l12) {
    dof <- dof + 2
  }
  if(fit.kmf) {
    dof <- dof + 2
  }
  chisqr <- mini$value
  par <- mini$par
  
  ## compute observables
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
  sr1 <- numeric()
  sr2 <- numeric()
  Sigma <- numeric()
  rssq <- numeric()
  
  r0mu.phys <- uniroot(fovermps.pion, c(0.001, 0.012),
                     tol=1.e-12, par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1., N=N)$root
  r0TwoB <- par[4]
  r0sqTwoBmu <- r0TwoB*r0mu.phys
  r0F <- par[3]
  fpi <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf)
  lmpssq <-  log( 0.1396^2 )
  l3 <- log((par[1]/fpi*0.1307)^2) - lmpssq
  l4 <- log((par[2]/fpi*0.1307)^2) - lmpssq
  if(fit.l12) {
    l1 <- log((par[7+2*N]/fpi*0.1307)^2) - lmpssq
    l2 <- log((par[8+2*N]/fpi*0.1307)^2) - lmpssq
  }
  F <- r0F/fpi*0.1307
  B0 <- r0TwoB/fpi*0.1307/2.
  r0 <- fpi/0.1307*0.1973
  mu.phys <-  r0mu.phys/fpi*0.1307
  Sigma <- (B0*F^2/2)^(1/3)
  rssq <- 12./(4*pi*F)^2*(l4-13./12.)*0.1973^2
  sr1 <- par[5+N]
  sr2 <- par[6+N]
  for(i in 1:N) {
    a[i] <- fpi/0.1307*0.1973/par[4+i]
    Zp[i] <- par[6+i+N]
  }
  
  boot.result <- NULL
  boots <- NULL
  ZPbootsamples <- NULL
  if(method != "no") {
    set.seed(seed)
    ZPbootsamples <- array(0., dim=c(boot.R,N))
    for(i in 1:N) {
      ZPbootsamples[,i] <- rnorm(boot.R, mean=ZPdata$ZP[i], sd=ZPdata$dZP[i])
    }

    boots <- array(0., dim=c(boot.R, (3*N+length(startvalues)+10)))
    for(s in 1:boot.R) {
      df <- list(data.frame(mu=data[[1]]$mu, mps=bootsamples[[1]][s, 1,],
                            dmps=data[[1]]$dmps,
                            fps=bootsamples[[1]][s, 2,],
                            dfps=data[[1]]$dfps, L=data[[1]]$L,
                            r0a=rnorm(length(data[[1]]$r0a), mean=data[[1]]$r0a, sd=data[[1]]$dr0a),
                            dr0a=data[[1]]$dr0a))
      ZPdf <- data.frame(ZP=ZPbootsamples[s,], dZP=ZPdata$dZP[1:N])
      if(N > 1) for(i in 2:N) {
        df[[i]] <- data.frame(mu=data[[i]]$mu, mps=bootsamples[[i]][s,1,],
                              dmps=data[[i]]$dmps,
                              fps=bootsamples[[i]][s,2,],
                              dfps=data[[i]]$dfps, L=data[[i]]$L,
                              r0a=rnorm(length(data[[i]]$r0a), mean=data[[i]]$r0a, sd=data[[i]]$dr0a),
                              dr0a=data[[i]]$dr0a)
      }
      boot.priors <- list(l1=rnorm(1, mean=priors$l1, sd=priors$dl1), dl1=priors$dl1,
                          l2=rnorm(1, mean=priors$l2, sd=priors$dl2), dl2=priors$dl2,
                          kM=rnorm(1, mean=priors$kM, sd=priors$dkM), dkM=priors$dkM,
                          kF=rnorm(1, mean=priors$kF, sd=priors$dkF), dkF=priors$dkF)
      
      mini.boot <- optim(par=par, fn=chisqr.piononly, method="BFGS", hessian=FALSE,
                         control=list(maxit=150, parscale=par, trace=0, REPORT=100),
                         data=df, ii=ii, priors=boot.priors,  r0exp=r0exp,
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
          ## a in fm for each beta
          boots[s,(i+N)] <- fpi/0.1307*0.1973/par.boot[4+i]
          ## 
          boots[s,(i+2*N)] <- -1.
        }
        lmpssqr <- log( 0.1396^2 )
        # l3 and l4
        boots[s,(1+3*N)] <- log((par.boot[1]/fpi*0.1307)^2) - lmpssqr
        boots[s,(2+3*N)] <- log((par.boot[2]/fpi*0.1307)^2) - lmpssqr
        if(fit.l12) {
          # l1 and l2
          boots[s,(3+3*N)] <- log((par.boot[7+2*N]/fpi*0.1307)^2) - lmpssqr
          boots[s,(4+3*N)] <- log((par.boot[8+2*N]/fpi*0.1307)^2) - lmpssqr
        }
        # f0 and B0
        boots[s,(5+3*N)] <- par.boot[3]/fpi*0.1307
        boots[s,(6+3*N)] <- par.boot[4]/fpi*0.1307/2.
        # r0
        boots[s,(7+3*N)] <- fpi/0.1307*0.1973
        # <r^2>_s
        boots[s,(8+3*N)] <- 12./(4*pi*boots[s,(5+3*N)])^2*(boots[s,(2+3*N)]-12./13.)*0.1973^2
        # sigma (B0*F^2/2)^(1/3)
        boots[s,(9+3*N)] <- (boots[s,(6+3*N)]*(boots[s,(5+3*N)])^2/2)^(1/3)
        for(i in 1:length(par.boot)) {
          boots[s,(3*N+9+i)] <- par.boot[i]
        }
        boots[s,(3*N+9+length(par.boot)+1)] <- mini.boot$value
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
  
  result <- list(par=par, result=list(mu.phys=mu.phys, F=F, a=a, l1=l1, l2=l2, l3=l3, l4=l4, r0=r0,
                            B0=B0, ZP=Zp, Sigma=Sigma, rssq=rssq, sr1=sr1, sr2=sr2,
                            chisqr=chisqr, dof=dof), fit=mini,
                 data=data, boot.result=boot.result, boots=boots, ZPdata=ZPdata, seed=seed,
                 fmcorrmatrix=corrmatrix, bootsamples=bootsamples, ZPbootsamples=ZPbootsamples,
                 method=method, ii=ii, fit.l12=fit.l12, boot.R=boot.R, fsmethod=fsmethod, r0exp=r0exp,
                 fit.asq=fit.asq, fit.kmf=fit.kmf, fit.nnlo=fit.nnlo, fit.mN=FALSE, fit.corr=fit.corr)
  attr(result, "class") <- c("pionChiPTfit", "list")  
  return(invisible(result))
}

chisqr.piononly <- function(par, data, ii, ZPdata, fsmethod="gl", a.guess, r0exp=2,
                            priors = list(l1=-0.4, dl1=0.6, l2=4.3, dl2=0.1, kM=0, dkM=10., kF=0., dkF=10.),
                            fit.nnlo=FALSE, fit.l12=FALSE, fit.asq=FALSE, fit.kmf=FALSE,
                            printit=FALSE, fit.mN=TRUE, cm=NULL) {

  fit.a <- -1.
  N <- length(ii)
  chisum <- 0.
  r0TwoB <- par[4]
  r0F <- par[3]
  fpi <- numeric()
  if(any(par[c(4:(4+N), (6+N+1):(6+2*N))] <= 0)) {
    return(invisible(NaN))
  }

  ## r0mu.phys is the renormalised quark mass tol=1.e-12,
  r0mu.phys <- try(uniroot(fovermps.pion,
                           interval=c(0.007, 0.01), tol=1.e-8,
                           par=par, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=-1.,
                           N=N)$root, silent=T)
  if(inherits(r0mu.phys, "try-error") || is.nan(r0mu.phys)) {
    ##return(invisible(NaN))
    r0mu.phys <- 0.009
  }
  r0fm <- getfps.pion(r0sqTwoBmu=r0TwoB*r0mu.phys, par=par, N=N, fit.nnlo=fit.nnlo,
                      fit.kmf=fit.kmf, fit.asq=FALSE)/0.1307*0.1973
  ## r0 Lambda = 0.62 from hep-lat/0411025
  ## alphas at 2 GeV
  ##alpha2GeV <- alphas(mu = 2.0, nl = 3, lam0 = 0.62/r0fm*0.1973)
  ## with lam0 = 0.2567068 GeV from arXiv:0804.3383 [hep-lat] 
  alpha2GeV <- 0.01914803
  for( i in 1:N) {
    if(fit.asq) {
      fit.a <- i
    }
    ij <- ii[[i]]
    a_fm <- r0fm/par[4+i]
    
    ## fit r0/a as a function of (a mu)^2 first
    ## r0 = r0(mu_pi) + Dr0*(mu_q^\gamma - (Z_P*r0mu_pi/r0)^\gamma)
    ## r0 <- par[4+i]  + (data[[i]]$mu^r0exp - (par[4+2*N+i]*r0mu.phys/par[4+i])^r0exp)*par[4+N+i]
    ##r0 <- par[4+i] + par[4+N+i]*data[[i]]$mu^r0exp
    r0 <- par[4+i] *(1 + par[5+N]*data[[i]]$mu*par[4+i]/par[6+N+i] + par[6+N]*(data[[i]]$mu*par[4+i]/par[6+N+i])^2 )

    chisum <- chisum + sum(((data[[i]]$r0a-r0)/data[[i]]$dr0)^2, na.rm=TRUE)

    ## now the rest
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[6+N+i]
    ##lres <- getfpsmps(r0sqTwoBmu,
    ##                  par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
    ## in principle first fit.asq=FALSE and do it after FS correction
    lres <- getfpsmps(r0sqTwoBmu,
                      par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=FALSE)
    mpssq <- lres$mps^2
    mpsV <- lres$mps/par[4+i]
    fpsV <- lres$fps/par[4+i]
    if(any(is.nan(mpssq))) {
      return(NaN)
    }

    ##mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
    ##fpsV <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)/par[4+i]
    ##mpsV <- sqrt(mpssq)/par[4+i]
    ##if(any(is.nan(mpssq)) || any(mpssq <= 0)) {
    ##  return(NaN)
    ##}


    if(fsmethod=="cdh" || fsmethod=="cdhnew") {
      if(fit.l12) {
        aLamb1=par[7+2*N]/par[4+i]
        aLamb2=par[8+2*N]/par[4+i]
      }
      else {
        aLamb1=exp(0.5*priors$l1)*(0.1396*a_fm/0.1973)
        aLamb2=exp(0.5*priors$l2)*(0.1396*a_fm/0.1973)
      }

      if(fsmethod=="cdhnew") {
        res <- cdhnew(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                      aLamb4=par[2]/par[4+i], ampiV=mpsV,
                      afpiV=fpsV, aF0=par[3]/par[4+i],
                      a2B0mu=r0sqTwoBmu/par[4+i]^2, L=data[[i]]$L[ij], rev=1,
                      printit=FALSE)
      }
      else {
        res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]/par[4+i],
                   aLamb4=par[2]/par[4+i], ampiV=mpsV, afpiV=fpsV,
                   aF0=fpsV, a_fm=a_fm, L=data[[i]]$L[ij], rev=1, printit=FALSE)
      }
      mpsV <- res$mpiFV
      fpsV <- res$fpiFV
    }
    else {
      r <-  mpssq/(4.0*pi*par[3])^2*g1( data[[i]]$L[ij]*mpsV )/par[4+i]

      mpsV <- mpsV*(1.+0.5*r)
      fpsV <- fpsV*(1.0-2.0*r)
    }
    ## for correct order of FS corrections and a^2 artifacts
    if(fit.asq) {
      npar <- length(par)
      fpsV <- fpsV + par[npar]*par[3]/par[4+i]^3
      mpsV <- sqrt(mpsV^2 + par[npar-1]*r0sqTwoBmu/par[4+i]^4)
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
    if(fit.l12 && i==1) {
      chisum <- chisum + ((priors$l1-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/priors$dl1)^2 +
        ((priors$l2-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/priors$dl2)^2
    }
    if(fit.kmf && i==1) {
      chisum <- chisum + ((priors$kM-par[7+2*N])/priors$dkM)^2
      + ((priors$kF-par[8+2*N])/priors$dkF)^2
    }
    ## Z_P: data is in RI' scheme at 1/a, fit par in MSbar at 2GeV
    ## alphas at 1/a
    ##alpha1ova <- alphas(mu = 1./a_fm*0.1973, nl = 3, lam0 = 0.62/r0fm*0.1973)
    alpha1ova <- alphas(mu = 1./a_fm*0.1973, nl = 3, lam0 = 0.2567068)
    chisum <- chisum + ((ZPdata$ZP[i]-zetazp(0.8178*par[6+N+i], alpha2GeV, alpha1ova))/ZPdata$dZP[i])^2
    if(printit) {
      cat("r0chiral", par[4+i], "\n")
      cat("datar0  ", data[[i]]$r0a, "\n")
      cat("modelr0 ", r0, "\n")
      cat("chir0   ", ((data[[i]]$r0a-r0)/data[[i]]$dr0), "\n")
      cat("sumr0   ", sum(((data[[i]]$r0a-r0)/data[[i]]$dr0)^2, na.rm=TRUE), "\n")
      cat("modelZP ", zetazp(0.8178*par[6+N+i], alpha2GeV, alpha1ova), "\n")
      cat("ZPdata  ", ZPdata$ZP[i], "\n")
      cat("chiZP   ", ((ZPdata$ZP[i]-zetazp(0.8178*par[6+N+i], alpha2GeV, alpha1ova))/ZPdata$dZP[i]), "\n")
      cat("modelm  ", mpsV, "\n")
      cat("datam   ", data[[i]]$mps[ij], "\n")
      cat("errm    ", data[[i]]$dmps[ij], "\n")
      cat("chim    ", ((data[[i]]$mps[ij]-mpsV)/data[[i]]$dmps[ij]), "\n")
      cat("summ    ", sum(((data[[i]]$mps[ij]-mpsV)/data[[i]]$dmps[ij])^2), "\n")
      cat("modelf  ", fpsV, "\n")
      cat("dataf   ", data[[i]]$fps[ij], "\n")
      cat("errf    ", data[[i]]$dfps[ij], "\n")
      cat("chif    ", ((data[[i]]$fps[ij]-fpsV)/data[[i]]$dfps[ij]), "\n")
      cat("sumf    ", sum(((data[[i]]$fps[ij]-fpsV)/data[[i]]$dfps[ij])^2), "\n")
      if(fit.l12 && i==1) {
        cat("l1,2    ", priors$l1, priors$l2, "\n")
        cat("chil    ", ((priors$l1-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/priors$dl1) + ((priors$l2-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/priors$dl2), "\n")
      }
      if(fit.kmf && i==1) {
        cat("kM,F    ",  priors$kM, priors$kF, "\n")
        cat("chik    ", ((priors$kM-par[7+2*N])/priors$dkM) + ((priors$kF-par[8+2*N])/priors$dkF), "\n")
      }
    }

  }
  ### ZP
  ### chisum <- chisum + sum(((par[(5+2*N):(4+3*N)]-ZPdata$ZP[1:N])/ZPdata$dZP[1:N])^2)
  if(printit) {
    cat("chisqr ", chisum, "\n\n")
  }
  return(invisible(chisum))
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
    rln1 <- log(r0sqTwoBmu/par[7+2*N]^2)
    rln2 <- log(r0sqTwoBmu/par[8+2*N]^2)
    # this line makes it work even if we do not fit k_M,F
    # however, fit.k=0. then...
    np <- ifelse((2*N+9)<=npar, 2*N+9, npar)
    return(r0sqTwoBmu*(1. + xi*rln3 + asq/par[4+fit.asq]^2*par[npar-1] +
                       17./2.*xi^2*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)^2
                                    +fit.k*8./17.*par[np]
                                    )
                       )
           )
  }
  return(r0sqTwoBmu*(1.0 + asq/par[4+fit.asq]^2*par[npar-1] + r0sqTwoBmu*(rln3/(4.0*pi*par[3])^2 )))
}

getfpsmps <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1.) {
  dl <- length(r0sqTwoBmu)
  res <- .Call("getbothc", r0sqTwoBmu, par, N, fit.nnlo, fit.kmf, fit.asq)
  return(invisible(list(mps=res[(dl+1):(2*dl)], fps=res[1:dl])))
}

getfpsmps.R <- function(r0sqTwoBmu, par, N, fit.nnlo=FALSE, fit.kmf=FALSE, fit.asq=-1.) {

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
    rln1 <- log(r0sqTwoBmu/par[2*N+7]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+8]^2)
    np <- ifelse((2*N+10)<=npar, 2*N+10, npar)
    fps <- par[3]*(1. - 2.*xi*rln4 +  asq/par[4+fit.asq]^2*par[npar] -
                   5*xi^2*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)^2
                           + fit.k*4./5.*par[np]
                           )
                   )
    np <- ifelse((2*N+9)<=npar, 2*N+9, npar)
    mpssq <- r0sqTwoBmu*(1. + xi*rln3 + asq/par[4+fit.asq]^2*par[npar-1] +
                         17./2.*xi^2*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)^2
                                      +fit.k*8./17.*par[np]
                                      )
                         )
            
  }
  else {
    fps <- par[3]*(1.0 + asq/par[4+fit.asq]^2*par[npar]  - 2.0*r0sqTwoBmu*(rln4)/(4.0*pi*par[3])^2 )
    mpssq <- r0sqTwoBmu*(1.0 + asq/par[4+fit.asq]^2*par[npar-1] + r0sqTwoBmu*(rln3/(4.0*pi*par[3])^2 ))
  }
  return(invisible(list(mps=sqrt(mpssq), fps=fps)))
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
    rln1 <- log(r0sqTwoBmu/par[2*N+7]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+8]^2)
    np <- ifelse((2*N+10)<=npar, 2*N+10, npar)
    return(par[3]*(1. - 2.*xi*rln4 +  asq/par[4+fit.asq]^2*par[npar] -
                   5*xi^2*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)^2
                           + fit.k*4./5.*par[np]
                           )
                   )
           )
  }
  return(par[3]*(1.0 + asq/par[4+fit.asq]^2*par[npar]  - 2.0*r0sqTwoBmu*(rln4)/(4.0*pi*par[3])^2 ))
}

anova.pionChiPTfit <- function(fit1, fit2) {
  Fvalue <- ((fit1$result$chisqr-fit2$result$chisqr)/(fit1$result$dof-fit2$result$dof))/(fit2$result$chisqr/fit2$result$dof)
  return(list(Fvalue=Fvalue, pvalue=1-pf(Fvalue, fit1$result$dof-fit2$result$dof, fit2$result$dof)))
}

average.pionChiPTfit <- function(list.fits, av.weight=TRUE) {
  e1 <- new.env()
  nl <- length(list.fits)

  ## first we copute the weights for each fit in list.fits
  ## we also build a list of all the fits -> fitlist
  weights <- rep(1., times=nl)
  i <- 1
  load(list.fits[i], e1)
  fit.cur <- get(ls(e1), e1)
  N <- length(fit.cur$data)
  boot.R <- fit.cur$boot.R
  if(av.weight) weights[i] <- 1-pgamma(fit.cur$result$chisqr/2, fit.cur$result$dof/2)
  fitlist <- list(fit.cur)
  rm(list=ls(e1), envir=e1)
  for(i in 2:nl) {
    load(list.fits[i], e1)
    fit.cur <- get(ls(e1), e1)
    if(boot.R > fit.cur$boot.R) {
      boot.R <- fit.cur$boot.R
    }
    if(av.weight) weights[i] <- 1-pgamma(fit.cur$result$chisqr/2, fit.cur$result$dof/2)
                                        #    weights[i] <- 1./pgamma(fit.cur$result$chisqr/2, fit.cur$result$dof/2)
    fitlist[[i]] <- fit.cur
    rm(list=ls(e1), envir=e1)
  }

  ## here we create the list of variables to be averaged
  ## their indices as well as their names
  ii <- c(1, (3*N+1):(3*N+7), (9+3*N+5+N), (9+3*N+6+N), 3*N+8, 
	2*N+10, 1, (9+3*N+5):(9+3*N+4+N), (N+1):(2*N), (9+3*N+6+N+1):(9+3*N+6+N+N), 
	NA, NA )
  nlist <- c("m_ud", "l3", "l4", "l1", "l2", "f0", "B0", "r0", "C1", "C2", "r^2_s", 
	"sigma", "fpif0", rep("r0/a", times=N), rep("a", times=N), rep("ZP", times=N), "Dm", "Df")

  ## bres will hold the results for each bootstrap sample
  bres <- array(0., dim=c(boot.R, length(ii)))
  tmp <- numeric(nl)
  for(j in 1:length(ii)) {
    if(nlist[j] == "Dm" || nlist[j] == "Df") {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          pidx <- 9+3*N+length(fitlist[[k]]$par)
	  if(nlist[j] == "Dm") pidx <- pidx-1
          tmp[k] <- NA
          if(fitlist[[k]]$fit.asq) {
            tmp[k] <- fitlist[[k]]$boots[i, pidx]
          }
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights, na.rm = TRUE)
      }
    }
    else if(nlist[j] == "sigma") {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          tmp[k] <- ((fitlist[[k]]$boots[i, 6+3*N]*fitlist[[k]]$boots[i, 5+3*N]^2)/2)^(1/3)
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights)
      }
    }
    else if(nlist[j] == "r^2_s") {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          tmp[k] <- 12./(4*pi*fitlist[[k]]$boots[i, 5+3*N])^2*(fitlist[[k]]$boots[i, ii[3]]-13./12.)*0.1973^2
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights)
      }
    }
    else if(nlist[j] == "B0") {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          tmp[k] <- fitlist[[k]]$boots[i,(6+3*N)]
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights)
      }
    }
    else if(nlist[j] == "fpif0") {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          tmp[k] <- 0.1307/fitlist[[k]]$boots[i,(5+3*N)]
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights)
      }
    }
    else {
      for(i in 1:boot.R) {
        for(k in 1:nl) {
          tmp[k] <- fitlist[[k]]$boots[i, ii[j]]
        }
        bres[i, j] <- weighted.median(x = tmp, w=weights)
      }
    }
  }

  jj <- c(1, 6, 7, 4, 5, 2, 9, 8, 5+N, 6+N, 12, 11, 2, 5:(4+N), 1:N, (6+N+1):(6+2*N), 2*N+11, 2*N+12)
  nr <- length(jj)

  res <- numeric(nr)
  histres <- array(0., dim=c(nr, nl))
  
  for(i in 1:nr) {
    if(nlist[i] == "fpif0") {
      for(j in 1:nl) {
        histres[i,j] <- 0.1307/fitlist[[j]]$result[[2]]
      }
    }
    else if(nlist[i] == "r^2_s") {
      for(j in 1:nl) {
        histres[i,j] <- 12./(4*pi*fitlist[[j]]$result[[2]])^2*(fitlist[[j]]$result[[7]]-13./12.)*0.1973^2
      }
    }
    else if(nlist[i] == "r0/a" || nlist[i] == "ZP" || nlist[i] == "C1" || nlist[i] == "C2") {
      for(j in 1:nl) {
        histres[i,j] <- fitlist[[j]]$par[jj[i]]
      }
    }
    else if(nlist[i] == "Dm" || nlist[i] == "Df") {
      for(j in 1:nl) {
        pidx <- length(fitlist[[j]]$par)
	if(nlist[i] == "Dm") pidx <- pidx-1
        histres[i,j] <- NA
        if(fitlist[[j]]$fit.asq) {
          histres[i,j] <- fitlist[[j]]$par[pidx]
	  cat(fitlist[[j]]$par[pidx], pidx, "\n")
        }
      }
    }
    else if(nlist[i] == "a") {
      for(j in 1:nl) {
        histres[i,j] <- fitlist[[j]]$result$a[jj[i]]
      }
    }
    else {
      for(j in 1:nl) {
        histres[i,j] <- fitlist[[j]]$result[[jj[i]]]
      }
    }
    res[i] <- weighted.median(histres[i,], w=weights, na.rm=TRUE)
  }

  kk <- c(1:length(ii))
  for(i in 1:length(ii)) {
    cat(nlist[kk[i]], "\t", signif(res[kk[i]], digits=4), "\t+-", signif(sd(bres[,kk[i]], na.rm=TRUE), digits=4), "\t+", signif(weighted.quantile(histres[kk[i],], w=weights, prob=c(0.8427), na.rm=TRUE)-res[kk[i]], digits=4), "\t-", signif(-weighted.quantile(histres[kk[i],], w=weights, prob=c(0.1573), na.rm=TRUE)+res[kk[i]], digits=4), "\t bias:", res[kk[i]]-mean(bres[,kk[i]], na.rm=TRUE), "\n")
    ##cat(nlist[kk[i]], "\t", res[kk[i]], "+-", sd(bres[,kk[i]], na.rm=TRUE), "+-", sqrt(cov.wt(data.frame(histres[kk[i],]), weights)$cov[1]), "bias:", res[kk[i]]-mean(bres[,kk[i]], na.rm=TRUE), "\n")
  }
  for(i in 1:length(ii)) {
    if(nlist[kk[i]] == "sigma" || nlist[kk[i]] == "B0" || nlist[kk[i]] == "f0" || nlist[kk[i]] == "m_ud") {
      printtab(res[kk[i]], sd(bres[,kk[i]], na.rm=TRUE), c=1000.)
    }
    else {
      printtab(res[kk[i]], sd(bres[,kk[i]], na.rm=TRUE))
    }
  }
  cat("---\n")
  for(i in 1:length(ii)) {
    if(nlist[kk[i]] == "sigma" || nlist[kk[i]] == "B0" || nlist[kk[i]] == "f0" || nlist[kk[i]] == "m_ud") {
      fac <- 1000.
    }
    else {
      fac <- 1.
    }
    printtab(res[kk[i]], sqrt(var(bres[,kk[i]], na.rm=TRUE) + max((weighted.quantile(histres[kk[i],], w=weights, prob=c(0.8427), na.rm=TRUE)-res[kk[i]]), (-weighted.quantile(histres[kk[i],], w=weights, prob=c(0.1573), na.rm=TRUE)+res[kk[i]]))^2), c=fac)
  }
  
  return(invisible(list(bootres=bres, res=res, histres=histres, weights=weights)))
}


predict.pion <- function(fit, r0, ZP, to=0.2) {
  N <- length(fit$data)
  npar <- length(fit$par)
  par <- fit$par
  xfit <- seq(from=0., to=to, length.out=500)

  fit.a <- -1
  if(fit$fit.asq) {
    fit.a <- 1
    par[4+1] <- r0
  }

  r0TwoB <- par[4]
  r0sqTwoBmu <- r0TwoB*xfit
  msq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit$fit.nnlo, fit.kmf=fit$fit.kmf, fit.asq=fit.a)
  return(data.frame(mu=xfit*ZP/r0, msq=msq/0.43^2*0.198^2, mpi=sqrt(msq)/0.43*0.198))
}
