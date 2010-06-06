## fit parameters
## 1:     f0K/f0
## 2:     fmK/f0
## 3:     b
## 4:     bK
## 5:     bmK
## 6-5+N: af0
## 6+N:   Dfpi
## 7+N:   DfK
## 8+N:   DfKm
##
## data is a double list



OSChiPTfit <- function(data, startvalues_, ii, bootsamples=NULL, cmatrix,
                        boot.R=100, debug=TRUE, fit.a=FALSE) {
  ## number of lattice spacings
  Nbeta <- length(data)
  ## number of points with unitary light mu
  Nens <- c()
  np <- 5 + Nbeta
  dof <- 0
  par <- startvalues_[1:np]
  for(i in 1:Nbeta) {
    Nens <- c(Nens, length(data[[i]]))
    
    for(k in 1:Nens[i]) {
      dof <- dof + length(ii[[i]][[k]]$uii) + length(ii[[i]][[k]]$lsii)
      np <- np + length(ii[[i]][[k]]$uii) + length(ii[[i]][[k]]$ssii)
      par <- c(par, data[[i]][[k]]$m[ii[[i]][[k]]$uii], data[[i]][[k]]$m[ii[[i]][[k]]$ssii])
    }
  }
  dof <- dof-5-Nbeta
  
  ## here we minimise chisqr
  mini <- optim(par=par, fn=chisqr.os, method="BFGS", hessian=FALSE,
                control=list(maxit=100, trace=debug, parscale=par, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, cmatrix=cmatrix, fit.a=fit.a)
  mini <- optim(par=mini$par, fn=chisqr.os, method="BFGS", hessian=FALSE,
                control=list(maxit=100, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, cmatrix=cmatrix, fit.a=fit.a)
  mini <- optim(par=mini$par, fn=chisqr.os, method="BFGS", hessian=FALSE,
                control=list(maxit=500, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, cmatrix=cmatrix, fit.a=fit.a)


  if(mini$convergence != 0) {
    warning("minimisation did not converge")
  }

  ## now the bootstrap
  if(!is.null(bootsamples)) {
    boots <- array(0., dim=c(boot.R, length(mini$par)+1))

    datan <- data
    for(s in 1:boot.R) {
      for(i in 1:Nbeta) {
        for(k in 1:Nens[i]) {
          datan[[i]][[k]]$m <- bootsamples[[i]][[k]][s,1,]
          ##datan[[i]][[k]]$f <- bootsamples[[i]][[k]][s,2,]
          datan[[i]][[k]]$f <- bootsamples[[i]][[k]][s,2,]*bootsamples[[i]][[k]][s,1,]/sinh(bootsamples[[i]][[k]][s,1,])
        }
      }
      ## minimise for each bootstrap sample
      bmini <- optim(par=mini$par, fn=chisqr.os, method="BFGS", hessian=FALSE,
                     control=list(maxit=100, trace=FALSE, parscale=mini$par, REPORT=50),
                     data=datan, ii=ii, Nens=Nens, Nbeta=Nbeta, cmatrix=cmatrix, fit.a=fit.a)
      bmini <- optim(par=bmini$par, fn=chisqr.os, method="BFGS", hessian=FALSE,
                     control=list(maxit=500, trace=debug, parscale=bmini$par, REPORT=50),
                     data=datan, ii=ii, Nens=Nens, Nbeta=Nbeta, cmatrix=cmatrix, fit.a=fit.a)
      boots[s, 1:length(mini$par)] <- bmini$par
      boots[s, 1+length(mini$par)] <- bmini$value
      if(bmini$convergence !=0 ) {
        warning("minimisation during bootstrap for s=", s, "did not converge")
      }
    }
  }
  else {
    boots <- NULL
  }
  
  result <- list(af0=mini$par[6:(5+Nbeta)], mini=mini, data=data, cmatrix=cmatrix, boot.R=boot.R,
                 bootsamples=bootsamples, boots=boots, 
                 dof=dof, ii=ii, startvalues=startvalues_)
  attr(result, "class") <- c("OSChiPTfit", "list")  
  return(invisible(result))
}

summary.OSChiPTfit <- function(fit) {
  Nbeta <- length(fit$data)
  par <- fit$mini$par
  bootres <- array(0, dim=c(fit$boot.R, 4, 2))
  cat("L / dof\t =", fit$mini$value, "/", fit$dof, "\n")
  cat("boot.R =", fit$boot.R, "\n")
  for(i in 1:Nbeta) {
    ampi <- uniroot(fpsovmps.os, interval = c(0.001, 0.12),
                    tol=1.e-12, par=par, i=i, N=Nbeta, fit.a=FALSE)$root
    afpi <- par[5+i] * getfpi.os(par=par, ampi^2, i=i, N=Nbeta, fit.a=FALSE)
    a <- ampi/0.135*0.198
    a <- afpi/0.13070*0.198
    ## neutral Kaon mass
    ## charged would be 0.49368 MeV
    amK <- 0.49765*a/0.198
    amDs = 1.9685*a/0.198
    amD = 1.8696*a/0.198
    afK <- par[5+i]*getfK.os(par=par, mpisq=ampi^2, msssq=(2*amK^2-ampi^2), i=i, N=Nbeta, fit.a=FALSE)
    Vus <- afpi/afK*0.97422*0.27599

    for(s in 1:fit$boot.R) {
      ## a mpi
      bootres[s, 1, i] <- uniroot(fpsovmps.os, interval = c(0.001, 0.12),
                                  tol=1.e-12, par=fit$boots[s,1:length(par)], i=i, N=Nbeta, fit.a=FALSE)$root
      ## a fpi
      bootres[s, 2, i] <- fit$boots[s, 5+i] * getfpi.os(par=fit$boots[s,1:length(par)], bootres[s, 1, i]^2, i=i, N=Nbeta, fit.a=FALSE)
      ## a
      bootres[s, 3, i] <- bootres[s, 1, i]/0.135*0.198
      bootres[s, 3, i] <- bootres[s, 2, i]/0.13070*0.198
      mKK <- 0.49765*bootres[s, 1, i]/0.135
      ## a fK
      bootres[s, 4, i] <- fit$boots[s, 5+i]*getfK.os(par=fit$boots[s,1:length(par)], mpisq=bootres[s, 1, i]^2, msssq=(2*mKK^2-bootres[s, 1, i]^2), i=i, N=Nbeta, fit.a=FALSE)
    }

    cat("*** beta value", i, "***\n")
    cat("ampi\t =", ampi, "+-", sd( bootres[, 1, i]), "\n")
    cat("amK\t =", amK, "+-", sd(0.49765*bootres[, 1, i]/0.135), "\n")
    ##cat("amDs\t =", amDs, "\n")
    ##cat("amD\t =", amD, "\n")
    cat("af0\t =", par[5+i], "+-", sd(fit$boots[,5+i]), "\n")
    cat("afpi\t =", afpi, "+-", sd(bootres[, 2, i]), "\n")
    cat("afK\t =", afK, "+-", sd(bootres[, 4, i]), "\n")
    cat("a\t =", a, "+-", sd(bootres[, 3, i]), "fm \n")
  }
  cat("*** independent results ***\n")
  cat("f0\t =", par[5+i]/a*0.198, "+-", sd(fit$boots[,5+i]/bootres[, 3, i])*0.198, "GeV\n")
  cat("fK\t =", afK/a*0.198, "+-", sd(bootres[, 4, i]/bootres[, 3, i])*0.198, "GeV\n")
  cat("fK/fpi\t =", afK / afpi, "+-", sd(bootres[, 4, i]/bootres[, 2, i]), "\n")
  cat("|Vus|\t =", Vus, "+-", sd(bootres[, 2, i]/bootres[, 4, i]*0.97422*0.27599), "\n")
  cat("l4\t =", par[3]/2. + 2*log(4*pi*par[5+i]/a*0.198/0.1396), "+-",
      sd(fit$boots[3]/2. + 2*log(4.*pi*fit$boots[5+i]/bootres[, 3, i]*0.198/0.1396)), "\n")
  return(invisible(list(bootres)))
}

chisqr.os <- function(par, data, ii, Nbeta, Nens, cmatrix, debug=FALSE, fit.a=FALSE) {

  chisum <- 0.
  npar <- 5+Nbeta
  for(i in 1:Nbeta) {
    for(k in 1:Nens[i]) {
      L <- data[[i]][[k]]$L[1]
      uij <- ii[[i]][[k]]$uii
      if(length(uij) != 1) {
        error("length(uij) != 1")
      }
      lsij <- ii[[i]][[k]]$lsii
      ssij <- ii[[i]][[k]]$ssii
      ## mpisq <- data[[i]][[k]]$m[uij]^2
      mpisq <- par[npar+1]^2
      ##rpi <-  mpisq/(4.0*pi*par[5+i])^2*g1( L*data[[i]][[k]]$m[uij] )
      rpi <-  mpisq/(4.0*pi*par[5+i])^2*g1( L*par[npar+1] )
      ##mpisq <-  par[npar+1]^2*(1.+0.5*rpi)^2
      fpi <-  par[5+i]*getfpi.os(par, mpisq, i, N=Nbeta, fit.a=fit.a)*(1.-2.*rpi)
      ##msssq <- data[[i]][[k]]$m[ssij]^2
      ##msssq <- par[(npar+2):(npar+length(uij)+length(ssij))]^2*(1.+0.5*rpi)^2
      msssq <- par[(npar+2):(npar+length(uij)+length(ssij))]^2
      nl <- length(msssq)
      ## mlsq <- rep(mpisq[1], times=nl)
      ##mlsq <- rep(par[npar+1]^2, times=nl)
      mlsq <- rep(mpisq, time=nl)
      ##rK <-  mlsq/(4.0*pi*par[5+i])^2*g1( L*sqrt(mlsq) )
      fK <- par[5+i]*getfK.os(par, mlsq, msssq, i, N=Nbeta, fit.a=fit.a)*(1-3*rpi/4.)

      x <- (c(fpi, fK, par[(npar+1):(npar+length(uij)+length(ssij))]*(1.+0.5*rpi)) -
            c(data[[i]][[k]]$f[uij], data[[i]][[k]]$fsinh[lsij], data[[i]][[k]]$m[uij],  data[[i]][[k]]$m[ssij]))
      chisum <- chisum + t(x) %*% cmatrix[[i]][[k]] %*% x

      npar <- npar + length(uij) + length(ssij)
      ##if(debug) {
      ##  cat("fpi:\n", fpi, "\n")
      ##  cat(data[[i]]$fsinh[uij], "\n")
      ##  cat("fK:\n", fK, "\n")
      ##  cat(data[[i]]$fsinh[lsij], "\n")
      ##}
    }
  }
  if(Nbeta>7) {
    chisum <- chisum + (par[6]/par[7] - 5.71/7.8)^2/(0.1)^2
  }
  return(chisum/2.)
}

getfpi.os <- function(par, mpisq, i, N, fit.a=TRUE) {
  xipi <- mpisq/(4*pi*par[5+i])^2
  if(fit.a) return(1-xipi*(2.*log(xipi) - par[3]) + par[6+N]*par[5+i]^2)
  return(1-xipi*(2.*log(xipi) - par[3]))
}

getfK.os <- function(par, mpisq, msssq, i, N, fit.a=TRUE) {
  xipi <- mpisq/(4.*pi*par[5+i])^2
  xiss <- msssq/(4.*pi*par[5+i])^2
  if(fit.a) return( (par[1] + par[2]*xiss) * (1 - xipi*(3./4.*log(xipi) - (par[4] + par[5]*xiss)) + (par[7+N] + par[8+N]*xiss)*par[5+i]^2) )
  else return( (par[1] + par[2]*xiss) * (1 - xipi*(3./4.*log(xipi) - (par[4] + par[5]*xiss))))
}

getfDssqrtmDs.os <- function(par, mpisq, msssq, mDs, af0, i, fit.a=TRUE) {
  xipi <- mpisq/(4.*pi*af0)^2
  xiss <- msssq/(4.*pi*af0)^2
  if(fit.a) {
    return( (par[1] + par[2]*xiss) *(1 + xipi*(par[3]+par[4]*xiss) + par[5]*mDs^2) + af0*par[6]/mDs +
           af0^2*(par[13] + par[14]*xiss))
  }
  return( (par[1] + par[2]*xiss) *(1 + xipi*(par[3]+par[4]*xiss) + par[5]*mDs^2) + af0*par[6]/mDs)
}

getRatio.os <- function(par, mpisq, msssq, mDs, af0, i, gc=0.61, fit.a=TRUE) {
  xipi <- mpisq/(4.*pi*af0)^2
  xiss <- msssq/(4.*pi*af0)^2
  if(fit.a) {
    return( (par[7] + par[8]*xiss) * (1 + 3*(1+3*gc^2)*xipi*log(xipi)/4. + xipi*(par[9]+par[10]*xiss) +
                                      par[11]*mDs^2 ) +af0*par[12]/mDs +af0^2*(par[15]+par[16]*xiss))
  }
  return( (par[7] + par[8]*xiss) * (1 + 3*(1+3*gc^2)*xipi*log(xipi)/4. + xipi*(par[9]+par[10]*xiss) +
                                    par[11]*mDs^2 ) +af0*par[12]/mDs)
}

fpsovmps.os <- function(x, par, i, N, fit.a=TRUE) {
  fpi <- par[5+i]*getfpi.os(par, x^2, i, N, fit.a)
  return(fpi / x - (130.7/135))
}


compindices <- function(Nens, Nstrange, Ncharm, Nlight) {
  Ntotal <- Nstrange+Ncharm+Nlight
  uii <- numeric(Nens)
  ssii <- numeric(Nstrange*Nens)
  lsii <- numeric(Nstrange*Nens)
  lcii <- numeric(Nens*Ncharm)
  scii <- numeric(Nens*Ncharm*Nstrange)

  N <- Ntotal
  uii[1] <- 1
  if(Nens > 1) {
    for(i in 2:Nens) {
      uii[i] <- uii[i-1] + (N*(N+1))/2
      N <- N-1
    }
  }
  
  Nl <- Nlight
  for(j in 1:Nens) {
    temp <- numeric(Nstrange)
    temp[1] <- uii[j]+Nl*(Nl+1)/2 + Nl*(Ntotal-Nlight)
    for(i in 2:Nstrange) {
      temp[i] <- temp[i-1] + (Ntotal - Nlight - i + 2 )
    }
    ssii[((j-1)*Nstrange+1):(j*Nstrange)] <- temp
    Nl <- Nl-1
  }

  Nl <- Nlight
  for(j in 1:Nens) {
    lsii[((j-1)*Nstrange+1):(j*Nstrange)] <- uii[j] + c(Nl:(Nl+Nstrange-1))
    Nl <- Nl-1
  }

  Nl <- Ntotal-Ncharm
  for(j in 1:Nens) {
    lcii[((j-1)*Ncharm+1):(j*Ncharm)] <- uii[j] + c(Nl:(Nl+Ncharm-1))
    Nl <- Nl-1
  }

  temp <- ssii + c(Nstrange:1)
  for(j in 1:(Nens*Nstrange)) {
    scii[((j-1)*Ncharm+1):(j*Ncharm)] <- c(temp[j]:(temp[j]+Ncharm-1))
  }
  return(invisible(list(uii=uii, lsii=lsii, ssii=ssii, lcii=lcii, scii=scii)))
}


OSChiPTfit2 <- function(data, startvalues_, ii, debug=TRUE, fit.a=FALSE) {
  ## number of lattice spacings
  Nbeta <- length(data)
  ## number of points with unitary light mu
  Nens <- c()
  np <- 5 + Nbeta
  dof <- -np
  for(i in 1:Nbeta) {
##    dof <- dof + length(ii[[uii[[i]]) + length(lsii[[i]])
    Nens <- c(Nens, length(data[[i]]))
  }

  mini <- optim(par=startvalues_, fn=chisqr.os2, method="BFGS", hessian=FALSE,
                control=list(maxit=100, trace=debug, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, fit.a=fit.a)
  mini <- optim(par=mini$par, fn=chisqr.os2, method="BFGS", hessian=FALSE,
                control=list(maxit=100, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, fit.a=fit.a)
  mini <- optim(par=mini$par, fn=chisqr.os2, method="BFGS", hessian=FALSE,
                control=list(maxit=500, trace=debug, parscale=mini$par, REPORT=50),
                data=data, ii=ii, Nens=Nens, Nbeta=Nbeta, fit.a=fit.a)

  ##cat("chisqr", chisqr.os(mini$par, data, uii, lsii, ssii, N, debug=TRUE), "\n")

  print(mini)
  i <- 1

  ampi.phys <- uniroot(fpsovmps.os, interval = c(0.001, 0.12),
                       tol=1.e-12, par=mini$par, i=i, N=Nbeta, fit.a=FALSE)$root
  afpi.phys <- mini$par[5+i] * getfpi.os(par=mini$par, ampi.phys^2, i=i, N=Nbeta, fit.a=FALSE)
  a <- ampi.phys/0.135*0.198*mini$par[6:(5+Nbeta)]/mini$par[5+i]
  ## neutral Kaon mass
  ## charged would be 0.49368 MeV
  amK.phys <- 0.49765*a/0.198
  amDs.phys = 1.9685*a/0.198
  amD.phys = 1.8696*a/0.198
  afK.phys <- mini$par[5+i]*getfK.os(par=mini$par, mpisq=ampi.phys^2, msssq=(2*amK.phys[i]^2-ampi.phys^2), i=i, N=Nbeta, fit.a=FALSE)
  
  miniD <- NULL

  para2 <- rep(1, times=16)
  ##chisqrD.os(para2, data, uii, ssii, lcii, scii, N, af0=mini$par[6:(5+N)], debug=TRUE)
  
  miniD <- optim(par=para2, fn=chisqrD.os, method="BFGS", hessian=FALSE,
                 control=list(maxit=150, trace=debug, REPORT=50),
                 data=data, ii=ii, af0=mini$par[6:(5+Nbeta)], Nbeta=Nbeta, Nens=Nens)
  miniD <- optim(par=miniD$par, fn=chisqrD.os, method="BFGS", hessian=FALSE,
                 control=list(maxit=500, trace=debug, parscale=miniD$par, REPORT=50),
                 data=data, ii=ii, af0=mini$par[6:(5+Nbeta)], Nbeta=Nbeta, Nens=Nens)
  print(miniD)
  afDs.phys <- mini$par[5+i]^(3/2) * getfDssqrtmDs.os(par=miniD$par, mpisq=ampi.phys^2, msssq=(2*amK.phys[i]^2-ampi.phys^2),
                                                    mDs=amDs.phys[i], af0=mini$par[5+i], i=i, fit.a=FALSE) / sqrt(amDs.phys[i])
  ratio.phys <- getRatio.os(par=miniD$par, mpisq=ampi.phys^2, msssq=(2*amK.phys[i]^2-ampi.phys^2),
                            mDs=amDs.phys[i], af0=mini$par[5+i], i=i, fit.a=FALSE)

  cat("ampi\t =", ampi.phys, "\n")
  cat("amK\t =", amK.phys, "\n")
  cat("amDs\t =", amDs.phys, "\n")
  cat("amD\t =", amD.phys, "\n")
  cat("afpi\t =", afpi.phys, "\n")
  cat("afK\t =", afK.phys, "\n")
  cat("a\t =", a, "fm \n")
  cat("f0\t =", mini$par[5+i]/a[i]*0.198, "GeV\n")
  cat("fK\t =", afK.phys/a[i]*0.198, "GeV\n")
  cat("fK/fpi\t =", afK.phys / afpi.phys, "\n")
  cat("l4\t =", mini$par[3]/2. + 2*log(4*pi*mini$par[5+i]/a[i]*0.198/0.1396), "\n")
  cat("chisqr/dof\t =", mini$value, "/", 1, "\n")
  cat("afDs\t =", afDs.phys, "\n")
  cat("fDs\t =", afDs.phys/a[i]*0.198, "GeV\n")
  cat("afD\t =", afDs.phys*sqrt(amDs.phys[i])*sqrt(amD.phys[i])/ratio.phys, "\n")
  cat("fD\t =", afDs.phys*sqrt(amDs.phys[i])*sqrt(amD.phys[i])/ratio.phys/a[i]*0.198, "GeV\n")
  cat("fDs/fD\t =", ratio.phys/sqrt(amDs.phys[i])/sqrt(amD.phys[i]), "\n")
  if(Nbeta == 2) {
    cat("f01/f02\t =", mini$par[6]/mini$par[7], "\n")
  }
}

chisqr.os2 <- function(par, data, ii, Nbeta, Nens, debug=FALSE, fit.a=FALSE) {

  chisum <- 0.
  for(i in 1:Nbeta) {
    for(k in 1:Nens[i]) {
      L <- data[[i]][[k]]$L[1]
      uij <- ii[[i]][[k]]$uii
      lsij <- ii[[i]][[k]]$lsii
      ssij <- ii[[i]][[k]]$ssii
      mpisq <- data[[i]][[k]]$m[uij]^2
      rpi <-  mpisq/(4.0*pi*par[5+i])^2*g1( L*data[[i]][[k]]$m[uij] )
      fpi <-  par[5+i]*getfpi.os(par, mpisq, i, N=Nbeta, fit.a=fit.a)*(1.-2.*rpi)
      msssq <- data[[i]][[k]]$m[ssij]^2
      mlsq <- numeric(length(msssq))
      nl <- length(msssq)/length(mpisq)
      for(j in 1:length(mpisq)) {
        mlsq[((j-1)*nl+1):(j*nl)] <- rep(mpisq[j], times=nl)
      }
      rK <-  mlsq/(4.0*pi*par[5+i])^2*g1( L*sqrt(mlsq) )
      fK <- par[5+i]*getfK.os(par, mlsq, msssq, i, N=Nbeta, fit.a=fit.a)*(1-3*rK/4.)
      chisum <- chisum + sum(( fpi - data[[i]][[k]]$fsinh[uij] )^2 / ( data[[i]][[k]]$df[uij] )^2) +
        sum(( fK - data[[i]][[k]]$fsinh[lsij] )^2 / ( data[[i]][[k]]$df[lsij] )^2)
      if(debug) {
        cat("fpi:\n", fpi, "\n")
        cat(data[[i]]$fsinh[uij], "\n")
        cat("fK:\n", fK, "\n")
        cat(data[[i]]$fsinh[lsij], "\n")
        
      }
    }
  }
  if(Nbeta>1) {
    chisum <- chisum + (par[6]/par[7] - 5.71/7.8)^2/(0.1)^2
  }
  return(chisum)
}


chisqrD.os <- function(par, data, ii, Nbeta, Nens, af0, gc=0.62, debug=FALSE) {
  chisum <- 0.
  for(i in 1:Nbeta) {
    for(k in 1:Nens[i]) {
      uij <- ii[[i]][[k]]$uii
      ##lsij <- lsii[[i]]
      ssij <- ii[[i]][[k]]$ssii
      lcij <- ii[[i]][[k]]$lcii
      scij <- ii[[i]][[k]]$scii
      Nstrange <- length(ssij)
      Ncharm <- length(scij)/Nstrange
      ## we do have Nens*Ncharm values for fD amd mD
      ## and Nens*Ncharm*Nstrange values for fDs and mDs
      mDDs <- numeric(Ncharm*Nstrange)
      fDDs <- numeric(Ncharm*Nstrange)
      mlsqDs <- numeric(Ncharm*Nstrange)
      msssqDs <- numeric(Ncharm*Nstrange)
      
      mpisq <- data[[i]][[k]]$m[uij]^2
      msssq <- data[[i]][[k]]$m[ssij]^2
      mDs <- data[[i]][[k]]$m[scij]
      fDs <- data[[i]][[k]]$fsinh[scij]
      mD <- data[[i]][[k]]$m[lcij]
      fD <- data[[i]][[k]]$fsinh[lcij]
      ## now we construct all vectors of the same length
      mlsqDs[1:(Ncharm*Nstrange)] <- rep(mpisq[1], times=Ncharm*Nstrange)
      for(j in 1:(Nstrange)) {
        msssqDs[((j-1)*Ncharm+1):(j*Ncharm)] <- rep(msssq[j], times=Ncharm)
      }
      for(j in 1:(Ncharm)) {
        mDDs[((j-1)*Nstrange+1):(j*Nstrange)] <- rep(mD[j], times=Nstrange)
        fDDs[((j-1)*Nstrange+1):(j*Nstrange)] <- rep(fD[j], times=Nstrange)
      }
      fDssqrtmDs <- af0[i]^(3/2) * getfDssqrtmDs.os(par, mpisq=mlsqDs, msssq=msssqDs, mDs=mDs, af0=af0[i], i=i)
      Ratio <- getRatio.os(par, mpisq=mlsqDs, msssq=msssqDs, mDs=mDs, af0=af0[i], i=i, gc=gc)
      if(debug) {
        cat("Nens", Nens[i], "\n")
        cat("Nstrange", Nstrange, "Ncharm", Ncharm, "\n")
        cat("mDs\t", mDs,"\n")
        cat("fDs\t", fDs,"\n")
        cat("mD\t", mDDs,"\n")
        cat("fD\t", fDDs,"\n")
        cat("mpi^2\t", mlsqDs,"\n")
        cat("mss^2\t", msssqDs,"\n")
        cat("fDssqrtmDs", fDssqrtmDs, "\n")
        cat("data      ", fDs*sqrt(mDs), "\n")
        cat("Ratio", Ratio, "\n")
        cat("data ", fDs*sqrt(mDs)/fDDs/sqrt(mDDs), "\n")
      }
      chisum <- chisum + sum( (fDs - fDssqrtmDs/sqrt(mDs))^2/data[[i]][[k]]$df[scij]) +
        sum( (Ratio - fDs*sqrt(mDs)/fDDs/sqrt(mDDs))^2 ) 
    }
  }
  return(chisum)
}
