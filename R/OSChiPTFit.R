## fit parameters
## 1:     f0K/f0
## 2:     fmK/f0
## 3:     b
## 4:     bK
## 5:     bmK
## 6-5+N: af0 
##
## data is a double list


OSChiPTfit <- function(data, startvalues_, bootsamples, uii, lsii, ssii, lcii, scii, debug=TRUE) {
  ## number of lattice spacings
  N <- length(data)
  ## number of points with unitary light mu

  np <- 5 + N
  dof <- -np
  for(i in 1:N) {
    dof <- dof + length(uii[[i]]) + length(lsii[[i]])
  }
  ##para <- c(1.6, 5., 1., 1., 1., 0.05)
  
  mini <- optim(par=startvalues_, fn=chisqr.os, method="BFGS", hessian=TRUE,
                control=list(maxit=150, trace=debug, REPORT=50),
                data=data, uii=uii, lsii=lsii, ssii=ssii, N=N)
  mini <- optim(par=mini$par, fn=chisqr.os, method="BFGS", hessian=TRUE,
                control=list(maxit=500, trace=debug, parscale=mini$par, REPORT=50),
                data=data, uii=uii, lsii=lsii, ssii=ssii, N=N)
  ##cat("chisqr", chisqr.os(mini$par, data, uii, lsii, ssii, N, debug=TRUE), "\n")

  print(mini)
  i <- 1

  ampi.phys <- uniroot(fpsovmps.os, interval = c(0.001, 0.12),
                       tol=1.e-12, par=mini$par, i=1)$root
  afpi.phys <- mini$par[6] * getfpi.os(par=mini$par, ampi.phys^2, i=1)
  a <- ampi.phys/0.135*0.198
  ## neutral Kaon mass
  ## charged would be 0.49368 MeV
  amK.phys <- 0.49765*a/0.198
  amDs.phys = 1.9685*a/0.198
  amD.phys = 1.8696*a/0.198
  afK.phys <- mini$par[6]*getfK.os(par=mini$par, mpisq=ampi.phys^2, msssq=(2*amK.phys^2-ampi.phys^2), i=1)
  
  miniD <- NULL
  if(!missing(lcii) && !missing(scii)) {
    para2 <- rep(1, times=12)
    ##chisqrD.os(para2, data, uii, ssii, lcii, scii, N, af0=mini$par[6:(5+N)], debug=TRUE)
    
    miniD <- optim(par=para2, fn=chisqrD.os, method="BFGS", hessian=TRUE,
                   control=list(maxit=150, trace=debug, REPORT=50),
                   data=data, uii=uii, ssii=ssii, lcii=lcii, scii=scii, af0=mini$par[6:(5+N)], N=N)
    miniD <- optim(par=miniD$par, fn=chisqrD.os, method="BFGS", hessian=TRUE,
                   control=list(maxit=500, trace=debug, parscale=miniD$par, REPORT=50),
                   data=data, uii=uii, ssii=ssii, lcii=lcii, scii=scii, af0=mini$par[6:(5+N)], N=N)
    print(miniD)
    afDs.phys <- mini$par[6]^(3/2) * getfDssqrtmDs.os(par=miniD$par, mpisq=ampi.phys^2, msssq=(2*amK.phys^2-ampi.phys^2),
                                                      mDs=amDs.phys, af0=mini$par[6], i=1) / sqrt(amDs.phys)
    ratio.phys <- getRatio.os(par=miniD$par, mpisq=ampi.phys^2, msssq=(2*amK.phys^2-ampi.phys^2),
                                                      mDs=amDs.phys, af0=mini$par[6], i=1)
  }
  cat("ampi\t =", ampi.phys, "\n")
  cat("amK\t =", amK.phys, "\n")
  cat("amDs\t =", amDs.phys, "\n")
  cat("amD\t =", amD.phys, "\n")
  cat("afK\t =", afK.phys, "\n")
  cat("a\t =", a, "fm \n")
  cat("f0\t =", mini$par[6]/a*0.198, "GeV\n")
  cat("fK\t =", afK.phys/a*0.198, "GeV\n")
  cat("fK/fpi\t =", afK.phys / afpi.phys, "\n")
  cat("l4\t =", mini$par[3]/2. + 2*log(4*pi*mini$par[6]/a*0.198/0.1396), "\n")
  cat("chisqr/dof\t =", mini$value, "/", dof, "\n")
  if(!missing(lcii) && !missing(scii)) {
    cat("afDs\t =", afDs.phys, "\n")
    cat("fDs\t =", afDs.phys/a*0.198, "GeV\n")
    cat("afD\t =", afDs.phys*sqrt(amDs.phys)*sqrt(amD.phys)/ratio.phys, "\n")
    cat("fD\t =", afDs.phys*sqrt(amDs.phys)*sqrt(amD.phys)/ratio.phys/a*0.198, "GeV\n")
    cat("fDs/fD\t =", ratio.phys/sqrt(amDs.phys)/sqrt(amD.phys), "\n")
  }
}



chisqr.os <- function(par, data, uii, lsii, ssii, N, debug=FALSE) {

  L <- 32
  chisum <- 0.
  for(i in 1:N) {
    uij <- uii[[i]]
    lsij <- lsii[[i]]
    ssij <- ssii[[i]]
    mpisq <- data[[i]]$m[uij]^2
    rpi <-  mpisq/(4.0*pi*par[5+i])^2*g1( L*data[[i]]$m[uij] )
    fpi <-  par[5+i]*getfpi.os(par, mpisq, i)*(1.-2.*rpi)
    msssq <- data[[i]]$m[ssij]^2
    mlsq <- numeric(length(msssq))
    nl <- length(msssq)/length(mpisq)
    for(j in 1:length(mpisq)) {
      mlsq[((j-1)*nl+1):(j*nl)] <- rep(mpisq[j], times=nl)
    }
    rK <-  mlsq/(4.0*pi*par[5+i])^2*g1( L*sqrt(mlsq) )
    fK <- par[5+i]*getfK.os(par, mlsq, msssq, i)*(1-3*rK/4.)
    chisum <- chisum + sum(( fpi - data[[i]]$f[uij] )^2 / ( data[[i]]$df[uij] )^2) +
      sum(( fK - data[[i]]$f[lsij] )^2 / ( data[[i]]$df[lsij] )^2)
    if(debug) {
      cat("fpi:\n", fpi, "\n")
      cat(data[[i]]$f[uij], "\n")
      cat("fK:\n", fK, "\n")
      cat(data[[i]]$f[lsij], "\n")
      
    }
  }
  return(chisum)
}

chisqrD.os <- function(par, data, uii, ssii, lcii, scii, N, af0, gc=0.62, debug=FALSE) {
  chisum <- 0.
  for(i in 1:N) {
    uij <- uii[[i]]
    ##lsij <- lsii[[i]]
    ssij <- ssii[[i]]
    lcij <- lcii[[i]]
    scij <- scii[[i]]
    Nens <- length(uij)
    Nstrange <- length(ssij)/Nens
    Ncharm <- length(scij)/Nstrange/Nens
    ## we do have Nens*Ncharm values for fD amd mD
    ## and Nens*Ncharm*Nstrange values for fDs and mDs
    mDDs <- numeric(Nens*Ncharm*Nstrange)
    fDDs <- numeric(Nens*Ncharm*Nstrange)
    mlsqDs <- numeric(Nens*Ncharm*Nstrange)
    msssqDs <- numeric(Nens*Ncharm*Nstrange)

    mpisq <- data[[i]]$m[uij]^2
    msssq <- data[[i]]$m[ssij]^2
    mDs <- data[[i]]$m[scij]
    fDs <- data[[i]]$f[scij]
    mD <- data[[i]]$m[lcij]
    fD <- data[[i]]$f[lcij]
    ## now we construct all vectors of the same length
    for(j in 1:Nens) {
      mlsqDs[((j-1)*Ncharm*Nstrange+1):(j*Ncharm*Nstrange)] <- rep(mpisq[j], times=Ncharm*Nstrange)
    }
    for(j in 1:(Nens*Nstrange)) {
      msssqDs[((j-1)*Ncharm+1):(j*Ncharm)] <- rep(msssq[j], times=Ncharm)
    }
    for(j in 1:(Nens*Ncharm)) {
      mDDs[((j-1)*Nstrange+1):(j*Nstrange)] <- rep(mD[j], times=Nstrange)
      fDDs[((j-1)*Nstrange+1):(j*Nstrange)] <- rep(fD[j], times=Nstrange)
    }
    fDssqrtmDs <- af0[i]^(3/2) * getfDssqrtmDs.os(par, mpisq=mlsqDs, msssq=msssqDs, mDs=mDs, af0=af0[i], i=i)
    Ratio <- getRatio.os(par, mpisq=mlsqDs, msssq=msssqDs, mDs=mDs, af0=af0[i], i=i, gc=gc)
    if(debug) {
      cat("Nens", Nens, "Nstrange", Nstrange, "Ncharm", Ncharm, "\n")
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
    chisum <- chisum + sum( (fDs*sqrt(mDs) - fDssqrtmDs)^2) +
      sum( (Ratio - fDs*sqrt(mDs)/fDDs/sqrt(mDDs))^2 ) 
  }
  return(chisum)
}

getfpi.os <- function(par, mpisq, i) {
  xipi <- mpisq/(4*pi*par[5+i])^2
  return(1-xipi*(2.*log(xipi) - par[3]))
}

getfK.os <- function(par, mpisq, msssq, i) {
  xipi <- mpisq/(4.*pi*par[5+i])^2
  xiss <- msssq/(4.*pi*par[5+i])^2
  return( (par[1] + par[2]*xiss) * (1 - xipi*(3./4.*log(xipi) - (par[4] + par[5]*xiss))) )
}

getfDssqrtmDs.os <- function(par, mpisq, msssq, mDs, af0, i) {
  xipi <- mpisq/(4.*pi*af0)^2
  xiss <- msssq/(4.*pi*af0)^2
  return( (par[1] + par[2]*xiss) *(1 + xipi*(par[3]+par[4]*xiss) + par[5]*mDs^2) + af0*par[6]/mDs )
}

getRatio.os <- function(par, mpisq, msssq, mDs, af0, i, gc=0.61) {
  xipi <- mpisq/(4.*pi*af0)^2
  xiss <- msssq/(4.*pi*af0)^2
  return( (par[7] + par[8]*xiss) * (1 + 3*(1+3*gc^2)*xipi*log(xipi)/4. + xipi*(par[9]+par[10]*xiss) + par[11]*mDs^2 ) +af0*par[12]/mDs )
}

fpsovmps.os <- function(x, par, i) {
  fpi <- par[5+i]*getfpi.os(par, x^2, i)
  return(fpi / x - (130.7/135))
}

