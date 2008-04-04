g1 <- function(x) {

  weights <- c(6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.)
  ex <- c(1:20)
  res <- x
  for( i in 1:length(x)) {
    sex <- x[i]*sqrt(ex)
    res[i] <- sum(4*weights*besselK(sex, 1)/(sex))
  }
  return(res)
}

fitit <- function(data, corr=NULL, twoB0=5.05, f0=0.0534, xln3=-1.92, xln4=-1.06,
                  fsmethod="gl", a.guess=0.0087) {
  
  Bmu <- twoB0*data$mu
  mpssq <- Bmu*(1.0+Bmu*(log(Bmu)-xln3)/(4.0*pi*f0)^2 )
  fps <- f0*(1.0-2.0*Bmu*(log(Bmu)-xln4)/(4.0*pi*f0)^2 )

  mpsv <- data$L*sqrt(mpssq)
  r <- mpssq/(4.0*pi*f0)^2*g1(mpsv)

  mps <- sqrt(mpssq)*(1.+0.5*r)
  fps <- fps*(1.0-2.0*r)
  par <- c(twoB0, f0, xln3, xln4)
  csq <- chisqr.comb(par, data=data)
  mini <- optim(par=par, fn=chisqr.comb, method="BFGS", control=list(maxit=150),
                hessian=TRUE, data=data, fsmethod=fsmethod,
                a.guess=a.guess)
  
  # compute some errors
  # invert Hessian
  Hinv <- solve(mini$hessian)
  # compute l3, l4, a and mu.phys
  mu.phys <- uniroot(fovermps, c(0.0001,0.004), tol = 1.e-12, par=mini$par)$root
  l3 <- mini$par[3] -
    log(mini$par[1]*mu.phys*(1.0+mini$par[1]*mu.phys*(log(mini$par[1]*mu.phys)-mini$par[3])/(4.0*pi*mini$par[2])^2 ))
  l4 <- mini$par[4] -
    log(mini$par[1]*mu.phys*(1.0+mini$par[1]*mu.phys*(log(mini$par[1]*mu.phys)-mini$par[3])/(4.0*pi*mini$par[2])^2 ))
  a <- (mini$par[2]*(1.0-2.0*mini$par[1]*mu.phys*(log(mini$par[1]*mu.phys)-mini$par[4])/(4.0*pi*mini$par[2])^2 ))/0.1307*0.198
  result <- data.frame(TwoB0=mini$par[1], dTwoB0=sqrt(2*Hinv[1,1]),
                       F0=mini$par[2], dF0=sqrt(2*Hinv[2,2]),
                       L3=mini$par[3], dL3=sqrt(2.*Hinv[3,3]),
                       L4=mini$par[4], dL4=sqrt(2.*Hinv[4,4]),
                       l3=l3, l4=l4, a=a, mu.phys=mu.phys)
  
  if(!missing(corr)) {
    Z <- array(0., dim=c(2,2,4))
    for(i in 1:4) {
      corr[1,1,i] <- corr[1,1,i]*data$dmps[i]^2
      corr[2,2,i] <- corr[2,2,i]*data$dfps[i]^2
      corr[2,1,i] <- corr[2,1,i]*data$dmps[i]*data$dfps[i]
      corr[1,2,i] <- corr[1,2,i]*data$dmps[i]*data$dfps[i]
      Z[,,i] <- solve(corr[,,i])
    }

    mini.cor <- optim(par=par, fn=chisqr.comb.corr, method="BFGS", control=list(type=3, maxit=150),
                      hessian=TRUE, data=data, cor.inv=Z)
    Hinv.cor <- solve(mini.cor$hessian)
    # compute l3, l4, a and mu.phys
    mu.phys.cor <- uniroot(fovermps, c(0.0001,0.004), tol = 1.e-12, par=mini.cor$par)$root
    l3.cor <- mini.cor$par[3] -
      log(mini.cor$par[1]*mu.phys.cor*(1.0+mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[3])/(4.0*pi*mini.cor$par[2])^2 ))
    l4.cor <- mini.cor$par[4] -
      log(mini.cor$par[1]*mu.phys.cor*(1.0+mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[3])/(4.0*pi*mini.cor$par[2])^2 ))
    a.cor <- (mini.cor$par[2]*(1.0-2.0*mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[4])/(4.0*pi*mini.cor$par[2])^2 ))/0.1307*0.198

    result.cor <- data.frame(TwoB0=mini.cor$par[1], dTwoB0=sqrt(2*Hinv.cor[1,1]),
                             F0=mini.cor$par[2], dF0=sqrt(2*Hinv.cor[2,2]),
                             L3=mini.cor$par[3], dL3=sqrt(2.*Hinv.cor[3,3]),
                             L4=mini.cor$par[4], dL4=sqrt(2.*Hinv.cor[4,4]),
                             l3=l3.cor, l4=l4.cor, a=a.cor, mu.phy=mu.phys.cor)

  }
  else {
    mini.cor <- NULL
    result.cor <- NULL
    chidiff.cor <- NULL
  }

  
  # Check with MINOS method
  # fix one \pm errors and refit the others
  # the compare Chi^2
  chidiff <- array(rep(0., times=8), dim=c(4,2))
  ind <- c(2:4)
  for(fix in 1:4) {
    if(fix!=1) {
      ind[(fix-1)]=(fix-1)
    }
    par <- mini$par[ind]
    value <- mini$par[fix]+sqrt(2*Hinv[fix,fix])
    mini2 <- optim(par=par, fn=chisqr.comb.fix, method="BFGS", control=list(type=1, maxit=150), hessian=FALSE, data=data, fix=fix, ind=ind, value=value)
    chidiff[fix,1] <- mini2$value-mini$value
    
    par <- mini$par[ind]
    value <- mini$par[fix]-sqrt(2*Hinv[fix,fix])
    mini2 <- optim(par=par, fn=chisqr.comb.fix, method="BFGS", control=list(type=1, maxit=150), hessian=FALSE, data=data, fix=fix, ind=ind, value=value)
    chidiff[fix,2] <- mini2$value-mini$value
  }

  if(!missing(corr)) {
    chidiff.cor <- array(rep(0., times=8), dim=c(4,2))
    ind <- c(2:4)
    for(fix in 1:4) {
      if(fix!=1) {
        ind[(fix-1)]=(fix-1)
      }
      par <- mini.cor$par[ind]
      value <- mini.cor$par[fix]+sqrt(2*Hinv.cor[fix,fix])
      mini.cor2 <- optim(par=par, fn=chisqr.comb.fix.corr, method="BFGS", control=list(type=1, maxit=150), 
                         hessian=FALSE, data=data, fix=fix, ind=ind, value=value, cor.inv=Z)
      chidiff.cor[fix,1] <- mini.cor2$value-mini.cor$value
    
      par <- mini.cor$par[ind]
      value <- mini.cor$par[fix]-sqrt(2*Hinv.cor[fix,fix])
      mini.cor2 <- optim(par=par, fn=chisqr.comb.fix.corr, method="BFGS", 
                         hessian=FALSE, data=data, fix=fix, ind=ind, value=value, cor.inv=Z)
      chidiff.cor[fix,2] <- mini.cor2$value-mini.cor$value
    }
  }

  xln3 <- mini$par[3]
  xln4 <- mini$par[4]
  f0 <- mini$par[2]
  twoB0 <- mini$par[1]

  Bmu <- twoB0*data$mu
  mpssq <- data$mps^2
  fps <- data$fps
  mpsv <- data$L*sqrt(mpssq)
  r <- mpssq/(4.0*pi*f0)^2*g1(mpsv)
  mps <- data$mps/(1.+0.5*r)
  fps <- fps/(1.0-2.0*r)
  df.corrected <- data.frame(mu=data$mu, mps=mps, dmps=data$dmps, fps=fps, dfps=data$dfps)
#  library(systemfit)

#  mps.formula <- mps/dmps ~ sqrt(TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 ))/dmps
#  fps.formula <- fps/dfps ~ F0*(1.0-2.0*TwoB*mu*(log(TwoB*mu)-L4)/(4.0*pi*F0)^2)/dfps
  
#  mps.formula2 <- mps/dmps ~ (sqrt(TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 ))
#                             *(1.+0.5*TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 )/
#                               (4.0*pi*f0)*g2(L*sqrt(TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 ))))
#                             )/dmps
#  fps.formula2 <- fps/dfps ~ (F0*(1.0-2.0*TwoB*mu*(log(TwoB*mu)-L4)/(4.0*pi*F0)^2)
#                             *(1.-2.*TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 )/
#                             (4.0*pi*f0)*g2(L*sqrt(TwoB*mu*(1.0+TwoB*mu*(log(TwoB*mu)-L3)/(4.0*pi*F0)^2 ))))
#                             )/dfps
#  start.values <- c(TwoB=5., F0=0.054, L3=-1.92, L4=-1.06)
#  model <- list(mps.formula, fps.formula)
#  model2 <- list(mps.formula2, fps.formula2)
#  model.ols <- nlsystemfit( "OLS", model, start.values, data=df.corrected)
#  print(model.ols)
#  model2.ols <- nlsystemfit( "SUR", model2, start.values, data=data)
#  return(invisible(list(result=result, fit=mini, combfit=model.ols)))
  return(invisible(list(result=result, fit=mini, chidiff=chidiff,
                        result.cor=result.cor, fit.cor=mini.cor,
                        chidiff.cor=chidiff.cor, fsmethod=fsmethod,
                        data=data)))
}

fitit.boot <- function(data, bootsamples, corr=NULL, twoB0=5.0, f0=0.0534, xln3=-1.92, xln4=-1.06,
                       fsmethod="gl", a.guess=0.0087) {
  boots <- data.frame(TwoB0 = rep(0., times=length(bootsamples[,1,1])),
                      F0 = rep(0., times=length(bootsamples[,1,1])),
                      L3 = rep(0., times=length(bootsamples[,1,1])),
                      L4 = rep(0., times=length(bootsamples[,1,1])),
                      mu = rep(0., times=length(bootsamples[,1,1])),
                      l3 = rep(0., times=length(bootsamples[,1,1])),
                      l4 = rep(0., times=length(bootsamples[,1,1])),
                      a  = rep(0., times=length(bootsamples[,1,1])),
                      F  = rep(0., times=length(bootsamples[,1,1]))
                      )

  par <- c(twoB0, f0, xln3, xln4)
  mini.first <- optim(par=par, fn=chisqr.comb, method="BFGS", 
                hessian=FALSE, data=data, fsmethod=fsmethod,
                a.guess=a.guess)
  
  # compute l3, l4, a and mu.phys
  mu.phys <- uniroot(fovermps, c(0.0001,0.001), tol = 1.e-12, par=mini.first$par)$root
  l3 <- mini.first$par[3] -
    log(mini.first$par[1]*mu.phys*(1.0+mini.first$par[1]*mu.phys*(log(mini.first$par[1]*mu.phys)-mini.first$par[3])/(4.0*pi*mini.first$par[2])^2 ))
  l4 <- mini.first$par[4] -
    log(mini.first$par[1]*mu.phys*(1.0+mini.first$par[1]*mu.phys*(log(mini.first$par[1]*mu.phys)-mini.first$par[3])/(4.0*pi*mini.first$par[2])^2 ))
  a <- (mini.first$par[2]*(1.0-2.0*mini.first$par[1]*mu.phys*(log(mini.first$par[1]*mu.phys)-mini.first$par[4])/(4.0*pi*mini.first$par[2])^2 ))/0.1307*0.1973
  F = mini.first$par[2]/a*0.1973
  result <- data.frame(TwoB0=mini.first$par[1], F0=mini.first$par[2], L3=mini.first$par[3], L4=mini.first$par[4],
                       l3=l3, l4=l4, a=a, mu=mu.phys, F=F)

  
  if(!missing(corr)) {
    Z <- array(0., dim=c(2,2,4))
    for(i in 1:4) {
      corr[1,1,i] <- corr[1,1,i]*data$dmps[i]^2
      corr[2,2,i] <- corr[2,2,i]*data$dfps[i]^2
      corr[2,1,i] <- corr[2,1,i]*data$dmps[i]*data$dfps[i]
      corr[1,2,i] <- corr[1,2,i]*data$dmps[i]*data$dfps[i]
      Z[,,i] <- solve(corr[,,i])
    }

    mini.cor <- optim(par=par, fn=chisqr.comb.corr, method="BFGS", control=list(type=3, maxit=150),
                      hessian=TRUE, data=data, cor.inv=Z)
    mu.phys.cor <- uniroot(fovermps, c(0.0001,0.001), tol = 1.e-12, par=mini.cor$par)$root
    l3.cor <- mini.cor$par[3] -
      log(mini.cor$par[1]*mu.phys.cor*(1.0+mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[3])/(4.0*pi*mini.cor$par[2])^2 ))
    l4.cor <- mini.cor$par[4] -
      log(mini.cor$par[1]*mu.phys.cor*(1.0+mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[3])/(4.0*pi*mini.cor$par[2])^2 ))
    a.cor <- (mini.cor$par[2]*(1.0-2.0*mini.cor$par[1]*mu.phys.cor*(log(mini.cor$par[1]*mu.phys.cor)-mini.cor$par[4])/(4.0*pi*mini.cor$par[2])^2 ))/0.1307*0.198
    F.cor = mini.cor$par[2]/a.cor*0.1973
    
    result.cor <- data.frame(TwoB0=mini.cor$par[1], F0=mini.cor$par[2], L3=mini.cor$par[3], L4=mini.cor$par[4],
                             l3=l3.cor, l4=l4.cor, a=a.cor, mu=mu.phys.cor, F=F.cor)

  }
  else {
    mini.cor <- NULL
    result.cor <- NULL
  }

  boots <- array(0., dim=c(length(bootsamples[,1,1]), 10))

  
  for(s in 1:length(bootsamples[,1,1])) {
    mps <- bootsamples[s,1,(1:length(data$mu))]
    fps <- bootsamples[s,2,(1:length(data$mu))]

    df <- data.frame(mu=data$mu, mps=mps, dmps=data$dmps, fps=fps, dfps=data$dfps, L=data$L)
    par <- mini.first$par
    mini <- optim(par=par, fn=chisqr.comb, method="BFGS", 
                  hessian=FALSE, data=df, fsmethod=fsmethod)
                                        # compute l3, l4, a and mu.phys
    boots[s,1] <- mini$par[1]
    boots[s,2] <- mini$par[2]
    boots[s,3] <- mini$par[3]
    boots[s,4] <- mini$par[4]
    boots[s,5] <- uniroot(fovermps, c(0.0001,0.001), tol = 1.e-12, par=mini$par)$root
    boots[s,6] <- mini$par[3] -
      log(mini$par[1]*boots[s,5]*(1.0+mini$par[1]*boots[s,5]*(log(mini$par[1]*boots[s,5])-mini$par[3])/(4.0*pi*mini$par[2])^2 ))
    boots[s,7] <- mini$par[4] -
      log(mini$par[1]*boots[s,5]*(1.0+mini$par[1]*boots[s,5]*(log(mini$par[1]*boots[s,5])-mini$par[3])/(4.0*pi*mini$par[2])^2 ))
    boots[s,8] <- (mini$par[2]*(1.0-2.0*mini$par[1]*boots[s,5]*(log(mini$par[1]*boots[s,5])-mini$par[4])/(4.0*pi*mini$par[2])^2 ))/0.1307*0.1973
    boots[s,9] = mini$value
    boots[s,10] = mini$par[2]/boots[s,8]*0.1973
    
  }
  boot.result <- data.frame(TwoB0=mean(boots[,1]), dTwoB0=sd(boots[,1]),
                            F0=mean(boots[,2]), dF0=sd(boots[,2]),
                            L3=mean(boots[,3]), dL3=sd(boots[,3]),
                            L4=mean(boots[,4]), dL4=sd(boots[,4]),
                            mu=mean(boots[,5]), dmu=sd(boots[,5]),
                            a=mean(boots[,8]), da=sd(boots[,8]),
                            l3=mean(boots[,6]), dl3=sd(boots[,6]),
                            l4=mean(boots[,7]), dl4=sd(boots[,7]),
                            F=mean(boots[,10]), dF=sd(boots[,10]))

  return(invisible(list(result=result, result.cor=result.cor,result.boot=boot.result,
                        t=boots, mini=mini.first, mini.cor=mini.cor, data=data,
                        fsmethod=fsmethod)))
}

correct.mf <- function(par, data) {
  fps <- data$fps
  mps <- data$mps
  mpsv <- data$L*mps
  r <-  mps^2/(4.0*pi*par[2])^2*g1(mpsv)
  
  mpsinf <- mps/(1.+0.5*r)
  fpsinf <- fps/(1.0-2.0*r)
  return(invisible(data.frame(mu=data$mu, mpsinf=mpsinf, fpsinf=fpsinf, mpsL=mps, fpsL=fps, dmps=data$dmps, dfps=data$dfps, L=data$L)))
}

prediction.mf <- function(TwoB, aF, aL3, aL4, x, L, fsmodel="gl", a_fm, predict.inf=FALSE) {
  if(predict.inf) rev <- -1
  else rev <- 1
  TwoBmu <- TwoB*x
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-log(aL3^2))/(4.0*pi*aF)^2 )
  fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu)-log(aL4^2))/(4.0*pi*aF)^2 )
  mps <- sqrt(mpssq) 
  mpsv <- L*sqrt(mpssq)
  if(fsmodel == "cdh") {
    aLamb1=sqrt(exp(-0.4+log((0.135*a_fm/0.1973)^2)))
    aLamb2=sqrt(exp(4.3+log((0.135*a_fm/0.1973)^2)))
    res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aL3,
               aLamb4=aL4, ampiV=mps, afpiV=fps,
               aF0=fps, a_fm=a_fm, L=data$L, rev=1)
    mpsL <- res$mpiFV
    fpsL <- res$fpiFV
  }
  else {
    r <-  mpssq/(4.0*pi*aF)^2*g1(mpsv)
  # r <-  TwoBmu/(4.0*pi*par[2])^2*g1(mpsv)
    mpsL <- sqrt(mpssq)*(1.+rev*0.5*r)
    fpsL <- fps*(1.0-rev*2.0*r)
  }
  return(invisible(data.frame(x=x, mps=mps, fps=fps, mpsL=mpsL, fpsL=fpsL, TwoBmu=TwoBmu)))
}

correct.datamf <- function(data, aF, aL3, aL4, fsmodel="gl", a_fm) {

  mpssq <- data$mps^2
  mpsv <- data$L*sqrt(mpssq)
  if(fsmodel == "cdh") {
    aLamb1=sqrt(exp(-0.4+log((0.135*a_fm/0.1973)^2)))
    aLamb2=sqrt(exp(4.3+log((0.135*a_fm/0.1973)^2)))
    res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=aL3,
               aLamb4=aL4, ampiV=data$mps, afpiV=data$fps,
               aF0=data$fps, a_fm=a_fm, L=data$L, rev=-1)
    mps <- res$mpiFV
    fps <- res$fpiFV    
  }
  else {
    r <- mpssq/(4.0*pi*aF)^2*g1(mpsv)
    mps <- data$mps/(1.+0.5*r)
    fps <- data$fps/(1.0-2.0*r)
  }
  df.corrected <- data.frame(mu=data$mu, mps=mps, dmps=data$dmps, fps=fps, dfps=data$dfps, L=data$L)
  return(invisible(df.corrected))
}

fovermps <- function(x, par) {
  TwoBmu <- par[1]*x
  mps <- sqrt( TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 ) )
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )
  return(fps/mps - (130.7/135))
}

fs.model <-  function(x, m, L, par) {
  return(m-x-par[1]*exp(-x*L))
}


chisqr.comb <- function(par, data, fsmethod="gl", model.par=c(9.08874565, -6.12895863),
                        a.guess=0.00078) {

  TwoBmu <- par[1]*data$mu
  if(par[1] < 0) {
    return(invisible(100000))
  }
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 )
  if(any(mpssq < 0)) {
    return(NaN)
  }
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )

  if(fsmethod == "cdh") {
    mu.phys <- try(uniroot(fovermps, c(0.0005,0.0009), tol = 1.e-12, par=par)$root, silent=T)
    if(inherits(mu.phys, "try-error") || is.nan(mu.phys)) a_fm <- a.guess
    else a_fm <- (par[2]*(1.0-2.0*par[1]*mu.phys*(log(par[1]*mu.phys)-par[4])/(4.0*pi*par[2])^2 ))/0.1307*0.1973
    aLamb1=sqrt(exp(-0.4+log((0.135*a_fm/0.1973)^2)))
    aLamb2=sqrt(exp(4.3+log((0.135*a_fm/0.1973)^2)))
    res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=sqrt(exp(par[3])),
               aLamb4=sqrt(exp(par[4])), ampiV=sqrt(mpssq), afpiV=fps,
               aF0=fps, a_fm=a_fm, L=data$L, rev=1)
    mps <- res$mpiFV
    fps <- res$fpiFV
  }
  else if(fsmethod == "model") {
    mps <- sqrt(mpssq)
    emps <- exp(-mps*data$L)/data$L^(3/2)
    mps <- mps + model.par[1]*emps
    fps <- fps + model.par[2]*emps
  }
  else {
    mpsv <- data$L*sqrt(mpssq)
    r <-  mpssq/(4.0*pi*par[2])^2*g1(mpsv)
    
    mps <- sqrt(mpssq)*(1.+0.5*r)
    fps <- fps*(1.0-2.0*r)
  }
  return(invisible(sum(((data$mps-mps)/data$dmps)^2) + sum(((data$fps-fps)/data$dfps)^2)))

}

chisqr.comb.corr <- function(par, data, cor.inv) {
  
  TwoBmu <- par[1]*data$mu
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 )
  if(any(mpssq < 0)) return(NaN)
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )
  
  mpsv <- data$L*sqrt(mpssq)
  r <-  mpssq/(4.0*pi*par[2])^2*g1(mpsv)
  
  mps <- sqrt(mpssq)*(1.+0.5*r)
  fps <- fps*(1.0-2.0*r)
  Sum <- 0.
  Sum <- sum((data$mps-mps)^2*cor.inv[1,1,]
             + (data$fps-fps)^2*cor.inv[2,2,]
             + (data$mps-mps)*(data$fps-fps)*cor.inv[2,1,])
  return(invisible(Sum))
}

chisqr.comb.fix <- function(par, data, ind, fix, value) {

  partmp <- rep(0., times=4)
  for(i in 1:length(par)) {
    partmp[ind[i]] <- par[i]
  }
  partmp[fix] <- value
  
  TwoBmu <- partmp[1]*data$mu
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-partmp[3])/(4.0*pi*partmp[2])^2 )
  if(any(mpssq < 0)) return(NaN)
  fps <- partmp[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-partmp[4])/(4.0*pi*partmp[2])^2 )
  
  mpsv <- data$L*sqrt(mpssq)
  r <-  mpssq/(4.0*pi*partmp[2])^2*g1(mpsv)
  
  mps <- sqrt(mpssq)*(1.+0.5*r)
  fps <- fps*(1.0-2.0*r)
  
  return(invisible(sum(((data$mps-mps)/data$dmps)^2) + sum(((data$fps-fps)/data$dfps)^2)))

}

chisqr.comb.fix.corr <- function(par, data, ind, fix, value, cor.inv) {

  partmp <- rep(0., times=4)
  for(i in 1:length(par)) {
    partmp[ind[i]] <- par[i]
  }
  partmp[fix] <- value
  
  TwoBmu <- partmp[1]*data$mu
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-partmp[3])/(4.0*pi*partmp[2])^2 )
  if(any(mpssq < 0)) return(NaN)
  fps <- partmp[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-partmp[4])/(4.0*pi*partmp[2])^2 )
  
  mpsv <- data$L*sqrt(mpssq)
  r <-  mpssq/(4.0*pi*partmp[2])^2*g1(mpsv)
  
  mps <- sqrt(mpssq)*(1.+0.5*r)
  fps <- fps*(1.0-2.0*r)

  Sum <- 0.
  Sum <- sum((data$mps-mps)^2*cor.inv[1,1,]
             + (data$fps-fps)^2*cor.inv[2,2,]
             + (data$mps-mps)*(data$fps-fps)*cor.inv[2,1,])
  return(invisible(Sum))
}

# data must be a list containing the N data.frames for the N lattice spacings
# startvalues must be the fit parameter vector
# (L3/F, L4/F, a_1F, 2a_1B, a_2F, 2a_2B, ...)

fit.Na <- function(data, startvalues, bootsamples, fsmethod="gl", a.guess) {
  if(missing(data) || missing(startvalues)) {
    stop("data and startvalues must be provided!\n")
  }
  if(missing(bootsamples)) {
    cat("bootstrap samples missing, will not compute error estimates!\n")
  }
  N <- length(data)
  if(fsmethod=="cdh") {
    if(missing(a.guess)) {
      cat("a.guess was missing! Guessing myself!\n")
      a.guess <- rep(0.1, times=N) 
    }
  }

  
  mini <- optim(par=startvalues, fn=chisqr.Na, method="BFGS", hessian=TRUE, data=data,
                fsmethod=fsmethod, a.guess=a.guess)
  dof <- 0
  for(i in 1:length(data)) {
    dof = dof + 2*length(data[[i]]$mu)
  }
  dof <- dof - length(startvalues)
  chisqr <- mini$value
  par <- mini$par
  # compute observables
  mu.phys <- c(0.)
  a <- c(0.)
  l3 <- c(0.)
  l4 <- c(0.)
  F <- c(0.)
  for(i in 1:length(data)) {
    mu.phys[i] <- uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par, indd=i)$root
    TwoaB <- par[(2+2*i-1)]
    TwoBmu <- TwoaB*mu.phys[i]
    aF <- par[(2+2*i)]
    a[i] <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
    l3[i] <- log(par[1]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    l4[i] <- log(par[2]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    F[i] = aF/a[i]*0.1973
  }
  
  boot.result <- NULL
  boots <- NULL
  if(!missing(bootsamples)) {

    boots <- array(0., dim=c(length(bootsamples[[1]][,1,1]), (2*length(data)+length(startvalues)+3)))

    for(s in 1:length(bootsamples[[1]][,1,1])) {
      mps <- bootsamples[[1]][s,1,(1:length(data[[1]]$mu))]
      fps <- bootsamples[[1]][s,2,(1:length(data[[1]]$mu))]
      df <- list(data.frame(mu=data[[1]]$mu, mps=mps, dmps=data[[1]]$dmps,
                         fps=fps, dfps=data[[1]]$dfps, L=data[[1]]$L))
      for(i in 2:length(data)) {
        mps <- bootsamples[[i]][s,1,(1:length(data[[i]]$mu))]
        fps <- bootsamples[[i]][s,2,(1:length(data[[i]]$mu))]
        df[[i]] <- data.frame(mu=data[[i]]$mu, mps=mps, dmps=data[[i]]$dmps,
                            fps=fps, dfps=data[[i]]$dfps, L=data[[i]]$L)
      }
      mini.boot <- optim(par=par, fn=chisqr.Na, method="BFGS", 
                         hessian=FALSE, data=df, fsmethod=fsmethod,
                         a.guess=a.guess)
      par.boot <- mini.boot$par
      TwoaB <- 0.
      TwoBmu <- 0.
      aF <- 0.
      for(i in 1:length(data)) {
        boots[s,i] <- uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par.boot, indd=i)$root
        TwoaB <- par.boot[(2+2*i-1)]
        TwoBmu <- TwoaB*boots[s,i]
        aF <- par.boot[(2+2*i)]
        boots[s,(i+length(data))] <-
          aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
      }
      boots[s,(1+2*length(data))] <- log(par.boot[1]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      boots[s,(2+2*length(data))] <- log(par.boot[2]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      boots[s,(3+2*length(data))] <- aF/boots[s,(2*length(data))]*0.1973
      for(i in 1:length(startvalues)) {
        boots[s,(2*length(data)+3+i)] <- par.boot[i]
      }
    }
    boot.result <- array(0., dim=c(length(boots[1,]), 2))
    for(i in 1:(length(boots[1,]))) {
      boot.result[i,1] <-  mean(boots[,i])
      boot.result[i,2] <-  sd(boots[,i])
    }
  }
  
  result <- list(par=par, mu.phys=mu.phys, F=F, a=a, l3=l3, l4=l4, chisqr=chisqr, dof=dof,
                 data=data, boot.result=boot.result, boots=boots)
  return(invisible(result))
}

chisqr.Na <- function(par, data, fsmethod="gl", a.guess) {

  N <- length(data)
  chisum <- 0.
  for( i in 1:N) {
    
    TwoaB <- par[(2+2*i-1)]
    if(any(TwoaB < 0)) {
      return(invisible(100000))
    }
    aF <- par[(2+2*i)]
    TwoBmu <- TwoaB*data[[i]]$mu
    mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 )
    if(any(mpssq < 0)) {
      return(NaN)
    }
    fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )

    if(fsmethod=="cdh") {
      mu.phys <- try(uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par, indd=i)$root, silent=T)
      if(inherits(mu.phys, "try-error") || is.nan(mu.phys)) a_fm <- a.guess[i]
      else {
        TwoBmu <- TwoaB*mu.phys
        a_fm <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
      }
      #a_fm <- a.guess[i]
      aLamb1=sqrt(exp(-0.4+log((0.135*a_fm/0.1973)^2)))
      aLamb2=sqrt(exp(4.3+log((0.135*a_fm/0.1973)^2)))
      res <- cdh(aLamb1=aLamb1, aLamb2=aLamb2, aLamb3=par[1]*aF,
                 aLamb4=par[2]*aF, ampiV=sqrt(mpssq), afpiV=fps,
                 aF0=fps, a_fm=a_fm, L=data[[i]]$L, rev=1, printit=F)
      mps <- res$mpiFV
      fps <- res$fpiFV
    }
    else {
      mpsv <- data[[i]]$L*sqrt(mpssq)
      r <-  mpssq/(4.0*pi*aF)^2*g1(mpsv)

      mps <- sqrt(mpssq)*(1.+0.5*r)
      fps <- fps*(1.0-2.0*r)
    }
    chisum <- chisum + (sum(((data[[i]]$mps-mps)/data[[i]]$dmps)^2) + sum(((data[[i]]$fps-fps)/data[[i]]$dfps)^2))
  }
  return(invisible(chisum))

}

fit.Na.cdh <- function(data, startvalues, bootsamples, a.guess) {
  if(missing(data) || missing(startvalues)) {
    stop("data and startvalues must be provided!\n")
  }
  if(missing(bootsamples)) {
    cat("bootstrap samples missing, will not compute error estimates!\n")
  }
  N <- length(data)
  np <- length(startvalues)

  
  mini <- optim(par=startvalues, fn=chisqr.Na.cdh, method="BFGS", hessian=TRUE, data=data,
                a.guess=a.guess)
  dof <- 0
  for(i in 1:length(data)) {
    dof = dof + 2*length(data[[i]]$mu)
  }
  dof <- dof - length(startvalues)
  chisqr <- mini$value
  par <- mini$par
  # compute observables
  mu.phys <- c(0.)
  a <- c(0.)
  l1 <- c(0.)
  l2 <- c(0.)
  l3 <- c(0.)
  l4 <- c(0.)
  F <- c(0.)
  for(i in 1:length(data)) {
    mu.phys[i] <- uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par, indd=i)$root
    TwoaB <- par[(2+2*i-1)]
    TwoBmu <- TwoaB*mu.phys[i]
    aF <- par[(2+2*i)]
    a[i] <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
    l1[i] <- log(par[np-1]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    l2[i] <- log(par[np]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    l3[i] <- log(par[1]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    l4[i] <- log(par[2]^2*aF^2) -
      log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 ) )
    F[i] = aF/a[i]*0.1973
  }
  
  boot.result <- NULL
  boots <- NULL
  if(!missing(bootsamples)) {

    boots <- array(0., dim=c(length(bootsamples[[1]][,1,1]), (2*N+np+5)))

    for(s in 1:length(bootsamples[[1]][,1,1])) {
      mps <- bootsamples[[1]][s,1,(1:length(data[[1]]$mu))]
      fps <- bootsamples[[1]][s,2,(1:length(data[[1]]$mu))]
      df <- list(data.frame(mu=data[[1]]$mu, mps=mps, dmps=data[[1]]$dmps,
                         fps=fps, dfps=data[[1]]$dfps, L=data[[1]]$L))
      for(i in 2:length(data)) {
        mps <- bootsamples[[i]][s,1,(1:length(data[[i]]$mu))]
        fps <- bootsamples[[i]][s,2,(1:length(data[[i]]$mu))]
        df[[i]] <- data.frame(mu=data[[i]]$mu, mps=mps, dmps=data[[i]]$dmps,
                            fps=fps, dfps=data[[i]]$dfps, L=data[[i]]$L)
      }
      mini.boot <- optim(par=par, fn=chisqr.Na.cdh, method="BFGS", 
                         hessian=FALSE, data=df,
                         a.guess=a.guess)
      par.boot <- mini.boot$par
      TwoaB <- 0.
      TwoBmu <- 0.
      aF <- 0.
      for(i in 1:N) {
        boots[s,i] <- uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par.boot, indd=i)$root
        TwoaB <- par.boot[(2+2*i-1)]
        TwoBmu <- TwoaB*boots[s,i]
        aF <- par.boot[(2+2*i)]
        boots[s,(i+length(data))] <-
          aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
      }
      # l1
      boots[s,(1+2*N)] <- log(par.boot[np-1]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      # l2
      boots[s,(2+2*N)] <- log(par.boot[np]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      # l3
      boots[s,(3+2*N)] <- log(par.boot[1]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      #l4
      boots[s,(4+2*N)] <- log(par.boot[2]^2*aF^2) -
        log( TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par.boot[1]^2))/(4.0*pi*aF)^2 ) )
      # f0 in GeV
      boots[s,(5+2*N)] <- aF/boots[s,(2*length(data))]*0.1973
      # fit parameter
      for(i in 1:length(startvalues)) {
        boots[s,(2*length(data)+5+i)] <- par.boot[i]
      }
    }
    boot.result <- array(0., dim=c(length(boots[1,]), 2))
    for(i in 1:(length(boots[1,]))) {
      boot.result[i,1] <-  mean(boots[,i])
      boot.result[i,2] <-  sd(boots[,i])
    }
  }
  else {
    bootsamples <- NULL
  }
  
  result <- list(par=par, mu.phys=mu.phys, F=F, a=a, l1=l1, l2=l2, l3=l3, l4=l4, chisqr=chisqr, dof=dof,
                 data=data, boot.result=boot.result, boots=boots, bootsamples=bootsamples)
  return(invisible(result))
}

chisqr.Na.cdh <- function(par, data, a.guess) {

  N <- length(data)
  np <- length(par)
  chisum <- 0.
  for( i in 1:N) {
    
    TwoaB <- par[(2+2*i-1)]
    if(any(TwoaB < 0)) {
      return(invisible(100000))
    }
    aF <- par[(2+2*i)]
    TwoBmu <- TwoaB*data[[i]]$mu
    mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 )
    if(any(mpssq < 0)) {
      return(NaN)
    }
    fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )

    mu.phys <- try(uniroot(fovermps.Na, c(0.0001, 0.004), tol=1.e-12, par=par, indd=i)$root, silent=T)
    if(inherits(mu.phys, "try-error") || is.nan(mu.phys)) a_fm <- a.guess[i]
    else {
      TwoBmu <- TwoaB*mu.phys
      a_fm <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )/0.1307*0.1973
    }
                                        #a_fm <- a.guess[i]
    res <- cdh(aLamb1=par[np-1]*aF, aLamb2=par[np]*aF, aLamb3=par[1]*aF,
               aLamb4=par[2]*aF, ampiV=sqrt(mpssq), afpiV=fps,
               aF0=fps, a_fm=a_fm, L=data[[i]]$L, rev=1, printit=F, parm=2.)
    mps <- res$mpiFV
    fps <- res$fpiFV
    
    chisum <- chisum + (sum(((data[[i]]$mps-mps)/data[[i]]$dmps)^2) + sum(((data[[i]]$fps-fps)/data[[i]]$dfps)^2))
  }
  return(invisible(chisum))

}

fovermps.Na <- function(x, par, indd) {
  i <- indd
  TwoaB <- par[(2+2*i-1)]
  aF <- par[(2+2*i)]
  TwoBmu <- TwoaB*x
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 )
  if(any(mpssq < 0)) return(NaN)
  fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )
  mps <- sqrt(mpssq)
  return(fps/mps - (130.7/135))
}

fovermps.Na.withr0 <- function(x, par, indd) {
  i <- indd
  r0TwoB <- par[4]
  r0F <- par[3]
  r0sqTwoBmu <- r0TwoB*x
  mpssq <- r0sqTwoBmu*(1.0+r0sqTwoBmu*(log(r0sqTwoBmu/par[1]^2))/(4.0*pi*r0F)^2 )
  if(any(mpssq < 0)) return(NaN)
  fpsV <- r0F*(1.0-2.0*r0sqTwoBmu*(log(r0sqTwoBmu/par[2]^2))/(4.0*pi*r0F)^2 )
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
# 6+2*N:          c_1
# 7+2*N:          g_A
# 7+2*N+1:        r0*L1
# 7+2*N+2:        r0*L2
# 7+2*N+3:        k_M
# 7+2*N+4:        k_F

fit.Na.withr0ZP <- function(data, startvalues, bootsamples, fsmethod="gl", a.guess,
                            r0bootsamples, r0data, ZPdata, ZPbootsamples,
                            fit.l12=FALSE, fit.asq=FALSE,
                            ii, boot.R=100, debug=FALSE) {
  if(missing(data) || missing(startvalues)) {
    stop("data and startvalues must be provided!\n")
  }
  if(missing(bootsamples)) {
    cat("bootstrap samples missing, will not compute error estimates!\n")
  }
  else if(missing(boot.R)) {
    boot.R <- length(bootsamples[[1]][,1,1])
  }
  N <- length(data)
  if(missing(ii)){
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
  if(length(startvalues) != np) stop("length of startvalues", length(startvalues),
             "must match number of fit parameters ", np, "!\n")
  if(fsmethod=="cdh") {
    if(missing(a.guess)) {
      cat("a.guess was missing! Guessing myself!\n")
      a.guess <- rep(0.1, times=N) 
    }
  }

  if(debug) {
    chisqr.Na.withr0ZP(par=startvalues, data=data, ii=ii,fsmethod=fsmethod,
                       a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                       fit.l12=fit.l12, fit.asq=fit.asq, printit=debug)
  }

  mini <- optim(par=startvalues, fn=chisqr.Na.withr0ZP, method="BFGS", hessian=FALSE,
                control=list(maxit=500), data=data, ii=ii,
                fsmethod=fsmethod, a.guess=a.guess, r0data=r0data, ZPdata=ZPdata,
                fit.l12=fit.l12, fit.asq=fit.asq)
  if(mini$convergence != 0) {
    warning("Attention: optim did not converge in initial run!\n Please adjust start values!")
  }
  if(debug) {
    chisqr.Na.withr0ZP(par=mini$par, data=data, ii=ii,fsmethod=fsmethod,
                     a.guess=a.guess, r0data=r0data, ZPdata=ZPdata, fit.l12=fit.l12, printit=debug)
  }
  dof <- 0
  for(i in 1:length(data)) {
    dof = dof + 2*length( ii[[i]] ) + length(na.omit(data[[i]]$mN[ii[[i]]]))
  }
  dof <- dof + length(r0data$r0) + length(ZPdata$ZP) - length(startvalues)
  chisqr <- mini$value
  par <- mini$par

  # compute observables
  mu.phys <- numeric()
  a <- numeric()
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
  for(i in 1:N) {
    mu.phys[i] <- uniroot(fovermps.Na.withr0, c(r0data$r0[i]*0.0001, r0data$r0[i]*0.004),
                          tol=1.e-12, par=par, indd=i)$root
    r0TwoB <- par[4]
    r0sqTwoBmu <- r0TwoB*mu.phys[i]
    r0F <- par[3]
    r0mN <- par[5+2*N]
    mpisq <- getmpssq(r0sqTwoBmu, par, N, fit.l12)
    fpi <- getfps(r0sqTwoBmu, par, N, fit.l12)
    mN <- getmN(r0sqTwoBmu, par, N)/fpi*0.130
    a[i] <- fpi/0.1307*0.1973/par[4+i]
    lmpssq <-  log( 0.1396^2 )
    l3 <- log((par[1]/fpi*0.130)^2) - lmpssq
    l4 <- log((par[2]/fpi*0.130)^2) - lmpssq
    if(fit.l12) {
      l1 <- log((par[8+2*N]/fpi*0.130)^2) - lmpssq
      l2 <- log((par[9+2*N]/fpi*0.130)^2) - lmpssq
    }
    F <- r0F/fpi*0.130
    mN0 <- r0mN/fpi*0.130
    c1 <- par[6+2*N]*fpi/0.130
    gA <- par[7+2*N]
    B0 <- r0TwoB/fpi*0.130/2.
    R0[i] <- par[4+i]
    Zp[i] <- par[4+i+N]
    Sigma <- (B0*F^2/2)^(1/3)
    rssq <- 12./(4*pi*F)^2*(l4-12./13.)*0.1973^2
  }
  
  boot.result <- NULL
  boots <- NULL
  if(!missing(bootsamples)) {
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

    boots <- array(0., dim=c(boot.R, (2*N+length(startvalues)+14)))
    for(s in 1:boot.R) {
      df <- list(data.frame(mu=data[[1]]$mu[ ii[[1]] ], mps=bootsamples[[1]][s, 1, ii[[1]] ],
                            dmps=data[[1]]$dmps[ ii[[1]] ],
                            fps=bootsamples[[1]][s, 2, ii[[1]] ],
                            dfps=data[[1]]$dfps[ ii[[1]] ], L=data[[1]]$L[ ii[[1]] ],
                            mN=rnorm(length(ii[[1]]), mean=data[[1]]$mN[ ii[[1]] ], sd=data[[1]]$dmN[ ii[[1]] ]),
                            dmN=data[[1]]$dmN[ ii[[1]] ] ))
      r0df <- data.frame(r0=r0bootsamples[s,], dr0=r0data$dr0)
      ZPdf <- data.frame(ZP=ZPbootsamples[s,], dZP=ZPdata$dZP)
      for(i in 2:N) {
        df[[i]] <- data.frame(mu=data[[i]]$mu[ ii[[i]] ], mps=bootsamples[[i]][s,1, ii[[i]] ],
                              dmps=data[[i]]$dmps[ ii[[i]] ],
                              fps=bootsamples[[i]][s,2, ii[[i]] ],
                              dfps=data[[i]]$dfps[ ii[[i]] ], L=data[[i]]$L[ ii[[i]] ],
                              mN=rnorm(length(ii[[i]]), mean=data[[i]]$mN[ ii[[i]] ], sd=data[[i]]$dmN[ ii[[i]] ]),
                              dmN=data[[i]]$dmN[ ii[[i]] ])
      }

      mini.boot <- optim(par=par, fn=chisqr.Na.withr0ZP, method="BFGS", hessian=FALSE,
                         control=list(maxit=500, trace=0), data=df, ii=ii,
                         fsmethod=fsmethod, a.guess=a.guess, r0data=r0df, ZPdata=ZPdf,
                         fit.l12=fit.l12, fit.asq=fit.asq)
      if(mini.boot$convergence != 0) {
        warning("optim did not converge for one bootsample")
        par.boot <- rep(NA, length(mini.boot$par))
        boots[s,] <- NA
      }
      else {
        par.boot <- mini.boot$par
        r0F <- numeric()
        mpssqr <- numeric()
        for(i in 1:N) {
          boots[s,i] <- uniroot(fovermps.Na.withr0, c(r0df$r0[i]*0.0001, r0df$r0[i]*0.004),
                                tol=1.e-12, par=par.boot, indd=i)$root
                                        #        r0TwoB <- par.boot[4]
          r0sqTwoBmu <- par.boot[4]*boots[s,i]
          r0F <- par.boot[3]
          mpssqr <- getmpssq(r0sqTwoBmu, par.boot, N, fit.l12)
          fpi <- getfps(r0sqTwoBmu, par.boot, N, fit.l12)
          boots[s,(i+N)] <- fpi/0.1307*0.1973/par.boot[4+i]
        }
        lmpssqr <- log( 0.1396^2 )
        boots[s,(1+2*N)] <- log((par.boot[1]/fpi*0.130)^2) - lmpssqr
        boots[s,(2+2*N)] <- log((par.boot[2]/fpi*0.130)^2) - lmpssqr
        if(fit.l12) {
          boots[s,(3+2*N)] <- log((par.boot[8+2*N]/fpi*0.130)^2) - lmpssqr
          boots[s,(4+2*N)] <- log((par.boot[9+2*N]/fpi*0.130)^2) - lmpssqr
        }
        boots[s,(5+2*N)] <- par.boot[3]/fpi*0.130
        boots[s,(6+2*N)] <- par.boot[4]/fpi*0.130/2.
        boots[s,(7+2*N)] <- par.boot[5+2*N]/fpi*0.130
        boots[s,(8+2*N)] <- getmN(r0sqTwoBmu, par.boot, N)/fpi*0.130
        boots[s,(9+2*N)] <- par.boot[6+2*N]*fpi/0.130
        boots[s,(10+2*N)] <- (boots[s,(6+2*N)]*boots[s,(5+2*N)]^2/2)^(1/3)
        boots[s,(11+2*N)] <- 12./(4*pi*boots[s,(5+2*N)])^2*(boots[s,(2+2*N)]-12./13.)*0.1973^2
        for(i in 1:length(par.boot)) {
          boots[s,(2*N+11+i)] <- par.boot[i]
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
  
  result <- list(par=par, result=list(mu.phys=mu.phys, F=F, a=a, l1=l1, l2=l2, l3=l3, l4=l4,
                            B0=B0, mN0=mN0, mN=mN, c1=c1, gA=gA,r0=R0, ZP=Zp, Sigma=Sigma, rssq=rssq,
                            chisqr=chisqr, dof=dof), fit=mini,
                 data=data, boot.result=boot.result, boots=boots, r0data=r0data, ZPdata=ZPdata,
                 bootsamples=bootsamples, r0bootsamples=r0bootsamples, ZPbootsamples=ZPbootsamples,
                 ii=ii, fit.l12=fit.l12, boot.R=boot.R, fsmethod=fsmethod, fit.asq=fit.asq)
  attr(result, "class") <- c("chiralfit", "list")  
  return(invisible(result))
}

chisqr.Na.withr0ZP <- function(par, data, ii, r0data, ZPdata, fsmethod="gl", a.guess,
                               fit.l12=FALSE, fit.asq=FALSE, printit=FALSE) {

  fit.a <- -1.
  N <- length(ii)
  chisum <- 0.
  for( i in 1:N) {
    if(fit.asq) {
      fit.a <- 4+i
    }
    ij <- ii[[i]]
    r0TwoB <- par[4]
    if(any(r0TwoB <= 0)) {
      return(invisible(NaN))
    }
#    r0 <- par[4+i]
#    ZP <- par[4+N+i]
    r0F <- par[3]
    r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+N+i]
    mpssq <- getmpssq(r0sqTwoBmu, par, N, fit.l12, fit.asq=fit.a)
    fpsV <- getfps(r0sqTwoBmu, par, N, fit.l12, fit.asq=fit.a)
    if(FALSE) {
      cat("r0         ", r0, "\n")
      cat("ZP         ", ZP, "\n")
      cat("mu         ", mu[ij], "\n")
      cat("noFSm      ", sqrt(mpssq)/r0, "\n")
      cat("noFSf      ", fpsV/r0, "\n")
      cat("r0sqTwoBmu ", r0sqTwoBmu, "\n")
      cat("r0F        ", par[3], "\n")
      cat("par        ", par, "\n")
      cat("N          ", N, "\n")
      cat("fit.l12    ", fit.l12, "\n\n")
    }
    if(any(is.nan(mpssq)) || any(mpssq <= 0)) {
      return(NaN)
    }
    if(fsmethod=="cdh") {
      r0mu.phys <- try(uniroot(fovermps.Na.withr0,
                               c(r0*0.0001, r0*0.004),
                               tol=1.e-12, par=par, indd=i)$root, silent=T)
      if(inherits(r0mu.phys, "try-error") || is.nan(r0mu.phys)) a_fm <- a.guess[i]
      else {
        a_fm <- getfps(r0sqTwoBmu=r0TwoB*r0mu.phys, par, fit.l12)/0.1307*0.1973/par[4+i]
      }
      if(fit.l12) {
        aLamb1=par[8+2*N]/par[4+i]
        aLamb2=par[9+2*N]/par[4+i]
      }
      else {
#        aLamb1=sqrt(exp(-0.4+log((0.1396*a_fm/0.1973)^2)))
        aLamb1=sqrt(exp(-0.4)*(0.1396*a_fm/0.1973)^2)
#        aLamb2=sqrt(exp(4.3+log((0.1396*a_fm/0.1973)^2)))
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
    chisum <- chisum + (sum(((data[[i]]$mps[ij]-mpsV)/(data[[i]]$dmps[ij]))^2) +
                        sum(((data[[i]]$fps[ij]-fpsV)/(data[[i]]$dfps[ij]))^2))
    if(fit.l12 && i==1) {
      chisum <- chisum + ((-0.4-log((aLamb1)^2/(0.1396*a_fm/0.1973)^2))/0.6)^2 +
        ((4.3-log((aLamb2)^2/(0.1396*a_fm/0.1973)^2))/0.1)^2
#      cat(log((aLamb1)^2/(0.1396*a_fm/0.1973)^2), log((aLamb2)^2/(0.1396*a_fm/0.1973)^2), "\n")
    }
    if(printit) {
      cat("r0model ", par[4+i], "\n")
      cat("r0data  ", r0data$r0[i], "\n")
      cat("chir0   ", (r0data$r0[i]-par[4+i])/r0data$dr0[i], "\n")
      cat("modelZP ", par[4+N+i], "\n")
      cat("ZPdata  ", ZPdata$ZP[i], "\n")
      cat("chiZP   ", (ZPdata$ZP[i]-par[4+N+i])/ZPdata$dZP[i], "\n")
      cat("modelm  ", mpsV, "\n")
      cat("datam   ", data[[i]]$mps[ij], "\n")
      cat("errm    ", data[[i]]$dmps[ij], "\n")
      cat("chim    ", ((data[[i]]$mps[ij]-mpsV)/data[[i]]$dmps[ij]), "\n")
      cat("modelf  ", fpsV, "\n")
      cat("dataf   ", data[[i]]$fps[ij], "\n")
      cat("errf    ", data[[i]]$dfps[ij], "\n")
      cat("chif    ", ((data[[i]]$fps[ij]-fpsV)/data[[i]]$dfps), "\n")
    }

# the nucleon
    mN <- getmN(r0sqTwoBmu=r0sqTwoBmu, par, N, fit.asq=fit.a)/par[4+i]
    chisum <- chisum + sum(((data[[i]]$mN-mN)/data[[i]]$dmN)^2, na.rm=TRUE)
  }
  chisum <- chisum + sum(((par[(5):(4+N)]-r0data$r0)/r0data$dr0)^2) +
    sum(((par[(5+N):(4+2*N)]-ZPdata$ZP)/ZPdata$dZP)^2)
  if(printit) {
    cat("chisqr ", chisum, "\n\n")
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
    mpssq <- getmpssq(r0sqTwoBmu, par, N, fit.l12)
    fpsV <- getfps(r0sqTwoBmu, par, N, fit.l12)
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

getmpssq <- function(r0sqTwoBmu, par, N, fit.l12=FALSE, fit.asq=-1) {

  npar <- length(par)
  xi <- r0sqTwoBmu/(4.0*pi*par[3])^2
  rln3 <- log(r0sqTwoBmu/par[1]^2)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
    fit.asq <- 1
  }
  if(fit.l12) {
    rln1 <- log(r0sqTwoBmu/par[2*N+8]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+9]^2)
    return(r0sqTwoBmu*(1. + xi*rln3 + asq/par[fit.asq]^2*par[npar-2] +
                       17./2.*xi^2*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)^2
#                                    +8./17.*par[N+6]
                                    )
                       )
           )
  }
  return(r0sqTwoBmu*(1.0 + asq/par[fit.asq]^2*par[npar-2] + r0sqTwoBmu*(rln3/(4.0*pi*par[3])^2 )))
}

getfps <- function(r0sqTwoBmu, par, N, fit.l12=FALSE, fit.asq=-1.) {

  npar <- length(par)
  xi <- r0sqTwoBmu/(4.0*pi*par[3])^2
  rln3 <- log(r0sqTwoBmu/par[1]^2)
  rln4 <- log(r0sqTwoBmu/par[2]^2)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
    fit.asq <- 1
  }
  if(fit.l12) {
    rln1 <- log(r0sqTwoBmu/par[2*N+8]^2)
    rln2 <- log(r0sqTwoBmu/par[2*N+9]^2)
    return(par[3]*(1. - 2.*xi*rln4 +  asq/par[fit.asq]^2*par[npar-1] -
                   5*xi^2*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)^2
                                        #                        + 4./5.*par[N+7]
                           )
                   )
           )
  }
  return(par[3]*(1.0 + asq/par[fit.asq]^2*par[npar-1]  - 2.0*r0sqTwoBmu*(rln4)/(4.0*pi*par[3])^2 ))
}

getmN <- function(r0sqTwoBmu, par, N, fit.asq=-1.) {
  npar <- length(par)
  asq <- 1.
  if(fit.asq == -1) {
    asq <- 0.
    fit.asq <- 1
  }
  return(par[5+2*N]*(1 + asq/par[fit.asq]^2*par[npar]) - 4.*par[6+2*N]*r0sqTwoBmu - 6.*par[7+2*N]^2/(32*pi*par[3]^2)*r0sqTwoBmu^(3/2))
}
