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
                  fsmethod="gl") {
  
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
                hessian=TRUE, data=data, fsmethod=fsmethod)
  
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
                       fsmethod="gl") {
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
                hessian=FALSE, data=data, fsmethod=fsmethod)
  
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

prediction.mf <- function(par, x, L) {
  TwoBmu <- par[1]*x
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 )
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )
  mps <- sqrt(mpssq) 
  mpsv <- L*sqrt(mpssq)
  # mpsv <- L*sqrt(TwoBmu)
  r <-  mpssq/(4.0*pi*par[2])^2*g1(mpsv)
  # r <-  TwoBmu/(4.0*pi*par[2])^2*g1(mpsv)
  mpsL <- sqrt(mpssq)*(1.+0.5*r)
  fpsL <- fps*(1.0-2.0*r)
  return(invisible(data.frame(x=x, mps=mps, fps=fps, mpsL=mpsL, fpsL=fpsL)))
}

correct.mf <- function(data, f0) {

  mpssq <- data$mps^2
  mpsv <- data$L*sqrt(mpssq)
  r <- mpssq/(4.0*pi*f0)^2*g1(mpsv)
  mps <- data$mps/(1.+0.5*r)
  fps <- data$fps/(1.0-2.0*r)
  df.corrected <- data.frame(mu=data$mu, mps=mps, dmps=data$dmps, fps=fps, dfps=data$dfps, L=data$L, f0=f0)
  return(invisible(df.corrected))
}

fovermps <- function(x, par) {
  TwoBmu <- par[1]*x
  mps <- sqrt( TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 ) )
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )
  return(fps/mps - (130.7/139.6))
}

fs.model <-  function(x, m, L, par) {
  return(m-x-par[1]*exp(-x*L))
}


chisqr.comb <- function(par, data, fsmethod="gl", model.par=c(9.08874565, -6.12895863)) {

  TwoBmu <- par[1]*data$mu
  mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu)-par[3])/(4.0*pi*par[2])^2 )
  fps <- par[2]*(1.0-2.0*TwoBmu*(log(TwoBmu)-par[4])/(4.0*pi*par[2])^2 )

  if(fsmethod == "cdh") {
    mu.phys <- try(uniroot(fovermps, c(0.0005,0.0009), tol = 1.e-12, par=par)$root, silent=T)
    if(inherits(mu.phys, "try-error")) mu.phys <- 0.00078
    a_fm <- (par[2]*(1.0-2.0*par[1]*mu.phys*(log(par[1]*mu.phys)-par[4])/(4.0*pi*par[2])^2 ))/0.1307*0.198
    res <- cdh(aLamb3=sqrt(exp(par[3])) , aLamb4=sqrt(exp(par[4])), ampiV=sqrt(mpssq), afpiV=fps,
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

fit.Na <- function(data, startvalues, bootsamples) {
  N <- length(data)
  
  mini <- optim(par=startvalues, fn=chisqr.Na, method="BFGS", hessian=TRUE, data=data)
#  print(mini)
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
                         hessian=FALSE, data=df)
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
  else {
    boot.result=NULL
  }
  
  result <- list(par=par, mu.phys=mu.phys, F=F, a=a, l3=l3, l4=l4, chisqr=chisqr, dof=dof,
                 data=data, boot.result=boot.result, boots=boots)
  return(invisible(result))
}

chisqr.Na <- function(par, data) {

  N <- length(data)
  chisum <- 0.
  for( i in 1:N) {
    
    TwoaB <- par[(2+2*i-1)]
    aF <- par[(2+2*i)]
    TwoBmu <- TwoaB*data[[i]]$mu
    mpssq <- TwoBmu*(1.0+TwoBmu*(log(TwoBmu/aF^2)-log(par[1]^2))/(4.0*pi*aF)^2 )
    fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )
  
    mpsv <- data[[i]]$L*sqrt(mpssq)
    r <-  mpssq/(4.0*pi*aF)^2*g1(mpsv)
  
    mps <- sqrt(mpssq)*(1.+0.5*r)
    fps <- fps*(1.0-2.0*r)
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
  fps <- aF*(1.0-2.0*TwoBmu*(log(TwoBmu/aF^2)-log(par[2]^2))/(4.0*pi*aF)^2 )
  mps <- sqrt(mpssq)
  return(fps/mps - (130.7/139.6))
}
