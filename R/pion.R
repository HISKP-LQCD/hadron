pion.boot <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, pl=FALSE,
                      boot.R=99, boot.l=10, tsboot.sim="geom", skip=0) {
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  library(boot)
  Time <-  2*max(cmicor$V3)
  Thalf <- max(cmicor$V3)
  T1 <- max(cmicor$V3)+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- max(cmicor$V1)
  Skip <- (skip*(T1)*nrObs*4+1)
  Z <- array(cmicor$V4[(Skip):length(cmicor$V4)]+cmicor$V5[(Skip):length(cmicor$V5)], 
             dim=c(nrObs*(T1)*4,(length(cmicor$V4[(Skip):length(cmicor$V4)])/(nrObs*(T1)*4))))
  rm(cmicor)
  W <- Z[1:(4*(T1)),]
  for(i in 1:(T1)) {
    if(i!=1 && i!=(T1)) {
      W[i,] <- W[i,]/2.
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),])/4.
      W[(i+2*T1),] <- W[(i+T1),]
      W[(i+3*T1),] <- W[(i+3*T1),]/2.
    }
    else {
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),])/2.
      W[(i+2*T1),] <- W[(i+T1),]
    }
  }
  rm(Z)

  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5) 

  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)
  
  Cor <- rep(0., times=4*T1)
  E <- rep(0., times=4*T1)
  
  for(i in 1:(T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
  }
  for(i in (T1):(2*T1)) {
    Cor[i] = mean(W[(i),])
    Cor[i+T1] = Cor[i]
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
    E[i+T1] = E[i]
  }
  for(i in (3*T1):(4*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
  }
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=E[(t1p1):(t2p1)], e2=E[(t1p1+T1):(t2p1+T1)],
                   e3=E[(t1p1+2*T1):(t2p1+2*T1)], e4=E[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(c(1.,0.1,0.12), chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  fit.mass <- abs(pionfit$par[3])
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass^3)*exp(fit.mass*Time/4)
  
  fit.dof <- (t2-t1+1)*3-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }
  
  pionfit.boot <- boot(data=t(W), statistic=fit2by2mass.boot, R=boot.R, stype="i",
                       Time=Time, t1=t1, t2=t2, Err=E, kappa=kappa, mu=mu)

  df.cov <- data.frame(x=pionfit.boot$t[,1], y=pionfit.boot$t[,2])
  Covmf <- cov(df.cov)
  Cormf <- cor(df.cov)
  
  if(pl) {
    X11()
    plot(pionfit.boot)
  }
  pionfit.boot.ci <- boot.ci(pionfit.boot, type = c("norm", "basic", "perc", "stud"))
  pionfit.tsboot <- tsboot(t(W), statistic=fit2by2mass.tsboot, R=boot.R, l=boot.l, sim=tsboot.sim,
                           Time=Time, t1=t1, t2=t2, Err=E, kappa=kappa, mu=mu)

  df.covts <- data.frame(x=pionfit.tsboot$t[,1], y=pionfit.tsboot$t[,2])
  Covtsmf <- cov(df.cov)
  Cortsmf <- cor(df.cov)

  if(pl) {
    X11()
    plot(pionfit.tsboot)
  }
  pionfit.tsboot.ci <- boot.ci(pionfit.tsboot, type = c("norm", "basic", "perc", "stud"))
  res <- list(fitresult=pionfit, boot=pionfit.boot, boot.ci=pionfit.boot.ci,
              tsboot=pionfit.tsboot, tsboot.ci=pionfit.tsboot.ci,
              Covmf=Covmf, Cormf=Cormf, Covtsmf=Covtsmf, Cortsmf=Cortsmf,
              effmass=pion.eff, kappa=kappa, mu=mu, t1=t1, t2=t2,
              Time=Time, N=length(W[1,]), matrix.size=2, no.masses=1)
  attr(res, "class") <- c("cfit", "list")    
  return(invisible(res))
}

pion <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  Time <-  2*max(cmicor$V3)
  Thalf <- max(cmicor$V3)
  T1 <- max(cmicor$V3)+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- max(cmicor$V1)
  Skip <- (skip*(T1)*nrObs*4+1)
  Z <- array(cmicor$V4[(Skip):length(cmicor$V4)]+cmicor$V5[(Skip):length(cmicor$V5)], 
             dim=c(nrObs*(T1)*4,(length(cmicor$V4[(Skip):length(cmicor$V4)])/(nrObs*(T1)*4))))
  rm(cmicor)
  W <- Z[1:(4*(T1)),]
  for(i in 1:(T1)) {
    if(i!=1 && i!=(T1)) {
      W[i,] <- W[i,]/2.
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),])/4.
      W[(i+2*T1),] <- W[(i+T1),]
      W[(i+3*T1),] <- W[(i+3*T1),]/2.
    }
    else {
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),])/2.
      W[(i+2*T1),] <- W[(i+T1),]
    }
  }
  rm(Z)

  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5) 

  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)
  
  Cor <- rep(0., times=4*T1)
  E <- rep(0., times=4*T1)
  
  for(i in 1:(T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
  }
  for(i in (T1):(2*T1)) {
    Cor[i] = mean(W[(i),])
    Cor[i+T1] = Cor[i]
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
    E[i+T1] = E[i]
  }
  for(i in (3*T1):(4*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
  }
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=E[(t1p1):(t2p1)], e2=E[(t1p1+T1):(t2p1+T1)],
                   e3=E[(t1p1+2*T1):(t2p1+2*T1)], e4=E[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(c(1.,0.1,0.12), chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  fit.mass <- abs(pionfit$par[3])
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass^3)*exp(fit.mass*Time/4)
  
  fit.dof <- (t2-t1+1)*3-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }

  fit.uwerrm <- uwerrderived(f=fit2by2mass, data=W, S=S, pl=pl, Time=Time, t1=t1, t2=t2, Err=E)
  fit.uwerrf <- uwerrderived(f=fit2by2fpi, data=W, S=S, pl=pl, Time=Time, t1=t1, t2=t2, Err=E)
  
  Chi <- rep(0., times=4*T1)
  Fit <- rep(0., times=4*T1)

  for( i in t1p1:t2p1) {
    Fit[i] <- pionfit$par[1]^2*cosh(fit.mass*((i-1)-Thalf))
    Chi[i] <- (Fit[i]-Cor[i])/E[i]
  }
  for( i in (T1+t1p1):(T1+t2p1)) {
    Fit[i] <- pionfit$par[1]*pionfit$par[2]*cosh(fit.mass*((i-T1-1)-Thalf))
    Chi[i] <- (Fit[i]-Cor[i])/E[i]
  }
  for( i in (2*T1+t1p1):(2*T1+t2p1)) {
    Fit[i] <- pionfit$par[1]*pionfit$par[2]*cosh(fit.mass*((i-2*T1-1)-Thalf))
    Chi[i] <- (Fit[i]-Cor[i])/E[i]
  }
  for( i in (3*T1+t1p1):(3*T1+t2p1)) {
    Fit[i] <- pionfit$par[1]*pionfit$par[2]*cosh(fit.mass*((i-3*T1-1)-Thalf))
    Chi[i] <- (Fit[i]-Cor[i])/E[i]
  }

  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]), Time=Time,
              fitdata=data.frame(t=c(0:Thalf), Fit=Fit, Cor=Cor, Err=E, Chi=Chi),
              uwerrresultmps=fit.uwerrm, uwerrresultfps=fit.uwerrf, effmass=pion.eff,
              kappa=kappa, mu=mu, matrix.size=2, no.masses=1)
  attr(res, "class") <- c("cfit", "list")  
  return(invisible(res))

#  return(invisible(data.frame(t=c(0:Thalf), Fit=Fit, Cor=Cor, Err=E, Chi=Chi)))
}

pion2 <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, skip=0) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! time range must be specified!")
  }
  Time <-  2*max(cmicor$V3)
  Thalf <- max(cmicor$V3)
  T1 <- max(cmicor$V3)+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- max(cmicor$V1)
  #positive times

  Skip <- (skip*(T1)*nrObs*4+1)
  Z <- array(cmicor$V4[(Skip):length(cmicor$V4)], 
             dim=c(nrObs*(T1)*4,(length(cmicor$V4[(Skip):length(cmicor$V4)])/(nrObs*(T1)*4))))
  # negative times
  W <- array(cmicor$V5[(Skip):length(cmicor$V5)], 
             dim=c(nrObs*(T1)*4,(length(cmicor$V4[(Skip):length(cmicor$V4)])/(nrObs*(T1)*4))))

  rm(cmicor)
  W <- arrangeCor(T1=T1, W=W, Z=Z)
  rm(Z)

  pion.eff.ll <- effectivemass(from=(t1+1), to=(t2+1), Time, W[1:T1,] , pl=FALSE, S=1.5)
  pion.eff.lf <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(T1+1):(2*T1),] , pl=FALSE, S=1.5)
  pion.eff.ff <- effectivemass(from=(t1+1), to=(t2+1), Time, W[(3*T1+1):(4*T1),] , pl=FALSE, S=1.5) 

  pion.eff <- data.frame(t=pion.eff.ll$t, mll=pion.eff.ll$mass, dmll=pion.eff.ll$dmass,
                         mlf=pion.eff.lf$mass, dmlf=pion.eff.lf$dmass,
                         mff=pion.eff.ff$mass, dmff=pion.eff.ff$dmass)
  
  Cor <- rep(0., times=16*T1)
  E <- rep(1., times=16*T1)
  
  for(i in 1:(16*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
#    E[i] = stderr(W[(i),])
#    E[i] = 1.
  }

  df <- data.frame(x=c((t1):(t2)),
                   y=c(Cor[(t1p1):(t2p1)], Cor[(t1p1+T1):(t2p1+T1)],
                     Cor[(t1p1+2*T1):(t2p1+2*T1)], Cor[(t1p1+3*T1):(t2p1+3*T1)],
                     Cor[(t1p1+4*T1):(t2p1+4*T1)], Cor[(t1p1+5*T1):(t2p1+5*T1)],
                     Cor[(t1p1+6*T1):(t2p1+6*T1)], Cor[(t1p1+7*T1):(t2p1+7*T1)],
#                     Cor[(t1p1+8*T1):(t2p1+8*T1)], Cor[(t1p1+9*T1):(t2p1+9*T1)],
#                     Cor[(t1p1+10*T1):(t2p1+10*T1)], Cor[(t1p1+11*T1):(t2p1+11*T1)],
                     Cor[(t1p1+12*T1):(t2p1+12*T1)], Cor[(t1p1+13*T1):(t2p1+13*T1)],
                     Cor[(t1p1+14*T1):(t2p1+14*T1)], Cor[(t1p1+15*T1):(t2p1+15*T1)]),
                   e=c(E[(t1p1):(t2p1)], E[(t1p1+T1):(t2p1+T1)],
                     E[(t1p1+2*T1):(t2p1+2*T1)], E[(t1p1+3*T1):(t2p1+3*T1)],
                     E[(t1p1+4*T1):(t2p1+4*T1)], E[(t1p1+5*T1):(t2p1+5*T1)],
                     E[(t1p1+6*T1):(t2p1+6*T1)], E[(t1p1+7*T1):(t2p1+7*T1)],
#                     E[(t1p1+8*T1):(t2p1+8*T1)], E[(t1p1+9*T1):(t2p1+9*T1)],
#                     E[(t1p1+10*T1):(t2p1+10*T1)], E[(t1p1+11*T1):(t2p1+11*T1)],
                     E[(t1p1+12*T1):(t2p1+12*T1)], E[(t1p1+13*T1):(t2p1+13*T1)],
                     E[(t1p1+14*T1):(t2p1+14*T1)], E[(t1p1+15*T1):(t2p1+15*T1)])
                   )
  # method="BFGS"
  pionfit <- optim(c(.35, 0.45, 0., 0.7, 0.135), chisqr4, method="BFGS", Thalf=Thalf,
                   x=df$x[1:(t2-t1+1)], y=df$y, err=df$e, tr =(t2-t1+1))

  fit.mass <- abs(pionfit$par[(length(pionfit$par))])
  cat(fit.mass, " ", pionfit$par[1], " ", pionfit$par[2], " ", pionfit$par[3], " ", pionfit$par[4], " ", pionfit$par[5], " \n")
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass^3)*exp(fit.mass*Time/4)
  
  fit.dof <- 2*(t2-t1+1)*3+(t2-t1+1)*4-length(pionfit$par)
  fit.chisqr <- pionfit$value

  if(pl) {
    plot.effmass(m=fit.mass, ll=pion.eff.ll, lf=pion.eff.lf, ff=pion.eff.ff)
  }
  
  ChiPLPL <- rep(0., times=(t2-t1+1))
  ChiPLPF <- rep(0., times=(t2-t1+1))
  ChiPFPF <- rep(0., times=(t2-t1+1))
  ChiPLAL <- rep(0., times=(t2-t1+1))
  ChiPLAF <- rep(0., times=(t2-t1+1))
  ChiPFAL <- rep(0., times=(t2-t1+1))
  ChiPFAF <- rep(0., times=(t2-t1+1))
  ChiALAL <- rep(0., times=(t2-t1+1))
  ChiALAF <- rep(0., times=(t2-t1+1))
  ChiAFAF <- rep(0., times=(t2-t1+1))
  FitPLPL <- rep(0., times=(t2-t1+1))
  FitPLPF <- rep(0., times=(t2-t1+1))
  FitPFPF <- rep(0., times=(t2-t1+1))
  FitPLAL <- rep(0., times=(t2-t1+1))
  FitPLAF <- rep(0., times=(t2-t1+1))
  FitPFAL <- rep(0., times=(t2-t1+1))
  FitPFAF <- rep(0., times=(t2-t1+1))
  FitALAL <- rep(0., times=(t2-t1+1))
  FitALAF <- rep(0., times=(t2-t1+1))
  FitAFAF <- rep(0., times=(t2-t1+1))

  b <- (t2-t1+1)
  for( i in 1:b) {
    FitPLPL[i] <- pionfit$par[1]^2*cosh(fit.mass*(df$x[i]-Thalf))
    ChiPLPL[i] <- (FitPLPL[i]-df$y[i])/df$e[i]
    FitPLPF[i] <- pionfit$par[1]*pionfit$par[2]*cosh(fit.mass*(df$x[i]-Thalf))
    ChiPLPF[i] <- (FitPLPF[i]-df$y[i+1*b])/df$e[i+1*b]
#    FitPFPL[i+2*b] <- pionfit$par[1]*pionfit$par[2]*cosh(fit.mass*(df$x[i]-Thalf))
#    ChiPFPL[i+2*b] <- (Fit[i+2*b]-df$y[i+2*b])/df$e[i+2*b]
    FitPFPF[i] <- pionfit$par[2]*pionfit$par[2]*cosh(fit.mass*(df$x[i]-Thalf))
    ChiPFPF[i] <- (FitPFPF[i]-df$y[i+3*b])/df$e[i+3*b]
    FitPLAL[i] <- -pionfit$par[1]*pionfit$par[3]*sinh(fit.mass*(df$x[i]-Thalf))
    ChiPLAL[i] <- (FitPLAL[i]-df$y[i+4*b])/df$e[i+4*b]
    FitPLAF[i] <- -pionfit$par[1]*pionfit$par[4]*sinh(fit.mass*(df$x[i]-Thalf))
    ChiPLAF[i] <- (FitPLAF[i]-df$y[i+5*b])/df$e[i+5*b]
    FitPFAL[i] <- -pionfit$par[2]*pionfit$par[3]*sinh(fit.mass*(df$x[i]-Thalf))
    ChiPFAL[i] <- (FitPFAL[i]-df$y[i+6*b])/df$e[i+6*b]
    FitPFAF[i] <- -pionfit$par[2]*pionfit$par[4]*sinh(fit.mass*(df$x[i]-Thalf))
    ChiPFAF[i] <- (FitPFPF[i]-df$y[i+7*b])/df$e[i+7*b]
    FitALAL[i] <- pionfit$par[3]*pionfit$par[3]*cosh(fit.mass*(df$x[i]-Thalf))
    ChiALAL[i] <- (FitALAL[i]-df$y[i+8*b])/df$e[i+8*b]
    FitALAF[i] <- pionfit$par[3]*pionfit$par[4]*cosh(fit.mass*(df$x[i]-Thalf))
    ChiALAF[i] <- (FitALAF[i]-df$y[i+9*b])/df$e[i+5*b]
#    FitaAFAL[i+10*b] <- pionfit$par[3]*pionfit$par[4]*cosh(fit.mass*(df$x[i]-Thalf))
#    ChiAFAL[i+10*b] <- (Fit[i+10*b]-df$y[i+10*b])/df$e[i+10*b]
    FitAFAF[i] <- pionfit$par[4]*pionfit$par[4]*cosh(fit.mass*(df$x[i]-Thalf))
    ChiAFAF[i] <- (FitAFAF[i]-df$y[i+11*b])/df$e[i+11*b]
  }


  fit.uwerrm <- uwerrderived(f=fit4by4pcac, data=rbind(W[(t1+1):(t2+1),], W[(T1+t1+1):(T1+t2+1),],
                                              W[(2*T1+t1+1):(2*T1+t2+1),], W[(3*T1+t1+1):(3*T1+t2+1),],
                                              W[(4*T1+t1+1):(4*T1+t2+1),], W[(5*T1+t1+1):(5*T1+t2+1),],
                                              W[(6*T1+t1+1):(6*T1+t2+1),], W[(7*T1+t1+1):(7*T1+t2+1),],
                                              W[(12*T1+t1+1):(12*T1+t2+1),], W[(13*T1+t1+1):(13*T1+t2+1),],
                                              W[(14*T1+t1+1):(14*T1+t2+1),], W[(15*T1+t1+1):(15*T1+t2+1),]),
                             S=S, pl=pl, Time=Time, t1=t1, t2=t2, Err=df$e)

  res <- list(fitresult=pionfit, t1=t1, t2=t2, Time=Time, N=length(W[1,]),
              fitdata=data.frame(t=c(1:b),
                FitPLPL=FitPLPL, CorPLPL=df$y[1:b], ErrPLPL=df$e[1:b], ChiPLPL=ChiPLPL,
                FitPLPL=FitPLPF, CorPLPF=df$y[(b+1):(2*b)], ErrPLPF=df$e[(b+1):(2*b)], ChiPLPF=ChiPLPF,
                FitPFPF=FitPFPF, CorPFPF=df$y[(3*b+1):(4*b)], ErrPFPF=df$e[(3*b+1):(4*b)], ChiPFPF=ChiPFPF,
                FitPLAL=FitPLAL, CorPLAL=df$y[(4*b+1):(5*b)], ErrPLAL=df$e[(4*b+1):(5*b)], ChiPLAL=ChiPLAL,
                FitPLAF=FitPLAF, CorPLAF=df$y[(5*b+1):(6*b)], ErrPLAF=df$e[(5*b+1):(6*b)], ChiPLAF=ChiPLAF,
                FitPFAL=FitPFAL, CorPFAL=df$y[(6*b+1):(7*b)], ErrPFAL=df$e[(6*b+1):(7*b)], ChiPFAL=ChiPFAL,
                FitPFAF=FitPFAF, CorPFAF=df$y[(7*b+1):(8*b)], ErrPFAF=df$e[(7*b+1):(8*b)], ChiPFAF=ChiPFAF,
                FitALAL=FitALAL, CorALAL=df$y[(8*b+1):(9*b)], ErrALAL=df$e[(8*b+1):(9*b)], ChiALAL=ChiALAL,
                FitALAL=FitALAF, CorALAF=df$y[(9*b+1):(10*b)], ErrALAF=df$e[(9*b+1):(10*b)], ChiALAF=ChiALAF,
                FitAFAF=FitAFAF, CorAFAF=df$y[(11*b+1):(12*b)], ErrAFAF=df$e[(11*b+1):(12*b)], ChiAFAF=ChiAFAF
                ),
              uwerrresultmpcac=fit.uwerrm, effmass=pion.eff, mu=mu, kappa=kappa,
              matrix.size=4, no.masses=1)
  attr(res, "class") <- c("cfit", "list")
 return(invisible(res))
}

pion3 <- function(cmicor, mu=0.1, kappa=0.156, t1, t2, S=1.5, pl=FALSE, nrmasses=1) {
  
  if(missing(cmicor)) {
    stop("Error! Data is missing!")
  }
  if(missing(t1) || missing(t2)) {
    stop("Error! time range must be specified!")
  }
  Time <-  2*max(cmicor$V3)
  Thalf <- max(cmicor$V3)
  T1 <- max(cmicor$V3)+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  nrObs <- max(cmicor$V1)
  #positive times
  W <- array(cmicor$V4, dim=c(nrObs*(T1)*4,length(cmicor$V4)/(nrObs*(T1)*4)))
  #negative times
  Z <- array(cmicor$V5, dim=c(nrObs*(T1)*4,length(cmicor$V4)/(nrObs*(T1)*4)))
  rm(cmicor)
  W <- arrangeCor(T1=T1, W=W, Z=Z)
  rm(Z)

  Cor <- rep(0., times=nrObs*T1)
  E <- rep(1., times=nrObs*T1)
  
  for(i in 1:(4*nrObs*T1)) {
    Cor[i] = mean(W[(i),])
    E[i] = uwerrprimary(W[(i),], pl=F)$dvalue
#    E[i] = stderr(W[(i),])
#    E[i] = 1.
  }

  ta <- 4
  tb <- 5 
  C1 <- get6x6matrix(Cor=Cor, t=ta, T1=T1)
  C2 <- get6x6matrix(Cor=Cor, t=tb, T1=T1)

  C3 <- (solve(C1) %*% C2)
  variational <- eigen(C3, symmetric=FALSE, only.values = FALSE, EISPACK=TRUE)
  print(variational)
  rm(C1, C2, C3)
  df <- data.frame(x=c((t1):(t2)),
                   y=c(Cor[(t1p1):(t2p1)], Cor[(t1p1+T1):(t2p1+T1)],
                     Cor[(t1p1+2*T1):(t2p1+2*T1)], Cor[(t1p1+3*T1):(t2p1+3*T1)],
                     Cor[(t1p1+4*T1):(t2p1+4*T1)], Cor[(t1p1+5*T1):(t2p1+5*T1)],
                     Cor[(t1p1+6*T1):(t2p1+6*T1)], Cor[(t1p1+7*T1):(t2p1+7*T1)],
#                     Cor[(t1p1+8*T1):(t2p1+8*T1)], Cor[(t1p1+9*T1):(t2p1+9*T1)],
#                     Cor[(t1p1+10*T1):(t2p1+10*T1)], Cor[(t1p1+11*T1):(t2p1+11*T1)],
                     Cor[(t1p1+12*T1):(t2p1+12*T1)], Cor[(t1p1+13*T1):(t2p1+13*T1)],
                     Cor[(t1p1+14*T1):(t2p1+14*T1)], Cor[(t1p1+15*T1):(t2p1+15*T1)],
                     Cor[(t1p1+16*T1):(t2p1+16*T1)], Cor[(t1p1+17*T1):(t2p1+17*T1)],
                     Cor[(t1p1+18*T1):(t2p1+18*T1)], Cor[(t1p1+19*T1):(t2p1+19*T1)],
                     Cor[(t1p1+20*T1):(t2p1+20*T1)], Cor[(t1p1+21*T1):(t2p1+21*T1)],
                     Cor[(t1p1+22*T1):(t2p1+22*T1)], Cor[(t1p1+23*T1):(t2p1+23*T1)],
#                     Cor[(t1p1+24*T1):(t2p1+24*T1)], Cor[(t1p1+25*T1):(t2p1+25*T1)],
#                     Cor[(t1p1+26*T1):(t2p1+26*T1)], Cor[(t1p1+27*T1):(t2p1+27*T1)],
                     Cor[(t1p1+28*T1):(t2p1+28*T1)], Cor[(t1p1+29*T1):(t2p1+29*T1)],
                     Cor[(t1p1+30*T1):(t2p1+30*T1)], Cor[(t1p1+31*T1):(t2p1+31*T1)],
#                     Cor[(t1p1+32*T1):(t2p1+32*T1)], Cor[(t1p1+33*T1):(t2p1+33*T1)],
#                     Cor[(t1p1+34*T1):(t2p1+34*T1)], Cor[(t1p1+35*T1):(t2p1+35*T1)],
                     ),
                   e=c(E[(t1p1):(t2p1)], E[(t1p1+T1):(t2p1+T1)],
                     E[(t1p1+2*T1):(t2p1+2*T1)], E[(t1p1+3*T1):(t2p1+3*T1)],
                     E[(t1p1+4*T1):(t2p1+4*T1)], E[(t1p1+5*T1):(t2p1+5*T1)],
                     E[(t1p1+6*T1):(t2p1+6*T1)], E[(t1p1+7*T1):(t2p1+7*T1)],
#                     E[(t1p1+8*T1):(t2p1+8*T1)], E[(t1p1+9*T1):(t2p1+9*T1)],
#                     E[(t1p1+10*T1):(t2p1+10*T1)], E[(t1p1+11*T1):(t2p1+11*T1)],
                     E[(t1p1+12*T1):(t2p1+12*T1)], E[(t1p1+13*T1):(t2p1+13*T1)],
                     E[(t1p1+14*T1):(t2p1+14*T1)], E[(t1p1+15*T1):(t2p1+15*T1)],
                     E[(t1p1+16*T1):(t2p1+16*T1)], E[(t1p1+17*T1):(t2p1+17*T1)],
                     E[(t1p1+18*T1):(t2p1+18*T1)], E[(t1p1+19*T1):(t2p1+19*T1)],
                     E[(t1p1+20*T1):(t2p1+20*T1)], E[(t1p1+21*T1):(t2p1+21*T1)],
                     E[(t1p1+22*T1):(t2p1+22*T1)], E[(t1p1+23*T1):(t2p1+23*T1)],
#                     E[(t1p1+24*T1):(t2p1+24*T1)], E[(t1p1+25*T1):(t2p1+25*T1)],
#                     E[(t1p1+26*T1):(t2p1+26*T1)], E[(t1p1+27*T1):(t2p1+27*T1)],
                     E[(t1p1+28*T1):(t2p1+28*T1)], E[(t1p1+29*T1):(t2p1+29*T1)],
                     E[(t1p1+30*T1):(t2p1+30*T1)], E[(t1p1+31*T1):(t2p1+31*T1)],
#                     E[(t1p1+32*T1):(t2p1+32*T1)], E[(t1p1+33*T1):(t2p1+33*T1)],
#                     E[(t1p1+34*T1):(t2p1+34*T1)], E[(t1p1+35*T1):(t2p1+35*T1)]
                     )
                   )
  # method="BFGS"
#  par <- rep(c(.35, 0.45, 0.1, 0.7, 0.1, 0.7, 0.135), times=nrmasses)
  par <- c(0.2211494, 0.2356294, -0.0002134285, 0.00733777, 0.02733800,0.05561932 , 0.1696234, 0.001, 0.00, 0.00, 0.00, 0.00, 0.00, 0.671399)
  pionfit <- optim(par[1:(nrmasses*7)], chisqr6nexp, method="CG", Thalf=Thalf,
                   x=df$x[1:(t2-t1+1)], y=df$y, err=df$e, tr =(t2-t1+1), n=nrmasses)
  print(pionfit)

  fit.mass <- rep(0., times=nrmasses)
  for(i in 1:nrmasses) {
    fit.mass[i] <- abs(pionfit$par[(i*length(pionfit$par)/nrmasses)])
  }
  fit.mass <- sort(fit.mass)
  
  fit.fpi <- 2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(fit.mass[1]^3)*exp(fit.mass[1]*Time/4)

  fit.dof <- 3*(t2-t1+1)*3+3*(t2-t1+1)*4-length(pionfit$par)
  fit.chisqr <- pionfit$value
  cat("mu     = ", mu, "\n")
  cat("kappa  = ", kappa, "\n")
  cat("nrmasses =", nrmasses, "\n")
  cat("Nr of measurements = ", length(W[1,]), "\n")
  cat("fitrange = ", t1, "-", t2, "\n")
  cat("mpi    = ", fit.mass[1], "\n")
  if(nrmasses > 1) {
    for(i in 2:nrmasses) {
      cat("m", i, "=", fit.mass[i], "\n")
    }
  }
  cat("fpi    = ", fit.fpi, "\n")
  cat("m_pcac = ", fit.mass[1]*pionfit$par[3]/pionfit$par[1]/2., "\n")
  cat("chi^2  = ", fit.chisqr, "\n")
  cat("dof    = ", fit.dof, "\n")
  cat("chi^2/dof = ", fit.chisqr/fit.dof, "\n")
  cat("P_L    = ", pionfit$par[1], "\n")
  cat("P_F    = ", pionfit$par[2], "\n")
  cat("A_L    = ", pionfit$par[3], "\n")
  cat("A_F    = ", pionfit$par[4], "\n")
  cat("4_L    = ", pionfit$par[5], "\n")
  cat("4_F    = ", pionfit$par[6], "\n")
  if(nrmasses > 1) {
    for(i in 2:nrmasses) {
      cat("P_L", i, "  = ", pionfit$par[1+(i-1)*7], "\n")
      cat("P_F", i, "  = ", pionfit$par[2+(i-1)*7], "\n")
      cat("A_L", i, "  = ", pionfit$par[3+(i-1)*7], "\n")
      cat("A_F", i, "  = ", pionfit$par[4+(i-1)*7], "\n")
      cat("4_L", i, "  = ", pionfit$par[5+(i-1)*7], "\n")
      cat("4_F", i, "  = ", pionfit$par[6+(i-1)*7], "\n")
    }
  }
  cat("Eigenvalues at times", ta-1, "-", tb-1, ":", sort(-log(abs(variational$values))), "\n")
  
  Chi <- rep(0., times=24*(t2-t1+1))
  Fit <- rep(0., times=24*(t2-t1+1))

  b <- (t2-t1+1)
  for( i in 1:(t2-t1+1)) {
    Fit[i] <- pionfit$par[1]^2*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[i] <- (Fit[i]-df$y[i])/df$e[i]
    Fit[(i+1*b)] <-  pionfit$par[1]*pionfit$par[2]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+1*b)] <-  (Fit[(i+1*b)]-df$y[(i+1*b)])/df$e[(i+1*b)]
    Fit[(i+2*b)] <-  pionfit$par[1]*pionfit$par[2]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+2*b)] <-  (Fit[(i+2*b)]-df$y[(i+2*b)])/df$e[(i+2*b)]
    Fit[(i+3*b)] <-  pionfit$par[2]*pionfit$par[2]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+3*b)] <-  (Fit[(i+3*b)]-df$y[(i+3*b)])/df$e[(i+3*b)]
    Fit[(i+4*b)] <- -pionfit$par[1]*pionfit$par[3]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+4*b)] <-  (Fit[(i+4*b)]-df$y[(i+4*b)])/df$e[(i+4*b)]
    Fit[(i+5*b)] <- -pionfit$par[1]*pionfit$par[4]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+5*b)] <-  (Fit[(i+5*b)]-df$y[(i+5*b)])/df$e[(i+5*b)]
    Fit[(i+6*b)] <- -pionfit$par[2]*pionfit$par[3]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+6*b)] <-  (Fit[(i+6*b)]-df$y[(i+6*b)])/df$e[(i+6*b)]
    Fit[(i+7*b)] <- -pionfit$par[2]*pionfit$par[4]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+7*b)] <-  (Fit[(i+7*b)]-df$y[(i+7*b)])/df$e[(i+7*b)]
    Fit[(i+8*b)] <-  pionfit$par[3]*pionfit$par[3]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+8*b)] <-  (Fit[(i+8*b)]-df$y[(i+8*b)])/df$e[(i+8*b)]
    Fit[(i+9*b)] <-  pionfit$par[3]*pionfit$par[4]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+9*b)] <-  (Fit[(i+9*b)]-df$y[(i+9*b)])/df$e[(i+5*b)]
    Fit[(i+10*b)] <-  pionfit$par[3]*pionfit$par[4]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+10*b)] <-  (Fit[(i+10*b)]-df$y[(i+10*b)])/df$e[(i+10*b)]
    Fit[(i+11*b)] <-  pionfit$par[4]*pionfit$par[4]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+11*b)] <-  (Fit[(i+11*b)]-df$y[(i+11*b)])/df$e[(i+11*b)]
    Fit[(i+12*b)] <-  pionfit$par[5]*pionfit$par[5]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+12*b)] <-  (Fit[(i+12*b)]-df$y[(i+12*b)])/df$e[(i+12*b)]
    Fit[(i+13*b)] <-  pionfit$par[5]*pionfit$par[6]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+13*b)] <-  (Fit[(i+13*b)]-df$y[(i+13*b)])/df$e[(i+13*b)]
    Fit[(i+14*b)] <-  pionfit$par[5]*pionfit$par[6]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+14*b)] <-  (Fit[(i+14*b)]-df$y[(i+14*b)])/df$e[(i+14*b)]
    Fit[(i+15*b)] <-  pionfit$par[6]*pionfit$par[6]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+15*b)] <-  (Fit[(i+15*b)]-df$y[(i+15*b)])/df$e[(i+15*b)]
    Fit[(i+16*b)] <- -pionfit$par[1]*pionfit$par[5]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+16*b)] <-  (Fit[(i+16*b)]-df$y[(i+16*b)])/df$e[(i+16*b)]
    Fit[(i+17*b)] <- -pionfit$par[1]*pionfit$par[6]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+17*b)] <-  (Fit[(i+17*b)]-df$y[(i+17*b)])/df$e[(i+17*b)]
    Fit[(i+18*b)] <- -pionfit$par[2]*pionfit$par[5]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+18*b)] <-  (Fit[(i+18*b)]-df$y[(i+18*b)])/df$e[(i+18*b)]
    Fit[(i+19*b)] <- -pionfit$par[2]*pionfit$par[6]*sinh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+19*b)] <-  (Fit[(i+19*b)]-df$y[(i+19*b)])/df$e[(i+19*b)]
    Fit[(i+20*b)] <-  pionfit$par[3]*pionfit$par[5]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+20*b)] <-  (Fit[(i+20*b)]-df$y[(i+20*b)])/df$e[(i+20*b)]
    Fit[(i+21*b)] <-  pionfit$par[4]*pionfit$par[5]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+21*b)] <-  (Fit[(i+21*b)]-df$y[(i+21*b)])/df$e[(i+21*b)]
    Fit[(i+22*b)] <-  pionfit$par[3]*pionfit$par[6]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+22*b)] <-  (Fit[(i+22*b)]-df$y[(i+22*b)])/df$e[(i+22*b)]
    Fit[(i+23*b)] <-  pionfit$par[4]*pionfit$par[6]*cosh(pionfit$par[7]*(df$x[i]-Thalf))
    Chi[(i+23*b)] <-  (Fit[(i+23*b)]-df$y[(i+23*b)])/df$e[(i+23*b)]
  }

  cat("--- Autocorrelation analysis for m_pcac ---\n")
  cat("    Fitting ", 2*4*(t2-t1+1)*6, " times \n")
  fit.uwerrm <- uwerrderived(f=fit6by6pcac, data=rbind(W[(t1+1):(t2+1),], W[(T1+t1+1):(T1+t2+1),],
                                              W[(2*T1+t1+1):(2*T1+t2+1),], W[(3*T1+t1+1):(3*T1+t2+1),],
                                              W[(4*T1+t1+1):(4*T1+t2+1),], W[(5*T1+t1+1):(5*T1+t2+1),],
                                              W[(6*T1+t1+1):(6*T1+t2+1),], W[(7*T1+t1+1):(7*T1+t2+1),],
                                              W[(12*T1+t1+1):(12*T1+t2+1),], W[(13*T1+t1+1):(13*T1+t2+1),],
                                              W[(14*T1+t1+1):(14*T1+t2+1),], W[(15*T1+t1+1):(15*T1+t2+1),],
                                              W[(16*T1+t1+1):(16*T1+t2+1),], W[(17*T1+t1+1):(17*T1+t2+1),],
                                              W[(18*T1+t1+1):(18*T1+t2+1),], W[(19*T1+t1+1):(19*T1+t2+1),],
                                              W[(20*T1+t1+1):(20*T1+t2+1),], W[(21*T1+t1+1):(21*T1+t2+1),],
                                              W[(22*T1+t1+1):(22*T1+t2+1),], W[(23*T1+t1+1):(23*T1+t2+1),],
                                              W[(28*T1+t1+1):(28*T1+t2+1),], W[(29*T1+t1+1):(29*T1+t2+1),],
                                              W[(30*T1+t1+1):(30*T1+t2+1),], W[(31*T1+t1+1):(31*T1+t2+1),]
                                              ),
                             S=S, pl=pl, Time=Time, t1=t1, t2=t2, Err=df$e, par=pionfit$par, n=nrmasses)
  cat("\nS      = ", S, "\n")
  cat("mpcac  = ", fit.uwerrm$value, "\n")
  cat("dmpcac = ", fit.uwerrm$dvalue, "\n")
  cat("ddmpcac= ", fit.uwerrm$ddvalue, "\n")
  cat("tauint = ", fit.uwerrm$tauint, "\n")
  cat("dtauint= ", fit.uwerrm$dtauint, "\n")
  cat("Wopt   = ", fit.uwerrm$Wopt, "\n")
  res <- list(fitresult=pionfit, t1=t1, t2=t2, N=length(W[1,]),
                        fitdata=data.frame(t=df$x, Fit=Fit, Cor=df$y, Err=df$e, Chi=Chi),
                        uwerrresultmpcac=fit.uwerrm, kappa=kappa, mu=mu, Time=Time,
              matrix.size=6, no.masses=1)
  attr(res, "class") <- c("cfit", "list")
  return(invisible(res))
}

arrangeCor <- function(T1, W, Z) {
#  dim(W2) <- dim(W)
  for(i in 1:(T1)) {
    if(i!=1 && i!=(T1)) {
      # PP for LL, (LF + FL)/2, FF -> symmetric cosh
      W[i,] <- (W[i,] + Z[i,])/2.
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/4.
      W[(i+2*T1),] <- W[(i+T1),]
      W[(i+3*T1),] <- (W[(i+3*T1),] + Z[(i+3*T1),])/2.
      # PA for LL, LF, FL, FF -> antisymmetric -sinh
      W[(i+4*T1),] <- (W[(i+4*T1),] - Z[(i+4*T1),])/2.
      W[(i+5*T1),] <- (W[(i+5*T1),] - Z[(i+5*T1),])/2.
      W[(i+6*T1),] <- (W[(i+6*T1),] - Z[(i+6*T1),])/2.
      W[(i+7*T1),] <- (W[(i+7*T1),] - Z[(i+7*T1),])/2.
      # use again PA -> -sinh
      W[(i+8*T1),] <- W[(i+4*T1),]
      W[(i+9*T1),] <- W[(i+5*T1),]
      W[(i+10*T1),] <- W[(i+6*T1),]
      W[(i+11*T1),] <- W[(i+7*T1),]
      # AA for LL, (LF + FL)/2, FF -> symmetric cosh
      W[(i+12*T1),] <- (W[(i+12*T1),] + Z[(i+12*T1),])/2.
      W[(i+13*T1),] <- (W[(i+13*T1),] + W[(i+14*T1),] + Z[(i+13*T1),] + Z[(i+14*T1),])/4.
      W[(i+14*T1),] <-  W[(i+13*T1),]
      W[(i+15*T1),] <- (W[(i+15*T1),] + Z[(i+15*T1),])/2.
      # 44 for LL, LF + FL/2, FF -> -cosh
      W[(i+16*T1),] <- -(W[(i+16*T1),] + Z[(i+16*T1),])/2.
      W[(i+17*T1),] <- -(W[(i+17*T1),] + W[(i+17*T1),] + Z[(i+18*T1),] + Z[(i+18*T1),])/4.
      W[(i+18*T1),] <-   W[(i+17*T1),]
      W[(i+19*T1),] <- -(W[(i+19*T1),] + Z[(i+19*T1),])/2.
      # P4 for LL, LF, FL, FF -> antisymmetric -sinh
      W[(i+20*T1),] <- (W[(i+20*T1),] - Z[(i+20*T1),])/2.
      W[(i+21*T1),] <- (W[(i+21*T1),] - Z[(i+21*T1),])/2.
      W[(i+22*T1),] <- (W[(i+22*T1),] - Z[(i+22*T1),])/2.
      W[(i+23*T1),] <- (W[(i+23*T1),] - Z[(i+23*T1),])/2.
      # 4P use again P4 -> -sinh
      W[(i+24*T1),] <- W[(i+20*T1),]
      W[(i+25*T1),] <- W[(i+21*T1),]
      W[(i+26*T1),] <- W[(i+22*T1),]
      W[(i+27*T1),] <- W[(i+23*T1),]
      # A4 use 4A -> cosh
      W[(i+28*T1),] <- (W[(i+32*T1),] + Z[(i+32*T1),])/2.
      W[(i+29*T1),] <- (W[(i+34*T1),] + Z[(i+34*T1),])/2.
      W[(i+30*T1),] <- (W[(i+33*T1),] + Z[(i+33*T1),])/2.
      W[(i+31*T1),] <- (W[(i+35*T1),] + Z[(i+35*T1),])/2.
      # 4A for LL, LF, FL, FF -> cosh
      W[(i+32*T1),] <- W[(i+28*T1),]
      W[(i+33*T1),] <- W[(i+29*T1),]
      W[(i+34*T1),] <- W[(i+30*T1),]
      W[(i+35*T1),] <- W[(i+31*T1),]
    }
    else {
      # PP for LL, (LF + FL)/2, FF -> symmetric cosh
      W[i,] <- (W[i,] + Z[i,])
      W[(i+T1),] <- (W[(i+T1),] + W[(i+2*T1),] + Z[(i+T1),] + Z[(i+2*T1),])/2.
      W[(i+2*T1),] <- W[(i+T1),]
      W[(i+3*T1),] <- (W[(i+3*T1),] + Z[(i+3*T1),])
      # PA for LL, LF + FL/2, FF -> antisymmetric -sinh
      W[(i+4*T1),] <- (W[(i+4*T1),] - Z[(i+4*T1),])
      W[(i+5*T1),] <- (W[(i+5*T1),] - Z[(i+5*T1),])
      W[(i+6*T1),] <- (W[(i+6*T1),] - Z[(i+6*T1),])
      W[(i+7*T1),] <- (W[(i+7*T1),] - Z[(i+7*T1),])
      # use again PA -> -sinh
      W[(i+8*T1),] <-  W[(i+4*T1),]
      W[(i+9*T1),] <-  W[(i+5*T1),]
      W[(i+10*T1),] <- W[(i+6*T1),]
      W[(i+11*T1),] <- W[(i+7*T1),]
      # AA for LL, (LF + FL)/2, FF -> symmetric cosh
      W[(i+12*T1),] <- (W[(i+12*T1),] + Z[(i+12*T1),])
      W[(i+13*T1),] <- (W[(i+13*T1),] + W[(i+14*T1),] + Z[(i+13*T1),] + Z[(i+14*T1),])/2.
      W[(i+14*T1),] <-  W[(i+13*T1),]
      W[(i+15*T1),] <- (W[(i+15*T1),] + Z[(i+15*T1),])
      # 44 for LL, LF + FL/2, FF -> -cosh
      W[(i+16*T1),] <- -(W[(i+16*T1),] + Z[(i+16*T1),])
      W[(i+17*T1),] <- -(W[(i+17*T1),] + W[(i+17*T1),] + Z[(i+18*T1),] + Z[(i+18*T1),])/2.
      W[(i+18*T1),] <-   W[(i+17*T1),]
      W[(i+19*T1),] <- -(W[(i+19*T1),] + Z[(i+19*T1),])
      # P4 for LL, LF, FL, FF -> antisymmetric -sinh
      W[(i+20*T1),] <- (W[(i+20*T1),] - Z[(i+20*T1),])
      W[(i+21*T1),] <- (W[(i+21*T1),] - Z[(i+21*T1),])
      W[(i+22*T1),] <- (W[(i+22*T1),] - Z[(i+22*T1),])
      W[(i+23*T1),] <- (W[(i+23*T1),] - Z[(i+23*T1),])
      # 4P use again P4 -> -sinh
      W[(i+24*T1),] <- W[(i+20*T1),]
      W[(i+25*T1),] <- W[(i+21*T1),]
      W[(i+26*T1),] <- W[(i+22*T1),]
      W[(i+27*T1),] <- W[(i+23*T1),]
      # A4 use 4A -> -sinh (swaped sign here...!)
      W[(i+28*T1),] <- (W[(i+32*T1),] + Z[(i+32*T1),])
      W[(i+29*T1),] <- (W[(i+34*T1),] + Z[(i+34*T1),])
      W[(i+30*T1),] <- (W[(i+33*T1),] + Z[(i+33*T1),])
      W[(i+31*T1),] <- (W[(i+35*T1),] + Z[(i+35*T1),])
      # 4A for LL, LF, FL, FF -> antisymmetric -sinh
      W[(i+32*T1),] <- W[(i+28*T1),]
      W[(i+33*T1),] <- W[(i+29*T1),]
      W[(i+34*T1),] <- W[(i+30*T1),]
      W[(i+35*T1),] <- W[(i+31*T1),]
    }
  }
  return(invisible(W))
}

get6x6matrix <- function(Cor, T1, t) {
  C1 <- array(0., dim=c(6,6))
  C1[1,1] = Cor[t]
  C1[1,2] = Cor[(t+T1)]
  C1[2,1] = Cor[(t+2*T1)]
  C1[2,2] = Cor[(t+3*T1)]
  C1[1,3] = Cor[(t+4*T1)]
  C1[3,1] = C1[1,3]
  C1[1,4] = Cor[(t+5*T1)]
  C1[4,1] = C1[1,4]
  C1[2,3] = Cor[(t+6*T1)]
  C1[3,2] = C1[2,3]
  C1[2,4] = Cor[(t+7*T1)]
  C1[4,2] = C1[2,4]
  C1[3,3] = Cor[(t+12*T1)]
  C1[3,4] = Cor[(t+13*T1)]
  C1[4,3] = C1[3,4]
  C1[4,4] = Cor[(t+15*T1)]
  C1[5,5] = Cor[(t+16*T1)]
  C1[5,6] = Cor[(t+17*T1)]
  C1[6,5] = C1[5,6]
  C1[6,6] = Cor[(t+19*T1)]
  C1[1,5] = Cor[(t+20*T1)]
  C1[5,1] = C1[1,5]
  C1[1,6] = Cor[(t+21*T1)]
  C1[6,1] = C1[1,6]
  C1[2,5] = Cor[(t+22*T1)]
  C1[5,2] = C1[2,5]
  C1[2,6] = Cor[(t+23*T1)]
  C1[6,2] = C1[2,6]
  C1[3,5] = Cor[(t+28*T1)]
  C1[5,3] = C1[3,5]
  C1[3,6] = Cor[(t+29*T1)]
  C1[6,3] = C1[3,6]
  C1[4,5] = Cor[(t+30*T1)]
  C1[5,4] = C1[4,5]
  C1[4,6] = Cor[(t+31*T1)]
  C1[6,4] = C1[4,6]
  return(invisible(C1))
}

fit2by2mass <- function(Cor, Err, t1, t2, Time, par=c(1.,0.1,0.12)) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=Err[(t1p1):(t2p1)], e2=Err[(t1p1+T1):(t2p1+T1)],
                   e3=Err[(t1p1+2*T1):(t2p1+2*T1)], e4=Err[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(par, chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  rm(df)
  return(abs(pionfit$par[3]))
}

fit2by2mass.boot <- function(Z, d, Err, t1, t2, Time, kappa, mu) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  W <- t(Z[d,])

  Cor <- rep(0., times=4*T1)
  
  for(i in 1:(T1)) {
    Cor[i] = mean(W[(i),])
  }
  for(i in (T1):(2*T1)) {
    Cor[i] = mean(W[(i),])
    Cor[i+T1] = Cor[i]
  }
  for(i in (3*T1):(4*T1)) {
    Cor[i] = mean(W[(i),])
  }
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=Err[(t1p1):(t2p1)], e2=Err[(t1p1+T1):(t2p1+T1)],
                   e3=Err[(t1p1+2*T1):(t2p1+2*T1)], e4=Err[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(c(1.,0.1,0.12), chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  rm(df)
  return(c(abs(pionfit$par[3]),
           2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(abs(pionfit$par[3])^3)*exp(abs(pionfit$par[3])*Time/4),
           pionfit$value, pionfit$par[1], pionfit$par[2]))
}

fit2by2mass.tsboot <- function(Z, Err, t1, t2, Time, kappa, mu) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  W <- t(Z)

  Cor <- rep(0., times=4*T1)
  
  for(i in 1:(T1)) {
    Cor[i] = mean(W[(i),])
  }
  for(i in (T1):(2*T1)) {
    Cor[i] = mean(W[(i),])
    Cor[i+T1] = Cor[i]
  }
  for(i in (3*T1):(4*T1)) {
    Cor[i] = mean(W[(i),])
  }
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=Err[(t1p1):(t2p1)], e2=Err[(t1p1+T1):(t2p1+T1)],
                   e3=Err[(t1p1+2*T1):(t2p1+2*T1)], e4=Err[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(c(1.,0.1,0.12), chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  rm(df)
  return(c(abs(pionfit$par[3]),
           2*kappa*2*mu/sqrt(2)*abs(pionfit$par[1])/sqrt(abs(pionfit$par[3])^3)*exp(abs(pionfit$par[3])*Time/4),
           pionfit$value, pionfit$par[1], pionfit$par[2]))
}

fit2by2fpi <- function(Cor, Err, t1, t2, Time) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
#  print(Err)
#  print(Cor)
  df <- data.frame(x=c((t1):(t2)), y1=Cor[(t1p1):(t2p1)], y2=Cor[(t1p1+T1):(t2p1+T1)],
                   y3=Cor[(t1p1+2*T1):(t2p1+2*T1)], y4=Cor[(t1p1+3*T1):(t2p1+3*T1)],
                   e1=Err[(t1p1):(t2p1)], e2=Err[(t1p1+T1):(t2p1+T1)],
                   e3=Err[(t1p1+2*T1):(t2p1+2*T1)], e4=Err[(t1p1+3*T1):(t2p1+3*T1)])
  
  pionfit <- optim(c(1.,0.1,0.12), chisqr2, method="BFGS", Thalf=Thalf,
                   x=df$x, y1=df$y1, y2=df$y2, y3=df$y3, y4=df$y4,
                   e1=df$e1, e2=df$e2, e3=df$e3, e4=df$e4)
  rm(df)
  return(abs(pionfit$par[1])/sqrt(abs(pionfit$par[3])^3)*exp(abs(pionfit$par[3])*Time/4))
}

fit4by4pcac <- function(Cor, Err, t1, t2, Time) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  df <- data.frame(x=c((t1):(t2)),
                   y=Cor,
                   e=Err
                   )
  pionfit <- optim(c(.35, 0.45, 0., 0.7, 0.135), chisqr4, method="BFGS", Thalf=Thalf,
                   x=df$x[1:(t2-t1+1)], y=df$y, e=df$e, tr =(t2-t1+1))
  cat(".")
  rm(df)
  return(abs(pionfit$par[5])*pionfit$par[3]/pionfit$par[1]/2.)
}

fit6by6pcac <- function(Cor, Err, t1, t2, Time, par, n) {
  Thalf <- Time/2
  T1 <- Thalf+1
  t1p1 <- (t1+1)
  t2p1 <- (t2+1)
  df <- data.frame(x=c((t1):(t2)),
                   y=Cor,
                   e=Err
                   )
  pionfit <- optim(par, chisqr6nexp, method="BFGS", Thalf=Thalf,
                   x=df$x[1:(t2-t1+1)], y=df$y, e=df$e, tr =(t2-t1+1), n=n)
  cat(".")
  rm(df)
  return(abs(pionfit$par[7])*pionfit$par[3]/pionfit$par[1]/2.)
}

chisqr2 <- function(par, Thalf, x, y1, y2, y3, y4, e1, e2, e3, e4) {
  return((sum((y1-par[1]*par[1]*cosh(abs(par[3])*(x[1:(length(y1))]-Thalf)))^2/e1^2)
          + sum((y2-par[1]*par[2]*cosh(abs(par[3])*(x[1:(length(y1))]-Thalf)))^2/e2^2)
          + sum((y4-par[2]*par[2]*cosh(abs(par[3])*(x[1:(length(y1))]-Thalf)))^2/e4^2)))
}

chisqr4 <- function(par, Thalf, x, y, err, tr) {

  # PP
  Sumall <- sum(((y[1:tr]                     - par[1]*par[1]*cosh(abs(par[5])*(x-Thalf)))/err[1:tr])^2)
  Sumall = Sumall + sum(((y[(tr+1):(2*tr)]    - par[1]*par[2]*cosh(abs(par[5])*(x-Thalf)))/err[(tr+1):(2*tr)])^2)
  Sumall = Sumall + sum(((y[(3*tr+1):(4*tr)]  - par[2]*par[2]*cosh(abs(par[5])*(x-Thalf)))/err[(3*tr+1):(4*tr)])^2)
  # PA
  Sumall = Sumall + sum(((y[(4*tr+1):(5*tr)]  - par[1]*par[3]*sinh(abs(par[5])*(x-Thalf)))/err[(4*tr+1):(5*tr)])^2)
  Sumall = Sumall + sum(((y[(5*tr+1):(6*tr)]  - par[1]*par[4]*sinh(abs(par[5])*(x-Thalf)))/err[(5*tr+1):(6*tr)])^2)
  Sumall = Sumall + sum(((y[(6*tr+1):(7*tr)]  - par[2]*par[3]*sinh(abs(par[5])*(x-Thalf)))/err[(6*tr+1):(7*tr)])^2)
  Sumall = Sumall + sum(((y[(7*tr+1):(8*tr)]  - par[2]*par[4]*sinh(abs(par[5])*(x-Thalf)))/err[(7*tr+1):(8*tr)])^2)
  # AA
  Sumall = Sumall + sum(((y[(8*tr+1):(9*tr)]  - par[3]*par[3]*cosh(abs(par[5])*(x-Thalf)))/err[(8*tr+1):(9*tr)])^2)
  Sumall = Sumall + sum(((y[(9*tr+1):(10*tr)] - par[3]*par[4]*cosh(abs(par[5])*(x-Thalf)))/err[(9*tr+1):(10*tr)])^2)
  Sumall = Sumall + sum(((y[(11*tr+1):(12*tr)]- par[4]*par[4]*cosh(abs(par[5])*(x-Thalf)))/err[(11*tr+1):(12*tr)])^2)
  return(Sumall)
}

fncosh <- function(par, Thalf, x, p1, p2, n, l) {
  s <- rep(0., times=length(x))
  for(i in 0:(n-1)) {
    s = s + par[(p1+i*l)]*par[(p2+i*l)]*cosh(abs(par[(l+i*l)])*(x-Thalf))
  }
  return(s)
}

fnsinh <- function(par, Thalf, x, p1, p2, n, l) {
  s <- rep(0., times=length(x))
  for(i in 0:(n-1)) {
    s = s + par[(p1+i*l)]*par[(p2+i*l)]*sinh(abs(par[(l+i*l)])*(x-Thalf))
  }
  return(s)
}

chisqr6nexp <- function(par, Thalf, x, y, err, tr, n=1) {

  if(length(par) != n*7) {
    stop("wrong number of fit parameters!")
  }
  Sumall <- 0.
  for(i in 0:(n-1)) {
    # PP
    Sumall = Sumall+sum(((y[1:tr]           -fncosh(par, Thalf, x, 1, 1, n, 7))/err[1:tr])^2)
    Sumall = Sumall+sum(((y[(tr+1):(2*tr)]  -fncosh(par, Thalf, x, 1, 2, n, 7))/err[(tr+1):(2*tr)])^2)
    Sumall = Sumall+sum(((y[(3*tr+1):(4*tr)]-fncosh(par, Thalf, x, 2, 2, n, 7))/err[(3*tr+1):(4*tr)])^2)
    # PA
    Sumall = Sumall+sum(((y[(4*tr+1):(5*tr)] -fnsinh(par, Thalf, x, 1, 3, n, 7))/err[(4*tr+1):(5*tr)])^2)
    Sumall = Sumall+sum(((y[(5*tr+1):(6*tr)] -fnsinh(par, Thalf, x, 1, 4, n, 7))/err[(5*tr+1):(6*tr)])^2)
    Sumall = Sumall+sum(((y[(6*tr+1):(7*tr)] -fnsinh(par, Thalf, x, 2, 3, n, 7))/err[(6*tr+1):(7*tr)])^2)
    Sumall = Sumall+sum(((y[(7*tr+1):(8*tr)] -fnsinh(par, Thalf, x, 2, 4, n, 7))/err[(7*tr+1):(8*tr)])^2)
    # AA
    Sumall = Sumall+sum(((y[(8*tr+1):(9*tr)]  -fncosh(par, Thalf, x, 3, 3, n, 7))/err[(8*tr+1):(9*tr)])^2)
    Sumall = Sumall+sum(((y[(9*tr+1):(10*tr)] -fncosh(par, Thalf, x, 3, 4, n, 7))/err[(9*tr+1):(10*tr)])^2)
    Sumall = Sumall+sum(((y[(11*tr+1):(12*tr)]-fncosh(par, Thalf, x, 4, 4, n, 7))/err[(11*tr+1):(12*tr)])^2)
    # 44
    Sumall = Sumall+sum(((y[(12*tr+1):(13*tr)] -fncosh(par, Thalf, x, 5, 5, n, 7))/err[(12*tr+1):(13*tr)])^2)
    Sumall = Sumall+sum(((y[(13*tr+1):(14*tr)] -fncosh(par, Thalf, x, 5, 6, n, 7))/err[(13*tr+1):(14*tr)])^2)
    Sumall = Sumall+sum(((y[(15*tr+1):(16*tr)] -fncosh(par, Thalf, x, 6, 6, n, 7))/err[(15*tr+1):(16*tr)])^2)
    # P4
    Sumall = Sumall+sum(((y[(16*tr+1):(17*tr)] -fnsinh(par, Thalf, x, 1, 5, n, 7))/err[(16*tr+1):(17*tr)])^2)
    Sumall = Sumall+sum(((y[(17*tr+1):(18*tr)] -fnsinh(par, Thalf, x, 1, 6, n, 7))/err[(17*tr+1):(18*tr)])^2)
    Sumall = Sumall+sum(((y[(18*tr+1):(19*tr)] -fnsinh(par, Thalf, x, 2, 5, n, 7))/err[(18*tr+1):(19*tr)])^2)
    Sumall = Sumall+sum(((y[(19*tr+1):(20*tr)] -fnsinh(par, Thalf, x, 2, 6, n, 7))/err[(19*tr+1):(20*tr)])^2)
    # 4A
    Sumall = Sumall+sum(((y[(20*tr+1):(21*tr)] -fncosh(par, Thalf, x, 3, 5, n, 7))/err[(20*tr+1):(21*tr)])^2)
    Sumall = Sumall+sum(((y[(21*tr+1):(22*tr)] -fncosh(par, Thalf, x, 3, 6, n, 7))/err[(21*tr+1):(22*tr)])^2)
    Sumall = Sumall+sum(((y[(22*tr+1):(23*tr)] -fncosh(par, Thalf, x, 4, 5, n, 7))/err[(22*tr+1):(23*tr)])^2)
    Sumall = Sumall+sum(((y[(23*tr+1):(24*tr)] -fncosh(par, Thalf, x, 4, 6, n, 7))/err[(23*tr+1):(24*tr)])^2)
  }
  return(Sumall)
}


